import pathlib

import sklearn.preprocessing
import torch
import numpy as np
import os
import re
import random
import scipy
import csv
import pandas

import matplotlib.pyplot as plt

def get_params_from_text(text) -> np.ndarray: # copied from POMtool tests
    # Matches to "number.number"
    params = re.findall("[0-9]{1,}\\.[0-9]{1,}", text)

    return np.array(params, dtype=np.float32)



class Dataset3_all():
    def __init__(self, dir_path: str):
        self.dir_path = dir_path
        self.var_mean_file = "C:/Users/imspr/thesis_2025/Thesis_emulator/05-20-50000_22/mean_and_var.csv"
        self.dirs = []
        self.POMtool_output_files = []
        self.nan_files = []
        self.mat_files = scipy.io.loadmat(dir_path + "/all_data.mat")


        self.normalize_factors = {}
        self.norm_flag = False
        self.biomarkers = {}
        self.mean_array = []
        self.variance_array = []
        self.mean_amplitude_array = []
        self.get_normalization_scales()
        self.data = {}
        self.data_params = {}
        self.params_dict = {}
        self.get_params_dict()
        self.filter_mat_keys()
        self.data_size = len(self.mat_files.keys())





    def get_normalization_scales(self):  # TODO get this from elsewhere
        if self.norm_flag:
            return
        if pathlib.Path(self.var_mean_file).exists():
            arr = pandas.read_csv(self.var_mean_file).to_numpy()
            self.mean_array = arr[0,:]
            self.variance_array = arr[1,:]
            return

        self.norm_flag = True
        mean_array = []
        variance_array = []

        for key in self.mat_files.keys():
            if key[:4] != "cell": #there are a few extra keys from matlab
                continue
            cell_mat_files = self.mat_files[key][0][0][-1]


            idx = 0
            for signal in cell_mat_files.T:
                mean_signal = np.mean(signal)
                variance_signal = np.var(signal)

                try:
                    variance_array[idx] = (variance_array[idx] + variance_signal)
                    mean_array[idx] = (mean_array[idx] + mean_signal)
                except:
                    variance_array.append(variance_signal)
                    mean_array.append(mean_signal)

                idx +=1

        self.variance_array = np.array(variance_array) / len(self.mat_files.keys())
        self.mean_array = np.array(mean_array) / len(self.mat_files.keys())
        self.save_mean_var()


    def save_mean_var(self):
        file_name = self.var_mean_file
        combined = np.vstack((self.mean_array, self.variance_array))

        df = pandas.DataFrame(combined)
        df.to_csv(file_name, index=False)


    def normalize_data(self, input, file): # normalize the conducatnces as well?
        arr = np.zeros((input.shape[0], input.shape[1]))
        idxs = [0, 1, 2, 3, 4]
        for idx in range(input.shape[0]):
            i = idxs[idx]
            mean_input = self.mean_array[i]
            variance_input = self.variance_array[i]

           # mean_input = np.mean(input[idx,:])
          #  variance_input = np.var(input[idx, :])


            normalized = ((input[idx,:] - mean_input)/ np.sqrt(variance_input))
           # normalized = ((input[idx, :] - mean_input) / (
            #            variance_input - mean_input))


            diff = np.diff(normalized,n=1)*100
            diff_x0 = np.append(diff, [normalized[0]])
       #     dct = scipy.fft.dct(normalized, type=2, norm='ortho')

           # arr[idx, :] = np.append(diff, [normalized[0]])

            arr[idx, :] = normalized

        return torch.tensor(arr).float()


    def get_biomarkers(self):
        file = self.dir_path + "/biomarkers.csv"
        biomarkers = pandas.read_csv(file)
        biomarkers = biomarkers.to_dict()
        self.biomarkers = biomarkers.pop("directory", None)


    def get_params_dict(self):
        file = self.dir_path + "/simulation_manifest.csv"
        temp = pandas.read_csv(file, header=None)
        temp = temp.to_numpy()

        for arr in temp:

            values = get_params_from_text(str(arr[0:]))
            self.params_dict[str(arr[0])] = torch.tensor(values).float()


    def __len__(self):
        return self.data_size

    def filter_mat_keys(self):
        new_dict = self.mat_files.copy()

        for key in self.mat_files.keys():
            if key[:4] != "cell":
                new_dict.pop(key)
                continue
            cell_mat_files = self.mat_files[key][0][0][-1]
            time = self.mat_files[key][0][0][-2]
            ind = 802.5
            ind_10s_start = np.argmin(time < ind)
            ind_10s_end = np.argmin(time < ind + 1.5)
            ind_10s_high = np.argmin(time < ind + 0.3)
            idx = 0
            for signal in cell_mat_files.T:
                if idx == 0:
                    part = signal[ind_10s_start:ind_10s_end] + 10
                    part_ = signal[ind_10s_start:ind_10s_end]
                    part_start = signal[ind_10s_start:ind_10s_high]
                    maximium_signal = np.max(signal)
                    maximium_start = np.max(part_start)

                    diff_sum = np.sum(np.diff(part_))

                    if part_[0] > -0.052 or maximium_signal < 0 or maximium_signal*0.9 >= maximium_start:
                       # print(f"pop: {key}")
                        new_dict.pop(key)
                        idx +=1
                        continue


                idx += 1
        self.mat_files = new_dict

    def __getitem__(self, idx):
        # TODO save the already used ones? or consider it
        # TODO clean this code and make it nice
        # TODO optimization <- use gpu? Use torch interp
        #TODO check if matlab gives error or the output is garbage?
        #TODO version control

        keys = list(self.mat_files.keys())
        dir = keys[idx]
     #   dir = random.choice(keys)
        file = dir + "/res.mat"

        params = self.params_dict[dir] # 2 was the maximium

        if dir in self.data:
            input_data = self.data[dir]
            return input_data, params, self.normalize_factors

        vals = self.mat_files[dir][0][0][-1]
        time = self.mat_files[dir][0][0][-2]

        ind = 802.4 #TODO clean Dataset so that this feature can be used (cond to AP, random here can not be there :D)
        ind_10s_start = np.argmin(time < ind)
        ind_10s_end = np.argmin(time < ind+1.5)

        time_10s = time[ind_10s_start:ind_10s_end].flatten()
        fs = 1000
        time_new = np.linspace(np.min(time_10s), np.max(time_10s),num=int(2 * fs))  # resampling to constant frequancy of 1/fs
        idx = 0
        idxs = [0,1, 2, 3, 4]
        input_data = []
        for signal in vals.T:
            # Choosing random 10s from the sample
            if idx not in idxs:
                idx += 1
                continue
            signal_10s = signal[ind_10s_start:ind_10s_end]


            # Resampling, but matching the times since the ODE does not give equal number of points for each AP
            signal_new = np.interp(time_new, time_10s, signal_10s)
            input_data.append(signal_new)

            idx +=1


        input_data = np.array(input_data)
        input_data = self.normalize_data(input_data, file)
        self.data[dir] = input_data



        self.time_new = time_new
        return input_data, params, self.normalize_factors