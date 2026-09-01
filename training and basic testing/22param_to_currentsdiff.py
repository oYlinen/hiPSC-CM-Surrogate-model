import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt

from torch.utils.data import random_split, DataLoader
from datasets3 import Dataset3_all
from networks3 import param_to_feature4_20
skip_training = False
import pandas

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device, flush=True)

mean_var_path = "/scratch/svc_td_cbig/olli/pom/run/matlab_50k_cntr/train_and_val2/mean_and_var.csv"
def load_train_data():

    train_path = "/scratch/svc_td_cbig/olli/pom/run/matlab_50k_cntr/train_and_val2"
    train_dataset = Dataset3_all(train_path)
    batch_size = 100
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    print(f"train data set length: {len(train_dataset)}", flush=True)
    print(f"train Batch size = {batch_size}", flush=True)
    return train_loader

def load_validation_data():
    test_path = "/scratch/svc_td_cbig/olli/pom/run/matlab_5k_cntr/validation2"
    test_dataset = Dataset3_all(test_path)
    batch_size = 100
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=True)
    print(f"test data set length: {len(test_dataset)}", flush=True)
    print(f"test Batch size = {batch_size}", flush=True)
    return test_loader

def load_testing_data():
    test_path = "/scratch/svc_td_cbig/olli/pom/run/matlab_10k_cntr_testing/testing"
    test_dataset = Dataset3_all(test_path)
    batch_size = 100
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=True)
    print(f"test data set length: {len(test_dataset)}", flush=True)
    print(f"test Batch size = {batch_size}", flush=True)
    return test_loader


def loss_and_optimizer(model):
    loss = nn.L1Loss()
    lr = 0.001
    print(f"lr = {lr}", flush=True)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    return loss, optimizer


model_20_param = param_to_feature4_20()
model = model_20_param
criterion, optimizer = loss_and_optimizer(model)


model = model.to(device)
criterion = criterion.to(device)

def plotting_function(time, arr, params, orig_arr):
    for idx in range(arr.shape[0]):
        current1 = arr[idx, 0, :]
        current2 = arr[idx, 1, :]
        current3 = arr[idx, 2, :]
        current4 = arr[idx, 3, :]
        current5 = arr[idx, 4, :]

        original1 = orig_arr[idx, 0, :]
        original2 = orig_arr[idx, 1, :]
        original3 = orig_arr[idx, 2, :]
        original4 = orig_arr[idx, 3, :]
        original5 = orig_arr[idx, 4, :]




        fig, axs = plt.subplots(2, 3)
        axs[0, 0].plot(time, original1, label="Original",  linewidth=3)
        axs[0, 0].plot(time, current1, label="Estimated", linewidth=3)
        axs[0, 0].set_title("Vm", fontsize=30)
        axs[0, 0].set_xlabel("Time (s)")
        axs[0, 0].set_ylabel("V")
        axs[0, 0].legend()

        axs[1, 0].plot(time, original3, label="Original", linewidth=3)
        axs[1, 0].plot(time, current3, label="Estimated", linewidth=3)
        axs[1, 0].set_title("Force (AT)", fontsize=30)
        axs[1, 0].set_ylabel("mN / mm2")
        axs[1, 0].set_xlabel("Time (s)")
        axs[1, 0].legend()

        axs[0, 1].plot(time, original2, label="Original", linewidth=3)
        axs[0, 1].plot(time, current2, label="Estimated", linewidth=3)
        axs[0, 1].set_title("Cai", fontsize=30)
        axs[0, 1].set_ylabel("mMol")
        axs[0, 1].set_xlabel("Time (s)")
        axs[0, 1].legend()

        axs[1, 1].plot(time, original4, label="Original", linewidth=3)
        axs[1, 1].plot(time, current4, label="Estimated", linewidth=3)
        axs[1, 1].set_title("ICaL", fontsize=30)
        axs[1, 1].set_ylabel("pA/pF")
        axs[1, 1].set_xlabel("Time (s)")
        axs[1, 1].legend()

        axs[0, 2].plot(time, original5, label="Original", linewidth=3)
        axs[0, 2].plot(time, current5, label="Estimated", linewidth=3)
        axs[0, 2].set_title("IKr", fontsize=30)
        axs[0, 2].set_ylabel("pA/pF")
        axs[0, 2].set_xlabel("Time (s)")
        axs[0, 2].legend()


        fig.tight_layout()
      #  fig.suptitle(str(params[idx]), fontsize=5)
        plt.show()
        plt.plot(time, original1, label="Original", linewidth=2)
        plt.plot(time, current1, label="Estimated", linewidth=2)
        plt.suptitle("Vm", fontsize=20)
        plt.xlabel("Time (s)", fontsize=15)
        plt.ylabel("V", fontsize=15)
        plt.legend()
        plt.show()





def test(model, criterion, test_loader, plot=False):

    model.eval()
    test_loss = []
    maximium = -np.inf
    with torch.no_grad():
       for input, params, _ in test_loader:
            input = input.to(device)  # Move the input to 'device' (CPU or GPU)

            params = params.to(device)
            features = model(params)

            if torch.max(input) > maximium:

                maximium = torch.max(input)
            loss = criterion(features, input)

            test_loss.append(loss.item())


            if plot:
                time = np.linspace(0, 1.5, 2000)
                pred = model.get_original_inputs(features, mean_var_path)
                orig = model.get_original_inputs(input, mean_var_path)
                plotting_function(time, pred, params, orig)

    test_psnr = 20 * torch.log10(maximium) - 10 * np.log10(np.mean(test_loss))

    return np.mean(test_loss), np.mean(test_psnr)


def train(model, optimizer, criterion, train_loader, test_loader, num_epochs=10):

    loss_list = []
    epoch_list = []
    epoch_val_loss = []
    epoch_val_psnr = []
    psnr_list = []

    model.train()
    for epoch in range(num_epochs):
        epoch_loss = []
        epoch_psnr = []
        model.train()

        for input, params, _ in train_loader:

            input = input.to(device)
            params = params.to(device)

            optimizer.zero_grad()
            pred = model(params)

            loss = criterion(pred, input)
            loss *= 1e3
            loss.backward()
            optimizer.step()

            epoch_loss.append(loss.item())
            psnr = (20 * torch.log10(torch.max(input)) - 10 * np.log10(loss.item() /1e3))
            epoch_psnr.append(psnr.detach().cpu())

        loss_list.append(np.mean(epoch_loss))
        psnr_list.append(np.mean(epoch_psnr))
        epoch_list.append(epoch+1)

        test_loss, test_psnr = test(model, criterion, test_loader)
        epoch_val_loss.append(test_loss)
        epoch_val_psnr.append(test_psnr)
        print(f'Train Epoch [{epoch + 1}/{num_epochs}], Train Loss: {np.mean(epoch_loss):.7f}, Train psnr: {np.mean(epoch_psnr):.7f},  Val Loss: {test_loss:.7f}, Val psnr: {test_psnr:.7f}', flush=True)



num_epochs = 1000
skip_training = True


if not skip_training:
    train_loader = load_train_data()
    validation_loader = load_validation_data()

    train_loss = train(model, optimizer, criterion, train_loader, validation_loader, num_epochs)
    torch.save(model.state_dict(), 'temp.pth')
    print("Your trained model is saved successfully!")
    test(model, criterion, validation_loader, True)
else:

    model.load_state_dict(torch.load("1000_1e4_10017_1302.pth", map_location=device))


    print("Loaded weights from your saved model successfully!")

    test_loader = load_testing_data()
    print(test(model, criterion, test_loader))
    print(test(model, criterion, test_loader,True))



















