import numpy as np
import scipy.fft
import torch
import torch.nn as nn
import pandas
import scipy


class param_to_feature4_20(nn.Module):
    def __init__(self):
        super().__init__()

        self.layer1 = nn.Sequential(
            nn.Linear(in_features=22, out_features=500),
            nn.Tanh(),
        )
        self.layer2 = nn.Sequential(
            nn.Linear(in_features=500, out_features=500),
            nn.Tanh(),
        )
        self.layer3 = nn.Sequential(
            nn.Linear(in_features=500, out_features=500),
            nn.Tanh(),
        )
        self.layer4 = nn.Sequential(
            nn.Linear(in_features=500, out_features=1000),
            nn.Tanh(),
        )
        self.layer5 = nn.Sequential(
            nn.Linear(in_features=1000, out_features=500),
            nn.Tanh(),
        )
        self.layer6 = nn.Sequential(
            nn.Linear(in_features=500, out_features=500),
            nn.Tanh(),
        )
        self.layer7 = nn.Sequential(
            nn.Linear(in_features=500, out_features=250),
            nn.Tanh(),
        )
        #self.rnn = nn.LSTM(input_size=20, hidden_size=50, num_layers=1, batch_first=True)


        self.dec1 = nn.Sequential(
            nn.ConvTranspose1d(in_channels=1, out_channels=256, kernel_size=4, padding=1, stride=2),
            nn.Tanh(),
        )

        self.dec2 = nn.Sequential(
            nn.ConvTranspose1d(in_channels=256, out_channels=128, kernel_size=4, padding=1, stride=2),
            nn.Tanh(),
        )
        self.dec3 = nn.Sequential(
            nn.ConvTranspose1d(in_channels=128, out_channels=5, kernel_size=4, padding=1, stride=2),
            nn.Tanh(),
        )


        self.output =  nn.Sequential(
            nn.Linear(2000, 2000),
        )

        self.linear_x0 = nn.Sequential(
            nn.Linear(in_features=20, out_features=500),
            nn.Tanh(),
            nn.Linear(in_features=500, out_features=100),
            nn.Tanh(),
            nn.Linear(in_features=100, out_features=4),
        )
        self.confidence = nn.Sequential(
            nn.Linear(in_features=250, out_features=200),
            nn.Tanh(),
            nn.Linear(in_features=200, out_features=50),
            nn.Tanh(),
            nn.Linear(in_features=50, out_features=4),
            nn.Tanh(),
        )

        self.mean_array = np.ndarray([])
        self.variance_array = np.ndarray([])
        self.norm = nn.BatchNorm1d(num_features=20)





    def forward(self, x):
        x = (x-0.5)/(2-0.5)

        layer1 = self.layer1(x)

        layer2 = self.layer2(layer1)
        layer3 = self.layer3(layer2)
        layer4 = self.layer4(layer3)
        layer5 = self.layer5(layer4)
        layer6 = self.layer6(layer5)
        layer7 = self.layer7(layer6)

        dim_addition = layer7[:, None, :] # exand the by 1 dimension to fit to the decoder
        decode1 = self.dec1(dim_addition)

        decode2 = self.dec2(decode1)
        decode3 = self.dec3(decode2)

        output = self.output(decode3)

  #      x0 = self.linear_x0(x)
  #      x0 = x0[:,:, None]
  #      output = torch.cat((output, x0), dim=2)

      #  conf = self.confidence(layer3)
      #  conf = conf[:,:, None]
      #  output = torch.cat((output, conf), dim=2)


        return output
    def get_original_inputs(self, output, path): #TODO think about where this belongs
        if not self.mean_array.any():
            arr = pandas.read_csv(path).to_numpy()
            self.mean_array = arr[0, [0, 1, 2, 3, 4]]
            self.variance_array = arr[1, [0, 1, 2, 3, 4]]

      #  output = np.dstack((output.cpu().detach().numpy(), np.zeros((output.shape[0], output.shape[1], 1500))))

        #output = scipy.fft.idct(output.cpu().detach().numpy(), norm='ortho')
        output = output.cpu().detach().numpy()
        output0_0 = output[:, :, -1]
        output0_0 = output0_0[..., np.newaxis]

       # output0 = np.dstack((output0_0, output[:, :, :-1].cpu().detach().numpy()/100)).cumsum(axis=2)
       # output0_unnormalized = 2*output0 * np.sqrt(self.variance_array[...,np.newaxis]) + self.mean_array[...,np.newaxis]  - np.sqrt(self.variance_array[...,np.newaxis])
        output0_unnormalized = output * np.sqrt(self.variance_array[...,np.newaxis]) + self.mean_array[...,np.newaxis]

        return output0_unnormalized