import numpy as np
import matplotlib.pyplot as plt
from networks3 import param_to_feature4_20
import torch
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device, flush=True)
mean_var_path = "mean_and_var.csv"


emulator = param_to_feature4_20()
emulator.to(device)
emulator.load_state_dict(torch.load("1000_1e4_10017_1302.pth", map_location=device))



params = torch.ones((1,22))
print(params.shape)


output_emu = emulator(params)
output_emu_original_units = emulator.get_original_inputs(output_emu, mean_var_path)[0]
AP_emu = output_emu_original_units[0,:]
CaT_emu = output_emu_original_units[1,:]
AT_emu = output_emu_original_units[2,:]
ICaL_emu = output_emu_original_units[3,:]
IKr_emu = output_emu_original_units[4,:]
time = np.linspace(0, 1500, 2000)



fig, axs = plt.subplots(2, 3)

axs[0, 0].plot(time, AP_emu, linewidth=3)
axs[0, 0].set_title("AP", fontsize=20)
axs[0, 0].set_xlabel("Time (s)", fontsize=15)
axs[0, 0].set_ylabel("mV", fontsize=15)

axs[1, 0].plot(time, AT_emu, linewidth=3)
axs[1, 0].set_title("Force (AT)", fontsize=20)
axs[1, 0].set_ylabel("mN / mm2", fontsize=15)
axs[1, 0].set_xlabel("Time (s)", fontsize=15)

axs[0, 1].plot(time, CaT_emu, linewidth=3)
axs[0, 1].set_title("CaT", fontsize=20)
axs[0, 1].set_ylabel("uMol", fontsize=15)
axs[0, 1].set_xlabel("Time (s)", fontsize=15)

axs[1, 1].plot(time, ICaL_emu, linewidth=3)
axs[1, 1].set_title("ICaL", fontsize=20)
axs[1, 1].set_ylabel("pA/pF", fontsize=15)
axs[1, 1].set_xlabel("Time (s)", fontsize=15)


axs[0, 2].plot(time, IKr_emu, linewidth=3)
axs[0, 2].set_title("IKr", fontsize=20)
axs[0, 2].set_ylabel("pA/pF", fontsize=15)
axs[0, 2].set_xlabel("Time (s)", fontsize=15)

axs[1, 2].axis('off')

fig.tight_layout()
plt.show()