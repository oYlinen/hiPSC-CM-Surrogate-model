import numpy as np
import matplotlib.pyplot as plt
from networks3 import surrogate_model_cnn
import torch
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device, flush=True)
mean_var_path = "mean_and_var.csv"


surrogate = surrogate_model_cnn()
surrogate.to(device)
surrogate.load_state_dict(torch.load("1000_1e4_10017_1302.pth", map_location=device))

gNa_multiplier = 1 #1
g_f_multiplier = 1 #2
g_CaL_multiplier = 1 #3
g_to_multiplier = 1 #4
g_Ks_multiplier = 1 #5
g_Kr_multiplier = 1 #6
g_K1_multiplier = 1 #7
kNaCa_multiplier = 1 #8
PNaK_multiplier = 1 #9
g_PCa_multiplier = 1 #10
GNaLmax_multiplier = 1 #11
RyRtauadapt_multiplier = 1 #12
RyRtauact_multiplier = 1 #13
RyRtauinact_multiplier = 1 #14
g_irel_max_multiplier = 1 #15
tau_m_multiplier = 1 #16
tau_h_j_multiplier = 1 #17
tau_d_multiplier = 1 #18
tau_f1_2_multiplier = 1 #19
RyRchalf_multiplier = 1 #20
F1_multiplier = 1 #21
F2_multiplier = 1 #22


params = [[gNa_multiplier, g_f_multiplier, g_CaL_multiplier, g_to_multiplier, g_Ks_multiplier, g_Kr_multiplier, g_K1_multiplier,
          kNaCa_multiplier, PNaK_multiplier, g_PCa_multiplier, GNaLmax_multiplier, RyRtauadapt_multiplier, RyRtauact_multiplier,
          RyRtauinact_multiplier, g_irel_max_multiplier, tau_m_multiplier, tau_h_j_multiplier, tau_d_multiplier, tau_f1_2_multiplier,
          RyRchalf_multiplier, F1_multiplier, F2_multiplier]]

params = torch.tensor(params)
print(params.shape)


output_sur = surrogate.forward(params)
output_sur_original_units = surrogate.get_original_inputs(output_sur, mean_var_path)[0]
AP_sur = output_sur_original_units[0,:]
CaT_sur = output_sur_original_units[1,:]
AT_sur = output_sur_original_units[2,:]
ICaL_sur = output_sur_original_units[3,:]
IKr_sur = output_sur_original_units[4,:]
time = np.linspace(0, 1500, 2000)



fig, axs = plt.subplots(2, 3)

axs[0, 0].plot(time, 1000*AP_sur, linewidth=3)
axs[0, 0].set_title("AP", fontsize=20)
axs[0, 0].set_xlabel("Time (s)", fontsize=15)
axs[0, 0].set_ylabel("mV", fontsize=15)

axs[1, 0].plot(time, AT_sur, linewidth=3)
axs[1, 0].set_title("Force (AT)", fontsize=20)
axs[1, 0].set_ylabel("mN / mm2", fontsize=15)
axs[1, 0].set_xlabel("Time (s)", fontsize=15)

axs[0, 1].plot(time, 1000*CaT_sur, linewidth=3)
axs[0, 1].set_title("CaT", fontsize=20)
axs[0, 1].set_ylabel("uMol", fontsize=15)
axs[0, 1].set_xlabel("Time (s)", fontsize=15)

axs[1, 1].plot(time, ICaL_sur, linewidth=3)
axs[1, 1].set_title("ICaL", fontsize=20)
axs[1, 1].set_ylabel("pA/pF", fontsize=15)
axs[1, 1].set_xlabel("Time (s)", fontsize=15)


axs[0, 2].plot(time, IKr_sur, linewidth=3)
axs[0, 2].set_title("IKr", fontsize=20)
axs[0, 2].set_ylabel("pA/pF", fontsize=15)
axs[0, 2].set_xlabel("Time (s)", fontsize=15)

axs[1, 2].axis('off')

fig.tight_layout()
plt.show()