import numpy as np
import matplotlib.pyplot as plt
from networks3 import param_to_feature4_20
import torch
import os

from matplotlib.widgets import Slider, Button, RadioButtons


os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device, flush=True)
mean_var_path = "mean_and_var.csv"


emulator = param_to_feature4_20()
emulator.to(device)
emulator.load_state_dict(torch.load("1000_1e4_10017_1302.pth", map_location=device))

input_params = ["gNa", "g_f", "g_CaL", "g_to", "g_Ks", "g_Kr", "g_K1", "kNaCa", "PNaK", "g_PCa", "GNaLmax", "RyRtauadapt",
                "RyRtauact", "RyRtauinact", "g_irel_max", "tau_m", "tau_h_j", "tau_d", "tau_f1_2", "RyRchalf", "F1", "F2"]

print(len(input_params))
params0 = torch.ones((1, len(input_params)))
print(params0.shape)


output_emu = emulator(params0)
output_emu_original_units = emulator.get_original_inputs(output_emu, mean_var_path)[0]
AP_emu0 = 1000*output_emu_original_units[0,:]
CaT_emu0 = 1000*output_emu_original_units[1,:]
AT_emu0 = output_emu_original_units[2,:]
ICaL_emu0 = output_emu_original_units[3,:]
IKr_emu0 = output_emu_original_units[4,:]
time = np.linspace(0, 1500, 2000)



fig, axs = plt.subplots(2, 3, figsize=(10, 20))

line0, = axs[0, 0].plot(time, AP_emu0, linewidth=3)
axs[0, 0].set_title("AP", fontsize=20)
axs[0, 0].set_xlabel("Time (ms)", fontsize=15)
axs[0, 0].set_ylabel("mV", fontsize=15)

line1, = axs[1, 0].plot(time, AT_emu0, linewidth=3)
axs[1, 0].set_title("Force (AT)", fontsize=20)
axs[1, 0].set_ylabel("mN / mm2", fontsize=15)
axs[1, 0].set_xlabel("Time (ms)", fontsize=15)

line2, = axs[0, 1].plot(time, CaT_emu0, linewidth=3)
axs[0, 1].set_title("CaT", fontsize=20)
axs[0, 1].set_ylabel("uMol", fontsize=15)
axs[0, 1].set_xlabel("Time (ms)", fontsize=15)

line3, = axs[1, 1].plot(time, ICaL_emu0, linewidth=3)
axs[1, 1].set_title("ICaL", fontsize=20)
axs[1, 1].set_ylabel("pA/pF", fontsize=15)
axs[1, 1].set_xlabel("Time (ms)", fontsize=15)


line4, = axs[0, 2].plot(time, IKr_emu0, linewidth=3)
axs[0, 2].set_title("IKr", fontsize=20)
axs[0, 2].set_ylabel("pA/pF", fontsize=15)
axs[0, 2].set_xlabel("Time (ms)", fontsize=15)

axs[1, 2].axis('off')

init_multiplier = 1

color = 'lightgoldenrodyellow'
sliders = []

y_cord = 0

for input_name in input_params:
    slider_ax = fig.add_axes([0.70, y_cord, 0.20, 0.03], facecolor=color)
    slider = Slider(slider_ax, input_name, 0.5, 2.0, valinit=init_multiplier)
    sliders.append(slider)
    y_cord += 0.02



def sliders_on_changed(val):
    params = params0
    index = 0
    for slider in sliders:
        params[0, index] = slider.val
        index +=1

    output_emu = emulator(params)
    output_emu_original_units = emulator.get_original_inputs(output_emu, mean_var_path)[0]
    AP_emu = 1000*output_emu_original_units[0, :]
    CaT_emu = 1000*output_emu_original_units[1, :]
    AT_emu = output_emu_original_units[2, :]
    ICaL_emu = output_emu_original_units[3, :]
    IKr_emu = output_emu_original_units[4, :]

    line0.set_ydata(AP_emu)
    line1.set_ydata(AT_emu)
    line2.set_ydata(CaT_emu)
    line3.set_ydata(ICaL_emu)
    line4.set_ydata(IKr_emu)

    axs[0, 0].set_ylim(min(np.min(AP_emu0), np.min(AP_emu))*1.1, max(np.max(AP_emu0), np.max(AP_emu))*1.1)
    axs[1, 0].set_ylim(min(np.min(AT_emu0), np.min(AT_emu), 0), max(np.max(AT_emu0), np.max(AT_emu))*1.1)
    axs[0, 1].set_ylim(min(np.min(CaT_emu0), np.min(CaT_emu), 0), max(np.max(CaT_emu0), np.max(CaT_emu))*1.1)
    axs[1, 1].set_ylim(min(np.min(ICaL_emu), np.min(ICaL_emu))*1.1, max(np.max(ICaL_emu), np.max(ICaL_emu))+0.05)
    axs[0, 2].set_ylim(min(np.min(IKr_emu), np.min(IKr_emu), 0), max(np.max(IKr_emu), np.max(IKr_emu))*1.1)



    fig.canvas.draw_idle()

for slider in sliders:
    slider.on_changed(sliders_on_changed)


reset_button_ax = fig.add_axes([0.05, 0.015, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=color, hovercolor='0.9')
def reset_button_on_clicked(mouse_event):
    for slider in sliders:
        slider.reset()
reset_button.on_clicked(reset_button_on_clicked)




random_button_ax = fig.add_axes([0.25, 0.015, 0.1, 0.04])
random_button = Button(random_button_ax, 'Random', color=color, hovercolor='0.9')
def random_button_on_clicked(mouse_event):
    params = params0
    index = 0
    for slider in sliders:
        random_value = np.random.random()*1.5+0.5
        params[0, index] = random_value
        slider.set_val(random_value)
        index += 1

    output_emu = emulator(params)
    output_emu_original_units = emulator.get_original_inputs(output_emu, mean_var_path)[0]
    AP_emu = 1000*output_emu_original_units[0, :]
    CaT_emu = 1000*output_emu_original_units[1, :]
    AT_emu = output_emu_original_units[2, :]
    ICaL_emu = output_emu_original_units[3, :]
    IKr_emu = output_emu_original_units[4, :]

    line0.set_ydata(AP_emu)
    line1.set_ydata(AT_emu)
    line2.set_ydata(CaT_emu)
    line3.set_ydata(ICaL_emu)
    line4.set_ydata(IKr_emu)

    axs[0, 0].set_ylim(min(np.min(AP_emu0), np.min(AP_emu))*1.1, max(np.max(AP_emu0), np.max(AP_emu))*1.1)
    axs[1, 0].set_ylim(min(np.min(AT_emu0), np.min(AT_emu), 0), max(np.max(AT_emu0), np.max(AT_emu))*1.1)
    axs[0, 1].set_ylim(min(np.min(CaT_emu0), np.min(CaT_emu), 0), max(np.max(CaT_emu0), np.max(CaT_emu))*1.1)
    axs[1, 1].set_ylim(min(np.min(ICaL_emu), np.min(ICaL_emu))*1.1, max(np.max(ICaL_emu), np.max(ICaL_emu))+0.05)
    axs[0, 2].set_ylim(min(np.min(IKr_emu), np.min(IKr_emu), 0), max(np.max(IKr_emu), np.max(IKr_emu))*1.1)



    fig.canvas.draw_idle()


random_button.on_clicked(random_button_on_clicked)


plt.show()






































