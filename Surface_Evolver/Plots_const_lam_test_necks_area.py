import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Parameters
grid_size = 50
h_zero_min = -0.3
h_zero_max = 0.6
Necks=np.array([1, 2, 3, 4])
lambda_ = 50

interval_number = 2

if interval_number == 0:
    h_zero_min, h_zero_max = -0.1, 0.4
elif interval_number == 1:
    h_zero_min, h_zero_max = -0.2, 0.2
elif interval_number == 2:
    h_zero_min, h_zero_max = -0.1, 0.1
elif interval_number == 3:
    h_zero_min, h_zero_max = -0.025, 0.012
elif interval_number == 4:
    h_zero_min, h_zero_max = -0.075, 0.012


'''
grid_size_array=np.array([50,50,50,50,50,50,50])
if lambda_ == 2:
    grid_size_array=np.array([50,50,50,50,50,50,50]) 
if lambda_ == 5:
    grid_size_array=np.array([50,50,50,50,50,36,50]) 
if lambda_ == 15:
    grid_size_array=np.array([50,50,50,50,50,50,50]) 
if lambda_ == 50:
    grid_size_array=np.array([50,50,50,50,50,50,50]) 

'''
kappa_param = 25

min_array = np.empty(len(Necks))

for i,Neck_nr in enumerate(Necks):
    sigma_array = kappa_param/lambda_**2
    

    #grid_size = grid_size_array[i]
    h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
    csv_file = f"Surface_Energy_csv/{Neck_nr}_energy_{grid_size}_1_{lambda_}_{interval_number}.csv"
    df = pd.read_csv(csv_file, header=None)
    surface_energy = df.to_numpy()

    area = surface_energy / sigma_array

    arg_dist_min = np.nanargmin(area)
    print(Neck_nr,':')
    print(f'Nr_min: {arg_dist_min}')
    print(f'min: {area[arg_dist_min] }')
    print(f'h_zero: {h_zero_array[arg_dist_min]}')
    line, = plt.plot(h_zero_array, area, marker='.', label = f"Neck {Neck_nr}")
    plt.vlines(h_zero_array[arg_dist_min], 0, area[arg_dist_min], 
                  ls='--', color=line.get_color())

    min_array[i] = h_zero_array[arg_dist_min]

print(min_array)
plt.xlabel(r"$H_0 \:/\: \mathrm{nm^{-1}}$")
plt.ylabel(r"Area $\:/\: \mathrm{nm²}$")
plt.ylim(bottom = 1500)
plt.grid()
plt.legend()
plt.savefig(f"Plots/H_0_zoom_{interval_number}_test_necks_area.pdf")

#np.savetxt('Min_H_0.csv',min_array)

