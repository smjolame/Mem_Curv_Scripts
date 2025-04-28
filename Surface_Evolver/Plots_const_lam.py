import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Parameters
grid_size = 50
h_zero_min = -0.3
h_zero_max = 0.6
Necks=np.array([331, 334, 335, 363, 574])
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


min_array = np.empty(len(Necks))

for i,Neck_nr in enumerate(Necks):

    #grid_size = grid_size_array[i]
    h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
    csv_file = f"Distances_csv/{Neck_nr}_distances_{grid_size}_1_{lambda_}_{interval_number}.csv"
    df = pd.read_csv(csv_file, header=None)
    grid_data = df.to_numpy()

    distances = grid_data[:,0]
    arg_dist_min = np.nanargmin(distances)
    print(Neck_nr,':')
    print(f'Nr_min: {arg_dist_min}')
    print(f'min: {distances[arg_dist_min] }')
    print(f'h_zero: {h_zero_array[arg_dist_min]}')
    if interval_number==3:
        line, = plt.plot(h_zero_array, distances, marker='.', label = rf"Neck {Neck_nr}, $H_0^\text{{est}}={h_zero_array[arg_dist_min]:.3f} \, \mathrm{{nm^{{-1}}}}$")
    else:
        line, = plt.plot(h_zero_array, distances, marker='.', label = f"Neck {Neck_nr}")
    plt.vlines(h_zero_array[arg_dist_min], 0, distances[arg_dist_min], 
                  ls='--', color=line.get_color())

    min_array[i] = h_zero_array[arg_dist_min]
print(min_array)
plt.axvline(-0.025, ls=':', alpha = 0.8,color='grey', lw='3')
plt.axvline(0.012, ls=':', alpha = 0.8,color='grey', lw='3')
plt.xlabel(r"$H_0 \:/\: \mathrm{nm^{-1}}$")
plt.ylabel(r"Integrated Distance $\:/\: \mathrm{nm}$")
plt.ylim(bottom = 0)
plt.grid()
plt.legend()
plt.savefig(f"Plots/H_0_zoom_{interval_number}.pdf")

#np.savetxt('Min_H_0.csv',min_array)

plt.clf()
for i,Neck_nr in enumerate(Necks):

    #grid_size = grid_size_array[i]
    h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
    csv_file = f"Energy_csv/{Neck_nr}_energy_{grid_size}_1_{lambda_}_{interval_number}.csv"
    df = pd.read_csv(csv_file, header=None)
    grid_data = df.to_numpy()

    energy = grid_data[:,0]
    vmin = 0
    vmax = 15000
    energy = np.where((energy >= vmin) & (energy <= vmax), energy, np.nan)
    valid = ~np.isnan(energy)
    arg_dist_min = np.nanargmin(energy[:26])
    print(Neck_nr,':')
    print(f'Nr_min: {arg_dist_min}')
    print(f'min: {energy[arg_dist_min] }')
    print(f'h_zero: {h_zero_array[arg_dist_min]}')
    line, = plt.plot(h_zero_array[valid], energy[valid], marker='.', label = f"Neck {Neck_nr}")
    plt.vlines(h_zero_array[arg_dist_min], 0, energy[arg_dist_min], 
                  ls='--', color=line.get_color())

    min_array[i] = h_zero_array[arg_dist_min]
print(min_array)

plt.xlabel(r"$H_0 \:/\: \mathrm{nm^{-1}}$")
plt.ylabel(r"Energy $\:/\: \mathrm{k_BT}$")
plt.ylim(bottom = 0)
plt.grid()
plt.legend()
plt.savefig(f"Plots/H_0_zoom_{interval_number}_energy.pdf")

#[-0.01428571 -0.01428571 -0.01428571  0.00204082 -0.01020408]
#[-0.01518367 -0.01216327 -0.01518367  0.00369388 -0.00914286]