import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Parameters
grid_size = 50
h_zero_min = -0.3
h_zero_max = 0.6
Necks=np.array([1, 2, 3, 4])
lambda_ = 10

interval_number = 4

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
    line, = plt.plot(h_zero_array, distances, marker='.', label = f"Neck {Neck_nr}")
    plt.vlines(h_zero_array[arg_dist_min], 0, distances[arg_dist_min], 
                  ls='--', color=line.get_color())

    min_array[i] = h_zero_array[arg_dist_min]

print(min_array)
plt.xlabel(r"$H_0 \:/\: \mathrm{nm^{-1}}$")
plt.ylabel(r"Integrated Distance $\:/\: \mathrm{nm}$")
plt.ylim(bottom = 0)
plt.grid()
plt.legend()
plt.savefig(f"Plots/H_0_zoom_{interval_number}_test_necks.pdf")

#np.savetxt('Min_H_0.csv',min_array)

plt.clf()
for i,Neck_nr in enumerate(Necks):

    #grid_size = grid_size_array[i]
    h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
    csv_file = f"Energy_csv/{Neck_nr}_energy_{grid_size}_1_{lambda_}_{interval_number}.csv"
    df = pd.read_csv(csv_file, header=None)
    grid_data = df.to_numpy()

    energy = grid_data[:,0]
    vmin = 200
    vmax = 1000
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
plt.ylim(bottom = 400)
plt.grid()
plt.legend()
plt.savefig(f"Plots/H_0_zoom_{interval_number}_energy_test_necks.pdf")




# 1,2,3,4:  
# lambda = 10: [-0.01995918 -0.01108163 -0.04836735 -0.05191837]
# lambda = 50: [-0.01285714 -0.00930612 -0.06079592 -0.05191837]