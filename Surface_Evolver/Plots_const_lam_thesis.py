import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Parameters
grid_size = 50
h_zero_min = -0.3
h_zero_max = 0.6
Necks=np.array([334, 574])
lambda_array = np.array([2,5,15,50])

#0: h_zero[-0.1, 0.4] -> for the region between bifork
#1: h_zero[-0.2,0.2] -> for energy minimum
#2: h_zero[-0.1,0.1]

Neck_nr = Necks[0]


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)


for ax, Neck_nr in zip(axes, Necks):
    for i, lambda_ in enumerate(lambda_array):
        h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
        csv_file = f"Distances_csv/{Neck_nr}_distances_{grid_size}_1_{lambda_}.csv"
        df = pd.read_csv(csv_file, header=None)
        grid_data = df.to_numpy()


        distances = grid_data[:, 0]

        valid = ~np.isnan(distances)
        arg_dist_min = np.nanargmin(distances)
        
        # Plot line and get color
        line, = ax.plot(h_zero_array[valid], distances[valid],marker='.', label = rf'$\lambda = {lambda_} \, \mathrm{{nm}}$')
        
        print(Neck_nr, ':')
        print(f'Nr_min: {arg_dist_min}')
        print(f'min: {distances[arg_dist_min]}')
        print(f'h_zero: {h_zero_array[arg_dist_min]}')

    ax.set_title(rf"Neck {Neck_nr}")
    ax.set_xlabel(r"$H_0 \:/\: \mathrm{nm^{-1}}$")
    ax.set_ylim(bottom=0)
    ax.grid()
    ax.legend()

axes[0].set_ylabel(r"Integrated Distance $\:/\: \mathrm{nm}$")
plt.tight_layout()
plt.savefig("Plots/distance_const_lam_smooth.pdf", format='pdf')


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)


for ax, Neck_nr in zip(axes, Necks):
    for i, lambda_ in enumerate(lambda_array):
        h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
        csv_file = f"Energy_csv/{Neck_nr}_energy_{grid_size}_1_{lambda_}.csv"
        df = pd.read_csv(csv_file, header=None)
        grid_data = df.to_numpy()

        energy = grid_data[:, 0]


        vmin = 0
        vmax = 15000
        energy = np.where((energy >= vmin) & (energy <= vmax), energy, np.nan)
        valid = ~np.isnan(energy)
        arg_dist_min = np.nanargmin(energy)
        
        # Plot line and get color
        line, = ax.plot(h_zero_array[valid], energy[valid],marker='.', label=rf"$\lambda =$ {lambda_}")
        ax.vlines(h_zero_array[arg_dist_min], 0, energy[arg_dist_min], 
                  ls='--', color=line.get_color())
        
        print(Neck_nr, ':')
        print(f'Nr_min: {arg_dist_min}')
        print(f'min: {energy[arg_dist_min]}')
        print(f'h_zero: {h_zero_array[arg_dist_min]}')

    ax.set_title(rf"Neck {Neck_nr}")
    ax.set_xlabel(r"$H_0 \:/\: \mathrm{nm^{-1}}$")
    ax.set_ylim(bottom=0)
    ax.grid()
    ax.legend()

axes[0].set_ylabel(r"Energy $\:/\: \mathrm{k_BT}$")
plt.tight_layout()
plt.savefig("Plots/energy_const_lam_smooth.pdf", format='pdf')



fig = plt.figure()

Neck_nr = 1


for i, lambda_ in enumerate(lambda_array):
    h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
    csv_file = f"Distances_csv/{Neck_nr}_distances_{grid_size}_1_{lambda_}.csv"
    df = pd.read_csv(csv_file, header=None)
    grid_data = df.to_numpy()


    distances = grid_data[:, 0]

    valid = ~np.isnan(distances)
    arg_dist_min = np.nanargmin(distances)
    
    # Plot line and get color
    line, = plt.plot(h_zero_array[valid], distances[valid], marker='.', label = rf'$\lambda = {lambda_} \, \mathrm{{nm}}$')
    
    print(Neck_nr, ':')
    print(f'Nr_min: {arg_dist_min}')
    print(f'min: {distances[arg_dist_min]}')
    print(f'h_zero: {h_zero_array[arg_dist_min]}')

plt.title(rf"Neck {Neck_nr}")
plt.xlabel(r"$H_0 \:/\: \mathrm{nm^{-1}}$")
plt.ylim(bottom=0)
plt.grid()
plt.legend()

plt.ylabel(r"Integrated Distance $\:/\: \mathrm{nm}$")
plt.tight_layout()
plt.savefig("Plots/distance_const_lam_Neck1.pdf", format='pdf')