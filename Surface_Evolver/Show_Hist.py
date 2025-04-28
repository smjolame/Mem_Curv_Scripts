import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker
# Use LaTeX for rendering math expressions
#mpl.rcParams['text.usetex'] = True
#
## Set font properties
#mpl.rcParams['font.family'] = 'serif'  # Serif for the main text
#mpl.rcParams['font.serif'] = ['Times New Roman']  # This can be customized to another serif font
mpl.rcParams['font.size'] = 17  # Set the font size
#
## Set math font to match "Latin Modern Math" used in your LaTeX document
#mpl.rcParams['mathtext.fontset'] = 'custom'
#mpl.rcParams['mathtext.rm'] = 'serif'  # Use serif for regular math text (this matches the font used for text)
#mpl.rcParams['mathtext.it'] = 'serif:italic'  # Italic for math symbols
#mpl.rcParams['mathtext.bf'] = 'serif:bold'  # Bold math symbols

# Parameters
grid_size = 50
h_zero_min = -0.3
h_zero_max = 0.6
lambda_min = 1
lambda_max = 50
kappa_param = 25 
Neck_nr = 574

point_interest = {
    1: [[9,5],[20,26],[34,34],[34,14],[34,35]],  
    5: [[9,5],[20,26],[34,34],[34,16],[34,35]],  
    334: [[9,5],[20,26],[34,33],[34,16],[34,34]],  
    574: [[9,5],[20,24],[34,29],[34,16],[34,30]]   
}


# Create the arrays for h_zero and lambda
h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
#h_zero_array = np.logspace(np.log10(h_zero_min), np.log10(h_zero_max), grid_size)
lambda_array = lambda_min * np.power(lambda_max / lambda_min, np.linspace(0, 1, grid_size))
sigma_array = kappa_param/lambda_array**2
lambda_target = np.array([2, 5, 15, 50])  # The λ values where you want vertical lines
lambda_indices = [np.abs(lambda_array - val).argmin() for val in lambda_target]

csv_file_dist = f"Distances_csv/{Neck_nr}_distances_{grid_size}_2_smooth.csv"
csv_file_energy = f"Energy_csv/{Neck_nr}_energy_{grid_size}_2_smooth.csv"
csv_file_curv_energy = f"Curv_Energy_csv/{Neck_nr}_energy_{grid_size}_2_smooth.csv"
csv_file_surface_energy = f"Surface_Energy_csv/{Neck_nr}_energy_{grid_size}_2_smooth.csv"

#csv_file_dist = f"Distances_csv/{Neck_nr}_distances_{grid_size}_2.csv"
#csv_file_energy = f"Energy_csv/{Neck_nr}_energy_{grid_size}_2.csv"
#csv_file_curv_energy = f"Curv_Energy_csv/{Neck_nr}_energy_{grid_size}_2.csv"
#csv_file_surface_energy = f"Surface_Energy_csv/{Neck_nr}_energy_{grid_size}_2.csv"
#csv_file = f"{Neck_nr}_distances_{grid_size}.csv"
# Read the CSV file (assuming it has no header and each value corresponds to a cell in the grid)
df_dist = pd.read_csv(csv_file_dist, header=None)
df_energy = pd.read_csv(csv_file_energy, header=None)
df_curv_energy = pd.read_csv(csv_file_curv_energy, header=None)
df_surface_energy = pd.read_csv(csv_file_surface_energy, header=None)
# Convert the dataframe to a numpy array (if necessary)
grid_dist = df_dist.to_numpy()
grid_energy = df_energy.to_numpy()
grid_curv_energy = df_curv_energy.to_numpy()
grid_surface_energy = df_surface_energy.to_numpy()


grid_area = grid_surface_energy / sigma_array


vmin = 0
vmax = 14000
masked_data = np.where((grid_energy >= vmin) & (grid_energy <= vmax), grid_energy, np.nan)
vmax = 7000
masked_data_curv = np.where((grid_curv_energy >= vmin) & (grid_curv_energy <= vmax), grid_curv_energy, np.nan)
masked_data_surface = np.where((grid_surface_energy >= vmin) & (grid_surface_energy <= vmax), grid_surface_energy, np.nan)
vmax = 14000
masked_data_area = np.where((grid_area >= vmin) & (grid_area <= vmax), grid_area, np.nan)


vmax = 0.6
masked_data_dist = np.where((grid_dist >= vmin) & (grid_dist <= vmax), grid_dist, np.nan)
############################
############################ Distance Plot
############################
fig, ax = plt.subplots(figsize=(16, 11))

# Plot the heatmap
cax_dist = ax.imshow(grid_dist, origin='lower', aspect='auto')

# Set tick positions
ax.set_xticks(np.arange(grid_size,step=2))
ax.set_yticks(np.arange(grid_size,step=2))

# Set tick labels
ax.set_xticklabels([f"{val:.2f}" for val in lambda_array[::2]], rotation=-45)
ax.set_yticklabels([f"{val:.2f}" for val in h_zero_array[::2]])

# Create a secondary x-axis on top
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.set_xticklabels([f"{val:.2f}" for val in sigma_array[::2]], rotation=45)

# Set axis labels
ax.set_xlabel(r'$\lambda \:/\: \mathrm{nm}$')
ax_top.set_xlabel(r'$\sigma \:/\: \mathrm{k_BT/nm²}$')
ax.set_ylabel(r'$H_0 \:/\: \mathrm{nm^{-1}}$')

# Add vertical lines at lambda positions
for idx in lambda_indices:
    ax.axvline(idx, color='white', linestyle='--', linewidth=1.5)
# Add colorbar
plt.colorbar(cax_dist, label=r"Integrated Distance $\:/\: \mathrm{nm}$")
plt.title(rf"Neck {Neck_nr}")

letters = ["a","b","c","d","e"]
for i, (x, y) in enumerate(point_interest[Neck_nr]):
    ax.annotate(
        letters[i],  # Label
        xy=(x, y),  # Point where arrow points
        xytext=(x - 2, y + 2),  # Offset position for label
        color='white',
        fontsize=15,
        fontweight='bold',
        ha='center',
        va='center',
        arrowprops=dict(arrowstyle="->", color='white', lw=1.5)  # Arrow properties
    )

plt.savefig(f"Plots/{Neck_nr}_distances_{grid_size}_2_smooth.pdf", format='pdf')


############################
############################ Energy Plot
############################

point_interest = {
    1: [[9,5],[20,26],[34,33],[49,17]],  
    5: [[9,5],[20,26],[34,33],[49,17]],  
    334: [[9,5],[20,26],[34,33],[49,14]],  
    574: [[9,5],[20,24],[34,29],[49,15]] 
}

fig, ax = plt.subplots(figsize=(16, 11))

# Plot the heatmap
cax_dist = ax.imshow(masked_data, origin='lower', aspect='auto')

# Set tick positions
ax.set_xticks(np.arange(grid_size,step = 2))
ax.set_yticks(np.arange(grid_size,step = 2))

# Set tick labels
ax.set_xticklabels([f"{val:.2f}" for val in lambda_array[::2]], rotation=-45)
ax.set_yticklabels([f"{val:.2f}" for val in h_zero_array[::2]])

# Create a secondary x-axis on top
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.set_xticklabels([f"{val:.2f}" for val in sigma_array[::2]], rotation=45)

# Set axis labels
ax.set_xlabel(r'$\lambda \:/\: \mathrm{nm}$')
ax_top.set_xlabel(r'$\sigma \:/\: \mathrm{k_BT/nm²}$')
ax.set_ylabel(r'$H_0 \:/\: \mathrm{nm^{-1}}$')

for idx in lambda_indices:
    ax.axvline(idx, color='white', linestyle='--', linewidth=1.5)

# Add colorbar
plt.colorbar(cax_dist, label=r"Energy $ \:/\: \mathrm{k_BT}$")
plt.title(rf"Neck {Neck_nr}")

letters = ["a","b","c","m"]
for i, (x, y) in enumerate(point_interest[Neck_nr]):
    ax.annotate(
        letters[i],  # Label
        xy=(x, y),  # Point where arrow points
        xytext=(x - 2, y + 2),  # Offset position for label
        color='white',
        fontsize=15,
        fontweight='bold',
        ha='center',
        va='center',
        arrowprops=dict(arrowstyle="->", color='white', lw=1.5)  # Arrow properties
    )


plt.savefig(f"Plots/{Neck_nr}_energy_{grid_size}_2_smooth.pdf", format='pdf')


############################
############################ Curv Energy Plot
############################

fig, ax = plt.subplots(figsize=(16, 11))

# Plot the heatmap
cax_dist = ax.imshow(masked_data_curv, origin='lower', aspect='auto')

# Set tick positions
ax.set_xticks(np.arange(grid_size,step = 2))
ax.set_yticks(np.arange(grid_size,step = 2))

# Set tick labels
ax.set_xticklabels([f"{val:.2f}" for val in lambda_array[::2]], rotation=-45)
ax.set_yticklabels([f"{val:.2f}" for val in h_zero_array[::2]])

# Create a secondary x-axis on top
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.set_xticklabels([f"{val:.2f}" for val in sigma_array[::2]], rotation=45)

# Set axis labels
ax.set_xlabel(r'$\lambda \:/\: \mathrm{nm}$')
ax_top.set_xlabel(r'$\sigma \:/\: \mathrm{k_BT/nm²}$')
ax.set_ylabel(r'$H_0 \:/\: \mathrm{nm^{-1}}$')

for idx in lambda_indices:
    ax.axvline(idx, color='white', linestyle='--', linewidth=1.5)

# Add colorbar
plt.colorbar(cax_dist, label=r"Curvature Energy $ \:/\: \mathrm{k_BT}$")
plt.title(rf"Neck {Neck_nr}")


plt.savefig(f"Plots/{Neck_nr}_curv_energy_{grid_size}_2_smooth.pdf", format='pdf')

############################
############################ surface Energy Plot
############################

fig, ax = plt.subplots(figsize=(16, 11))

# Plot the heatmap
cax_dist = ax.imshow(masked_data_surface, origin='lower', aspect='auto')

# Set tick positions
ax.set_xticks(np.arange(grid_size,step = 2))
ax.set_yticks(np.arange(grid_size,step = 2))

# Set tick labels
ax.set_xticklabels([f"{val:.2f}" for val in lambda_array[::2]], rotation=-45)
ax.set_yticklabels([f"{val:.2f}" for val in h_zero_array[::2]])

# Create a secondary x-axis on top
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.set_xticklabels([f"{val:.2f}" for val in sigma_array[::2]], rotation=45)

# Set axis labels
ax.set_xlabel(r'$\lambda \:/\: \mathrm{nm}$')
ax_top.set_xlabel(r'$\sigma \:/\: \mathrm{k_BT/nm²}$')
ax.set_ylabel(r'$H_0 \:/\: \mathrm{nm^{-1}}$')

for idx in lambda_indices:
    ax.axvline(idx, color='white', linestyle='--', linewidth=1.5)

# Add colorbar
plt.colorbar(cax_dist, label=r"Surface Energy $ \:/\: \mathrm{k_BT}$")
plt.title(rf"Neck {Neck_nr}")


plt.savefig(f"Plots/{Neck_nr}_surface_energy_{grid_size}_2_smooth.pdf", format='pdf')



############################
############################ surface Area Plot
############################

fig, ax = plt.subplots(figsize=(16, 11))

# Plot the heatmap
cax_dist = ax.imshow(masked_data_area, origin='lower', aspect='auto')

# Set tick positions
ax.set_xticks(np.arange(grid_size,step = 2))
ax.set_yticks(np.arange(grid_size,step = 2))

# Set tick labels
ax.set_xticklabels([f"{val:.2f}" for val in lambda_array[::2]], rotation=-45)
ax.set_yticklabels([f"{val:.2f}" for val in h_zero_array[::2]])

# Create a secondary x-axis on top
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.set_xticklabels([f"{val:.2f}" for val in sigma_array[::2]], rotation=45)

# Set axis labels
ax.set_xlabel(r'$\lambda \:/\: \mathrm{nm}$')
ax_top.set_xlabel(r'$\sigma \:/\: \mathrm{k_BT/nm²}$')
ax.set_ylabel(r'$H_0 \:/\: \mathrm{nm^{-1}}$')

for idx in lambda_indices:
    ax.axvline(idx, color='white', linestyle='--', linewidth=1.5)

# Add colorbar
plt.colorbar(cax_dist, label=r"Neck Area $ \:/\: \mathrm{nm²}$")
plt.title(rf"Neck {Neck_nr}")


plt.savefig(f"Plots/{Neck_nr}_area_{grid_size}_2_smooth.pdf", format='pdf')

############################
############################ Energy densitiy Plot
############################

fig, ax = plt.subplots(figsize=(16, 11))

# Plot the heatmap
cax_dist = ax.imshow(masked_data/masked_data_area, origin='lower', aspect='auto')

# Set tick positions
ax.set_xticks(np.arange(grid_size,step = 2))
ax.set_yticks(np.arange(grid_size,step = 2))

# Set tick labels
ax.set_xticklabels([f"{val:.2f}" for val in lambda_array[::2]], rotation=-45)
ax.set_yticklabels([f"{val:.2f}" for val in h_zero_array[::2]])

# Create a secondary x-axis on top
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.set_xticklabels([f"{val:.2f}" for val in sigma_array[::2]], rotation=45)

# Set axis labels
ax.set_xlabel(r'$\lambda \:/\: \mathrm{nm}$')
ax_top.set_xlabel(r'$\sigma \:/\: \mathrm{k_BT/nm²}$')
ax.set_ylabel(r'$H_0 \:/\: \mathrm{nm^{-1}}$')

for idx in lambda_indices:
    ax.axvline(idx, color='white', linestyle='--', linewidth=1.5)

# Add colorbar
plt.colorbar(cax_dist, label=r"Curvature Energy $ \:/\: \mathrm{k_BT}$")
plt.title(rf"Neck {Neck_nr}")


plt.savefig(f"Plots/{Neck_nr}_density_energy_{grid_size}_2_smooth.pdf", format='pdf')


min_index = np.unravel_index(np.nanargmin(grid_dist), grid_dist.shape)
h_i = min_index[0]
l_i = min_index[1]

h_zero = h_zero_min + h_i * (h_zero_max - h_zero_min) / (grid_size - 1)
ex_l = l_i/ (grid_size - 1) 

lambda_ = lambda_min * (lambda_max / lambda_min)**ex_l


sigma = kappa_param/(lambda_**2) 
print(f"The indices of the smallest value are: {min_index}")
print(f'H_zero = {h_zero}')
print(f'lambda = {lambda_}')
print(f'sigma = {sigma}')
print(f"The smallest value is: {grid_dist[min_index]}")


min_index = np.unravel_index(np.nanargmin(grid_energy), grid_energy.shape)
h_i = min_index[0]
l_i = min_index[1]

h_zero = h_zero_min + h_i * (h_zero_max - h_zero_min) / (grid_size - 1)
ex_l = l_i/ (grid_size - 1) 

lambda_ = lambda_min * (lambda_max / lambda_min)**ex_l


sigma = kappa_param/(lambda_**2) 
print(f"The indices of the smallest value are: {min_index}")
print(f'H_zero = {h_zero}')
print(f'lambda = {lambda_}')
print(f'sigma = {sigma}')
print(f"The smallest energy value is: {grid_energy[min_index]}")


min_index = np.unravel_index(np.nanargmin(grid_curv_energy), grid_curv_energy.shape)
h_i = min_index[0]
l_i = min_index[1]

h_zero = h_zero_min + h_i * (h_zero_max - h_zero_min) / (grid_size - 1)
ex_l = l_i/ (grid_size - 1) 

lambda_ = lambda_min * (lambda_max / lambda_min)**ex_l


sigma = kappa_param/(lambda_**2) 
print(f"The indices of the smallest value are: {min_index}")
print(f'H_zero = {h_zero}')
print(f'lambda = {lambda_}')
print(f'sigma = {sigma}')
print(f"The smallest curv energy value is: {grid_curv_energy[min_index]}")


min_index = np.unravel_index(np.nanargmin(grid_surface_energy), grid_surface_energy.shape)
h_i = min_index[0]
l_i = min_index[1]

h_zero = h_zero_min + h_i * (h_zero_max - h_zero_min) / (grid_size - 1)
ex_l = l_i/ (grid_size - 1) 

lambda_ = lambda_min * (lambda_max / lambda_min)**ex_l


sigma = kappa_param/(lambda_**2) 
print(f"The indices of the smallest value are: {min_index}")
print(f'H_zero = {h_zero}')
print(f'lambda = {lambda_}')
print(f'sigma = {sigma}')
print(f"The smallest surface energy value is: {grid_surface_energy[min_index]}")



