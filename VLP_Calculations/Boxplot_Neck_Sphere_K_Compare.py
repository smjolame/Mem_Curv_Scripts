import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os
import matplotlib as mpl
# Use LaTeX for rendering math expressions
mpl.rcParams['text.usetex'] = True

# Set font properties
mpl.rcParams['font.family'] = 'serif'  # Serif for the main text
mpl.rcParams['font.serif'] = ['Times New Roman']  # This can be customized to another serif font
mpl.rcParams['font.size'] = 16  # Set the font size

# Set math font to match "Latin Modern Math" used in your LaTeX document
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'serif'  # Use serif for regular math text (this matches the font used for text)
mpl.rcParams['mathtext.it'] = 'serif:italic'  # Italic for math symbols
mpl.rcParams['mathtext.bf'] = 'serif:bold'  # Bold math symbols


K_array_Neck_Sample = []
K_array_Sphere_Sample = []

Sample = 3
Sample_path = f'Results/Sample_{Sample}'
## sample loop:
while os.path.exists(Sample_path):
    K_array_Neck = []
    K_array_Sphere = []
    if Sample == 3: TS = 2
    elif Sample == 4: TS = 3
    elif Sample == 5: TS = 7
    TS_path = f'Results/Sample_{Sample}/TS00{TS}'
    ## TS loop:
    while os.path.exists(TS_path):

        VLP = 1
        file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'

        while os.path.exists(file_path):
            print('Sample: ', Sample, 'TS: ', TS, 'VLP: ', VLP)
            pv_mesh = pv.read(file_path)
            K = pv_mesh.cell_data['Gaussian_Curvature_smoothed']
            mask_neck = pv_mesh.cell_data['Neck'] == 1
            mask_sphere = pv_mesh.cell_data['Sphere'] == 1
            Neck_K = K[mask_neck]
            Sphere_K = K[mask_sphere]
            Neck_K = Neck_K[~np.isnan(Neck_K)]
            Sphere_K = Sphere_K[~np.isnan(Sphere_K)]
            K_array_Neck.append(np.median(Neck_K))
            K_array_Sphere.append(np.mean(Sphere_K))
            #K_array_Neck.extend(Neck_K)
            #K_array_Sphere.extend(Sphere_K)


            VLP +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'

        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'
    K_array_Neck_Sample.append(K_array_Neck)
    K_array_Sphere_Sample.append(K_array_Sphere)

    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'



################# Fancy Boxplot

# Number of arrays
n1 = len(K_array_Neck_Sample)

# Set up positions for the boxes
spacing = 1  # Reduce spacing between boxplot groups
shift = 0.2   # Reduce offset between paired boxes

positions1 = np.arange(1, n1 + 1) * spacing - shift  # Shift left
positions2 = np.arange(1, n1 + 1) * spacing + shift  # Shift right

xtick_positions = np.arange(1, n1 + 1) * spacing  
# Plotting
plt.figure(figsize=(8, 6))

# Boxplot for the first array of arrays (array1)
bp1 = plt.boxplot(K_array_Neck_Sample, positions=positions1, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'), label = 'Neck',showfliers=False)

# Boxplot for the second array of arrays (array2)
bp2 = plt.boxplot(K_array_Sphere_Sample, positions=positions2, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='red', color='red'),
                  capprops=dict(color='red'),
                  whiskerprops=dict(color='red'),
                  flierprops=dict(markerfacecolor='red', markeredgecolor='red'),
                  medianprops=dict(color='black'), label = 'Sphere',showfliers=False)

# Add scatter points for No Neck data 
for i, data in enumerate(K_array_Neck_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions1[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Add scatter points for With Neck data 
for i, data in enumerate(K_array_Sphere_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions2[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Customizing the plot
plt.ylabel(r'$K \:/\: \mathrm{nm^{-2}}$')
plt.title('Gaussian Curvature per VLP of the neck and of the sphere')
plt.xticks(xtick_positions, labels=[r'E$\Delta$H2',r'E knockout',r'E wildtype'])
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('Plots/Neck_Sphere/Boxplot_Neck_Sphere_K_compare.pdf')
plt.clf()



##### T Test
from scipy.stats import ttest_ind
from itertools import combinations
# Perform pairwise t-tests
pairwise_results = {}
for i, j in combinations(range(len(K_array_Neck_Sample)), 2):
    t_stat, p_value = ttest_ind(K_array_Neck_Sample[i], K_array_Neck_Sample[j], equal_var=True)
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('Neck:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")


    # Perform pairwise t-tests
pairwise_results = {}
for i, j in combinations(range(len(K_array_Sphere_Sample)), 2):
    t_stat, p_value = ttest_ind(K_array_Sphere_Sample[i], K_array_Sphere_Sample[j], equal_var=True)
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('Sphere:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")