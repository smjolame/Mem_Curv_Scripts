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

diameter_array_Sample = [] 
diameter_Neck_array_Sample = []
diameter_noNeck_array_Sample = []


sphere_energy_array_Sample = []
sphere_energy_Neck_array_Sample = []
sphere_energy_noNeck_array_Sample = []

H_Neck_array_Sample = []
H_noNeck_array_Sample = []

H_Neck_array_Sample_var = []
H_noNeck_array_Sample_var = []

H_dia_Neck_array_Sample = []
H_dia_noNeck_array_Sample = []

K_Neck_array_Sample = []
K_noNeck_array_Sample = []



neck_energy_array_Sample = []



Sample = 3
Sample_path = f'Results/Sample_{Sample}'
## sample loop:
while os.path.exists(Sample_path):
    diameter_Neck_array_TS = [] 
    diameter_noNeck_array_TS = [] 

    sphere_energy_Neck_array_TS = []
    sphere_energy_noNeck_array_TS = []

    H_Neck_array_TS = []
    H_noNeck_array_TS = []

    H_Neck_array_TS_var = []
    H_noNeck_array_TS_var = []

    H_dia_Neck_array_TS = []
    H_dia_noNeck_array_TS = []

    K_Neck_array_TS = []
    K_noNeck_array_TS = []

    neck_energy_array_TS = []

    if Sample == 3: TS = 2
    elif Sample == 4: TS = 3
    elif Sample == 5: TS = 7
    TS_path = f'Results/Sample_{Sample}/TS00{TS}'
    ## TS loop:
    while os.path.exists(TS_path):
        VLP = 1
        file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'
        Neck = True

        diameter_Neck_array = [] 
        diameter_noNeck_array = [] 
        
        sphere_energy_Neck_array = []
        sphere_energy_noNeck_array = []

        H_Neck_array = []
        H_noNeck_array = []

        H_Neck_array_var = []
        H_noNeck_array_var = []

        H_dia_Neck_array = []
        H_dia_noNeck_array = []

        K_Neck_array = []
        K_noNeck_array = []

        #neck data
        neck_energy_array = [] 
        while os.path.exists(file_path):
            print('Sample: ', Sample, 'TS: ', TS, 'VLP: ', VLP)
            pv_mesh_neck = pv.read(file_path)
            neck_energy = pv_mesh_neck.field_data['Energy_in_Neck_per_square_nm_kBT'][0]
            neck_energy_array.append(neck_energy)
            diameter = pv_mesh_neck.field_data['Radius_Sphere_fit'][0] * 2
            diameter_Neck_array.append(diameter)
            sphere_energy = pv_mesh_neck.field_data['Energy_in_Sphere_per_square_nm_kBT'][0]
            sphere_energy_Neck_array.append(sphere_energy)
            
            mask_sphere = pv_mesh_neck.cell_data['Sphere'] == 1

            H = pv_mesh_neck.cell_data['Mean_Curvature_smoothed']
            Sphere_H = H[mask_sphere]
            Sphere_H = Sphere_H[~np.isnan(Sphere_H)]


            H_Neck_array.append(np.median(Sphere_H))
            H_Neck_array_var.append(np.std(Sphere_H, ddof=1))

            K = pv_mesh_neck.cell_data['Gaussian_Curvature_smoothed']
            Sphere_K = K[mask_sphere]
            Sphere_K = Sphere_K[~np.isnan(Sphere_K)]
            K_Neck_array.append(np.mean(Sphere_K))

            

            R = -diameter/2
            H_dia = 1/R

            H_dia_Neck_array.append(H_dia)

            VLP +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'
        
        file_path = f'Results/Sample_{Sample}/TS00{TS}/noNeck/VLP_{VLP}.vtk'
        Neck = False
               
        while os.path.exists(file_path):
            print('Sample: ', Sample, 'TS: ', TS, 'VLP: ', VLP)
            pv_mesh_noNeck = pv.read(file_path)
            diameter = pv_mesh_noNeck.field_data['Radius_Sphere_fit'][0] * 2
            diameter_noNeck_array.append(diameter)
            sphere_energy = pv_mesh_noNeck.field_data['Energy_in_Sphere_per_square_nm_kBT'][0]
            sphere_energy_noNeck_array.append(sphere_energy)


            mask_sphere = pv_mesh_noNeck.cell_data['Sphere'] == 1

            H = pv_mesh_noNeck.cell_data['Mean_Curvature_smoothed']
            Sphere_H = H[mask_sphere]
            Sphere_H = Sphere_H[~np.isnan(Sphere_H)]
            H_noNeck_array.append(np.median(Sphere_H))
            H_noNeck_array_var.append(np.std(Sphere_H, ddof=1))
            
            K = pv_mesh_noNeck.cell_data['Gaussian_Curvature_smoothed']
            Sphere_K = K[mask_sphere]
            Sphere_K = Sphere_K[~np.isnan(Sphere_K)]
            K_noNeck_array.append(np.mean(Sphere_K))

            R = -diameter/2
            H_dia = 1/R

            H_dia_noNeck_array.append(H_dia)
            
            VLP +=1
            file_path= f'Results/Sample_{Sample}/TS00{TS}/noNeck/VLP_{VLP}.vtk'
        

        sphere_energy_Neck_array_TS.extend(sphere_energy_Neck_array)
        diameter_Neck_array_TS.extend(diameter_Neck_array)
        H_Neck_array_TS.extend(H_Neck_array)
        H_Neck_array_TS_var.extend(H_Neck_array_var)
        K_Neck_array_TS.extend(K_Neck_array)
        H_dia_Neck_array_TS.extend(H_dia_Neck_array)

        sphere_energy_noNeck_array_TS.extend(sphere_energy_noNeck_array)
        diameter_noNeck_array_TS.extend(diameter_noNeck_array)
        H_noNeck_array_TS.extend(H_noNeck_array)
        H_noNeck_array_TS_var.extend(H_noNeck_array_var)
        K_noNeck_array_TS.extend(K_noNeck_array)
        H_dia_noNeck_array_TS.extend(H_dia_noNeck_array)
    
        neck_energy_array_TS.extend(neck_energy_array)
        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'

    diameter_array_Sample.append(np.append(diameter_Neck_array_TS,diameter_noNeck_array_TS))
    diameter_Neck_array_Sample.append(diameter_Neck_array_TS)
    diameter_noNeck_array_Sample.append(diameter_noNeck_array_TS)

    sphere_energy_array_Sample.append(np.append(sphere_energy_Neck_array_TS,sphere_energy_noNeck_array_TS))
    sphere_energy_Neck_array_Sample.append(sphere_energy_Neck_array_TS)
    sphere_energy_noNeck_array_Sample.append(sphere_energy_noNeck_array_TS)

    H_Neck_array_Sample.append(H_Neck_array_TS)
    H_noNeck_array_Sample.append(H_noNeck_array_TS)

    H_Neck_array_Sample_var.append(H_Neck_array_TS_var)
    H_noNeck_array_Sample_var.append(H_noNeck_array_TS_var)

    H_dia_Neck_array_Sample.append(H_dia_Neck_array_TS)
    H_dia_noNeck_array_Sample.append(H_dia_noNeck_array_TS)

    K_Neck_array_Sample.append(K_Neck_array_TS)
    K_noNeck_array_Sample.append(K_noNeck_array_TS)

    neck_energy_array_Sample.append(neck_energy_array_TS)

    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'




Sample_array = np.arange(3,6)


# Boxplots
plt.figure(figsize=(10, 6))

plt.boxplot(diameter_array_Sample, positions=Sample_array, widths=0.3)

plt.xlabel(r'Sample')
plt.xticks(Sample_array)  # Set x-ticks to match list indices
plt.ylabel(r"Diameter $ \:/\: \mathrm{nm}$")
plt.title(r'Diameter of spheres')
plt.grid(True)
plt.savefig('Plots/Boxplot_diameter.pdf')
plt.close()




################# Fancy Boxplot

# Number of arrays
n1 = len(diameter_array_Sample)

# Set up positions for the boxes
spacing = 1  # Reduce spacing between boxplot groups
shift = 0.2   # Reduce offset between paired boxes

positions1 = np.arange(1, n1 + 1) * spacing - shift  # Shift left
positions2 = np.arange(1, n1 + 1) * spacing + shift  # Shift right

xtick_positions = np.arange(1, n1 + 1) * spacing  
# Plotting
plt.figure(figsize=(8, 6))

# Boxplot for the first array of arrays (array1)
bp1 = plt.boxplot(diameter_noNeck_array_Sample, positions=positions1, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'), label = 'post-scission',showfliers=False)

# Boxplot for the second array of arrays (array2)
bp2 = plt.boxplot(diameter_Neck_array_Sample, positions=positions2, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='red', color='red'),
                  capprops=dict(color='red'),
                  whiskerprops=dict(color='red'),
                  flierprops=dict(markerfacecolor='red', markeredgecolor='red'),
                  medianprops=dict(color='black'), label = 'late bud',showfliers=False)

# Add scatter points for post scission data 
for i, data in enumerate(diameter_noNeck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions1[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Add scatter points for late bud data 
for i, data in enumerate(diameter_Neck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions2[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)


# Customizing the plot
plt.ylabel(r"Diameter $ \:/\: \mathrm{nm}$")
plt.title(r'Diameter of spheres')
plt.xticks(xtick_positions, labels=[r'E$\Delta$H2',r'E knockout',r'E wildtype'])
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('Plots/noNeck_Neck/Boxplot_diameter_compare.pdf')




# Energy
plt.figure(figsize=(8, 6))

# Boxplot for the first array of arrays (array1)
bp1 = plt.boxplot(sphere_energy_noNeck_array_Sample, positions=positions1, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'), label = 'post-scission',showfliers=False)

# Boxplot for the second array of arrays (array2)
bp2 = plt.boxplot(sphere_energy_Neck_array_Sample, positions=positions2, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='red', color='red'),
                  capprops=dict(color='red'),
                  whiskerprops=dict(color='red'),
                  flierprops=dict(markerfacecolor='red', markeredgecolor='red'),
                  medianprops=dict(color='black'), label = 'late bud',showfliers=False)

# Add scatter points for post scission data 
for i, data in enumerate(sphere_energy_noNeck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions1[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Add scatter points for late bud data 
for i, data in enumerate(sphere_energy_Neck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions2[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Customizing the plot
plt.ylabel(r"Sphere Energy $ \:/\: \mathrm{k_BT/nm^2}$")
plt.title(r'Sphere energy')
plt.xticks(xtick_positions, labels=[r'E$\Delta$H2',r'E knockout',r'E wildtype'])
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('Plots/noNeck_Neck/Boxplot_Sphere_Energy_compare.pdf')




# H
plt.figure(figsize=(8, 6))

all_arrays = [H_noNeck_array_Sample, H_Neck_array_Sample, H_dia_noNeck_array_Sample, H_dia_Neck_array_Sample]



# Boxplot for the first array of arrays (array1)
bp1 = plt.boxplot(H_noNeck_array_Sample, positions=positions1, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'), label = 'post-scission',showfliers=False)

# Boxplot for the second array of arrays (array2)
bp2 = plt.boxplot(H_Neck_array_Sample, positions=positions2, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='red', color='red'),
                  capprops=dict(color='red'),
                  whiskerprops=dict(color='red'),
                  flierprops=dict(markerfacecolor='red', markeredgecolor='red'),
                  medianprops=dict(color='black'), label = 'late bud',showfliers=False)

# Add scatter points for post scission data 
for i, data in enumerate(H_noNeck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions1[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Add scatter points for late bud data 
for i, data in enumerate(H_Neck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions2[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)


all_values = [value for sublist in all_arrays for arr in sublist for value in arr]
H_min = min(all_values) - 0.001
H_max = max(all_values) + 0.001


# Customizing the plot
plt.ylabel(r'$H \:/\: \mathrm{nm^{-1}}$')
plt.title(r'Mean Curvature')
plt.xticks(xtick_positions, labels=[r'E$\Delta$H2',r'E knockout',r'E wildtype'])
plt.ylim(H_min,H_max)
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('Plots/noNeck_Neck/Boxplot_H_compare.pdf')


# H_dia
plt.figure(figsize=(8, 6))

# Boxplot for the first array of arrays (array1)
bp1 = plt.boxplot(H_dia_noNeck_array_Sample, positions=positions1, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'), label = 'post-scission',showfliers=False)

# Boxplot for the second array of arrays (array2)
bp2 = plt.boxplot(H_dia_Neck_array_Sample, positions=positions2, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='red', color='red'),
                  capprops=dict(color='red'),
                  whiskerprops=dict(color='red'),
                  flierprops=dict(markerfacecolor='red', markeredgecolor='red'),
                  medianprops=dict(color='black'), label = 'late bud',showfliers=False)

# Add scatter points for post scission data 
for i, data in enumerate(H_dia_noNeck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions1[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Add scatter points for late bud data 
for i, data in enumerate(H_dia_Neck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions2[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Customizing the plot
plt.ylabel(r'$H \:/\: \mathrm{nm^{-1}}$')
plt.title(r'Mean Curvature of the fitted spheres')
plt.xticks(xtick_positions, labels=[r'E$\Delta$H2',r'E knockout',r'E wildtype'])
plt.ylim(H_min,H_max)
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('Plots/noNeck_Neck/Boxplot_H_dia_compare.pdf')


# K
plt.figure(figsize=(8, 6))

# Boxplot for the first array of arrays (array1)
bp1 = plt.boxplot(K_noNeck_array_Sample, positions=positions1, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'), label = 'post-scission',showfliers=False)

# Boxplot for the second array of arrays (array2)
bp2 = plt.boxplot(K_Neck_array_Sample, positions=positions2, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='red', color='red'),
                  capprops=dict(color='red'),
                  whiskerprops=dict(color='red'),
                  flierprops=dict(markerfacecolor='red', markeredgecolor='red'),
                  medianprops=dict(color='black'), label = 'late bud',showfliers=False)

# Add scatter points for post scission data 
for i, data in enumerate(K_noNeck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions1[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Add scatter points for late bud data 
for i, data in enumerate(K_Neck_array_Sample):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions2[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Customizing the plot
plt.ylabel(r'$K \:/\: \mathrm{nm^{-2}}$')
plt.title(r'Gaussian Curvature')
plt.xticks(xtick_positions, labels=[r'E$\Delta$H2',r'E knockout',r'E wildtype'])
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('Plots/noNeck_Neck/Boxplot_K_compare.pdf')


# H_dia_var
all_values = [value for sublist in [H_noNeck_array_Sample_var,H_Neck_array_Sample_var] for arr in sublist for value in arr]
H_min = min(all_values) - 0.001
H_max = max(all_values) + 0.001
plt.figure(figsize=(8, 6))

# Boxplot for the first array of arrays (array1)
bp1 = plt.boxplot(H_noNeck_array_Sample_var, positions=positions1, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'), label = 'post-scission',showfliers=False)

# Boxplot for the second array of arrays (array2)
bp2 = plt.boxplot(H_Neck_array_Sample_var, positions=positions2, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='red', color='red'),
                  capprops=dict(color='red'),
                  whiskerprops=dict(color='red'),
                  flierprops=dict(markerfacecolor='red', markeredgecolor='red'),
                  medianprops=dict(color='black'), label = 'late bud',showfliers=False)

# Add scatter points for post scission data 
for i, data in enumerate(H_noNeck_array_Sample_var):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions1[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Add scatter points for late bud data 
for i, data in enumerate(H_Neck_array_Sample_var):
    x_jitter = np.random.normal(0, 0.05, size=len(data))  # Add small horizontal jitter
    plt.scatter(np.full_like(data, positions2[i]) + x_jitter, data, color='black', alpha=0.7, s=20, zorder=4)

# Customizing the plot
plt.ylabel(r'$H \:/\: \mathrm{nm^{-1}}$')
plt.title(r'Standard deviation of Mean Curvature ')
plt.xticks(xtick_positions, labels=[r'E$\Delta$H2',r'E knockout',r'E wildtype'])
plt.ylim(H_min,H_max)
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig('Plots/noNeck_Neck/Boxplot_H_dia_compare_var.pdf')



##### T Test
from scipy.stats import ttest_ind
from itertools import combinations
# Perform pairwise t-tests

print('Diameter')
pairwise_results = {}
for i, j in combinations(range(len(diameter_noNeck_array_Sample)), 2):
    t_stat, p_value = ttest_ind(diameter_noNeck_array_Sample[i], diameter_noNeck_array_Sample[j], equal_var=True)
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('post-scission:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")


    # Perform pairwise t-tests
pairwise_results = {}
for i, j in combinations(range(len(diameter_Neck_array_Sample)), 2):
    t_stat, p_value = ttest_ind(diameter_Neck_array_Sample[i], diameter_Neck_array_Sample[j], equal_var=True)
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('late bud:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")

# Perform pairwise t-tests


pairwise_results = {}
for i in range(len(diameter_noNeck_array_Sample)):
    t_stat, p_value = ttest_ind(diameter_noNeck_array_Sample[i], diameter_Neck_array_Sample[i], equal_var=True)
    pairwise_results[f"Sample {i+1}"] = (t_stat, p_value)

# Print results
print('between categories:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")


##################################################################################################################

print('H')
pairwise_results = {}
for i, j in combinations(range(len(H_noNeck_array_Sample)), 2):
    t_stat, p_value = ttest_ind(H_noNeck_array_Sample[i], H_noNeck_array_Sample[j], equal_var=True)
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('post-scission:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")


    # Perform pairwise t-tests
pairwise_results = {}
for i, j in combinations(range(len(H_Neck_array_Sample)), 2):
    t_stat, p_value = ttest_ind(H_Neck_array_Sample[i], H_Neck_array_Sample[j], equal_var=True)
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('late bud:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")

pairwise_results = {}
for i in range(len(diameter_noNeck_array_Sample)):
    t_stat, p_value = ttest_ind(H_noNeck_array_Sample[i], H_Neck_array_Sample[i], equal_var=True)
    pairwise_results[f"Sample {i+1}"] = (t_stat, p_value)

# Print results
print('between categories:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")


##################################################################################################################

print('K')
pairwise_results = {}
for i, j in combinations(range(len(K_noNeck_array_Sample)), 2):
    t_stat, p_value = ttest_ind(K_noNeck_array_Sample[i], K_noNeck_array_Sample[j], equal_var=True)
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('post-scission:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")


    # Perform pairwise t-tests
pairwise_results = {}
for i, j in combinations(range(len(K_Neck_array_Sample)), 2):
    t_stat, p_value = ttest_ind(K_Neck_array_Sample[i], K_Neck_array_Sample[j], equal_var=True)
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('late bud:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")

pairwise_results = {}
for i in range(len(diameter_noNeck_array_Sample)):
    t_stat, p_value = ttest_ind(K_noNeck_array_Sample[i], K_Neck_array_Sample[i], equal_var=True)
    pairwise_results[f"Sample {i+1}"] = (t_stat, p_value)

# Print results
print('between categories:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")

##################################################################################################################

print('H_var')
pairwise_results = {}
for i, j in combinations(range(len(H_noNeck_array_Sample_var)), 2):
    t_stat, p_value = ttest_ind(H_noNeck_array_Sample_var[i], H_noNeck_array_Sample_var[j], equal_var=True, alternative='two-sided')
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('post-scission:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")


    # Perform pairwise t-tests
pairwise_results = {}
for i, j in combinations(range(len(H_Neck_array_Sample_var)), 2):
    t_stat, p_value = ttest_ind(H_Neck_array_Sample_var[i], H_Neck_array_Sample_var[j], equal_var=True, alternative='two-sided')
    pairwise_results[f"Sample {i+1} vs Sample {j+1}"] = (t_stat, p_value)

# Print results
print('late bud:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")


pairwise_results = {}
for i in range(len(diameter_noNeck_array_Sample)):
    t_stat, p_value = ttest_ind(H_noNeck_array_Sample_var[i], H_Neck_array_Sample_var[i], equal_var=True)
    pairwise_results[f"Sample {i+1}"] = (t_stat, p_value)

# Print results
print('between categories:')
for pair, (t_stat, p_value) in pairwise_results.items():
    print(f"{pair}: t-statistic = {t_stat:.3f}, p-value = {p_value:.3f}")