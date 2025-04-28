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


Thickness_array_Neck_Sample = []
Thickness_array_Sphere_Sample = []

Sample = 3
Sample_path = f'Results/Sample_{Sample}'
## sample loop:
while os.path.exists(Sample_path):
    Thickness_array_Neck = []
    Thickness_array_Sphere = []
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
            Thickness = pv_mesh.cell_data['Thickness']
            mask_neck = pv_mesh.cell_data['Neck'] == 1
            mask_sphere = pv_mesh.cell_data['Sphere'] == 1
            Neck_Thickness = Thickness[mask_neck]
            Sphere_Thickness = Thickness[mask_sphere]
            Neck_Thickness = Neck_Thickness[~np.isnan(Neck_Thickness)]
            Sphere_Thickness = Sphere_Thickness[~np.isnan(Sphere_Thickness)]
            Thickness_array_Neck.extend(Neck_Thickness)
            Thickness_array_Sphere.extend(Sphere_Thickness)


            VLP +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'

        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'
    Thickness_array_Neck_Sample.append(Thickness_array_Neck)
    Thickness_array_Sphere_Sample.append(Thickness_array_Sphere)

    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'



################# Fancy Boxplot

# Number of arrays
n1 = len(Thickness_array_Neck_Sample)

# Set up positions for the boxes
spacing = 1  # Reduce spacing between boxplot groups
shift = 0.2   # Reduce offset between paired boxes

positions1 = np.arange(1, n1 + 1) * spacing - shift  # Shift left
positions2 = np.arange(1, n1 + 1) * spacing + shift  # Shift right

xtick_positions = np.arange(1, n1 + 1) * spacing  
# Plotting
plt.figure(figsize=(8, 6))

# Boxplot for the first array of arrays (array1)
bp1 = plt.boxplot(Thickness_array_Neck_Sample, positions=positions1, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'), label = 'Neck',showfliers=False)

# Boxplot for the second array of arrays (array2)
bp2 = plt.boxplot(Thickness_array_Sphere_Sample, positions=positions2, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='red', color='red'),
                  capprops=dict(color='red'),
                  whiskerprops=dict(color='red'),
                  flierprops=dict(markerfacecolor='red', markeredgecolor='red'),
                  medianprops=dict(color='black'), label = 'Sphere',showfliers=False)



# Customizing the plot
plt.ylabel(r"Thickness $ \:/\: \mathrm{nm}$")
plt.title('Thickness of neck and sphere')
plt.xticks(xtick_positions, labels=[r'E$\Delta$H2',r'E knockout',r'E wild-type'])
plt.grid(True)
plt.legend()
plt.savefig('Plots/Neck_Sphere/Boxplot_Neck_Sphere_Thickness_compare.pdf')
plt.clf()