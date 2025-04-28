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

Samp = [r'E$\Delta$H2',r'E knockout',r'E wild-type']
Thickness_array_Sample = []

Sample = 3
Sample_path = f'Results/Sample_{Sample}'
## sample loop:
while os.path.exists(Sample_path):
    Thickness_array = []
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
            Thickness_list = [Neck_Thickness, Sphere_Thickness]
            Thickness_array.append(Thickness_list)


            VLP +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'

        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'
    Thickness_array_Sample.append(Thickness_array)

    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'



# Define the number of samples
num_samples = len(Thickness_array_Sample)

# Calculate the maximum number of pairs across all samples
max_pairs = max(len(sample) for sample in Thickness_array_Sample)


# Create a figure with subplots (one for each sample)
fig, axs = plt.subplots(num_samples, 1, figsize=(8, 8), sharex=True)

# Adjust the space between subplots
fig.subplots_adjust(hspace=0.5)

# Loop through each sample and plot its boxplots in a separate subplot
for i, (ax, sample) in enumerate(zip(axs, Thickness_array_Sample)):
    # Generate x positions for the boxplots, one pair at a time
    x_positions = np.arange(len(sample)) * 2  # One boxplot per pair

    for j, (arr1, arr2) in enumerate(sample):
        # Create boxplots for each pair of arrays
        ax.boxplot([arr1, arr2], positions=[x_positions[j], x_positions[j] + 0.5], widths=0.3, showfliers=False)

    # Set y-label for each sample region (e.g., Sample 3, Sample 4, Sample 5)
    ax.set_ylabel(r"Thickness $ \:/\: \mathrm{nm}$")

    # Optionally, add a grid for clarity
    ax.grid(True)
    # Set the X-ticks at the center of each pair of boxplots
    pair_centers = np.arange(max_pairs) * 2 + 0.25
    ax.set_xticks(pair_centers)
    ax.set_xticklabels([f'{j+1}' for j in range(max_pairs)])
    ax.set_title(Samp[i])

# Set x-axis label at the bottom
axs[-1].set_xlabel('VLP')

# Add a title
fig.suptitle('Comparison between neck- (left) and sphere- (right) thickness within each sample')

plt.savefig('Plots/Neck_Sphere/Boxplot_Neck_Sphere_Thickness.pdf', format = 'pdf')