import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os
import pandas as pd
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
H_array_Sample = []
H_array_Neck_Sample = []
H_array_Sphere_Sample = []
Sample = 3
Sample_path = f'Results/Sample_{Sample}'
## sample loop:
while os.path.exists(Sample_path):
    H_array = []
    H_array_Neck = []
    H_array_Sphere = []
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
            H = pv_mesh.cell_data['Mean_Curvature_smoothed']
            mask_neck = pv_mesh.cell_data['Neck'] == 1
            mask_sphere = pv_mesh.cell_data['Sphere'] == 1
            Neck_H = H[mask_neck]
            Sphere_H = H[mask_sphere]
            Neck_H = Neck_H[~np.isnan(Neck_H)]
            Sphere_H = Sphere_H[~np.isnan(Sphere_H)]
            H_array_Neck.append(Neck_H)
            H_array_Sphere.append(Sphere_H)
            H_list = [Neck_H, Sphere_H]
            H_array.append(H_list)


            VLP +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'

        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'
    H_array_Sample.append(H_array)
    H_array_Neck_Sample.append(H_array_Neck)
    H_array_Sphere_Sample.append(H_array_Sphere)

    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'


categories = ['E_deltah', 'E_ko', 'E_wt']  # Names for the first dimension
### Neck
# Find the maximum length across all VLPs and categories
max_length = max(
    max(len(vlp_data) for vlp_data in category_data)
    for category_data in H_array_Neck_Sample
)

# Pad the data arrays with NaN to make them the same length
padded_data = []
for category_data in H_array_Neck_Sample:
    padded_category = []
    for vlp_data in category_data:
        padded_category.append(np.pad(vlp_data, (0, max_length - len(vlp_data)), constant_values=np.nan))
    padded_data.append(padded_category)

# Now create columns in the format "category_VLP"
columns = {}
for cat_idx, category in enumerate(categories):
    for vlp_idx, vlp_data in enumerate(padded_data[cat_idx], start=1):
        column_name = f"{category}_{vlp_idx}"
        columns[column_name] = vlp_data

# Create the DataFrame and save to CSV
df = pd.DataFrame(columns)

# Save to CSV
df.to_csv('H_Neck_Arrays.csv', index=False)

### Sphere
# Find the maximum length across all VLPs and categories
max_length = max(
    max(len(vlp_data) for vlp_data in category_data)
    for category_data in H_array_Sphere_Sample
)

# Pad the data arrays with NaN to make them the same length
padded_data = []
for category_data in H_array_Sphere_Sample:
    padded_category = []
    for vlp_data in category_data:
        padded_category.append(np.pad(vlp_data, (0, max_length - len(vlp_data)), constant_values=np.nan))
    padded_data.append(padded_category)

# Now create columns in the format "category_VLP"
columns = {}
for cat_idx, category in enumerate(categories):
    for vlp_idx, vlp_data in enumerate(padded_data[cat_idx], start=1):
        column_name = f"{category}_{vlp_idx}"
        columns[column_name] = vlp_data

# Create the DataFrame and save to CSV
df = pd.DataFrame(columns)

# Save to CSV
df.to_csv('H_Sphere_Arrays.csv', index=False)


# Define the number of samples
num_samples = len(H_array_Sample)

# Calculate the maximum number of pairs across all samples
max_pairs = max(len(sample) for sample in H_array_Sample)


# Create a figure with subplots (one for each sample)
fig, axs = plt.subplots(num_samples, 1, figsize=(8, 8), sharex=True)

# Adjust the space between subplots
fig.subplots_adjust(hspace=0.5)

# Loop through each sample and plot its boxplots in a separate subplot
for i, (ax, sample) in enumerate(zip(axs, H_array_Sample)):
    # Generate x positions for the boxplots, one pair at a time
    x_positions = np.arange(len(sample)) * 2  # One boxplot per pair

    for j, (arr1, arr2) in enumerate(sample):
        # Create boxplots for each pair of arrays
        ax.boxplot([arr1, arr2], positions=[x_positions[j], x_positions[j] + 0.5], widths=0.3, showfliers=False)

    # Set y-label for each sample region (e.g., Sample 3, Sample 4, Sample 5)
    ax.set_ylabel(r'$H \:/\: \mathrm{nm^{-1}}$')

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
fig.suptitle('Comparison between neck- (left) and sphere- (right) mean curvature within each sample')

plt.savefig('Plots/Neck_Sphere/Boxplot_Neck_Sphere_H.pdf', format = 'pdf')