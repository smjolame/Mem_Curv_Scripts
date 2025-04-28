import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
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

# Initialize dictionaries to store arrays for each property
data_sample = {
    "H": {},
    "K": {},
    "Thickness": {},
    "Energy": {}
}

# List of properties and corresponding cell data keys
properties = {
    "H": "Mean_Curvature_smoothed",
    "K": "Gaussian_Curvature_smoothed",
    "Thickness": "Thickness",
    "Energy": "Energy_per_square_nm_kBT"
}

nice_necks = [[3, 3, 1], [3, 3, 2], [3, 3, 4], [3, 3, 5], [3, 6, 3], [4, 5, 2], [5, 7, 1], [5, 7, 3], [5, 7, 4]]

# Function to extract and clean data
def extract_data(pv_mesh, mask_neck, mask_sphere, prop_name):
    prop_data = pv_mesh.cell_data[prop_name]
    Neck_data = prop_data[mask_neck][~np.isnan(prop_data[mask_neck])]
    Sphere_data = prop_data[mask_sphere][~np.isnan(prop_data[mask_sphere])]
    return [Neck_data, Sphere_data]

# Main loop through nice_necks
for Sample, TS, VLP in nice_necks:
    file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'
    pv_mesh = pv.read(file_path)

    print('Sample:', Sample, 'TS:', TS, 'VLP:', VLP)
    
    # Create masks for neck and sphere
    mask_neck = pv_mesh.cell_data['Neck'] == 1
    mask_sphere = pv_mesh.cell_data['Sphere'] == 1
    # Extract only the Neck region
    neck_mesh = pv_mesh.extract_cells(mask_neck)

    # Save the new mesh
    neck_mesh.save(f'Results/Nice_Necks/Neck{Sample}{TS}{VLP}.vtk')

    # Initialize sample entries if they don't exist
    for prop in properties:
        if Sample not in data_sample[prop]:
            data_sample[prop][Sample] = []
    # Extract and store data for each property
    for prop, prop_key in properties.items():
        data_sample[prop][Sample].append(extract_data(pv_mesh, mask_neck, mask_sphere, prop_key))
# Sort the data by sample and convert dictionaries to lists
sorted_data_sample = {prop: [data_sample[prop][s] for s in sorted(data_sample[prop].keys())] for prop in data_sample}

# Now sorted_data_sample contains lists of the cleaned arrays for each sample and property




def plot_boxplots(data_sample, ylabel, title, output_file):
    num_samples = len(data_sample)
    max_pairs = max(len(sample) for sample in data_sample)

    fig, axs = plt.subplots(num_samples, 1, figsize=(8, 8), sharex=True)
    fig.subplots_adjust(hspace=0.5)

    for i, (ax, sample) in enumerate(zip(axs, data_sample)):
        x_positions = np.arange(len(sample)) * 2
        for j, (arr1, arr2) in enumerate(sample):
            ax.boxplot([arr1, arr2], positions=[x_positions[j], x_positions[j] + 0.5], widths=0.3, showfliers=False)
        ax.set_ylabel(ylabel)
        ax.grid(True)
        pair_centers = np.arange(max_pairs) * 2 + 0.25
        ax.set_xticks(pair_centers)
        ax.set_xticklabels([f'{j+1}' for j in range(max_pairs)])
        ax.set_title(Samp[i])
    
    axs[-1].set_xlabel('VLP')
    fig.suptitle(title)
    plt.savefig(output_file, format='pdf')

# Plot for each property
plot_boxplots(sorted_data_sample['H'], r'$H \:/\: \mathrm{nm^{-1}}$', 'Comparison between neck(left)- and sphere(right)- mean curvature', 'Plots/Nice_Necks/Boxplot_Neck_Sphere_H_Curv.pdf')
plot_boxplots(sorted_data_sample['K'], r'$K \:/\: \mathrm{nm^{-2}}$', 'Comparison between neck(left)- and sphere(right)- Gaussian curvature', 'Plots/Nice_Necks/Boxplot_Neck_Sphere_K_Curv.pdf')
plot_boxplots(sorted_data_sample['Energy'], r"Thickness $ \:/\: \mathrm{k_BT/nm^2}$", 'Comparison between neck(left)- and sphere(right)- energy', 'Plots/Nice_Necks/Boxplot_Neck_Sphere_Energy.pdf')
plot_boxplots(sorted_data_sample['Thickness'], r"Thickness $ \:/\: \mathrm{nm}$", 'Comparison between neck(left)- and sphere(right)- thickness', 'Plots/Nice_Necks/Boxplot_Neck_Sphere_Thickness.pdf')

