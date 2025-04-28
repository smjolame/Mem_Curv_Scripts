import pyvista as pv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 17  # Set the font size
kappa = 25
sigma = 0.25
grid_size = 50
h_zero_min = -0.3
h_zero_max = 0.6
lambda_min = 1
lambda_max = 50
h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
lambda_array = lambda_min * np.power(lambda_max / lambda_min, np.linspace(0, 1, grid_size))
sigma_array = kappa/lambda_array**2

#nice_necks = [[3, 3, 1], [3, 3, 2], [3, 3, 4], [3, 3, 5], [3, 6, 3], [5, 7, 1], [5, 7, 3], [5, 7, 4]]
nice_necks = [334, 574]
for Neck_nr in nice_necks:


    file_path=f'Necks/Neck{Neck_nr}.vtk'
    mesh = pv.read(file_path)

    H_smoothed = mesh.cell_data['Mean_Curvature_smoothed']
    H = mesh.cell_data['Mean_Curvature']
    # Energy Calculations
    def Surface_Energy_per_Cell(Cell_area, sigma=1):
        return sigma*Cell_area
    def Curv_Energy_per_Cell(H,Cell_area,H_0=0):
        return kappa*(H-H_0)**2*Cell_area
    

    Cell_areas = mesh.compute_cell_sizes(length=False, area=True, volume=False)['Area']
    Surface_Energy_per_cell = Surface_Energy_per_Cell(Cell_areas) 
    Curv_Energy_per_cell = Curv_Energy_per_Cell(H,Cell_areas)
    Energy_per_Cell = Surface_Energy_per_cell + Curv_Energy_per_cell
    Total_Energy = np.sum(Energy_per_Cell)
    Total_Curv_Energy = np.sum(Curv_Energy_per_cell)
    Total_Surface_Energy = np.sum(Surface_Energy_per_cell)
    Average_Energy = Total_Energy/(np.sum(Cell_areas))

    print(f'Neck{Neck_nr}')
    print('Total Energy:',Total_Energy)
    print('Total Curv Energy:',Total_Curv_Energy)
    print('Total Surface Energy:',Total_Surface_Energy)
    print('Average_Energy:',Average_Energy)
    print('Total Area:', (np.sum(Cell_areas)))


    csv_file_energy = f"Energy_csv/{Neck_nr}_energy_{grid_size}_2_smooth.csv"
    df_energy = pd.read_csv(csv_file_energy, header=None)
    grid_energy = df_energy.to_numpy()
        # Prepare a DataFrame to store energies
    energy_matrix = np.zeros((grid_size, grid_size))

    for i, H_0 in enumerate(h_zero_array):
        for j, sigma in enumerate(sigma_array):
            Surface_Energy = Surface_Energy_per_Cell(Cell_areas, sigma)
            Curv_Energy = Curv_Energy_per_Cell(H, Cell_areas, H_0)
            Energy_per_Cell = Surface_Energy + Curv_Energy
            Total_Energy = np.sum(Energy_per_Cell)
            
            # Store the total energy in the matrix
            energy_matrix[i, j] = Total_Energy

    
    vmin = 0
    vmax = 15000
    masked_data = np.where((energy_matrix >= vmin) & (energy_matrix <= vmax), energy_matrix, np.nan)

    diff = energy_matrix - grid_energy
    vmin = -2000
    vmax = 2000
    masked_data = np.where((diff >= vmin) & (diff <= vmax), diff, np.nan)
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
    # Add colorbar
    plt.colorbar(cax_dist, label=r"Energy $ \:/\: \mathrm{k_BT}$")
    plt.title(rf"Neck {Neck_nr}")


    #plt.show()
    plt.savefig(f"Plots/Energy/{Neck_nr}_energy_difference_{grid_size}_2.pdf", format='pdf')
    plt.clf()
