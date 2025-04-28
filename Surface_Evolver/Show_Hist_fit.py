import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
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
    1: [[9,5],[20,26],[34,34],[34,35],[34,16]],  
    334: [[9,5],[20,26],[34,33],[34,34],[34,16]],  
    574: [[9,5],[20,24],[34,29],[34,30],[34,16]]   
}


# Create the arrays for h_zero and lambda
h_zero_array = np.linspace(h_zero_min, h_zero_max, grid_size)
#h_zero_array = np.logspace(np.log10(h_zero_min), np.log10(h_zero_max), grid_size)
lambda_array = lambda_min * np.power(lambda_max / lambda_min, np.linspace(0, 1, grid_size))
sigma_array = kappa_param/lambda_array**2
lambda_target = np.array([2, 5, 15, 50])  # The λ values where you want vertical lines
lambda_indices = [np.abs(lambda_array - val).argmin() for val in lambda_target]

csv_file_dist = f"Distances_csv/{Neck_nr}_distances_{grid_size}_2.csv"
csv_file_energy = f"Energy_csv/{Neck_nr}_energy_{grid_size}_2.csv"
csv_file_curv_energy = f"Curv_Energy_csv/{Neck_nr}_energy_{grid_size}_2.csv"
csv_file_surface_energy = f"Surface_Energy_csv/{Neck_nr}_energy_{grid_size}_2.csv"
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
vmax = 14000
masked_data_curv = np.where((grid_curv_energy >= vmin) & (grid_curv_energy <= vmax), grid_curv_energy, np.nan)
masked_data_surface = np.where((grid_surface_energy >= vmin) & (grid_surface_energy <= vmax), grid_surface_energy, np.nan)
vmax = 14000
masked_data_area = np.where((grid_area >= vmin) & (grid_area <= vmax), grid_area, np.nan)


vmax = 0.6
masked_data_dist = np.where((grid_dist >= vmin) & (grid_dist <= vmax), grid_dist, np.nan)


# Define h and l values
l_values = np.logspace(0, np.log10(50), 50)  # Logarithmic spacing from 1 to 50
h_values = np.linspace(-0.3, 0.6, 50)  # Linear spacing from -0.3 to 0.6


# Extract valid (h, l) pairs where grid_dist is not NaN
h_list, l_list = [], []
for i in range(50):  # Loop over h values (rows)
    for j in range(50):  # Loop over l values (columns)
        if not np.isnan(masked_data_dist[i, j]):  # Only take valid data points
            h_list.append(h_values[i])  # Store h value
            l_list.append(l_values[j])  # Store l value

# Convert lists to NumPy arrays
h_array = np.array(h_list)
l_array = np.array(l_list)


# Step 1: Extract valid (non-NaN) values
valid_mask = ~np.isnan(masked_data_dist)  # Mask for non-NaN values

grid_values = masked_data_dist[valid_mask]  # Valid values (to use as weights)

# Step 2: Define Weights (Higher weight for smaller values)
weights = 1 / (grid_values + 1e-6)  # Example: inverse weighting

# Define alternative fitting functions
def log_parabola(h, a, b, c):
    return np.exp(a * (h-b)**2 + c)



# Fit each function
popt_log_parabola, _ = curve_fit(log_parabola, h_array, l_array, sigma=weights, absolute_sigma=False)

# Generate fitted curves
h_fit = np.linspace(h_values.min(), h_values.max(), 100)
l_fit_log_parabola = log_parabola(h_fit, *popt_log_parabola)


# Plot the data
plt.scatter(h_array, l_array, color='blue', alpha=0.5, label="Valid Data")

# Plot fitted curves
plt.plot(h_fit, l_fit_log_parabola, 'r-', label="Log-Parabola")


# Labels and title
plt.xlabel("h")
plt.ylabel("l")
plt.yscale('log')
plt.ylim(0,50)

plt.legend()
plt.title("Testing Different Fits")

# Show plot
plt.show()