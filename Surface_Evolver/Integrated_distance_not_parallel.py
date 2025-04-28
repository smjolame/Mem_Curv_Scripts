import pyvista as pv
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import sys
import os

def extract_energy(obj_filename):
    energies = {"# Total Energy:": None, "# Curvature Energy:": None, "# Surface Energy:": None}
    
    with open(obj_filename, 'r') as file:
        for line in file:
            for key in energies:
                if line.startswith(key):
                    energies[key] = float(line.split(":")[1].strip())
                    if all(energies.values()):  # Stop early if all values are found
                        return tuple(energies.values())
    
    return tuple(energies.values())

# Input arguments
if len(sys.argv) != 4:
    print("Usage: python process_row.py <row_index> <Neck_nr> <grid_sice>")
    sys.exit(1)

i = int(sys.argv[1])  # Row index (i)
Neck_nr = sys.argv[2]  # Neck_nr
n = int(sys.argv[3])  # Number of rows/columns


csv_file=f"Distances_csv/{Neck_nr}_distances_{n}_2.csv"
csv_file_energy=f"Energy_csv/{Neck_nr}_energy_{n}_2.csv"
csv_file_curv_energy=f"Curv_Energy_csv/{Neck_nr}_energy_{n}_2.csv"
csv_file_surface_energy=f"Surface_Energy_csv/{Neck_nr}_energy_{n}_2.csv"



# Parameters
#data_file = f"Necks/Neck{Neck_nr}.stl"
data_file = f"Necks/Neck{Neck_nr}.obj"

# Define file paths and corresponding arrays
files = {
    "distance": csv_file,
    "energy": csv_file_energy,
    "curvature_energy": csv_file_curv_energy,
    "surface_energy": csv_file_surface_energy,
}

# Initialize arrays
arrays = {
    key: pd.read_csv(file, header=None).to_numpy() if os.path.exists(file) else np.full((n, n), np.nan)
    for key, file in files.items()
}

# Unpack if needed
distance_array = arrays["distance"]
energy_array = arrays["energy"]
curv_energy_array = arrays["curvature_energy"]
surface_energy_array = arrays["surface_energy"]

# Precompute triangle areas in data_mesh
data_mesh = pv.read(data_file)
triangle_areas = data_mesh.compute_cell_sizes()["Area"]

# Extract the faces from the data_mesh and compute centroids
faces = data_mesh.faces.reshape(-1, 4)[:, 1:]
vertices_data = data_mesh.points[faces]
centroids = np.mean(vertices_data, axis=1)

# Function to process each (i, j) pair
def process_pair(j):
    print(f"Processing (i={i}, j={j})...")

    # Load the min_mesh for this (i, j) pair
    min_file = f"Min_Necks_temp/Neck{Neck_nr}_{i}{j}.obj"
    if not os.path.exists(min_file):
        distance_array[i, j] = np.nan
        energy_array[i, j] = np.nan
        curv_energy_array[i, j] = np.nan
        surface_energy_array[i, j] = np.nan
    else:

        min_mesh = pv.read(min_file)
        total_energy, curv_energy, surface_energy = extract_energy(min_file)

        # Use cKDTree for faster closest-point search on min_mesh
        tree = cKDTree(min_mesh.points)

        # Find the closest points for all centroids in one go
        distances, _ = tree.query(centroids)

        # Compute weighted distances (area * distance)
        weighted_distances = distances * triangle_areas

        # Normalize the weighted distance sum by the total area
        total_area = np.sum(triangle_areas)
        normalized_distance = np.sum(weighted_distances) / total_area

        # Store the result in the distance_array
        distance_array[i, j] = normalized_distance
        energy_array[i, j] = total_energy
        curv_energy_array[i, j] = curv_energy
        surface_energy_array[i, j] = surface_energy

# Process all j indices sequentially for the given row i
for j in range(n):
    process_pair(j)

# Save the updated distance_array back to the CSV file
pd.DataFrame(distance_array).to_csv(csv_file, header=False, index=False)
pd.DataFrame(energy_array).to_csv(csv_file_energy, header=False, index=False)
pd.DataFrame(curv_energy_array).to_csv(csv_file_curv_energy, header=False, index=False)
pd.DataFrame(surface_energy_array).to_csv(csv_file_surface_energy, header=False, index=False)