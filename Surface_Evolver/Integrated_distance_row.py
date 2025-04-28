import pyvista as pv
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from concurrent.futures import ThreadPoolExecutor
import sys
import os

# Input arguments
if len(sys.argv) != 4:
    print("Usage: python process_row.py <row_index> <Neck_nr> <grid_sice>")
    sys.exit(1)

i = int(sys.argv[1])  # Row index (i)
Neck_nr = sys.argv[2]  # Neck_nr
n = int(sys.argv[3]) # rows

csv_file=f"Distances_csv/{Neck_nr}_distances_{n}.csv"
#csv_file=f"{Neck_nr}_distances_{n}.csv"
# Parameters
data_file = f"Necks/Neck{Neck_nr}.stl"


# Initialize or read the distance array
if os.path.exists(csv_file):
    distance_array = pd.read_csv(csv_file, header=None).to_numpy()
else:
    distance_array = np.full((n, n), np.nan)  # Initialize with NaN

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
    min_mesh = pv.read(min_file)

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

# Process all j indices in parallel for the given row i
with ThreadPoolExecutor() as executor:
    executor.map(process_pair, range(n))

# Save the updated distance_array back to the CSV file
pd.DataFrame(distance_array).to_csv(csv_file, header=False, index=False)

