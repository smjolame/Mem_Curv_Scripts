
import pyvista as pv
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

Neck_nr = 1
n = 50
lambda_ = 15
h_it = 23
# Parameters
point_interest_1 = [[2,5],[5,26],[15,34],[15,35], [15,14], [50,17]]
point_interest_334 = [[2,5],[5,26],[15,33],[15,34], [15,16], [50,14]]
point_interest_574 = [[2,5],[5,24],[15,29],[15,30], [15,16], [50,15]]

#point_interest_334 = [[2,5],[5,26],[15,33],[50,14]]
#point_interest_574 = [[2,5],[5,24],[15,29],[50,15]]

camera_positions = {
    1:[(56.57942499916897, -5.118316843652349, 36.20041279082784),
 (1.2597963428679804, -1.052099541567668, 4.207064064312306),
 (-0.4902529496568355, 0.1212481034198313, 0.863105406523367)],
    334: [(27.87791807759006, -66.9083845574688, 4.408938902763034),
 (-3.779945559637189, 1.832140739028579, -0.17730919400205936),
 (0.8424716702382711, 0.4112149206908405, 0.34805714164079987)],

    574: [(3.7852716849960375, -42.15629786253971, -30.738206879903267),
 (3.7150864805906227, 1.646715604709498, -1.2647225747167061),
 (-0.9806612072594487, -0.1103369653614953, 0.1616457566756923)], 

}

Neck_nr = 334

all_distances = []

# First loop: Gather min and max values
for (lambda_, h_it) in point_interest_334:
    data_file = f"Necks/Neck{Neck_nr}.stl"
    min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it}.obj"

    data_mesh = pv.read(data_file)
    min_mesh = pv.read(min_file)

    # Compute triangle centroids and areas
    data_mesh["Area"] = data_mesh.compute_cell_sizes()["Area"]
    faces = data_mesh.faces.reshape(-1, 4)[:, 1:]
    centroids = np.mean(data_mesh.points[faces], axis=1)

    # Compute distances
    tree = cKDTree(min_mesh.points)
    distances, _ = tree.query(centroids)
    weighted_distances = distances * data_mesh["Area"]
    normalized_distances = weighted_distances / np.sum(data_mesh["Area"])

    all_distances.extend(weighted_distances)

# Determine global min and max
vmin, vmax = min(all_distances), max(all_distances)
print(vmin)
print(vmax)
for (lambda_,h_it) in point_interest_334:
    
    data_file = f"Necks/Neck{Neck_nr}.stl"
    min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it}.obj"

    data_mesh = pv.read(data_file)
    min_mesh = pv.read(min_file)

    # Compute triangle centroids and areas
    data_mesh["Area"] = data_mesh.compute_cell_sizes()["Area"]
    faces = data_mesh.faces.reshape(-1, 4)[:, 1:]  # Extract triangle vertex indices
    centroids = np.mean(data_mesh.points[faces], axis=1)


    # Compute distances from centroids to nearest min_mesh points
    tree = cKDTree(min_mesh.points)
    distances, _ = tree.query(centroids)
    weighted_distances = distances * data_mesh["Area"]

    # Normalize weighted distances
    data_mesh["distances"] = weighted_distances

    # Distance visualization
    plotter = pv.Plotter(image_scale=2, off_screen=True)
    plotter.add_mesh(data_mesh, scalars="distances", cmap="viridis", clim=[vmin, vmax],  show_edges=True)
    plotter.add_mesh(min_mesh, opacity=0.7)
    camera_pos = camera_positions.get(Neck_nr)
    if camera_pos:
        plotter.camera_position = camera_pos
    #plotter.show()

    # Print the current camera position after adjustment
    #print(plotter.camera_position)





    # Manually set the camera position (replace with your printed values)

    # Display the plot
    
    plotter.screenshot(f"Plots/Screenshots_of_interest/Neck_{Neck_nr}_{lambda_}_{h_it}.png")


Neck_nr = 574

all_distances_1 = []

# First loop: Gather min and max values
for (lambda_, h_it) in point_interest_574:
    data_file = f"Necks/Neck{Neck_nr}.stl"
    min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it}.obj"

    data_mesh = pv.read(data_file)
    min_mesh = pv.read(min_file)

    # Compute triangle centroids and areas
    data_mesh["Area"] = data_mesh.compute_cell_sizes()["Area"]
    faces = data_mesh.faces.reshape(-1, 4)[:, 1:]
    centroids = np.mean(data_mesh.points[faces], axis=1)

    # Compute distances
    tree = cKDTree(min_mesh.points)
    distances, _ = tree.query(centroids)
    weighted_distances = distances * data_mesh["Area"]
    normalized_distances = weighted_distances / np.sum(data_mesh["Area"])

    all_distances_1.extend(weighted_distances)

# Determine global min and max
vmin, vmax = min(all_distances_1), max(all_distances_1)
print(vmin)
print(vmax)
for (lambda_,h_it) in point_interest_574:

    data_file = f"Necks/Neck{Neck_nr}.stl"
    data_mesh = pv.read(data_file)
    triangle_areas = data_mesh.compute_cell_sizes()["Area"]

    # Extract the faces from the data_mesh and compute centroids
    faces = data_mesh.faces.reshape(-1, 4)[:, 1:]
    vertices_data = data_mesh.points[faces]
    centroids = np.mean(vertices_data, axis=1)

    min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it}.obj"
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

    data_mesh['distances'] = weighted_distances

    # Create a PyVista plotter
    plotter = pv.Plotter(image_scale=2, off_screen=True)

    # Add the data mesh and map the 'distance' array as the color
    plotter.add_mesh(data_mesh, scalars='distances', cmap='viridis',clim=[vmin, vmax], show_edges=True)
    plotter.add_mesh(min_mesh, opacity=0.7)
    # Add a color bar
    # Show the plot interactively and adjust orientation manually
    #plotter.show()

    # Print the current camera position after adjustment
    #print(plotter.camera_position)



    camera_pos = camera_positions.get(Neck_nr)
    if camera_pos:
        plotter.camera_position = camera_pos

    # Manually set the camera position (replace with your printed values)

    # Display the plot
    
    plotter.screenshot(f"Plots/Screenshots_of_interest/Neck_{Neck_nr}_{lambda_}_{h_it}.png")


Neck_nr = 1

all_distances = []

# First loop: Gather min and max values
for (lambda_, h_it) in point_interest_574:
    data_file = f"Necks/Neck{Neck_nr}.obj"
    min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it}.obj"

    data_mesh = pv.read(data_file)
    min_mesh = pv.read(min_file)

    # Compute triangle centroids and areas
    data_mesh["Area"] = data_mesh.compute_cell_sizes()["Area"]
    faces = data_mesh.faces.reshape(-1, 4)[:, 1:]
    centroids = np.mean(data_mesh.points[faces], axis=1)

    # Compute distances
    tree = cKDTree(min_mesh.points)
    distances, _ = tree.query(centroids)
    weighted_distances = distances * data_mesh["Area"]
    normalized_distances = weighted_distances / np.sum(data_mesh["Area"])

    all_distances.extend(weighted_distances)

# Determine global min and max
vmin, vmax = min(all_distances), max(all_distances)
print(vmin)
print(vmax)
for (lambda_,h_it) in point_interest_1:

    data_file = f"Necks/Neck{Neck_nr}.obj"
    data_mesh = pv.read(data_file)
    triangle_areas = data_mesh.compute_cell_sizes()["Area"]

    # Extract the faces from the data_mesh and compute centroids
    faces = data_mesh.faces.reshape(-1, 4)[:, 1:]
    vertices_data = data_mesh.points[faces]
    centroids = np.mean(vertices_data, axis=1)

    min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it}.obj"
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

    data_mesh['distances'] = weighted_distances

    # Create a PyVista plotter
    plotter = pv.Plotter(image_scale=2, off_screen=True)

    # Add the data mesh and map the 'distance' array as the color
    plotter.add_mesh(data_mesh, scalars='distances', cmap='viridis',clim=[vmin, vmax], show_edges=True)
    plotter.add_mesh(min_mesh, opacity=0.7)
    # Add a color bar
    # Show the plot interactively and adjust orientation manually
    #plotter.show()

    # Print the current camera position after adjustment
    #print(plotter.camera_position)



    camera_pos = camera_positions.get(Neck_nr)
    if camera_pos:
        plotter.camera_position = camera_pos

    # Manually set the camera position (replace with your printed values)

    # Display the plot
    
    plotter.screenshot(f"Plots/Screenshots_of_interest/Neck_{Neck_nr}_{lambda_}_{h_it}.png")