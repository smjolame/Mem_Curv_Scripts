
import pyvista as pv
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

letters = ["a","b","c","d","e","m"]
n = 50
lambda_ = 15
h_it = 23
# Parameters
point_interest_1 = [[2,5],[5,26],[15,34], [15,14],[15,35], [50,17]]
point_interest_334 = [[2,5],[5,26],[15,33], [15,16],[15,34], [50,14]]
point_interest_574 = [[2,5],[5,24],[15,29], [15,16],[15,30], [50,15]]

#point_interest_334 = [[2,5],[5,26],[15,33],[50,18]]
#point_interest_574 = [[2,5],[5,24],[15,29],[50,14]]



camera_positions = {
    334: [(27.87791807759006, -66.9083845574688, 4.408938902763034),
 (-3.779945559637189, 1.832140739028579, -0.17730919400205936),
 (0.8424716702382711, 0.4112149206908405, 0.34805714164079987)],

    574: [(3.7852716849960375, -42.15629786253971, -30.738206879903267),
 (3.7150864805906227, 1.646715604709498, -1.2647225747167061),
 (-0.9806612072594487, -0.1103369653614953, 0.1616457566756923)], 

}
plotter = pv.Plotter(shape=(len(point_interest_334), 2), off_screen=True) 
Neck_nr = 334
for i, (lambda_, h_it) in enumerate(point_interest_334):
    letter = letters[i]

    
    data_file = f"Necks/Neck{Neck_nr}.stl"
    min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it}.obj"

    data_mesh = pv.read(data_file)
    min_mesh = pv.read(min_file)

    # Compute triangle centroids and areas
    data_mesh["Area"] = data_mesh.compute_cell_sizes()["Area"]
    faces = data_mesh.faces.reshape(-1, 4)[:, 1:]  # Extract triangle vertex indices
    centroids = np.mean(data_mesh.points[faces], axis=1)

    # Perform PCA to align principal axes
    pca = PCA(n_components=3)
    pca.fit(centroids)
    R = np.vstack(pca.components_)  # Rotation matrix

    # Rotate both meshes
    data_mesh.points[:] = data_mesh.points @ R.T
    min_mesh.points[:] = min_mesh.points @ R.T

    # Generate slices
    xz_slice_data = data_mesh.slice(normal=(0, 1, 0))  # XZ slice
    yz_slice_data = data_mesh.slice(normal=(1, 0, 0))  # YZ slice
    xz_slice_min = min_mesh.slice(normal=(0, 1, 0))
    yz_slice_min = min_mesh.slice(normal=(1, 0, 0))

    # Compute distances from centroids to nearest min_mesh points
    tree = cKDTree(min_mesh.points)
    distances, _ = tree.query(centroids)
    weighted_distances = distances * data_mesh["Area"]

    data_mesh["distances"] = weighted_distances


    # XZ slice
    plotter.subplot(i, 0)
    plotter.add_mesh(xz_slice_data, color="blue", label="neck data",line_width=3)
    plotter.add_mesh(xz_slice_min, color="red", label="min. surface",line_width=3)
    plotter.add_legend(size=(0.3,0.3), loc='center left', face='none')
    plotter.camera_position = [(0, -35, 0), (0, 0, 0), (0, 0, 1)]  # Side view
    plotter.add_text(rf"Slice ({letter})", position="upper_edge", font_size=12, color="black")
    plotter.add_axes(shaft_length=1, ylabel='',label_size=(0.5, 0.5))

    # YZ slice
    plotter.subplot(i, 1)
    plotter.add_mesh(yz_slice_data, color="blue", label="neck data",line_width=3)
    plotter.add_mesh(yz_slice_min, color="red", label="min. surface",line_width=3)
    plotter.add_legend(size=(0.3,0.3), loc='center left', face='none')
    plotter.camera_position = [(35, 0, 0), (0, 0, 0), (0, 0, 1)]  # Front view
    plotter.add_text(rf"Slice ({letter})", position="upper_edge", font_size=12, color="black")
    plotter.add_axes(shaft_length=1, xlabel='',label_size=(0.5, 0.5))
#plotter.show()
plotter.window_size = [1650, 1188]
plotter.save_graphic(f"Plots/Screenshots/Slices/Neck_{Neck_nr}_slices.pdf")


plotter = pv.Plotter(shape=(len(point_interest_574), 2), off_screen=True) 
Neck_nr = 574
for i, (lambda_, h_it) in enumerate(point_interest_574):
    letter = letters[i]
    
    data_file = f"Necks/Neck{Neck_nr}.stl"
    min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it}.obj"

    data_mesh = pv.read(data_file)
    min_mesh = pv.read(min_file)

    # Compute triangle centroids and areas
    data_mesh["Area"] = data_mesh.compute_cell_sizes()["Area"]
    faces = data_mesh.faces.reshape(-1, 4)[:, 1:]  # Extract triangle vertex indices
    centroids = np.mean(data_mesh.points[faces], axis=1)

    # Perform PCA to align principal axes
    pca = PCA(n_components=3)
    pca.fit(centroids)
    R = np.vstack(pca.components_)  # Rotation matrix

    # Rotate both meshes
    data_mesh.points[:] = data_mesh.points @ R.T
    min_mesh.points[:] = min_mesh.points @ R.T

    # Generate slices
    xz_slice_data = data_mesh.slice(normal=(0, 1, 0))  # XZ slice
    yz_slice_data = data_mesh.slice(normal=(1, 0, 0))  # YZ slice
    xz_slice_min = min_mesh.slice(normal=(0, 1, 0))
    yz_slice_min = min_mesh.slice(normal=(1, 0, 0))

    # XZ slice
    plotter.subplot(i, 0)
    plotter.add_mesh(xz_slice_data, color="blue", label="neck data",line_width=3)
    plotter.add_mesh(xz_slice_min, color="red", label="min. surface",line_width=3)
    plotter.add_legend(size=(0.3,0.3), loc='center left', face='none')
    #plotter.camera_position = [(0, -35, 0), (0, 0, 0), (0, 0, 1)]  # Side view
    plotter.camera_position='xz'
    plotter.camera.zoom(2)
    plotter.add_text(rf"Slice ({letter})", position="upper_edge", font_size=12, color="black")
    plotter.add_axes(shaft_length=1, ylabel='',label_size=(0.5, 0.5))

    # YZ slice
    plotter.subplot(i, 1)
    plotter.add_mesh(yz_slice_data, color="blue", label="neck data",line_width=3)
    plotter.add_mesh(yz_slice_min, color="red", label="min. surface",line_width=3)
    plotter.add_legend(size=(0.3,0.3), loc='center left', face='none')
    #plotter.camera_position = [(35, 0, 0), (0, 0, 0), (0, 0, 1)]  # Front view
    plotter.camera_position='yz'
    plotter.camera.zoom(2)
    plotter.add_text(rf"Slice ({letter})", position="upper_edge", font_size=12, color="black")
    plotter.add_axes(shaft_length=1, xlabel='',label_size=(0.5, 0.5))

plotter.window_size = [1650, 1188]
plotter.save_graphic(f"Plots/Screenshots/Slices/Neck_{Neck_nr}_slices.pdf")

camera = pv.Camera()
plotter = pv.Plotter(shape=(len(point_interest_1), 2), off_screen=True) 
Neck_nr = 1
for i, (lambda_, h_it) in enumerate(point_interest_1):
    letter = letters[i]
    
    data_file = f"Necks/Neck{Neck_nr}.obj"
    min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it}.obj"

    data_mesh = pv.read(data_file)
    min_mesh = pv.read(min_file)

    # Compute triangle centroids and areas
    data_mesh["Area"] = data_mesh.compute_cell_sizes()["Area"]
    faces = data_mesh.faces.reshape(-1, 4)[:, 1:]  # Extract triangle vertex indices
    centroids = np.mean(data_mesh.points[faces], axis=1)

    # Perform PCA to align principal axes
    pca = PCA(n_components=3)
    pca.fit(centroids)
    R = np.vstack(pca.components_)  # Rotation matrix

    # Rotate both meshes
    data_mesh.points[:] = data_mesh.points @ R.T
    min_mesh.points[:] = min_mesh.points @ R.T

    # Generate slices
    xz_slice_data = data_mesh.slice(normal=(0, 1, 0))  # XZ slice
    yz_slice_data = data_mesh.slice(normal=(1, 0, 0))  # YZ slice
    xz_slice_min = min_mesh.slice(normal=(0, 1, 0))
    yz_slice_min = min_mesh.slice(normal=(1, 0, 0))

    # XZ slice
    plotter.subplot(i, 0)
    plotter.add_mesh(xz_slice_data, color="blue", label="neck data",line_width=3)
    plotter.add_mesh(xz_slice_min, color="red", label="min. surface",line_width=3)
    plotter.add_legend(size=(0.3,0.3), loc='center left', face='none')
    #plotter.camera_position = [(0, -40, 0), (0, 0, 8), (0, 0, 1)]  # Side view
    plotter.camera_position='xz'
    plotter.camera.zoom(2)
    plotter.add_text(rf"Slice ({letter})", position="upper_edge", font_size=12, color="black")
    plotter.add_axes(shaft_length=1, ylabel='',label_size=(0.5, 0.5))

    # YZ slice
    plotter.subplot(i, 1)
    plotter.add_mesh(yz_slice_data, color="blue", label="neck data",line_width=3)
    plotter.add_mesh(yz_slice_min, color="red", label="min. surface",line_width=3)
    plotter.add_legend(size=(0.3,0.3), loc='center left', face='none')
    #plotter.camera_position = [(40, 0, 0), (0, 0, 8), (0, 0, 1)]  # Front view
    plotter.camera_position='yz'
    plotter.camera.zoom(2)
    plotter.add_text(rf"Slice ({letter})", position="upper_edge", font_size=12, color="black")
    plotter.add_axes(shaft_length=1, xlabel='',label_size=(0.5, 0.5))

plotter.window_size = [1650, 1188]
plotter.save_graphic(f"Plots/Screenshots/Slices/Neck_{Neck_nr}_slices.pdf")