import pyvista as pv
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os
# Parameters
grid_size = 50
h_zero_min = -0.3
h_zero_max = 0.6
Necks=np.array([331, 332, 334, 335, 363, 571, 574])
lambda_ = 50


interval_number = 2

if interval_number == 0:
    h_zero_min, h_zero_max = -0.1, 0.4
elif interval_number == 1:
    h_zero_min, h_zero_max = -0.2, 0.2
elif interval_number == 2:
    h_zero_min, h_zero_max = -0.1, 0.1

h_array = np.linspace(h_zero_min,h_zero_max,50)
h_it = np.arange(25,49,2)


# Parameters



for Neck_nr in Necks:
    plotter = pv.Plotter(shape=(6, 2), off_screen=True) 
    for i, (j, k) in enumerate(np.ndindex(6, 2)):
        
        
        
        data_file = f"Necks/Neck{Neck_nr}.stl"
        min_file = f"Min_Necks/Min_Necks_temp/Neck{Neck_nr}_{lambda_}_{h_it[i]}_{interval_number}.obj"
        if os.path.exists(min_file):
            data_mesh = pv.read(data_file)
            min_mesh = pv.read(min_file)

            # Compute triangle centroids and areas
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
            xz_slice_min = min_mesh.slice(normal=(0, 1, 0))

            # XZ slice
        
            plotter.subplot(j, k)
            plotter.add_mesh(xz_slice_data, color="blue", label="neck data",line_width=3)
            plotter.add_mesh(xz_slice_min, color="red", label="min. surface",line_width=3)
            plotter.add_legend(size=(0.3,0.3), loc='center left', face='none')
            plotter.camera_position = [(0, -40, 0), (0, 0, 0), (0, 0, 1)]  # Side view
            plotter.add_text(rf"H_0: %.3f" % h_array[h_it[i]], position="upper_left", font_size=8, color="black")

    #plotter.show()
    plotter.window_size = [1650, 1188]
    plotter.save_graphic(f"Plots/Screenshots/Slices/Neck_{Neck_nr}_{interval_number}_slices.pdf")
