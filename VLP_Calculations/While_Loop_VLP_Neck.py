import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os
from sklearn.decomposition import PCA

from sklearn.cluster import DBSCAN
from scipy.optimize import minimize
r = '10'



pca = PCA(n_components=2)

Sample = 3
Sample_path = f'Results/Sample_{Sample}'
## sample loop:
while os.path.exists(Sample_path):
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

            data = pv_mesh.cell_centers().points
            gaussian_curv = pv_mesh.cell_data['Gaussian_Curvature_smoothed'] 


            curv_mask = gaussian_curv <= -0.0005
            data_curv = data[curv_mask]

            dbscan = DBSCAN(eps=3, min_samples=50).fit(data_curv)
            labels = dbscan.labels_

            # Find unique clusters and their sizes (excluding noise)
            unique_labels, counts = np.unique(labels[labels != -1], return_counts=True)

            # Sort clusters by size in descending order
            sorted_indices = np.argsort(counts)[::-1]
            sorted_labels = unique_labels[sorted_indices]
            sorted_counts = counts[sorted_indices]

            # The largest cluster label and size
            largest_cluster_label = sorted_labels[0]
            largest_cluster_size = sorted_counts[0]

            neck_mask = labels == largest_cluster_label
            # Retrieve the points belonging to the largest cluster
            data_neck = data_curv[neck_mask]

            mask_data_to_neck = np.zeros(len(curv_mask), dtype=bool)

            indices_in_neck = np.where(neck_mask)[0]

            for idx in indices_in_neck:
                mask_data_to_neck[np.where(curv_mask)[0][idx]] = True

            data_neck_vtk = np.zeros(len(data))
            data_neck_vtk[mask_data_to_neck] = 1

            Cell_areas = pv_mesh.compute_cell_sizes(length=False, area=True, volume=False)['Area']
            Cell_areas_neck = Cell_areas[mask_data_to_neck]
            Neck_Area = np.sum(Cell_areas_neck)

            Energy_per_Cell_kBT = pv_mesh.cell_data['Energy_per_Cell_kBT']
            Energy_per_Area_kBT = pv_mesh.cell_data['Energy_per_square_nm_kBT']
            Energy_per_Cell_kBT_neck = Energy_per_Cell_kBT[mask_data_to_neck]
            Neck_Energy = np.sum(Energy_per_Cell_kBT_neck)



            pv_mesh.cell_data['Neck'] = data_neck_vtk

            pv_mesh.field_data['Neck_Area_square_nm'] = Neck_Area
            pv_mesh.field_data['Energy_in_Neck_per_square_nm_kBT'] = Neck_Energy/Neck_Area




            pv_mesh.save(f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk')

            VLP +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'
        
        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'


    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'


