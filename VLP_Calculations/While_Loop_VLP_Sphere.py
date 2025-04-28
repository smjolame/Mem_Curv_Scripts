import numpy as np
import pyvista as pv
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from Own_Functions import set_axes_equal
from scipy.optimize import minimize
import os


r = '10'

def objective(params, points):
    Cx, Cy, Cz, R = params
    distances = np.sqrt((points[:, 0] - Cx)**2 + (points[:, 1] - Cy)**2 + (points[:, 2] - Cz)**2)
    return np.sum((distances - R)**2)

def fit_sphere(points):
    # Initial guess for the center and radius
    initial_guess = np.mean(points, axis=0).tolist() + [np.mean(np.linalg.norm(points - np.mean(points, axis=0), axis=1))]

    
    # Optimize the objective function
    result = minimize(objective, initial_guess, args=(points,), method='L-BFGS-B')
    
    Cx, Cy, Cz, R = result.x
    return (Cx, Cy, Cz), R

def ellipsoid_objective(params, points):
    Cx, Cy, Cz, a, b, c = params
    distances = ((points[:, 0] - Cx) / a) ** 2 + ((points[:, 1] - Cy) / b) ** 2 + ((points[:, 2] - Cz) / c) ** 2
    return np.sum((distances - 1) ** 2)

def fit_ellipsoid(points):
    # Initial guess for the center and axis lengths
    initial_guess = np.mean(points, axis=0).tolist() + [np.mean(np.linalg.norm(points - np.mean(points, axis=0), axis=1)) for i in range(3)]
    
    # Optimize the objective function
    result = minimize(ellipsoid_objective, initial_guess, args=(points,), method='L-BFGS-B')
    
    Cx, Cy, Cz, a, b, c = result.x
    return (Cx, Cy, Cz), (a, b, c)

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
        Neck = True

        while os.path.exists(file_path):
            print('Sample: ', Sample, 'TS: ', TS, 'VLP: ', VLP)
            pv_mesh = pv.read(file_path)

            data = pv_mesh.cell_centers().points
            mean_curv = pv_mesh.cell_data['Mean_Curvature_smoothed'] 


            curv_mask = mean_curv <= -0.015
            data_curv = data[curv_mask]

            dbscan = DBSCAN(eps=3, min_samples=60).fit(data_curv)
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

            sphere_mask = labels == largest_cluster_label
            # Retrieve the points belonging to the largest cluster
            data_sphere = data_curv[sphere_mask]

            '''
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(*data_curv.T, c=labels, cmap='viridis', alpha=0.6)
            ax.scatter(*data_sphere.T, marker='x',linewidth=2, c = 'r')
            #plt.show()
            plt.close()
            '''

            center, radius = fit_sphere(data_sphere)
            #center_el, axes_lengths = fit_ellipsoid(data_sphere)

            #a = axes_lengths[2]
            #b = axes_lengths[0]
            #ecc = np.sqrt(1-(b/a)**2)

            #a, b, c = axes_lengths
            '''
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(*data_sphere[::20].T, marker='x',linewidth=2,alpha=0.5 ,c = 'c')
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x = radius*np.cos(u)*np.sin(v) + center[0]
            y = radius*np.sin(u)*np.sin(v) + center[1]
            z = radius*np.cos(v) + center[2]
            ax.plot_wireframe(x, y, z, color="b", alpha = 0.3, label = 'Sphere')
            '''
            '''    
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x = a * np.outer(np.cos(u), np.sin(v)) + center_el[0]
            y = b * np.outer(np.sin(u), np.sin(v)) + center_el[1]
            z = c * np.outer(np.ones_like(u), np.cos(v)) + center_el[2]
            ax.plot_wireframe(x, y, z, color="r", alpha = 0.3, label = 'Ellipse')
        '''
            #set_axes_equal(ax)
            #plt.legend()
            #plt.show()
            #plt.close()

            mask_data_to_sphere = np.zeros(len(curv_mask), dtype=bool)

            indices_in_sphere = np.where(sphere_mask)[0]

            for idx in indices_in_sphere:
                mask_data_to_sphere[np.where(curv_mask)[0][idx]] = True


            ### Data for fit:
            data_sphere_vtk = np.zeros(len(data))
            data_sphere_vtk[mask_data_to_sphere] = 1

            ## Areas:
            Cell_areas = pv_mesh.compute_cell_sizes(length=False, area=True, volume=False)['Area']
            Cell_areas_sphere = Cell_areas[mask_data_to_sphere]
            Area_Sphere = np.sum(Cell_areas_sphere)
            Area_Sphere_fit = np.pi*4*radius**2


            ### Energy:
            Energy_per_Cell_kBT = pv_mesh.cell_data['Energy_per_Cell_kBT']
            Energy_per_Area_kBT = pv_mesh.cell_data['Energy_per_square_nm_kBT']
            Energy_per_Cell_kBT_Sphere = Energy_per_Cell_kBT[mask_data_to_sphere]
            Sphere_Energy = np.sum(Energy_per_Cell_kBT_Sphere)  # Gesamtenergy in dem Datenbereich
            Sphere_Energy_Fit = Sphere_Energy/Area_Sphere*Area_Sphere_fit #Gesamtenergie pro Fläche in dem Datenbereich mal Fläche der gefitteten ganzen Kugel


            pv_mesh.cell_data['Sphere'] = data_sphere_vtk

            pv_mesh.field_data['Radius_Sphere_fit'] = radius
            #pv_mesh.field_data['Eccentricity_Ellipse_fit'] = ecc
            pv_mesh.field_data['Area_Sphere_fit'] = Area_Sphere_fit
            pv_mesh.field_data['Energy_in_Sphere_extrapolated_per_square_nm_fit_kBT'] = Sphere_Energy_Fit/Area_Sphere_fit
            pv_mesh.field_data['Energy_in_Sphere_per_square_nm_kBT'] = Sphere_Energy/Area_Sphere


            if Neck == True:
                pv_mesh.save(f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk')
            else:
                pv_mesh.save(f'Results/Sample_{Sample}/TS00{TS}/noNeck/VLP_{VLP}.vtk')


            VLP +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'

            if not os.path.exists(file_path):
                Neck = False

            if Neck == False:
                file_path = f'Results/Sample_{Sample}/TS00{TS}/noNeck/VLP_{VLP}.vtk'
        
        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'


    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'