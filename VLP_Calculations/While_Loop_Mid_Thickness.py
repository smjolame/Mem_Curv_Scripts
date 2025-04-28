import numpy as np
from sklearn.neighbors import KDTree
from sklearn.decomposition import PCA
from stl import mesh
from multiprocessing import Pool
import sys
import time
import pyvista as pv#
import os
r = '10'


n = 2500
sys.setrecursionlimit(n)


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
        file_path_stl = f'EM_Data/Sample_{Sample}/TS00{TS}/Neck/VLP{VLP}.stl'
        Neck = True

        while os.path.exists(file_path):
            pv_mesh = pv.read(file_path)
            stl_mesh = mesh.Mesh.from_file(file_path_stl)


            # load data
            data = stl_mesh.centroids * 1e9
            vectors = stl_mesh.vectors * 1e9
            normals_mesh = stl_mesh.normals * 1e9

            normals_mesh = normals_mesh/np.linalg.norm(normals_mesh, axis=1, keepdims=True)

            #### number of data points
            N_data = len(data)
            print('# Points: ', N_data)
            print('Sample: ', Sample, 'TS: ', TS, 'VLP: ', VLP)


            ##### global variables:
            tree_ind_smooth = None
            dist_smooth = None
            r_smooth = None


            #### functions



            def build_kdtree(data, metric='minkowski'):
                return KDTree(data, metric=metric)

            def find_neighbours(triangles, data, k=10):
                tree = build_kdtree(data)
                NN_ind = tree.query(data, return_distance=False, k=k)
                
                triangles_set = [set(map(tuple, tri)) for tri in triangles]
                
                def shared_vertex_count(tri_set1, tri_set2):
                    return len(tri_set1 & tri_set2)
                
                neighbour_list = []
                for i, tri_set in enumerate(triangles_set):
                    neighbours = [
                        ind for ind in NN_ind[i][1:]  # Skip the first one as it's the triangle itself
                        if shared_vertex_count(tri_set, triangles_set[ind]) >= 2
                    ]
                    neighbour_list.append(neighbours)
                
                return neighbour_list


                

            def dfs(neighbours, tri_id, current, visited):
                visited.add(current)
                for neighbor_id in neighbours[current]:
                    if neighbor_id not in visited and neighbor_id in tri_id:
                        dfs(neighbours, tri_id, neighbor_id, visited)




            def find_connected_components(neighbours, tri_id):
                visited = set()

                start_id = tri_id[0]
                tri_id_set = set(tri_id)
                dfs(neighbours, tri_id_set, start_id, visited)
                        
                return list(visited)


            def initialize_smoothing(data, r_smooth):
                global tree_ind_smooth, dist_smooth
                tree_smooth = KDTree(data, metric='minkowski')
                tree_ind_smooth, dist_smooth = tree_smooth.query_radius(data, r=r_smooth, sort_results=True, return_distance=True)



            def smoothing_mesh_normal(ind):

                ind_NN = tree_ind_smooth[ind]  # indices of the NN
                dist_NN = dist_smooth[ind]  # distances of NN region

                ind_components = find_connected_components(neighbours, ind_NN)
                normals_mesh_single = normals_mesh[ind_components]

                mask = np.isin(ind_NN, ind_components)
                dist_NN_single = dist_NN[mask]
                sigma = r_smooth / 2
                w_i_array = np.exp(-(dist_NN_single) ** 2 / (2 * sigma ** 2))
                w_i_array /= np.sum(w_i_array)  # Normalize the weights
                normal_mesh_smoothed = np.dot(w_i_array, normals_mesh_single)

                return normal_mesh_smoothed


            def smoothing_mesh_normals(data, normals_mesh, r=2):
                global r_smooth
                r_smooth = r
                # Initialize the KDTree and neighbors globally
                initialize_smoothing(data, r_smooth)

                N_data = len(data)
                normals_mesh_smoothed = np.empty_like(normals_mesh)

                # Use multiprocessing Pool to parallelize the computation
                with Pool() as pool:
                    results = pool.map(smoothing_mesh_normal, range(N_data))

                # Store results in the normals_mesh_smoothed array
                for i, normal_mesh_smoothed in enumerate(results):
                    normals_mesh_smoothed[i] = normal_mesh_smoothed

                return normals_mesh_smoothed




            ###### used radius
            r = 10
            #################################################################


            ########### approach with KDTree:

            # KD Tree
            ######################################
            start = time.time()

            tree = KDTree(data, metric='minkowski')
            tree_ind , dist = tree.query_radius(data,r=r, sort_results=True, return_distance=True)
            #dist, tree_ind = tree.query(data, sort_results=True, return_distance=True, k=2500) #k points including itself
            from sys import getsizeof
            print('Size of dist in Mb:', getsizeof(dist.base)/1000000)


            end = time.time()
            print('Time KDTree:', end - start)
            ######################################

            # Neighbour Calculation
            ######################################
            start = time.time()

            neighbours = find_neighbours(vectors, data)

            end = time.time()
            print('Time for Neighbour Calc.:', end - start)
            ######################################

            # smoothing of of normals:
            #########

            normals_mesh_smoothed = smoothing_mesh_normals(data, normals_mesh, r=2)



            def calc_mid_point_radius(ind):
                ind_NN = tree_ind[ind]
                data_NN = data[ind_NN]
                normals_mesh_NN = normals_mesh_smoothed[ind_NN]

                # calculate mesh neighbors in r-environment (so divide the different membranes) and use only that unconnected mesh
                ind_components = find_connected_components(neighbours, ind_NN)
                mask_not_connected = np.isin(ind_NN, ind_components, invert=True)
                data_not_connected = data_NN[mask_not_connected]
                closest_not_connected_ind = np.argmax(mask_not_connected)
                
                # in case the ind point has no unconnected mesh in r environment. Use splitting via normals:
                if len(data_not_connected) == 0:
                    scalar_products = np.dot(normals_mesh_NN, normals_mesh_smoothed[ind])
                    mask_single_membrane_sc = scalar_products < 0
                    closest_not_connected_ind = np.argmax(mask_single_membrane_sc)
                    if closest_not_connected_ind == 0:
                        #print(ind)
                        return np.nan

                # get the difference vector of the ind point and the closest not conected neighbor
                diff_vector = data_NN[closest_not_connected_ind] - data_NN[0]
                thickness = np.linalg.norm(diff_vector)

                return thickness 




            start = time.time()
            thicknesses = np.empty(N_data)

            # Use multiprocessing Pool to parallelize the computation
            with Pool() as pool:
                results = pool.map(calc_mid_point_radius, range(N_data))

            # Store results 
            for i, thickness in enumerate(results):
                thicknesses[i] = thickness


            pv_mesh.cell_data['Thickness'] = thicknesses

            #mid_surface_pcd_radius.save(f'mid_surface_Budding{Budd}_{r}_radius_Curv.vtk')
            if Neck == True:
                pv_mesh.save(f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk')
            else:
                pv_mesh.save(f'Results/Sample_{Sample}/TS00{TS}/noNeck/VLP_{VLP}.vtk')

            end = time.time()

            print('Time mid surface radius:', end - start)
            VLP +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'
            file_path_stl = f'EM_Data/Sample_{Sample}/TS00{TS}/Neck/VLP{VLP}.stl'

            if not os.path.exists(file_path):
                Neck = False

            if Neck == False:
                file_path = f'Results/Sample_{Sample}/TS00{TS}/noNeck/VLP_{VLP}.vtk'
                file_path_stl = f'EM_Data/Sample_{Sample}/TS00{TS}/noNeck/VLP{VLP}.stl'
        
        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'


    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'

