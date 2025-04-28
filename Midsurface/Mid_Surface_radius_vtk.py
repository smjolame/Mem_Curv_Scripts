import numpy as np
from sklearn.neighbors import KDTree
from sklearn.decomposition import PCA
from stl import mesh
from multiprocessing import Pool
import sys
import time
import pyvista as pv

n = 2500
sys.setrecursionlimit(n)

TS = '002'
Budd = '1'

file_path = f'TS{TS}/mid_surface_TS{TS}_Budding{Budd}.vtk'



pv_mesh = pv.read(file_path)
stl_mesh = mesh.Mesh.from_file(f'TS{TS}/Budding{Budd}.stl')


# load data
data = stl_mesh.centroids * 1e9
vectors = stl_mesh.vectors * 1e9
normals_mesh = stl_mesh.normals * 1e9

H = pv_mesh.cell_data['Mean_Curvature_smoothed']
K = pv_mesh.cell_data['Gaussian_Curvature_smoothed']
normals_mesh = normals_mesh/np.linalg.norm(normals_mesh, axis=1, keepdims=True)

#### number of data points
N_data = len(data)
print('# Points: ', N_data)


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
    H_NN = H[ind_NN]
    K_NN = K[ind_NN]
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
            print(ind)
            return np.nan, np.nan , np.nan, np.nan

    # get the difference vector of the ind point and the closest not conected neighbor
    diff_vector = data_NN[closest_not_connected_ind] - data_NN[0]
    thickness = np.linalg.norm(diff_vector)
    H_mid = H_NN[0]
    K_mid = K_NN[0]

    # if the normal points in a simuilar direction as the difference vector -> flip the difference vector
    if np.dot(diff_vector, normals_mesh_smoothed[ind]) > 0:
        diff_vector *=-1

    mid_point = data_NN[0] + diff_vector/2
    return mid_point, thickness , H_mid, K_mid




start = time.time()
mid_surface = np.empty_like(data)
thicknesses = np.empty(N_data)
H_mids = np.empty(N_data)
K_mids = np.empty(N_data)

# Use multiprocessing Pool to parallelize the computation
with Pool() as pool:
    results = pool.map(calc_mid_point_radius, range(N_data))

# Store results 
for i, (mid_point, thickness, H_mid, K_mid) in enumerate(results):
    mid_surface[i] = mid_point
    thicknesses[i] = thickness
    H_mids[i] = H_mid
    K_mids[i] = K_mid


mid_surface_pcd_radius = pv.PolyData(mid_surface)
mid_surface_pcd_radius.point_data['Thickness'] = thicknesses
mid_surface_pcd_radius.point_data['Mean_Curvature'] = H_mids
mid_surface_pcd_radius.point_data['Gaussian_Curvature'] = K_mids
mid_surface_pcd_radius.point_data['Normals'] = normals_mesh_smoothed
mid_surface_pcd_radius.save(f'mid_surface_Budding{Budd}_{r}_radius_Curv.vtk')

end = time.time()

print('Time mid surface radius:', end - start)


