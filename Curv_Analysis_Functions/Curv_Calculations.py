import numpy as np
from sklearn.neighbors import KDTree
from sklearn.decomposition import PCA
from scipy.spatial.transform import Rotation as R
from scipy.optimize import curve_fit 
from scipy import constants 
from sympy import symbols, diff, simplify, lambdify
from stl import mesh
from multiprocessing import Pool
import sys
import time
import pyvista as pv
import os


n = 2500
sys.setrecursionlimit(n)

Temp = 300
k_B = constants.k
kappa = 20*k_B*Temp
kappa_bar = -1/3*kappa

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



Sample = 3
TS = 3
VLP = 3
file_path = f'VLP2.stl'



your_mesh = mesh.Mesh.from_file(file_path)


data = your_mesh.centroids * 1e9
vectors = your_mesh.vectors * 1e9
normals_mesh = your_mesh.normals * 1e9
normals_mesh = normals_mesh/np.linalg.norm(normals_mesh, axis=1, keepdims=True)




###### used radius
r = 10
#################################################################

#### number of data points
N_data = len(data)
print('# Points: ', N_data)


##### global variables:
tree_ind_smooth = None
dist_smooth = None
r_smooth = None

########### with KDTree:

# KD Tree
######################################
start = time.time()

tree = KDTree(data, metric='minkowski')
tree_ind , dist = tree.query_radius(data,r=r, sort_results=True, return_distance=True)
#dist, tree_ind = tree.query(data, sort_results=True, return_distance=True, k=2500) #k points including itself



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

normals_mesh_smoothed = smoothing_mesh_normals(data, normals_mesh, r=1.5)


##### Calcualtion of H and K
#### fit function

def func(xy, a, b, c, d, e, f): 
    x, y = xy 
    return a*x**2 + b*y**2 + c*x*y + d*x + e*y + f

x_sym, y_sym, z_sym = symbols('x y z')
a_sym, b_sym, c_sym, d_sym, e_sym, f_sym= symbols('a b c d e f')
func_symbol = a_sym*x_sym**2 + b_sym*y_sym**2 + c_sym*x_sym*y_sym + d_sym*x_sym + e_sym*y_sym + f_sym
f_x =  diff(func_symbol, x_sym)
f_y =  diff(func_symbol, y_sym)
f_yy = diff(func_symbol, y_sym, y_sym)
f_xx = diff(func_symbol, x_sym, x_sym)
f_xy = diff(func_symbol, x_sym, y_sym)

H_symbol = (f_xx*(1+f_y**2) - 2*f_xy*f_x*f_y + f_yy*(1+f_x**2)) / (2*(1+f_x**2+f_y**2)**(3/2))
H_symbol = simplify(H_symbol)
H_func = lambdify([[x_sym,y_sym],a_sym,b_sym,c_sym,d_sym,e_sym,f_sym], H_symbol)

K_symbol = (f_xx*f_yy-f_xy**2)/(1+f_x**2+f_y**2)**2
K_symbol = simplify(K_symbol)
K_func = lambdify([[x_sym,y_sym],a_sym,b_sym,c_sym,d_sym,e_sym,f_sym], K_symbol)

'''
def func(xy, a, b, c, d, e, f, g, h, i, j): 
    x, y = xy 
    return g*x**3 + h*y**3 + a*x**2 + b*y**2 + c*x*y + d*x**2*y + e*y**2*x + f + i*x + j*y

x_sym, y_sym, z_sym = symbols('x y z')
a_sym, b_sym, c_sym, d_sym, e_sym, f_sym, g_sym, h_sym, i_sym, j_sym= symbols('a b c d e f g h i j')
func_symbol = g_sym*x_sym**3 + h_sym*y_sym**3 + a_sym*x_sym**2 + b_sym*y_sym**2 + c_sym*x_sym*y_sym + d_sym*x_sym**2*y_sym + e_sym*y_sym**2*x_sym + f_sym + i_sym*x_sym + j_sym*y_sym
f_x =  diff(func_symbol, x_sym)
f_y =  diff(func_symbol, y_sym)
f_yy = diff(func_symbol, y_sym, y_sym)
f_xx = diff(func_symbol, x_sym, x_sym)
f_xy = diff(func_symbol, x_sym, y_sym)

H_symbol = (f_xx*(1+f_y**2) - 2*f_xy*f_x*f_y + f_yy*(1+f_x**2)) / (2*(1+f_x**2+f_y**2)**(3/2))
H_symbol = simplify(H_symbol)
H_func = lambdify([[x_sym,y_sym],a_sym,b_sym,c_sym,d_sym,e_sym,f_sym,g_sym,h_sym, i_sym, j_sym], H_symbol)

K_symbol = (f_xx*f_yy-f_xy**2)/(1+f_x**2+f_y**2)**2
K_symbol = simplify(K_symbol)
K_func = lambdify([[x_sym,y_sym],a_sym,b_sym,c_sym,d_sym,e_sym,f_sym,g_sym,h_sym, i_sym, j_sym], K_symbol)

'''
############################# Calc. Normals with PCA.-> Rotate data -> Fit:
z_axis = np.array([0,0,1])

start = time.time()


def calc_curvatures(ind):
    ind_NN = tree_ind[ind]  # indices of the NN
    dist_NN = dist[ind]  # distances of NN region

    # calculate mesh neighbors in r-environment (so divide the different membranes) and use only that connected mesh
    ind_components = find_connected_components(neighbours, ind_NN)
    data_NN_single = data[ind_components]
    

    # avoid overhangs via calculating the saclarproduct between the considered normal and all other normals of the mesh
    normals_mesh_NN_single = normals_mesh_smoothed[ind_components]
    scalar_products = np.dot(normals_mesh_NN_single, normals_mesh_smoothed[ind])
    mask_single_membrane_sc = scalar_products > 0

    data_single_sc = data_NN_single[mask_single_membrane_sc]
    data_single_sc = data_single_sc - data[ind]


    if len(data_single_sc) < 20:
        H_value = np.nan
        K_value = np.nan
        return H_value, K_value

    

    normal = normals_mesh_smoothed[ind]
    

    rot, _ = R.align_vectors(z_axis, normal)
    data_local_coordinates = rot.apply(data_single_sc - data[ind])
    normal_rot = rot.apply(normal)
    if np.dot(normal_rot,z_axis) < 0:
        print('warning, rot failed: ind', ind,'Normal:', normal, 'Norm:', np.linalg.norm(normal), 'sc:', np.dot(normal_rot,z_axis))
        data_rot = data_single_sc - data[ind]
        data_rot = data_rot * np.array([1, 1, -1])

    x = data_local_coordinates[:, 0]
    y = data_local_coordinates[:, 1]
    z = data_local_coordinates[:, 2]

    # get the distances of the edited data points
    mask = np.isin(ind_NN, ind_components)
    dist_NN_single_sc = dist_NN[mask][mask_single_membrane_sc]

    #### calculating the weights:
    sigma_value = r  # use the Radius as distance parameter
    weights = np.exp(-(dist_NN_single_sc) ** 2 / (2 * sigma_value ** 2))
    sigma = 1.0 / weights

    popt, pcov = curve_fit(func, (x, y), z, sigma=sigma)
    X = x[0]
    Y = y[0]
    H_value = H_func((X, Y), *popt)
    K_value = K_func((X, Y), *popt)
    
    return H_value, K_value

def smooth_curvatures(ind):
    ind_NN = tree_ind[ind]  # indices of the NN
    dist_NN = dist[ind]  # distances of NN region
    normals_mesh_NN = normals_mesh_smoothed[ind_NN]  # mesh normals in NN region

    # split membrane via comparison between direction of pivot mesh normal and all NN mesh normals
    scalar_products = np.dot(normals_mesh_NN, normals_mesh_smoothed[ind])
    mask_single_membrane_sc = scalar_products > 0
    ind_NN_single = ind_NN[mask_single_membrane_sc]
    dist_NN_single = dist_NN[mask_single_membrane_sc]

    sigma = r_smooth/2
    w_i_array = np.exp(-(dist_NN_single)**2/(2*sigma**2))
    H_smoothed = np.sum(w_i_array*H[ind_NN_single])/(np.sum(w_i_array))
    K_smoothed = np.sum(w_i_array*K[ind_NN_single])/(np.sum(w_i_array))


    
    return H_smoothed, K_smoothed


# Initialize arrays to store results
H = np.empty(N_data)
K = np.empty(N_data)

# Use multiprocessing Pool to parallelize the computation
with Pool() as pool:
    results = pool.map(calc_curvatures, range(N_data))

# Store results in H and K arrays
for i, (H_value, K_value) in enumerate(results):
    H[i] = H_value
    K[i] = K_value

end = time.time()

print('Time Curvatue:', end - start)






# smoothing
####################################################


start = time.time()
r_smooth = 2

# KD Tree
######################################


tree = KDTree(data, metric='minkowski')
tree_ind , dist = tree.query_radius(data,r=r_smooth, sort_results=True, return_distance=True)

######################################


# Initialize arrays to store results
H_smoothed = np.empty(N_data)
K_smoothed = np.empty(N_data)

# Use multiprocessing Pool to parallelize the computation
with Pool() as pool:
    results = pool.map(smooth_curvatures, range(N_data))

# Store results in H and K arrays
for i, (H_value, K_value) in enumerate(results):
    H_smoothed[i] = H_value
    K_smoothed[i] = K_value

end = time.time()

print('Time Curvatue smoothed:', end - start)

points = vectors.reshape(-1, 3)
faces = np.arange(points.shape[0]).reshape(-1, 3)
faces = np.hstack([np.full((faces.shape[0], 1), 3), faces])

# Create the PyVista mesh
pv_mesh = pv.PolyData(points, faces)

# Energy Calculations
def Hamiltonian_per_Cell(H,K,Cell_area):
    return (2*kappa*H**2+kappa_bar*K)*Cell_area

Cell_areas = pv_mesh.compute_cell_sizes(length=False, area=True, volume=False)['Area']
    
Energy_per_Cell = Hamiltonian_per_Cell(H_smoothed,K_smoothed,Cell_areas)
Energy_per_Area = Energy_per_Cell / Cell_areas
Energy_per_Cell_kBT = Energy_per_Cell/(k_B*Temp)
Energy_per_Area_kBT = Energy_per_Area/(k_B*Temp)

pv_mesh.cell_data['Energy_per_Cell_kBT'] = Energy_per_Cell_kBT
pv_mesh.cell_data['Energy_per_square_nm_kBT'] = Energy_per_Area_kBT





##### not smoothed
#data_curvatures = np.column_stack((data, H, K))

#np.savetxt(f'Curvatures_Results/TS00{TS}_Obj_{OBJ}_r_{r}_Curvatures.csv', data_curvatures)


# Add the curvature data to the mesh
pv_mesh.cell_data['Mean_Curvature'] = H
pv_mesh.cell_data['Gaussian_Curvature'] = K

##### smoothed
#data_curvatures = np.column_stack((data, H_smoothed, K_smoothed))

#np.savetxt(f'Curvatures_Results/TS00{TS}_Obj_{OBJ}_r_{r}_Curvatures_smoothed.csv', data_curvatures)

# Add the curvature data to the mesh
pv_mesh.cell_data['Mean_Curvature_smoothed'] = H_smoothed
pv_mesh.cell_data['Gaussian_Curvature_smoothed'] = K_smoothed

pv_mesh.cell_data['Normals'] = normals_mesh_smoothed
            


pv_mesh.save(f'VLP2.vtk')

