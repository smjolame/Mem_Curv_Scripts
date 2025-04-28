import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import Normalize
from sklearn.neighbors import KDTree
from sklearn.decomposition import PCA
from scipy.spatial.transform import Rotation as R
from scipy.optimize import curve_fit 
from sympy import *
from scipy.stats import skewnorm
from Own_Functions import set_axes_equal
import matplotlib as mpl
import pyvista as pv
from multiprocessing import Pool
import sys
import time

# Use LaTeX for rendering math expressions
mpl.rcParams['text.usetex'] = True

# Set font properties
mpl.rcParams['font.family'] = 'serif'  # Serif for the main text
mpl.rcParams['font.serif'] = ['Times New Roman']  # This can be customized to another serif font
mpl.rcParams['font.size'] = 18  # Set the font size

# Set math font to match "Latin Modern Math" used in your LaTeX document
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'serif'  # Use serif for regular math text (this matches the font used for text)
mpl.rcParams['mathtext.it'] = 'serif:italic'  # Italic for math symbols
mpl.rcParams['mathtext.bf'] = 'serif:bold'  # Bold math symbols


n = 2500
sys.setrecursionlimit(n)

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


file_path = "noisy_torus_with_curvatures.vtk"
pv_mesh = pv.read(file_path)


data = pv_mesh.cell_centers().points 
faces_vtk = pv_mesh.faces.reshape((-1, 4)) 
vectors = pv_mesh.points[faces_vtk[:, 1:4]]
pv_mesh.compute_normals(cell_normals=True, inplace=True)  # Ensure cell normals are computed
normals_mesh = pv_mesh.cell_normals
normals_mesh /= np.linalg.norm(normals_mesh, axis=1, keepdims=True) 




K_analy = pv_mesh.cell_data['Gaussian_Curvature'] 
H_analy = pv_mesh.cell_data['Mean_Curvature'] 

###### used radius
r = 8
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



########################################################
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


############################# rotate data and fit:
z_axis = np.array([0,0,1])


H = np.empty(N_data)
K = np.empty(N_data)


############################# Calc. Normals with PCA.-> Rotate data -> Fit:
z_axis = np.array([0,0,1])

start = time.time()

def calc_curvatures(ind):
    ind_NN = tree_ind[ind]  # indices of the NN
    dist_NN = dist[ind]  # distances of NN region

    data_NN_single = data[ind_NN]

    data_single_sc = data_NN_single - data[ind]


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


    #### calculating the weights:
    sigma_value = r  # use the Radius as distance parameter
    weights = np.exp(-(dist_NN) ** 2 / (2 * sigma_value ** 2))
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


    sigma = r_smooth/2
    w_i_array = np.exp(-(dist_NN)**2/(2*sigma**2))
    H_smoothed = np.sum(w_i_array*H[ind_NN])/(np.sum(w_i_array))
    K_smoothed = np.sum(w_i_array*K[ind_NN])/(np.sum(w_i_array))


    
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

pv_mesh.cell_data['Gaussian_Curvature_pred'] = K
pv_mesh.cell_data['Mean_Curvature_pred'] = H


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

pv_mesh.cell_data['Mean_Curvature_smoothed'] = H_smoothed
pv_mesh.cell_data['Gaussian_Curvature_smoothed'] = K_smoothed

pv_mesh.cell_data['Normals_smoothed'] = normals_mesh_smoothed

pv_mesh.save("torus_curv_pred.vtk")
# set up a figure twice as wide as it is tall
fig = plt.figure(figsize=(13.5,10), layout = 'constrained',dpi = 300)
cmap = matplotlib.cm.viridis

# Calculate the IQR
q1 = np.percentile(H, 25)
q3 = np.percentile(H, 75)
iqr = q3 - q1

# Define the outlier thresholds
lower_bound = q1 - 1.5 * iqr
upper_bound = q3 + 1.5 * iqr

# Replace outliers with NaN
mask_H = ~((H < lower_bound) | (H > upper_bound))
H = H[mask_H]
H_analy_clean = H_analy[mask_H]

lose = 1

vmin = min(H.min(), H_analy_clean.min())
vmax = max(H.max(), H_analy_clean.max())
norm = Normalize(vmin=vmin, vmax=vmax)
colorbar = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

# =============
# First subplot
# =============
# set up the axes for the first plot
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
ax1.set_title(r'$H_{pred}$')
ax1.set_xticks([-20,0,20])
ax1.set_yticks([-20,0,20])
ax1.set_zticks([-20,0,20])

pred = ax1.scatter(*data[mask_H][::lose].T, c=H[::lose], cmap=cmap, linewidth=0.5, label = 'H_prediction', norm=norm)
set_axes_equal(ax1)
# ==============
# Second subplot
# ==============
# set up the axes for the second plot
ax2 = fig.add_subplot(2, 2, 2, projection='3d')
ax2.set_title(r'$H_{ana}$')
ax2.set_xticks([-20,0,20])
ax2.set_yticks([-20,0,20])
ax2.set_zticks([-20,0,20])
analy = ax2.scatter(*data[mask_H][::lose].T, c=H_analy_clean[::lose], cmap=cmap, linewidth=0.5, label = 'H_analy', norm=norm)
cbar = fig.colorbar(colorbar, ax=[ax1,ax2], label = r'$H$')
cbar.ax.ticklabel_format(scilimits=(-3, 3))


set_axes_equal(ax2)


ax3 = fig.add_subplot(2, 2, 3, projection='3d')
ax3.set_title(r'rel. error: $\left|(H_{ana}-H_{pred})/H_{ana}\right|$')
ax3.set_xticks([-20,0,20])
ax3.set_yticks([-20,0,20])
ax3.set_zticks([-20,0,20])
Diff_rel  = ax3.scatter(*data[mask_H][::lose].T, c=np.abs((H_analy_clean[::lose]-H[::lose])/H_analy_clean[::lose]), cmap='coolwarm', linewidth=0.5, label = 'H_analy-H_pred/(H_analy)')
cbar = fig.colorbar(Diff_rel, ax=ax3, label = 'rel. error')
cbar.ax.ticklabel_format(scilimits=(-3, 3))
set_axes_equal(ax3)


ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title(r'rel. error: $\left(H_{ana}-H_{pred}\right)/H_{ana}$')
ax4.hist((H_analy_clean-H)/H_analy_clean, bins = 50, density=True)
ax4.set_xlabel(r'rel. error')


plt.savefig('Plots/Torus_Curv_H.png', format='png')


# set up a figure twice as wide as it is tall
fig = plt.figure(figsize=(13.5,10), layout = 'constrained', dpi = 300)
cmap = matplotlib.cm.viridis


vmin = min(K.min(), K_analy.min())
vmax = max(K.max(), K_analy.max())
norm = Normalize(vmin=vmin, vmax=vmax)
colorbar = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

# =============
# First subplot
# =============
# set up the axes for the first plot
ax1 = fig.add_subplot(2, 2, 1, projection='3d')
ax1.set_title(r'$K_{pred}$')
ax1.set_xticks([-20,0,20])
ax1.set_yticks([-20,0,20])
ax1.set_zticks([-20,0,20])

pred = ax1.scatter(*data[::lose].T, c=K[::lose], cmap=cmap, linewidth=0.5, label = r'$K_{pred}$', norm=norm)
set_axes_equal(ax1)
# ==============
# Second subplot
# ==============
# set up the axes for the second plot
ax2 = fig.add_subplot(2, 2, 2, projection='3d')
ax2.set_title(r'$K_{ana}$')
ax2.set_xticks([-20,0,20])
ax2.set_yticks([-20,0,20])
ax2.set_zticks([-20,0,20])
analy = ax2.scatter(*data[::lose].T, c=K_analy[::lose], cmap=cmap, linewidth=0.5, label = r'$K_{ana}$', norm=norm)
cbar = fig.colorbar(colorbar, ax=[ax1,ax2],label = r'$K$')
cbar.ax.ticklabel_format(scilimits=(-3, 3))

set_axes_equal(ax2)
K_rel = (K_analy-K)/K_analy
print(np.median(K_rel))
mask = np.abs(K_rel) > 0.6
ax3 = fig.add_subplot(2, 2, 3, projection='3d')
ax3.set_title(r'rel. error: $\left|(K_{ana}-K_{pred})/K_{ana}\right|$')
ax3.set_xticks([-20,0,20])
ax3.set_yticks([-20,0,20])
ax3.set_zticks([-20,0,20])
Diff_rel  = ax3.scatter(*data[~mask][::lose].T, c=np.abs(K_rel[~mask][::lose]), cmap='coolwarm', linewidth=0.5, label = 'K_analy-K_pred/(K_analy)')
cbar = fig.colorbar(Diff_rel, ax=ax3, label = 'rel. error')
cbar.ax.ticklabel_format(scilimits=(-3, 3))
set_axes_equal(ax3)

print(np.median(K_rel[~mask]))
print(np.median((H_analy_clean-H)/H_analy_clean))


ax4 = fig.add_subplot(2, 2, 4)
ax4.set_title(r'rel. error: $\left(K_{ana}-K_{pred}\right)/K_{ana}$')
ax4.hist(K_rel[~mask], bins = 100, density=True)
ax4.set_xlim(-0.4,0.4)
ax4.set_xlabel(r'rel. error')

plt.savefig('Plots/Torus_Curv_K.png', format='png')

H_array = np.column_stack((H_analy_clean, H))
K_array = np.column_stack((K_analy, K))
print(H_array.shape)
print(K_array.shape)
n_bins = 40
fig, axs = plt.subplots(1, 2, tight_layout=True, figsize = (7,5) ,sharey=True)
axs[0].hist(H_array,bins = n_bins, density=False, histtype='step', stacked=False, label = ['ana', 'pred'])
axs[0].set_xlabel('mean curvature')
axs[0].set_ylabel('count of points')
axs[0].legend()
axs[1].hist(K_array,bins = n_bins, density=False, histtype='step', stacked=False, label = ['ana', 'pred'])
axs[1].set_xlabel('Gaussian curvature')
axs[1].set_ylabel('count of points')
axs[1].legend()

plt.savefig('Plots/Torus_Curv_Comparison.pdf', format='pdf')
plt.clf()


Radius = np.array([5,6,7,8,9,10,11,12,13])


fig_H_err, axs_H_err = plt.subplots(3, 3, sharey=True, tight_layout=True,figsize = (13.5,10))
fig_H_err.suptitle(r'$H_{ana}-H_{pred}$ with predicted normals', fontsize=10)
axs_H_err = axs_H_err.ravel()


fig_H_err_comparison, ax_H_err_comparison = plt.subplots(figsize = (13.5,10))


fig_H_calc_Norm_rel_err, axs_H_calc_Norm_rel_err = plt.subplots(3, 3, sharey=True, tight_layout=True,figsize = (13.5,10))
fig_H_calc_Norm_rel_err.suptitle(r'$(H_{ana}-H_{pred})/H_{ana}$ with predicted normals', fontsize=10)
axs_H_calc_Norm_rel_err = axs_H_calc_Norm_rel_err.ravel()


n_bins = 25
nr_hist = 0
for r in Radius:
    print(r)
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

    
    

    H_diff = H_analy-H
    H_rel = (H_analy-H)/H_analy
    

    # Calculate skewness using the moment-based formula
    skewness, loc, scale = skewnorm.fit(H_diff)

    # Generate a range of x values
    x = np.linspace(min(H_diff), max(H_diff), 1000)

    # Calculate the y values for a skewed Gaussian distribution
    y_skewed = skewnorm.pdf(x, skewness, loc=loc, scale=scale)

    # Plot the skewed Gaussian curve
    axs_H_err[nr_hist].plot(x, y_skewed, color='k', label='Skewed Gauß-Curve', ls='--')
    

    x = np.linspace(min(H_diff), max(H_diff), 1000)
    y_skewed = skewnorm.pdf(x, skewness, loc=loc, scale=scale)
    axs_H_err[nr_hist].hist(H_diff, bins='auto', density=True)
    axs_H_err[nr_hist].set_title(f'radius = {r}',fontsize='small')
    axs_H_err[nr_hist].set_xlim(left=-0.0005, right=0.02)
    axs_H_err[nr_hist].plot(x, y_skewed, color='k', label='pred. normals', ls='--')

    axs_H_calc_Norm_rel_err[nr_hist].hist(H_rel, bins='auto', density=True)
    axs_H_calc_Norm_rel_err[nr_hist].set_title(f'radius = {r}',fontsize='small')
    axs_H_calc_Norm_rel_err[nr_hist].set_xlim(left=-0.4, right=0.05)


    nr_hist += 1
    ax_H_err_comparison.plot(x,y_skewed, alpha = 0.6,label=f'R = {r}, std = {"{:.2e}".format(scale)}',linewidth=5.0)
    ax_H_err_comparison.set_xlabel(r'$H_{ana}-H_{pred}$')
    ax_H_err_comparison.set_ylabel(r'number of counts')

ax_H_err_comparison.legend()
fig_H_err_comparison.tight_layout()
ax_H_err_comparison.grid()
ax_H_err_comparison.set_xlim([-0.0005, 0.0205])
fig_H_err.savefig('Plots/Radii/fig_H_calc_Norm_err.pdf', format='pdf')
fig_H_err_comparison.savefig('Plots/Radii/fig_H_err_comparison.pdf',format='pdf')
fig_H_calc_Norm_rel_err.savefig('Plots/Radii/fig_H_calc_Norm_rel_err.pdf',format='pdf')
