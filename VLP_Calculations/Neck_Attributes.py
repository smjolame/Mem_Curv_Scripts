import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from Own_Functions import set_axes_equal
from sklearn.decomposition import PCA
from stl import mesh



neck_dia_array = []
sphere_dia_array = []
neck_dia_ratio_array = []

r = '10'
nice_necks = [[3, 3, 1], [3, 3, 2], [3, 3, 4], [3, 3, 5], [3, 6, 3], [5, 7, 1], [5, 7, 3], [5, 7, 4]]
for Sample, TS, VLP in nice_necks:
    file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'
    file_path_stl = f'Results/Nice_Necks/Neck{Sample}{TS}{VLP}.stl'
    pv_mesh = pv.read(file_path)
    stl_mesh = mesh.Mesh.from_file(file_path_stl)

    data = pv_mesh.cell_centers().points
    data_stl = stl_mesh.centroids
    data_stl = data_stl - np.mean(data_stl, axis=0)
    neck_data_mask = pv_mesh.cell_data['Neck'].astype(bool)
    neck_data =  data[neck_data_mask]
    neck_data = neck_data - np.mean(neck_data, axis=0)

    pca = PCA(n_components=3)

    neck_pca = pca.fit_transform(neck_data)
    eigenvalues = pca.explained_variance_

    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(*data_stl.T)
    #ax.quiver(*np.mean(neck_data, axis=0), *pca.components_[0], normalize = True, length = 8, color = 'r')
    #ax.quiver(*np.mean(neck_data, axis=0), *pca.components_[1], normalize = True, length = 8, color = 'b')
    #ax.quiver(*np.mean(neck_data, axis=0), *pca.components_[2], normalize = True, length = 8, color = 'k')
    #set_axes_equal(ax)
    ##plt.show()
    #plt.clf()



    neck_2d = np.column_stack((neck_pca[:,0],neck_pca[:,1]))


    threshold = 0.3

    horizontal = neck_2d[np.abs(neck_2d[:, 1]) <= threshold]
    x_array = horizontal[:,0]
    ind_x_left = np.where(x_array == np.max(x_array[x_array < 0]))[0][0]
    ind_x_right = np.where(x_array == np.min(x_array[x_array > 0]))[0][0]
    x_left = horizontal[ind_x_left]
    x_right = horizontal[ind_x_right]
    
    vertical = neck_2d[np.abs(neck_2d[:, 0]) <= threshold]
    y_array = vertical[:,1]
    ind_y_left = np.where(y_array == np.max(y_array[y_array < 0]))[0][0]
    ind_y_right = np.where(y_array == np.min(y_array[y_array > 0]))[0][0]
    y_left = vertical[ind_y_left]
    y_right = vertical[ind_y_right]

    x_diff = np.linalg.norm(x_left-x_right)
    y_diff = np.linalg.norm(y_left-y_right)
#
    #    
    #plt.figure(figsize=(8, 6))
    #plt.scatter(*neck_2d.T, c='blue', edgecolor='k', alpha=0.7)
    #plt.title(f'Neck {Sample}{TS}{VLP}')
    #plt.xlabel(r'x $ \:/\: \mathrm{nm}$')
    #plt.ylabel(r'y $ \:/\: \mathrm{nm}$')
    #plt.xlim(-25,25)
    #plt.ylim(-25,25)
    #plt.scatter(*x_left, color = 'red')
    #plt.scatter(*x_right, color = 'red')
    #plt.scatter(*y_left, color = 'deepskyblue')
    #plt.scatter(*y_right, color = 'deepskyblue')
    #plt.plot([x_left[0], x_right[0]], [x_left[1], x_right[1]], color='red', label=f'Diameter = {x_diff:.2f} $\mathrm{{nm}}$', lw = 4)
    #plt.plot([y_left[0], y_right[0]], [y_left[1], y_right[1]], color='deepskyblue', label=f'Diameter = {y_diff:.2f} $\mathrm{{nm}}$', lw = 4)
    #plt.grid(True)
    #plt.legend()
    ##plt.show()
    #plt.savefig(f'Plots/Neck_{Sample}_{TS}_{VLP}.pdf')
    #plt.clf()

    pv_mesh.field_data['Dia_1'] = x_diff
    pv_mesh.field_data['Dia_2'] = y_diff
    #pv_mesh.save(f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk')

    diameter = pv_mesh.field_data['Radius_Sphere_fit'][0] * 2

    if x_diff < y_diff:
        a = y_diff/2
        b = x_diff/2
    else:
        a = x_diff/2
        b = y_diff/2
    ecc = np.sqrt(1-b**2/a**2)
    neck_dia_array.append((x_diff+y_diff)/2)
    sphere_dia_array.append(diameter)
    neck_dia_ratio_array.append(ecc)

print(neck_dia_ratio_array)
print(neck_dia_array)

fig, ax = plt.subplots(figsize=(4, 6))  # Adjust width and height
plt.boxplot(neck_dia_ratio_array, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'),showfliers=False)


# Add jittered points
x_positions = np.ones_like(neck_dia_ratio_array)  # All points at x=1
x_jitter = np.random.normal(0, 0.05, size=len(neck_dia_ratio_array))  # Small random jitter
ax.scatter(x_positions[:-3] + x_jitter[:-3], neck_dia_ratio_array[:-3], color="black", alpha=1, s=10, zorder=10, label = r'E$\Delta$H2')
ax.scatter(x_positions[-3:] + x_jitter[-3:], neck_dia_ratio_array[-3:], color="grey", alpha=1, s=10, zorder=10, label = r'E wildtype', marker='x')
plt.legend()
# Remove x-label and add some styling
ax.set_xticks([])  # Remove x-axis ticks
ax.set_ylabel(r"Eccentricity", fontsize=12)  # Add y-axis label
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.grid()

plt.savefig('Plots/Boxplot_Eccentricity.pdf', format='pdf', bbox_inches='tight')


fig, ax = plt.subplots(figsize=(4, 6))  # Adjust width and height
plt.boxplot(neck_dia_array, widths=0.3, patch_artist=True,
                  boxprops=dict(facecolor='deepskyblue', color='deepskyblue'),
                  capprops=dict(color='deepskyblue'),
                  whiskerprops=dict(color='deepskyblue'),
                  flierprops=dict(markerfacecolor='deepskyblue', markeredgecolor='deepskyblue'),
                  medianprops=dict(color='black'),showfliers=False)


# Add jittered points
x_positions = np.ones_like(neck_dia_array)  # All points at x=1
x_jitter = np.random.normal(0, 0.05, size=len(neck_dia_array))  # Small random jitter
ax.scatter(x_positions[:-3] + x_jitter[:-3], neck_dia_array[:-3], color="black", alpha=1, s=10, zorder=10, label = r'E$\Delta$H2')
ax.scatter(x_positions[-3:] + x_jitter[-3:], neck_dia_array[-3:], color="grey", alpha=1, s=10, zorder=10, label = r'E wildtype', marker='x')
plt.legend()
# Remove x-label and add some styling
ax.set_xticks([])  # Remove x-axis ticks
ax.set_ylabel(r"Mean Diameter $ \:/\: \mathrm{{nm}}$", fontsize=12)  # Add y-axis label
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.grid()

plt.savefig('Plots/Boxplot_Diameter.pdf', format='pdf', bbox_inches='tight')


#n = 5

#plt.scatter(neck_dia_array[:n],sphere_dia_array[:n], marker = 'x', label = f'Sample 3')
#plt.scatter(neck_dia_array[n:],sphere_dia_array[n:], marker = 'x', label = f'Sample 5')
#plt.xlabel('Diameter of neck in nm')
#plt.ylabel('Diameter of VLP in nm')
#plt.grid()
#plt.legend()    
#plt.savefig('Plots/neck_dia_compare.png')
