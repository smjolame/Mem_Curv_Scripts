import numpy as np
import pyvista as pv


TS = '002'
r = 10
Budd = '1'

trans_array = [['002',1,3], ['002',2,3], ['02',1,-10], ['02',2,20], ['02',3,3], ['005',1,0],['005',2,-10],['07',1,0]]

file_path_radius = f'TS{TS}/mid_surface_TS{TS}_Budding{Budd}.vtk'
pv_mesh = pv.read(file_path_radius)
data = pv_mesh.cell_centers().points

print(pv_mesh)


z = np.array(data[:,2])
print('Nans_z:',np.count_nonzero(np.isnan(z)))

data = data[~np.isnan(z)]
z = z[~np.isnan(z)]

translation = 2

mask = (z > np.mean(z)+translation-np.std(z)/10) & (z < np.mean(z) + translation + np.std(z)/10)
thickness = pv_mesh.point_data['Thickness']
print('Nans_thick:',np.count_nonzero(np.isnan(thickness)))
mask_nan = ~np.isnan(thickness)
thickness = thickness[mask_nan]
H = pv_mesh.point_data['Mean_Curvature']
print('Nans_mean:',np.count_nonzero(np.isnan(H)))
H = H[mask_nan]
K = pv_mesh.point_data['Gaussian_Curvature']
print('Nans_gauss:',np.count_nonzero(np.isnan(K)))
K = K[mask_nan]
Normals = pv_mesh.point_data['Normals']
print('Nans_norm:',np.count_nonzero(np.isnan(Normals)))
Normals = Normals[mask_nan]

### only the slice:
H = H[mask]
K = K[mask]
Normals = Normals[mask]
thickness = thickness[mask]
data = data[mask]

pv_mesh = pv.PolyData(data)
pv_mesh.point_data['Thickness'] = thickness
pv_mesh.point_data['Mean_Curvature'] = H
pv_mesh.point_data['Gaussian_Curvature'] = K
pv_mesh.point_data['Normals'] = Normals
pv_mesh.save(f'Slices/TS{TS}_Budding{Budd}_Slice.vtk')
'''
plt.hist(thickness, bins=50)
plt.show()
'''