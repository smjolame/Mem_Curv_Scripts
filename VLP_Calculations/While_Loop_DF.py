import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os
from sklearn.decomposition import PCA
import pandas as pd
r = '10'
nice_necks = [[3,3,1],[3,3,2],[3,3,4],[3,3,5],[3,6,3],[5,7,1],[5,7,3],[5,7,4]]


rows = []

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
            mask_sphere = pv_mesh.cell_data['Sphere'] == 1
            H = pv_mesh.cell_data['Mean_Curvature_smoothed']
            H_sphere_std = np.std(H[mask_sphere])
            if Neck == True:
                neck_energy = pv_mesh.field_data['Energy_in_Neck_per_square_nm_kBT'][0]
                mask_neck = pv_mesh.cell_data['Neck'] == 1
                H_neck = np.mean(H[mask_neck])
                
            else:
                neck_energy = np.nan
            sphere_diameter = pv_mesh.field_data['Radius_Sphere_fit'][0] * 2
            sphere_energy = pv_mesh.field_data['Energy_in_Sphere_extrapolated_per_square_nm_fit_kBT'][0]

            #### Checking for nice necks:
            if [Sample, TS, VLP] in nice_necks:
                Neck_dia_1 = pv_mesh.field_data['Dia_1'][0]
                Neck_dia_2 = pv_mesh.field_data['Dia_2'][0]
            else:
                Neck_dia_1 = np.nan
                Neck_dia_2 = np.nan


            new_row = {'Sample':Sample, 'TS':TS,'VLP':VLP,'Neck':Neck ,'Neck Diameter 1':Neck_dia_1,'Neck Diameter 2':Neck_dia_2,'Neck Energy':neck_energy,'Sphere Diameter':sphere_diameter,'Sphere Energy':sphere_energy, 'Mean Curvature Neck': H_neck, 'Std of Mean Curvature Sphere': H_sphere_std}


            VLP +=1            
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'

            if not os.path.exists(file_path):
                Neck = False
            if Neck == False:
                file_path = f'Results/Sample_{Sample}/TS00{TS}/noNeck/VLP_{VLP}.vtk'
            rows.append(new_row)


            
        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'


    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'
df = pd.DataFrame(rows)
df.to_csv('CSV/Results.csv', index=False,na_rep='NaN')

