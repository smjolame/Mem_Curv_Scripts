import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import os
import pandas as pd




Sample = 3
Sample_path = f'Results/Sample_{Sample}'
## sample loop:
while os.path.exists(Sample_path):
    if Sample == 3: TS = 2
    elif Sample == 4: TS = 3
    elif Sample == 5: TS = 7
    TS_path = f'Results/Sample_{Sample}/TS00{TS}'
    print(f'Sample {Sample}:')
    VLP_Neck_count = 0
    VLP_noNeck_count = 0
    ## TS loop:
    while os.path.exists(TS_path):
        VLP = 1

        
        file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'

        while os.path.exists(file_path):
            VLP +=1
            VLP_Neck_count +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/Neck/VLP_{VLP}.vtk'

        
        file_path = f'Results/Sample_{Sample}/TS00{TS}/noNeck/VLP_{VLP}.vtk'

        while os.path.exists(file_path):
            VLP +=1
            VLP_noNeck_count +=1
            file_path = f'Results/Sample_{Sample}/TS00{TS}/noNeck/VLP_{VLP}.vtk'

        TS +=1
        TS_path = f'Results/Sample_{Sample}/TS00{TS}'
    print(VLP_Neck_count)
    print(VLP_noNeck_count)   


    Sample += 1
    Sample_path = f'Results/Sample_{Sample}'