import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

nice_necks = [[3, 3, 1], [3, 3, 4], [3, 3, 5], [3, 6, 3], [5, 7, 4]]

min_array = np.genfromtxt("Min_H_0.csv")
df_Results = pd.read_csv("Results.csv")


# Define column names that matter
columns_to_match = ["Sample", "TS", "VLP"]

# Convert nice_necks into a DataFrame
nice_necks_df = pd.DataFrame(nice_necks, columns=columns_to_match)

# Filter df_Results to get only matching rows
filtered_df = df_Results.merge(nice_necks_df, on=columns_to_match)

# Select only the "Sphere Diameter" column
sphere_diameters = filtered_df["Sphere Diameter"].to_numpy()
mean_curv = filtered_df["Mean Curvature Neck"].to_numpy()
sphere_radius = sphere_diameters/2

sphere_curv = 1/sphere_radius

print(sphere_diameters)
print(min_array)


Neck = np.arange(0,5)

plt.scatter(Neck,mean_curv, label ='Neck H')
plt.scatter(Neck,-min_array, label ='Neck H_0')
plt.scatter(Neck,sphere_curv, label='Sphere H')
plt.legend()
plt.show()