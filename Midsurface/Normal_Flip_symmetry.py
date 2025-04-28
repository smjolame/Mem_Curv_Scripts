import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import KDTree
import pandas as pd
from scipy.stats import linregress
import matplotlib as mpl

# Use LaTeX for rendering math expressions
mpl.rcParams['text.usetex'] = True

# Set font properties
mpl.rcParams['font.family'] = 'serif'  # Serif for the main text
mpl.rcParams['font.serif'] = ['Times New Roman']  # This can be customized to another serif font
mpl.rcParams['font.size'] = 14  # Set the font size

# Set math font to match "Latin Modern Math" used in your LaTeX document
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'serif'  # Use serif for regular math text (this matches the font used for text)
mpl.rcParams['mathtext.it'] = 'serif:italic'  # Italic for math symbols
mpl.rcParams['mathtext.bf'] = 'serif:bold'  # Bold math symbols

def flip_normals_against_mean(points, normals_):
    mean_point = np.mean(points, axis=0)
    shifted_points = points - mean_point
    sc = np.sum(shifted_points*normals_, axis=1)
    normals_[sc<0] *=-1
    return normals_

def calculate_local_density(points, radius):
    tree = KDTree(points)
    densities = np.array([len(tree.query_ball_point(point, radius)) for point in points])
    return densities

def find_startpoint_density_based(points):
    # Step 1: Calculate local densities
    densities = calculate_local_density(points, radius=15.0)  # Set radius according to your data scale
    
    # Step 2: get the point of smallest density
    low_density_point = points[np.argmin(densities)]
    
    return low_density_point

def sort_points_by_proximity(points, start_point):
    sorted_points = []
    remaining_points = points.tolist()  # Convert to list for easier manipulation
    
    # Start with the first point
    current_point = start_point
    sorted_points.append(current_point)
    
    while remaining_points:
        # Find the nearest neighbor to the current point
        distances = np.linalg.norm(np.array(remaining_points) - np.array(current_point), axis=1)
        nearest_idx = np.argmin(distances)
        
        # Move the nearest point to the sorted list
        current_point = remaining_points.pop(nearest_idx)
        sorted_points.append(current_point)
    
    return np.array(sorted_points)

def ensure_normals_consistency(points, normals_):
    consistent_normals = []
    valid_points = []
    
    # Start with the first normal as is
    consistent_normals.append(normals_[0])
    valid_points.append(points[0])
    
    # Go through each point and check the consistency with the previous normal
    for i in range(1, len(points)):
        previous_normal = consistent_normals[-1]
        current_normal = normals_[i]
        current_point = points[i]
        
        # Check the dot product between consecutive normals
        # if the sc is to small, the normal is not good
        #if np.abs(np.dot(previous_normal, current_normal)) < 0.09:
        #    print(f"Skipping normal at index {i} due to near-zero scalar product: {np.dot(previous_normal, current_normal)}")
        #    continue
        if np.dot(previous_normal, current_normal) < 0.15:
            # Flip the current normal if it's pointing in the opposite direction
            consistent_normals.append(-current_normal)
        else:
            consistent_normals.append(current_normal)
        # Only add the corresponding point if the normal was kept
        valid_points.append(current_point)
    
    return np.array(valid_points), np.array(consistent_normals)

def sort_and_ensure_normals(points, normals_):

    # Step 1: Find chain endpoints
    start_point= find_startpoint_density_based(points)
    
    # Step 2: Sort points by proximity
    sorted_points = sort_points_by_proximity(points, start_point)
    
    # Reorder normals to match the sorted points
    sorted_normals = normals_[np.argsort(np.linalg.norm(points - sorted_points[:, None], axis=2))[:,0]]
    
    # Step 3: Ensure normals are consistently oriented
    consistent_points, consistent_normals = ensure_normals_consistency(sorted_points, sorted_normals)
    
    return consistent_points, consistent_normals

def pixelize_pointcloud_with_scalars(data, scalars, normals_, pixel_size=0.5):
    # Define pixel size (grid resolution)

    # Determine the min bounds of the point cloud
    min_bound = np.min(data, axis=0) # results in a vector with x = smallest of the x; y = ...

    # Map points to pixel grid and assign scalar values
    pixel_indices = ((data - min_bound) / pixel_size).astype(int)

    # Find unique pixel indices and their corresponding counts
    unique_pixels, inverse_indices, pixel_counts = np.unique(pixel_indices, axis=0, return_inverse=True, return_counts=True)

    # Sum the scalar values for each pixel
    scalars_sums = np.zeros(len(unique_pixels))
    point_sums = np.zeros((len(unique_pixels), 2))
    np.add.at(scalars_sums, inverse_indices, scalars)
    np.add.at(point_sums, inverse_indices, data)


    # Compute the average scalar for each pixel
    scalars_ = scalars_sums / pixel_counts
    
    # Compute the centroid of each unique pixel
    centroids = point_sums / pixel_counts[:, np.newaxis]
    first_occurrence_indices = np.zeros(len(unique_pixels), dtype=int)
    for i, inv_idx in enumerate(inverse_indices):
        if pixel_counts[inv_idx] == 1 or first_occurrence_indices[inv_idx] == 0:
            first_occurrence_indices[inv_idx] = i

    normals_ = normals_[first_occurrence_indices]

    return centroids, scalars_, normals_

array = [['002','1',True],['002','2',True],['02','1',True],['02','2',True],['02','3',True], ['005','1',False],['005','2',False], ['07','1',True]]
#array = [['002','1',True],['002','2',False],['02','1',True],['02','2',True],['02','3',True], ['005','1',True],['005','2',True], ['07','1',False]]

fit_parameters = []

for TS, Budd, inside in array:
    file_path = f'Slices/TS{TS}_Budding{Budd}_Slice.vtk'
    pv_mesh_radius = pv.read(file_path)
    data = pv_mesh_radius.cell_centers().points
    valid_points_mask = ~np.isnan(data).any(axis=1)
    thickness = pv_mesh_radius.point_data['Thickness'][valid_points_mask]
    normals = pv_mesh_radius.point_data['Normals'][valid_points_mask]
    H = pv_mesh_radius.point_data['Mean_Curvature'][valid_points_mask]
    data = data[valid_points_mask]

    ### project on yx plane:

    z = data[:,2]
    mask = (z > np.mean(z)-np.std(z)/2) & (z < np.mean(z) + np.std(z)/2)
    ### only the slice:
    thickness = thickness[mask]
    normals = normals[mask]
    H = H[mask]
    data = data[mask]
    x = data[:,0]
    y = data[:,1]

    x_n = normals[:,0]
    y_n = normals[:,1]
    normals = np.column_stack((x_n,y_n))
    normals = normals / np.linalg.norm(normals, axis=1, keepdims=True)

    #plt.hist(thickness, bins=50)
    data = np.column_stack((x,y))

    # dbscan to get only the important part
    dbscan = DBSCAN(eps=6, min_samples=10).fit(data)
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

    db_mask = labels == largest_cluster_label
    # Retrieve the points belonging to the largest cluster
    data = data[db_mask]
    thickness = thickness[db_mask]
    normals = normals[db_mask]
    H = H[db_mask]



    
    pixel_size = 3  # Define pixel size (grid resolution)
    pixelized_points, pixelized_thickness , pixelized_normals = pixelize_pointcloud_with_scalars(data,thickness,normals, pixel_size)
    pixelized_normals = flip_normals_against_mean(pixelized_points, pixelized_normals)

    start_point = find_startpoint_density_based(pixelized_points)
    sorted_points = sort_points_by_proximity(pixelized_points, start_point)

    sorted_normals = pixelized_normals[np.argsort(np.linalg.norm(pixelized_points - sorted_points[:, None], axis=2))[:,0]]

    sorted_points , flipped_normals = sort_and_ensure_normals(pixelized_points,pixelized_normals)
    if inside == False:
        flipped_normals *= -1



    ######## template ready: flipping all normals and care for the curvatures

    def adjust_curvature_signs(points_, normals_, points_pattern ,normals_pattern, curvatures):
        kdtree = KDTree(points_pattern)
        _, nearest_indices = kdtree.query(points_)

        for i, nearest_idx in enumerate(nearest_indices):
        # Get the normal vector of the nearest template point
            template_normal = normals_pattern[nearest_idx]
        
            # Check the dot product between the dense normal and the template normal
            if np.dot(normals_[i], template_normal) < 0:
                # Flip the normal vector if they point in opposite directions
                normals_[i] *=-1
                curvatures[i] *=-1
        
        return curvatures, normals_
    
    n = 10

    H_adjusted , normals_adjusted = adjust_curvature_signs(data, normals, sorted_points ,flipped_normals, H)

    # Binning curvature into intervals of 0.01
    bins = np.arange(-0.07, 0.08, 0.01)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    # Create a dataframe to organize data
    df = pd.DataFrame({'curvature': H, 'thickness': thickness})

    # Bin curvature values
    df['curvature_bin'] = pd.cut(df['curvature'], bins=bins, include_lowest=True)

    # Group by bins and prepare data for boxplot
    grouped_data = [df[df['curvature_bin'] == bin]['thickness'] for bin in df['curvature_bin'].cat.categories]

    # Identify non-empty bins and their centers
    non_empty_bins = [bin_data for bin_data in grouped_data if not bin_data.empty]
    non_empty_bin_centers = [bin_centers[i] for i, bin_data in enumerate(grouped_data) if not bin_data.empty]

    # Sort bins by absolute distance to zero and take the first four on each side
    sorted_bins = sorted(zip(non_empty_bin_centers, non_empty_bins), key=lambda x: abs(x[0]))
    closest_bins = sorted_bins[:6]  # Get the closest 4 on each side of zero

    # Extract bin centers and statistics
    closest_bin_centers, closest_grouped_data = zip(*closest_bins)
    means = [bin_data.mean() for bin_data in closest_grouped_data]
    std_devs = [bin_data.std() for bin_data in closest_grouped_data]  # Use std deviation as uncertainty

    # Fit a line through the mean values
    slope, intercept, _, _, _ = linregress(closest_bin_centers, means)

    # Create the fitted line
    x_fit = np.linspace(min(closest_bin_centers), max(closest_bin_centers), 10)
    y_fit = slope * x_fit + intercept

    # Plot the boxplots
    plt.boxplot(closest_grouped_data, positions=closest_bin_centers, widths=0.01, showfliers=False)

    # Plot the fitted line
    plt.plot(x_fit, y_fit, color='green', label=f'linear fit')

    # Customize plot
    plt.title(rf'Fit of TS {TS}, Budd {Budd} for $H \in [-0.03,0.03]$')
    plt.legend()
    plt.grid()



    formatted_bins = [f"{x:.2f}" for x in bins]
    # Boxplot in the bottom subplot
    plt.xticks(bins, formatted_bins,rotation = 45)  # Set ticks on bin edges
    plt.xlabel(r'$H \:/\: \mathrm{nm^{-1}}$')
    plt.ylabel(r"Thickness $ \:/\: \mathrm{nm}$")
   
    plt.xlim([-0.05, 0.05])

    # Adjust layout
    plt.tight_layout()
    plt.savefig(f'Plots/Curv_thickness{TS}_{Budd}_Closest.png',format='png')
    plt.clf()


    fit_parameters.append([slope, intercept])


fit_parameters = np.array(fit_parameters)
for slope, intercept in fit_parameters:
    y_fit = slope * x_fit + intercept
    plt.plot(x_fit, y_fit, linewidth = 2)



plt.xlabel(r'$H \:/\: \mathrm{nm^{-1}}$')
plt.ylabel(r"Thickness $ \:/\: \mathrm{nm}$")
plt.grid()
plt.savefig(f'Plots/Curv_thickness_Linear.png',format='png')
plt.clf()
print(fit_parameters[:,0])

fig, ax = plt.subplots(figsize=(3, 5))  # Adjust width and height
ax.boxplot(fit_parameters[:, 0], vert=True, patch_artist=True, boxprops=dict(facecolor="lightgray"))


# Add jittered points
x_positions = np.ones_like(fit_parameters[:, 0])  # All points at x=1
x_jitter = np.random.normal(0, 0.05, size=len(fit_parameters[:, 0]))  # Small random jitter
ax.scatter(x_positions + x_jitter, fit_parameters[:, 0], color="black", alpha=0.7, s=10, zorder=10)

# Remove x-label and add some styling
ax.set_xticks([])  # Remove x-axis ticks
ax.set_ylabel("Slope", fontsize=12)  # Add y-axis label
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig('Plots/Boxplot_slopes.pdf', format='pdf', bbox_inches='tight')
plt.show()  # Show the plot for preview

