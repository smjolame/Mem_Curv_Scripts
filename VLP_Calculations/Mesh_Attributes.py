import pyvista as pv
from stl import mesh
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

Sample = '5'
TS = '7'
VLP = '1'


file_path = f'EM_Data/Sample_{Sample}/TS00{TS}/Neck/VLP{VLP}.stl'
your_mesh = mesh.Mesh.from_file(file_path)


import numpy as np
from stl import mesh

def face_area(v1, v2, v3):
    """Calculate the area of a triangle given its three vertices."""
    return 0.5 * np.linalg.norm(np.cross(v2 - v1, v3 - v1))

def calculate_mesh_density(stl_mesh):
    """
    Calculate number of vertices, number of faces, total surface area, and face density.
    
    stl_mesh: Mesh object loaded from an STL file.
    """
    # Get unique vertices and faces from the STL mesh
    vertices = np.unique(stl_mesh.vectors.reshape(-1, 3), axis=0)*1e9
    faces = stl_mesh.vectors*1e9
    
    Nᵥ = len(vertices)  # Number of unique vertices 
    Nᶠ = len(faces)     # Number of triangular faces
    
    # Calculate total surface area
    total_area = 0.0
    for face in faces:
        v1, v2, v3 = face
        total_area += face_area(v1, v2, v3)
    # Calculate average face area
    A_avg = total_area / Nᶠ

    # Calculate face density
    face_density = Nᶠ / total_area
    
    return {
        'Number of Vertices': Nᵥ,
        'Number of Faces': Nᶠ,
        'Total Surface Area (nm²)': total_area,  # Specify units
        'Average Face Area (nm²)': A_avg,
        'Face Density (faces/nm²)': face_density
    }

def estimate_voxel_size(stl_mesh, k=2, bins=50):
    """
    Estimate the voxel size from a point cloud.
    
    Parameters:
        points (ndarray): Nx3 array of point cloud data.
        k (int): Number of nearest neighbors to consider (default: 2 for closest neighbor).
        bins (int): Number of bins for the histogram.
        
    Returns:
        float: Estimated voxel size.
    """
    # Build KDTree for efficient neighbor search
    points = np.unique(stl_mesh.vectors.reshape(-1, 3), axis=0)*1e9
    tree = cKDTree(points)
    
    # Query the k nearest neighbors
    distances, _ = tree.query(points, k=k)
    
    # Extract the nearest neighbor distances (ignore self-distance at k=1)
    nn_distances = distances[:, 1]
    
    # Plot a histogram of the distances
    plt.hist(nn_distances, bins=bins, color='blue', alpha=0.7)
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.title('Nearest Neighbor Distance Histogram')
    plt.show()
    
    # Estimate the voxel size (mode of the histogram)
    voxel_size = np.median(nn_distances)  # Mode could also be considered here
    return voxel_size

# Load the STL file using numpy-stl
stl_mesh = mesh.Mesh.from_file(file_path)

# Calculate mesh density metrics
metrics = calculate_mesh_density(stl_mesh)


voxel_size = estimate_voxel_size(stl_mesh)

print(f"Estimated Voxel Size (nm): {voxel_size}")
# Print the results
for key, value in metrics.items():
    print(f'{key}: {value}')
