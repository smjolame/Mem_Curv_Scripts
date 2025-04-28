import numpy as np
import pyvista as pv

def torus_surface_area(R, r):
    """Calculate the surface area of a torus with major radius R and minor radius r."""
    return 4 * np.pi**2 * R * r

def estimate_required_faces(average_face_area, target_density):
    """Estimate the required number of faces based on the average face area and face density."""
    required_faces = int(np.ceil(target_density / average_face_area))
    return required_faces

def calculate_edge_length(average_face_area):
    """Calculate the edge length based on the average face area of triangles."""
    edge_length = np.sqrt(average_face_area * 4 / np.sqrt(3))
    return edge_length

def create_torus_points(ring_radius, tube_radius, num_points_ring, num_points_tube):
    """Generate points on the surface of a torus."""
    points = []
    
    for i in range(num_points_ring):
        theta = i * 2 * np.pi / num_points_ring
        for j in range(num_points_tube):
            phi = j * 2 * np.pi / num_points_tube
            x = (ring_radius + tube_radius * np.cos(phi)) * np.cos(theta)
            y = (ring_radius + tube_radius * np.cos(phi)) * np.sin(theta)
            z = tube_radius * np.sin(phi)
            points.append([x, y, z])

    return np.array(points)

def create_torus_faces(num_points_ring, num_points_tube):
    """Create triangular faces for the torus mesh."""
    faces = []
    
    for i in range(num_points_ring):
        for j in range(num_points_tube):
            next_i = (i + 1) % num_points_ring
            next_j = (j + 1) % num_points_tube
            
            # Triangle 1
            faces.append([3, 
                i * num_points_tube + j,
                next_i * num_points_tube + j,
                i * num_points_tube + next_j])
            # Triangle 2
            faces.append([3, 
                next_i * num_points_tube + j,
                next_i * num_points_tube + next_j,
                i * num_points_tube + next_j])
    
    return np.array(faces).flatten()

def calculate_mean_curvature(phi, R, r):
    """Analytical mean curvature formula."""
    return -(R + 2 * r * np.cos(phi)) / (2 * r * (R + r * np.cos(phi)))

def calculate_gaussian_curvature(theta, phi, R, r):
    """Gaussian curvature using fundamental forms."""
    # Partial derivatives
    dR_dtheta = np.array([
        -(R + r * np.cos(phi)) * np.sin(theta),
        (R + r * np.cos(phi)) * np.cos(theta),
        np.zeros_like(theta)
    ]).T
    dR_dphi = np.array([
        -r * np.sin(phi) * np.cos(theta),
        -r * np.sin(phi) * np.sin(theta),
        r * np.cos(phi)
    ]).T
    d2R_dtheta2 = np.array([
        -(R + r * np.cos(phi)) * np.cos(theta),
        -(R + r * np.cos(phi)) * np.sin(theta),
        np.zeros_like(theta)
    ]).T
    d2R_dphi2 = np.array([
        -r * np.cos(phi) * np.cos(theta),
        -r * np.cos(phi) * np.sin(theta),
        -r * np.sin(phi)
    ]).T

    # First fundamental form coefficients
    E = np.einsum('ij,ij->i', dR_dtheta, dR_dtheta)
    F = np.einsum('ij,ij->i', dR_dtheta, dR_dphi)
    G = np.einsum('ij,ij->i', dR_dphi, dR_dphi)

    # Normal vector
    normal = np.cross(dR_dtheta, dR_dphi)
    normal /= np.linalg.norm(normal, axis=1)[:, np.newaxis]

    # Second fundamental form coefficients
    L = np.einsum('ij,ij->i', d2R_dtheta2, normal)
    N = np.einsum('ij,ij->i', d2R_dphi2, normal)

    # Gaussian curvature
    K = (L * N - F**2) / (E * G - F**2)
    return K

def calculate_curvatures(torus_mesh, ring_radius, tube_radius):
    """Compute curvatures for the torus mesh."""
    face_centers = torus_mesh.cell_centers().points
    normals = torus_mesh.cell_normals

    # Convert Cartesian coordinates of face centers to toroidal parameters
    theta = np.arctan2(face_centers[:, 1], face_centers[:, 0])
    R = np.sqrt(face_centers[:, 0]**2 + face_centers[:, 1]**2) - ring_radius
    phi = np.arctan2(face_centers[:, 2], R)

    # Calculate mean and Gaussian curvatures
    mean_curvature = calculate_mean_curvature(phi, ring_radius, tube_radius)
    gaussian_curvature = calculate_gaussian_curvature(theta, phi, ring_radius, tube_radius)

    return mean_curvature, gaussian_curvature



# Parameters for the torus
ring_radius = 25  # Major radius
tube_radius = 10   # Minor radius
average_face_area = 0.3  # Average face area
target_density = 2.5  # Target density in faces per unit area

# Calculate the required number of faces
required_faces = estimate_required_faces(average_face_area, target_density)

# Calculate the edge length based on the average face area
edge_length = calculate_edge_length(average_face_area)

# Determine the number of points for the ring and tube
num_points_ring = int(np.ceil(2 * np.pi * ring_radius / edge_length))  # Based on the circumference of the ring
num_points_tube = int(np.ceil(2 * np.pi * tube_radius / edge_length))   # Based on the circumference of the tube

# Generate torus points and faces
points = create_torus_points(ring_radius, tube_radius, num_points_ring, num_points_tube)
faces = create_torus_faces(num_points_ring, num_points_tube)

# Create the PyVista mesh
torus_mesh = pv.PolyData(points, faces)
stl_filename = "torus.stl"
torus_mesh.save(stl_filename)
# Calculate curvatures for the torus
mean_curvature, gaussian_curvature = calculate_curvatures(torus_mesh, ring_radius, tube_radius)

# Add curvatures to the mesh as cell data
torus_mesh.cell_data["Mean_Curvature"] = mean_curvature
torus_mesh.cell_data["Gaussian_Curvature"] = gaussian_curvature

# Save the mesh as .vtk file with curvatures
vtk_filename = "torus_with_curvatures.vtk"
torus_mesh.save(vtk_filename)



# Adding noise to the torus mesh points
noise_level = 0.1  # Adjust the noise intensity as needed

# Add Gaussian noise to the vertex positions
noisy_points = points + np.random.normal(0, noise_level, points.shape)

# Create a new PyVista mesh with noisy points
noisy_torus_mesh = pv.PolyData(noisy_points, faces)

mean_curvature, gaussian_curvature = calculate_curvatures(noisy_torus_mesh, ring_radius, tube_radius)

# Preserve the curvature values from the original mesh
noisy_torus_mesh.cell_data["Mean_Curvature"] = mean_curvature
noisy_torus_mesh.cell_data["Gaussian_Curvature"] = gaussian_curvature
#noisy_torus_mesh.cell_data["Normals"] = torus_mesh.cell_normals

# Save the noisy mesh as a new file
noisy_vtk_filename = "noisy_torus_with_curvatures.vtk"
noisy_torus_mesh.save(noisy_vtk_filename)

print(f"Torus mesh with curvatures saved as {vtk_filename}")




