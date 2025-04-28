import pyvista as pv
from collections import defaultdict
import sys

Neck_nr = sys.argv[1]

def stl_to_fe(stl_file, fe_file):
    # Read the STL file
    mesh = pv.read(stl_file)

    vertices = mesh.points
    faces = mesh.faces.reshape(-1, 4)[:, 1:]  # Convert faces to vertex indices

    # Map edges with orientation
    edge_map = {}
    edge_count = 1
    edge_loops = defaultdict(list)

    def add_oriented_edge(v1, v2, face_idx, edge_count):
        # Ensure edge is oriented properly
        edge = (v1, v2)
        rev_edge = (v2, v1)

        if edge not in edge_map and rev_edge not in edge_map:
            edge_map[edge] = edge_count
            edge_count += 1

        # Assign the correct orientation
        oriented_edge = edge if edge in edge_map else rev_edge
        edge_loops[face_idx].append(edge_map[oriented_edge] if oriented_edge == edge else -edge_map[rev_edge])
        return edge_count
    for face_idx, face in enumerate(faces, start=1):
        edge_count = add_oriented_edge(face[0] + 1, face[1] + 1, face_idx, edge_count)
        edge_count = add_oriented_edge(face[1] + 1, face[2] + 1, face_idx, edge_count)
        edge_count = add_oriented_edge(face[2] + 1, face[0] + 1, face_idx, edge_count)



    # Write to the .fe file
    with open(fe_file, 'w') as fe:
        # Write vertices
        fe.write("vertices\n")
        for i, v in enumerate(vertices, start=1):
            fe.write(f"{i} {v[0]} {v[1]} {v[2]}\n")

        # Write edges
        fe.write("\nedges\n")
        for edge, idx in edge_map.items():
            fe.write(f"{idx} {edge[0]} {edge[1]}\n")

        # Write faces
        fe.write("\nfaces\n")
        for face_idx, edges in edge_loops.items():
            edge_str = " ".join(map(str, edges))
            fe.write(f"{face_idx} {edge_str}\n")

    print(f"Converted {stl_file} to {fe_file}")


# Example usage
#stl_to_fe(f"Necks/Neck{Neck_nr}.obj", f"fe/Neck{Neck_nr}.fe")
#stl_to_fe(f"elliptical_cylinder.stl", f"fe/elliptical_cylinder.fe")
stl_to_fe(f"cylinder.stl", f"fe/cylinder.fe")