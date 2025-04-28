import sys

Neck_nr = sys.argv[1]


def find_boundary(fe_file):
    vertices = {}
    edges = {}
    faces = []

    # Read the .fe file
    with open(fe_file, 'r') as fe:
        lines = fe.readlines()

        # Parsing vertices
        i = lines.index("vertices\n") + 1
        while lines[i].strip() != "edges":
            line = lines[i].strip()
            if line:  # Skip empty lines
                parts = line.split()
                if len(parts) >= 4:  # Ensure the line has enough data
                    vertex_id = int(parts[0])
                    x, y, z = map(float, parts[1:4])
                    vertices[vertex_id] = {"coords": (x, y, z)}
            i += 1

        # Parsing edges
        i = lines.index("edges\n") + 1
        while lines[i].strip() != "faces":
            line = lines[i].strip()
            if line:  # Skip empty lines
                parts = line.split()
                if len(parts) >= 3:  # Ensure the line has enough data
                    edge_id = int(parts[0])
                    vertex_1 = int(parts[1])
                    vertex_2 = int(parts[2])
                    if edge_id not in edges:
                        edges[edge_id] = {"vertices": (vertex_1, vertex_2), "faces": []}
                    else:
                        edges[edge_id]["vertices"] = (vertex_1, vertex_2)
            i += 1

        # Parsing faces
        i = lines.index("faces\n") + 1
        while i < len(lines):
            line = lines[i].strip()
            if line:  # Skip empty lines
                parts = line.split()
                if len(parts) >= 2:  # Ensure the line has enough data
                    face_id = int(parts[0])
                    edge_list = list(map(int, parts[1:]))
                    faces.append((face_id, edge_list))
                    for edge_id in edge_list:
                        # Handle negative edge IDs by considering their positive counterparts
                        if edge_id < 0:
                            edge_id = -edge_id  # Map negative edge to positive
                        if edge_id in edges:
                            edges[edge_id]["faces"].append(face_id)
                        else:
                            print(f"Warning: Edge {edge_id} not found in edges dictionary.")
            i += 1

    # Identify boundary edges (those with only one adjacent face)
    boundary_edges = [edge_id for edge_id, edge_data in edges.items() if len(edge_data["faces"]) == 1]
    
    # Identify boundary vertices (those that belong to a boundary edge)
    boundary_vertices_ids = set()
    for edge_id in boundary_edges:
        boundary_vertices_ids.update(edges[edge_id]["vertices"])

    # Print the boundary vertices
    print("Boundary vertices:")
    for vertex_id in boundary_vertices_ids:
        print(f"Vertex {vertex_id}: {vertices[vertex_id]['coords']}")

    return boundary_vertices_ids, boundary_edges

# Example usage
#fe_file = f"fe/Neck{Neck_nr}.fe"  # Replace with your .fe file path
fe_file = f"fe/cylinder.fe"  # Replace with your .fe file path

boundary_vertices_ids, boundary_edges_ids  = find_boundary(fe_file)




def append_keyword_to_vertices(fe_file, boundary_vertex_ids, boundary_edges_ids, keyword="fixed"):
    # Read the .fe file
    with open(fe_file, 'r') as fe:
        lines = fe.readlines()

    # Variables to store the updated lines and process the vertices section
    updated_lines = []
    vertices_section_started = False
    edges_section_started = False
    faces_section_started = False

    # Iterate through all lines in the .fe file
    for line in lines:
        stripped_line = line.strip()
        

        if stripped_line == "vertices":
            vertices_section_started = True
            updated_lines.append(line)  # Add the header for the vertices section
            continue
        if stripped_line == "edges":
            edges_section_started = True
            updated_lines.append(line)  # Add the header for the edges section
            continue
        if stripped_line == "faces":
            faces_section_started = True
            updated_lines.append(line)  # Add the header for the edges section
            continue
        
        # Process vertices section
        if vertices_section_started and not edges_section_started:
            parts = stripped_line.split()
            if len(parts) > 0:
                vertex_id = int(parts[0])

                # If the vertex is in the boundary list, append the custom keyword
                if vertex_id in boundary_vertex_ids:
                    parts.append(keyword)  # Append the keyword

                # Rebuild the vertex line, ensuring it's still one line per vertex
                updated_lines.append(" ".join(parts) + "\n")  # Keep each vertex on its own line
            else:
                updated_lines.append(line)  # For empty or malformed lines, just add them as is
        # Process edge section
        if edges_section_started and not faces_section_started:
            parts = stripped_line.split()
            if len(parts) > 0:
                edge_id = int(parts[0])

                # If the edge is in the boundary list, append the custom keyword
                if edge_id in boundary_edges_ids:
                    parts.append(keyword)  # Append the keyword

                # Rebuild the edge line, ensuring it's still one line per vertex
                updated_lines.append(" ".join(parts) + "\n")  # Keep each edge on its own line
            else:
                updated_lines.append(line)  # For empty or malformed lines, just add them as is
        # Just add other sections (edges, faces) without modifications
        elif faces_section_started:
            updated_lines.append(line)

    # Write the updated .fe file
    with open("fe/boundaries_" + fe_file.split('/')[-1], 'w') as updated_fe:
        updated_fe.writelines(updated_lines)



append_keyword_to_vertices(fe_file, boundary_vertices_ids,boundary_edges_ids, keyword="fixed")