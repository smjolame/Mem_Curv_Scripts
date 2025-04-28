import sys

Neck_nr = sys.argv[1]

def generate_code_section_top():
    code = f"""

PARAMETER Neck_nr = {Neck_nr}
quantity mean_curv energy method star_perp_sq_mean_curvature global 


// T := 300
// mesh in nm, energy in k_b T, kappa in k_b T, lambda in nm

"""
    return code

def generate_code_section_bottom():
    code = """



"""
    return code

#file_path = f"fe/boundaries_Neck{Neck_nr}.fe"

file_path = f"fe/boundaries_cylinder.fe"


def read_external_file(file_path):
    with open(file_path, "r") as file:
        return file.read()

def create_fe_file():
    # Open the .fe file in write mode
    #with open(f"code_Necks_fe/code_Neck{Neck_nr}.fe", "w") as fe_file:
    with open(f"code_Necks_fe/code_cylinder.fe", "w") as fe_file:
        
        # Writing some Python generated code
        fe_file.write(generate_code_section_top()) 
        fe_file.write(read_external_file(file_path))  
        fe_file.write(generate_code_section_bottom()) 

    
    print("The .fe file has been created and saved.")



create_fe_file()