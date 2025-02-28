import os
from pymol import cmd
import yaml

with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)


# Setting the path
receptor_file = "protein.pdbqt"  # Receptor file path
ligand_folder = os.path.join(config["step5_output_folder"],"docking_pose") 
output_folder = os.path.join(config["output_folder"], config["step6_output_folder"])  # output folder

# Create output folder (if it does not exist)
os.makedirs(output_folder, exist_ok=True)

# loaded receptor
cmd.load(receptor_file, "receptor")

# Iterate through each ligand file in the ligand folder
for ligand_file in os.listdir(ligand_folder):
    if ligand_file.endswith(".pdbqt"):
        ligand_path = os.path.join(ligand_folder, ligand_file)
        ligand_name = os.path.splitext(ligand_file)[0]  # Remove the extension
        
        # loaded ligand
        cmd.load(ligand_path, "ligand")
        
        # The chain name of the modified ligand is P
        cmd.alter("ligand", "chain='P'")
        cmd.sort()
        
        # Combining receptors and ligands
        merged_output = os.path.join(output_folder, f"{ligand_name}.pdb")
        cmd.save(merged_output, "receptor or ligand")
        
        # Delete the ligand object and prepare to process the next ligand
        cmd.delete("ligand")
        
        print(f"Merged: {ligand_file} -> {merged_output}")

# Deleting a receptor object
cmd.delete("receptor")

print("All consolidated tasks completed!")
