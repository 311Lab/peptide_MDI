import os
import yaml
import Bio.PDB
from PeptideBuilder import Geometry
import PeptideBuilder

# Configuring output directories and files
config_path = "config.yaml"
try:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
except Exception as e:
    print(f"Failed to read configuration file: {e}")
    exit(1)

# Read the input file path
input_file =  config["step2_output_file"]


# Configuring the output folder
output_folder = os.path.join(config["output_folder"], config["step3_output_folder"])

# If the output folder does not exist, create it
os.makedirs(output_folder, exist_ok=True)

OUT = Bio.PDB.PDBIO()

# Read the output file of Step 2 (de-emphasized peptide sequence).
try:
    with open(input_file, "r") as file:
        lines = file.readlines()
except FileNotFoundError:
    print(f"File not found: {input_file}")
    exit(1)
except Exception as e:
    print(f"Failed to read file: {e}")
    exit(1)

# Process peptide sequences one by one and generate .pdb files
for line in lines:
    peptide_sequence = line.strip()  # Remove line breaks and other whitespace characters
    if peptide_sequence:
        print(f"Processing peptide sequence: {peptide_sequence}")

        # Using PeptideBuilder to Generate Extended Structures of Peptides
        try:
            structure = PeptideBuilder.make_extended_structure(peptide_sequence)
        except Exception as e:
            print(f"Error generating extended structure: {peptide_sequence}, error: {e}")
            continue

        # Writes the structure to the specified output directory
        OUT.set_structure(structure)

        # Generate Output File Path
        pdb_filename = f"peptide_{peptide_sequence}.pdb"
        pdb_filepath = os.path.join(output_folder, pdb_filename)

        # Save the .pdb file
        try:
            OUT.save(pdb_filepath)
            print(f"{peptide_sequence} has successfully generated a .pdb file, save path: {pdb_filepath}")
        except Exception as e:
            print(f"Failed to save .pdb file: {pdb_filepath}, error: {e}")
            continue

print("Step3 Conversion to .pdb is complete.")
