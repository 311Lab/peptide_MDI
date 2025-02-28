import os
import subprocess
import yaml


with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)


# Configuration file path
config_gpf = "config.gpf"  # AutoGrid Configuration File
grid_log = "grid.log"  # Log files for AutoGrid
config_file = "config.fld"  # AutoDock-GPU Profile

# Number of runs and other parameters
nrun = 20
nev = 2500000
psize = 150

# Input and output paths
input_folder = config["step4_output_folder"]  # Folders where ligands are stored
output_folder = os.path.join(config["output_folder"], config["step5_output_folder"])

dlg_folder = os.path.join(output_folder, "dlg")
xml_folder = os.path.join(output_folder, "xml")
best_pose_folder = os.path.join(output_folder, "docking_pose")

# Create output folder (if it does not exist)
os.makedirs(dlg_folder, exist_ok=True)
os.makedirs(xml_folder, exist_ok=True)
os.makedirs(best_pose_folder, exist_ok=True)

# Step 1: Execute the AutoGrid command
print("Running AutoGrid...")
autogrid_cmd = ["autogrid4", "-p", config_gpf, "-l", grid_log]

try:
    subprocess.run(autogrid_cmd, check=True)
    print("AutoGrid execution completed.")
except subprocess.CalledProcessError as e:
    print(f"Error occurred while running AutoGrid: {e}")
    exit(1)

# Step 2: Execute AutoDock-GPU in Batch
print("Running AutoDock-GPU for each ligand...")
for ligand_file in os.listdir(input_folder):
    if ligand_file.endswith(".pdbqt"):
        # Get the base name of the ligand file (minus the extension)
        basename = os.path.splitext(ligand_file)[0]
        ligand_path = os.path.join(input_folder, ligand_file)

        # Output File Path
        output_log = os.path.join(dlg_folder, f"{basename}.dlg")
        output_xml = os.path.join(xml_folder, f"{basename}.xml")
        output_pdbqt = os.path.join(best_pose_folder, f"{basename}.pdbqt")

        # Build AutoDock-GPU command
        autogpu_cmd = [
            "autodock_gpu_128wi",
            "-lfile", ligand_path,
            "-ffile", config_file,
            "-nrun", str(nrun),
            "-resnam", basename,
            "-gbest", "1",
            "-nev", str(nev),
            "-psize", str(psize),
        ]

        # The print command is used for debugging
        print(f"Running command: {' '.join(autogpu_cmd)}")

        # Execute the AutoDock-GPU command
        try:
            subprocess.run(autogpu_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while processing {ligand_file}: {e}")
            continue

        # Moving .dlg files
        if os.path.exists(f"{basename}.dlg"):
            os.rename(f"{basename}.dlg", output_log)

        # Moving .xml files
        if os.path.exists(f"{basename}.xml"):
            os.rename(f"{basename}.xml", output_xml)

        # Moving the optimal configuration file
        if os.path.exists("best.pdbqt"):
            os.rename("best.pdbqt", output_pdbqt)

print("Batch docking is complete!")
