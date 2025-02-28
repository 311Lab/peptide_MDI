from openbabel import openbabel
import os
import yaml
import logging
from concurrent.futures import ThreadPoolExecutor
import sys

# Configuration logs, separation information and error logs
info_log_file = "process_info.log"
error_log_file = "process_error.log"

# Configuration Information Log
info_handler = logging.FileHandler(info_log_file, mode='w')
info_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
info_handler.setLevel(logging.INFO)

# Configuring the Error Log
error_handler = logging.FileHandler(error_log_file, mode='w')
error_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
error_handler.setLevel(logging.ERROR)

# Setting up the logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(info_handler)
logger.addHandler(error_handler)

# Redirect stderr to error_log_file to avoid writing to Nextflow .err file
sys.stderr = open(error_log_file, "a")

# Initialize Open Babel
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb", "pdbqt")
obConversion.SetOptions("K", openbabel.OBConversion.OUTOPTIONS)  # Ignore kekulize warnings

obBuilder = openbabel.OBBuilder()

# Load Configuration
config_path = "config.yaml"
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

input_folder = config["step3_output_folder"]
output_folder = os.path.join(config["output_folder"], config["step4_output_folder"])
os.makedirs(output_folder, exist_ok=True)

def ensure_3D_coordinates(mol, input_pdb):
    """
    Make sure the molecule has 3D coordinates, if not try to generate it.
    """
    if not mol.Has3D():
        logger.warning(f"{input_pdb} No 3D coordinates, trying to generate...")
        if not obBuilder.Build(mol):
            logger.error(f"{input_pdb} Failed to generate 3D coordinates.")
            return False
    return True

def convert_pdb_to_pdbqt(input_pdb, input_folder, output_folder):
    """
    Convert individual PDB files to PDBQT files.
    """
    mol = openbabel.OBMol()
    pdb_path = os.path.join(input_folder, input_pdb)
    logger.info(f"Trying to read a file: {pdb_path}")

    try:
        # Reading PDB files
        if not obConversion.ReadFile(mol, pdb_path):
            logger.error(f"Failed to read: {pdb_path}, please check the file format")
            return False

        # Adding hydrogen atoms and optimizing molecules
        mol.AddHydrogens()
        if not obBuilder.Build(mol):
            logger.error(f"Failed to build molecule: {input_pdb}")
            return False

        # Ensure 3D coordinates
        if not ensure_3D_coordinates(mol, input_pdb):
            return False

        # Write to PDBQT file
        pdbqt_path = os.path.join(output_folder, f"{os.path.splitext(input_pdb)[0]}.pdbqt")
        if not obConversion.WriteFile(mol, pdbqt_path):
            logger.error(f"Write failed: {pdbqt_path}")
            return False

        logger.info(f"{input_pdb} has been successfully converted to {pdbqt_path}.")
        return True

    except Exception as e:
        logger.error(f"Error processing file {input_pdb}: {e}")
        return False

def process_files(input_files, input_folder, output_folder):
    """
    Parallel processing of file conversions.
    """
    success_count = 0
    failed_files = []

    def process_file(input_pdb):
        nonlocal success_count
        if convert_pdb_to_pdbqt(input_pdb, input_folder, output_folder):
            success_count += 1
        else:
            failed_files.append(input_pdb)

    with ThreadPoolExecutor() as executor:
        executor.map(process_file, input_files)

    return success_count, failed_files

# Getting a list of PDB files
input_files = [f for f in os.listdir(input_folder) if f.endswith(".pdb")]

# carry out a conversion
logger.info(f"Starting processing {len(input_files)} PDB files...")
success_count, failed_files = process_files(input_files, input_folder, output_folder)

# Output result statistics
logger.info(f"Conversion Complete: Successful {success_count} files, Failed {len(failed_files)} files")
if failed_files:
    logger.warning("List of failed files.")
    for file in failed_files:
        logger.warning(f"  File: {file}")
