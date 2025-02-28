import os
import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
import logging
from time import time
import csv
import yaml
from multiprocessing import Pool, cpu_count

# Read configuration file
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

def init_logging(log_file="interface_analysis.log"):
    """
    Initialize the logging function.
    """
    logging.basicConfig(
        filename=log_file,
        filemode="w",
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO
    )
    logging.info("Logging initialized.")

def init_pyrosetta():
    """
    Initializes PyRosetta, allowing unknown residue types to be ignored.
    """
    logging.info("Initializing PyRosetta...")
    try:
        pyrosetta.init(extra_options="-ex1 -ex2aro -use_input_sc -ignore_zero_occupancy false -ignore_unrecognized_res -packing:linmem_ig 10")
        logging.info("PyRosetta initialized successfully.")
    except Exception as e:
        logging.error(f"Error initializing PyRosetta: {e}")
        raise

def prepare_pose(pdb_file):
    """
    Prepare the PDB's pose object and add OXT atoms for the C side.
    """
    try:
        pose = pyrosetta.pose_from_pdb(pdb_file)
        logging.info(f"OXT atom added to C-terminal of {pdb_file}")
        return pose
    except Exception as e:
        logging.error(f"Error preparing pose for {pdb_file}: {e}")
        raise

def get_chains_from_pdb(pdb_file, exclude_chain):
    """
    Extract all strand names from the PDB file, excluding the specified ligand strand names if they exist.
    """
    chains = set()
    try:
        with open(pdb_file, "r") as f:
            for line in f:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 22:
                    chain = line[21].strip()
                    # Exclude the specified chain
                    if chain and chain != exclude_chain:
                        chains.add(chain)
        logging.info(f"Chains identified in {pdb_file}: {chains}")
    except Exception as e:
        logging.error(f"Error reading chains from {pdb_file}: {e}")
        raise
    return list(chains)




def analyze_interface(pdb_file, ligand_chain, distance_threshold):
    """
    Analyze the interface characteristics of a single PDB file.
    """
    try:
        logging.info(f"Analyzing interface for {pdb_file}...")

        # Prepare Pose
        pose = prepare_pose(pdb_file)
        logging.info(f"Pose prepared successfully for {pdb_file}")

        # Confirmation of chain information in Pose
        pdbinfo = pose.pdb_info()
        pose_chains = set([pdbinfo.chain(i) for i in range(1, pose.total_residue() + 1)])
        if ligand_chain not in pose_chains:
            raise ValueError(f"Ligand chain '{ligand_chain}' not found in {pdb_file}. Pose contains chains: {pose_chains}")

        # Obtaining receptor chain information from raw files
        receptor_chains = get_chains_from_pdb(pdb_file, ligand_chain)
        if not receptor_chains:
            raise ValueError(f"No receptor chains found in {pdb_file}")

        logging.info(f"Receptor chains: {receptor_chains}, Ligand chain: {ligand_chain}")

        # Splicing receptor chains into contiguous strings
        chain_combination = f"{''.join(receptor_chains)}_{ligand_chain}"

        # Initializing InterfaceAnalyzerMover
        interface_analyzer = InterfaceAnalyzerMover(chain_combination, False)
        interface_analyzer.set_pack_input(False)
        interface_analyzer.set_pack_separated(False)
        interface_analyzer.set_compute_packstat(True)
        interface_analyzer.set_compute_interface_sc(True)
        interface_analyzer.set_compute_interface_delta_hbond_unsat(True)
        interface_analyzer.set_compute_separated_sasa(True)

        # Application Analysis
        interface_analyzer.apply(pose)

        # Extracted results
        metrics = {
            "description": os.path.basename(pdb_file),
            "interface_dG": interface_analyzer.get_interface_dG(),
            "dG_cross": interface_analyzer.get_crossterm_interface_energy(),
            "dG_cross_ratio": interface_analyzer.get_crossterm_interface_energy_ratio(),
            "dG_separated": interface_analyzer.get_separated_interface_energy(),
            "dG_separated_ratio": interface_analyzer.get_separated_interface_energy_ratio(),
            "dSASA_int": interface_analyzer.get_interface_delta_sasa(),
            "complexed_SASA": interface_analyzer.get_complexed_sasa(),
            "total_hbond_energy": interface_analyzer.get_total_Hbond_E(),
            "unsat_hbonds": interface_analyzer.get_interface_delta_hbond_unsat(),
            "interface_residues_count": interface_analyzer.get_num_interface_residues(),
            "side1_residues_count": interface_analyzer.get_side1_nres(),
            "side2_residues_count": interface_analyzer.get_side2_nres(),
            "interface_packstat": interface_analyzer.get_interface_packstat(),
            "side1_score": interface_analyzer.get_side1_score(),
            "side2_score": interface_analyzer.get_side2_score()
        }

        logging.info(f"Analysis completed for {pdb_file}: {metrics}")
        return metrics
    except Exception as e:
        logging.error(f"Error analyzing {pdb_file}: {e}")
        raise

def batch_analyze_interfaces(folder, ligand_chain, distance_threshold, output_csv):
    """
    Batch analyze all PDB files in a folder and output the results to CSV.
    """
    pdb_files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(".pdb")]
    total_files = len(pdb_files)
    logging.info(f"Total files to analyze: {total_files}")
    print(f"Total files to analyze: {total_files}")

    # Batch analyze all PDB files in a folder and output the results to CSV.
    with Pool(cpu_count()) as pool:
        results = pool.starmap(analyze_interface, [(pdb_file, ligand_chain, distance_threshold) for pdb_file in pdb_files])

    # Write CSV results
    if results:
        fieldnames = results[0].keys()
        with open(output_csv, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        logging.info(f"Analysis complete. Results saved to {output_csv}")
        print(f"\nAnalysis complete. Results saved to {output_csv}")

if __name__ == "__main__":
    # Get input folder and output folder from configuration file
    input_folder = config["step6_output_folder"]
    ligand_chain = "P"  # ligand chain name
    distance_threshold = 5.5  # Distance thresholds for interface analysis

    # Setting the output directory and CSV file path
    output_folder = os.path.join(config["output_folder"], config["step8_output_folder"])
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Output folder created: {output_folder}")

    # Define CSV file path
    output_csv = os.path.join(output_folder, "interface_analysis_results.csv")

    # Initialization Logging
    init_logging(log_file=os.path.join(output_folder, "interface_analysis.log"))

    try:
        # Initializing PyRosetta
        init_pyrosetta()

        # Batch analyze and save results to CSV
        batch_analyze_interfaces(input_folder, ligand_chain, distance_threshold, output_csv)
    except Exception as e:
        logging.critical(f"Critical error: {e}")
        print(f"Critical error: {e}")
