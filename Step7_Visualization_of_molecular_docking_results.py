import os
import subprocess
import yaml
from multiprocessing import Pool, cpu_count

# Read configuration file
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Root folder path
root_folder = config["step6_output_folder"]
output_root = os.path.join(config["output_folder"], config["step7_output_folder"])  # 结果存放目录

# Create results folder (if it does not exist)
os.makedirs(output_root, exist_ok=True)

# Log file path
log_file = os.path.join(output_root, "processing.log")


def process_pdb_file(pdb_file):
    """
    Process a single PDB file, run the PLIP command and record the results.
    """
    input_file = os.path.join(root_folder, pdb_file)
    pdb_name = os.path.splitext(pdb_file)[0]
    output_subfolder = os.path.join(output_root, f"{pdb_name}_results")
    os.makedirs(output_subfolder, exist_ok=True)

    cmd = [
        "plip",
        "-f", input_file,
        "--peptides", "P",
        "-t", "-x", "-y", "-p",
        "-o", output_subfolder
    ]

    try:
        subprocess.run(cmd, check=True)
        return f"Success: {pdb_file} processed. Results saved in {output_subfolder}"
    except subprocess.CalledProcessError as e:
        return f"Error: {pdb_file} failed with error: {e}"


def log_results(results):
    """
    Writes the results of the process to a log file.
    """
    with open(log_file, "w") as log:
        for result in results:
            log.write(result + "\n")
            print(result)
    print(f"Batch processing completed. Logs saved in {log_file}")


def main():
    # Finding PDB files in the root folder
    pdb_files = [f for f in os.listdir(root_folder) if f.endswith(".pdb")]

    if not pdb_files:
        print("No PDB files found in the root folder.")
        return

    # Using Multiple Processes for PDB Files
    with Pool(cpu_count()) as pool:
        results = pool.map(process_pdb_file, pdb_files)

    # Recording the results of processing
    log_results(results)


if __name__ == "__main__":
    # Check if the path exists
    if not os.path.exists(root_folder):
        raise FileNotFoundError(f"Root folder does not exist: {root_folder}")

    main()
