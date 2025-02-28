import os
import yaml
import json
import sys
from DeepDigest.the_main import DeepDigest

# Get current working directory
current_work_dir = os.getcwd()
print(f"Current working directory: {current_work_dir}")

# Receive JSON directory path
json_dir = sys.argv[2]  # The second command line argument is the JSON directory path
print(f"JSON directory: {json_dir}")

# Loading a YAML configuration file
config_path = sys.argv[1]  # The first command line argument is config.yaml
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

# Loading parameters
data_path = os.path.join(current_work_dir, os.path.basename(config["data_path"]))

# Ensure that the output folder path is out_put in the current working directory.
output_folder = os.path.join(current_work_dir, "out_put")
os.makedirs(output_folder, exist_ok=True)

# Ensure that the output folder path is out_put in the current working directory.
output_file = config["step1_output_file"]
res_path = os.path.join(output_folder, output_file)

# debug output
print(f"Input file path: {data_path}")
print(f"Output folder: {output_folder}")
print(f"Output file path: {res_path}")

# Dynamically loaded JSON file names
json_filename = config.get("protease", None) + ".json"  # Read JSON file name from configuration
json_file_path = os.path.join(json_dir, json_filename)

if not os.path.exists(json_file_path):
    raise FileNotFoundError(f"JSON file not found: {json_file_path}")

print(f"Using JSON file: {json_file_path}")

# Load JSON file content
with open(json_file_path, "r") as json_file:
    json_data = json.load(json_file)

# Output loaded JSON data
print(f"Loaded JSON data: {json_data}")

# Execute DeepDigest
DeepDigest(
    data_path=data_path,
    res_path=res_path,
    regular=config["regular"],
    protease=config["protease"],
    missed_cleavages=config["missed_cleavages"],
    min_len=config["min_len"],
    max_len=config["max_len"]
)

print("Step1 completed successfully.")
