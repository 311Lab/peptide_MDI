import pandas as pd
import yaml
import os
import sys

# Checking command line arguments
if len(sys.argv) < 3:
    raise ValueError("Usage: python Step2_Peptide_Recognition.py <config.yaml> <step1_file.txt>")

config_path = sys.argv[1]
step1_file = sys.argv[2]

# Load Configuration File
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

# Defining the output path
output_folder = config["output_folder"]
output_file = config["step2_output_file"]
res_path = os.path.join(output_folder, output_file)

# Make sure the output directory exists
os.makedirs(output_folder, exist_ok=True)

# Read the output file of Step1, specifying the tab as the separator.
try:
    df = pd.read_csv(step1_file, sep="\t")
    print(f"Successfully read file: {step1_file}")
except FileNotFoundError:
    raise FileNotFoundError(f"File not found: {step1_file}")
except Exception as e:
    raise RuntimeError(f"Error reading file: {e}")

# Remove duplicate values in 'Peptide sequence' column
if "Peptide sequence" not in df.columns:
    raise KeyError("The column 'Peptide sequence' does not exist in the input file.")

unique_peptides = df["Peptide sequence"].drop_duplicates()

# Save the de-duplicated sequence to a new TXT file
unique_peptides.to_csv(res_path, index=False, header=False, sep="\t")
print(f"The de-duplicated sequence has been saved to {res_path}")
