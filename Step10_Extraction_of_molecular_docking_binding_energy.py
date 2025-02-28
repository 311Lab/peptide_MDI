import os
import csv
import xml.etree.ElementTree as ET
import yaml

with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)


# Input and output paths
input_folder =  os.path.join(config["step5_output_folder"],"xml")  # Replace with the path to the folder where the .xml file is stored
output_folder  =  os.path.join(config["output_folder"], config["step10_output_folder"])  # Output CSV file name
os.makedirs(output_folder, exist_ok=True)
output_csv = os.path.join(output_folder, "output.csv")


# Initialize the result storage list
results = []

# Iterate through all .xml files in the folder
for file in os.listdir(input_folder):
    if file.endswith(".xml"):
        file_path = os.path.join(input_folder, file)

        # Extract the sequence (assuming the filename format is peptide_SEQUENCE.xml)
        sequence = os.path.splitext(file)[0].split("_")[-1]

        # Parsing XML files
        tree = ET.parse(file_path)
        root = tree.getroot()

        # Extract binding_energy for the best configuration (run node with rank=“1”)
        best_run = root.find(".//run[@rank='1']")
        if best_run is not None:
            binding_energy = best_run.attrib.get("binding_energy")
            results.append({"sequence": sequence, "binding_energy": binding_energy})
        else:
            results.append({"sequence": sequence, "binding_energy": "N/A"})  # If not found rank="1”

# Write the results to a CSV file
with open(output_csv, mode="w", newline="", encoding="utf-8") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=["sequence", "binding_energy"])
    writer.writeheader()
    writer.writerows(results)

print(f"Extraction completed! Results have been saved to {output_csv}")
