import os
import subprocess
import logging
import traceback
import csv
import re
import yaml
from multiprocessing import Pool, cpu_count

# Load Configuration File
config_path = "config.yaml"
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

# Log Configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s: %(message)s',
    filename='prodigy_peptide_analysis.log',
    filemode='w'
)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)

class PeptideProdigyAnalyzer:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.summary_csv_path = os.path.join(output_folder, 'prodigy_peptide_summary.csv')

        os.makedirs(output_folder, exist_ok=True)
        self.init_summary_csv()

    def init_summary_csv(self):
        """Initialize CSV file header"""
        headers = [
            'Filename',
            'Total_Chains',
            'Total_Residues',
            'Total_Contacts',
            'Charged_Charged_Contacts',
            'Charged_Polar_Contacts',
            'Charged_Apolar_Contacts',
            'Polar_Polar_Contacts',
            'Apolar_Polar_Contacts',
            'Apolar_Apolar_Contacts',
            'Apolar_NIS_Percentage',
            'Charged_NIS_Percentage',
            'Binding_Affinity',
            'Dissociation_Constant'
        ]

        with open(self.summary_csv_path, 'w', newline='') as csvfile:
            csv.writer(csvfile).writerow(headers)

    def get_protein_chains(self, pdb_file_path):
        """Get all protein chains except 'P'"""
        try:
            with open(pdb_file_path, 'r') as f:
                lines = f.readlines()

            # Extract all chains except 'P'
            protein_chains = sorted(set(
                line[21]
                for line in lines
                if line.startswith('ATOM') and
                   line[21].strip() and
                   line[21] != 'P'
            ))

            logging.info(f"Protein chain for file {os.path.basename(pdb_file_path)}: {protein_chains}")
            return protein_chains
        except Exception as e:
            logging.error(f"Failed to get protein chain information {pdb_file_path}: {e}")
            return []

    def parse_prodigy_output(self, output_file, filename):
        """Parsing Prodigy's complete output file"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()

            # Use regular expressions to extract all key information
            results = {
                'Filename': filename,
            }

            # Matching rules
            patterns = [
                ('Total_Chains', r'Parsed structure file .+ \((\d+) chains, (\d+) residues\)'),
                ('Total_Contacts', r'No\. of intermolecular contacts: (\d+)'),
                ('Charged_Charged_Contacts', r'No\. of charged-charged contacts: (\d+)'),
                ('Charged_Polar_Contacts', r'No\. of charged-polar contacts: (\d+)'),
                ('Charged_Apolar_Contacts', r'No\. of charged-apolar contacts: (\d+)'),
                ('Polar_Polar_Contacts', r'No\. of polar-polar contacts: (\d+)'),
                ('Apolar_Polar_Contacts', r'No\. of apolar-polar contacts: (\d+)'),
                ('Apolar_Apolar_Contacts', r'No\. of apolar-apolar contacts: (\d+)'),
                ('Apolar_NIS_Percentage', r'Percentage of apolar NIS residues: ([\d.]+)'),
                ('Charged_NIS_Percentage', r'Percentage of charged NIS residues: ([\d.]+)'),
                ('Binding_Affinity', r'Predicted binding affinity \(kcal\.mol-1\):\s*([-\d.]+)'),
                ('Dissociation_Constant', r'Predicted dissociation constant \(M\) at 25\.0°C:\s*([\de.-]+)')
            ]

            for key, pattern in patterns:
                match = re.search(pattern, content)
                if match:
                    # For the number of chains and residues, special treatment is required
                    if key == 'Total_Chains':
                        results['Total_Chains'] = match.group(1)
                        results['Total_Residues'] = match.group(2)
                    else:
                        results[key] = match.group(1)
                else:
                    results[key] = 'N/A'

            return results

        except Exception as e:
            logging.error(f"Parsing {output_file} failed: {e}")
            return None

    def run_prodigy_analysis(self, pdb_file_path, protein_chain, distance_threshold):
        """Analysis of peptide-protein interactions using Prodigy"""
        try:
            filename = os.path.splitext(os.path.basename(pdb_file_path))[0]
            output_file = os.path.join(self.output_folder, f"{filename}_peptide_affinity.out")

            # Fixed peptide chain is 'P' and protein_chain is the target protein chain
            command = [
                'prodigy',
                pdb_file_path,
                '--selection', 'P', protein_chain,
                '--distance-cutoff',  str(distance_threshold)
            ]

            with open(output_file, 'w') as outfile:
                result = subprocess.run(
                    command,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    text=True
                )

            if result.returncode == 0:
                logging.info(f"Successfully analyzed {filename} (peptide chain P with {protein_chain} chain)")
                parsed_results = self.parse_prodigy_output(output_file, filename)
                if parsed_results:
                    with open(self.summary_csv_path, 'a', newline='') as csvfile:
                        csv_writer = csv.DictWriter(csvfile, fieldnames=parsed_results.keys())
                        csv_writer.writerow(parsed_results)

                return True
            else:
                error_msg = result.stderr.strip()
                if "No contacts found for selection" in error_msg:
                    logging.warning(f"No contact found between P-chain and {protein_chain} in {filename}.")
                else:
                    logging.error(f"Analyzing {filename} failed: {error_msg}")
                return False

        except Exception as e:
            logging.error(f"Exception while analyzing {pdb_file_path}: {e}")
            logging.error(traceback.format_exc())
            return False

    def analyze_interface(self, pdb_file, distance_threshold):
        """Analyzing peptide-protein interactions in a single PDB file"""
        try:
            # Getting the protein chain
            protein_chains = self.get_protein_chains(pdb_file)

            # Try to analyze each protein chain
            for chain in protein_chains:
                if self.run_prodigy_analysis(pdb_file, chain, distance_threshold):
                    return True
            return False

        except Exception as e:
            logging.error(f"Exception while analyzing {pdb_file}: {e}")
            logging.error(traceback.format_exc())
            return False

    def process_pdb_files(self):
        """Batch processing of PDB files"""
        total_files = 0
        successful_analyses = 0
        failed_files = []

        pdb_files = [os.path.join(self.input_folder, pdb_file) for pdb_file in os.listdir(self.input_folder) if pdb_file.endswith('.pdb')]

        # Make sure the value of 'config[“Prodigy_distance_cutoff”]' is correct
        distance_threshold = config.get("Prodigy_distance_cutoff", 5.0)  # The default value is 5.0

        # Parallel processing with multiprocessing
        with Pool(cpu_count()) as pool:
            try:
                results = pool.starmap(self.analyze_interface, [(pdb_file, distance_threshold) for pdb_file in pdb_files])
            except Exception as e:
                logging.error(f"Exception occurred during parallel processing: {e}")
                results = []

        # Collection of analytical results
        for result in results:
            if result:
                successful_analyses += 1
            else:
                failed_files.append(result)

        # Print Analytics Summary
        logging.info("\n=== analysis and summary ===")
        logging.info(f"Total number of documents: {total_files}")
        logging.info(f"Success Analysis: {successful_analyses}")
        logging.info(f"Failed documents: {failed_files}")

def main():
    input_folder = config["step6_output_folder"]
    output_folder = os.path.join(config["output_folder"], config["step9_output_folder"])

    analyzer = PeptideProdigyAnalyzer(input_folder, output_folder)
    analyzer.process_pdb_files()

if __name__ == "__main__":
    main()
