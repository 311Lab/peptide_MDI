import os
import subprocess
import logging
import traceback
import csv
import re
import yaml

config_path = "config.yaml"
with open(config_path, "r") as f:
    config = yaml.safe_load(f)


# 日志配置
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
        """初始化CSV文件头"""
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
        """获取除'P'外的所有蛋白质链"""
        try:
            with open(pdb_file_path, 'r') as f:
                lines = f.readlines()

            # 提取除'P'外的所有链
            protein_chains = sorted(set(
                line[21]
                for line in lines
                if line.startswith('ATOM') and
                   line[21].strip() and
                   line[21] != 'P'
            ))

            logging.info(f"文件 {os.path.basename(pdb_file_path)} 的蛋白质链: {protein_chains}")
            return protein_chains
        except Exception as e:
            logging.error(f"获取蛋白质链信息失败 {pdb_file_path}: {e}")
            return []

    def parse_prodigy_output(self, output_file, filename):
        """解析Prodigy完整输出文件"""
        try:
            with open(output_file, 'r') as f:
                content = f.read()

            # 使用正则表达式提取所有关键信息
            results = {
                'Filename': filename,
            }

            # 匹配规则
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
                    # 对于链数和残基数，需要特殊处理
                    if key == 'Total_Chains':
                        results['Total_Chains'] = match.group(1)
                        results['Total_Residues'] = match.group(2)
                    else:
                        results[key] = match.group(1)
                else:
                    results[key] = 'N/A'

            return results

        except Exception as e:
            logging.error(f"解析 {output_file} 失败: {e}")
            return None

    def run_prodigy_analysis(self, pdb_file_path, protein_chain):
        """使用Prodigy分析肽-蛋白相互作用"""
        try:
            filename = os.path.splitext(os.path.basename(pdb_file_path))[0]
            output_file = os.path.join(self.output_folder, f"{filename}_peptide_affinity.out")

            # 固定肽链为'P'，protein_chain为目标蛋白链
            command = [
                'prodigy',
                pdb_file_path,
                '--selection', 'P', protein_chain,
                '--distance-cutoff',  str(config["Prodigy_distance_cutoff"])
            ]

            with open(output_file, 'w') as outfile:
                result = subprocess.run(
                    command,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    text=True
                )

            if result.returncode == 0:
                logging.info(f"成功分析 {filename} (肽链P 与 {protein_chain} 链)")
                parsed_results = self.parse_prodigy_output(output_file, filename)
                if parsed_results:
                    with open(self.summary_csv_path, 'a', newline='') as csvfile:
                        csv_writer = csv.DictWriter(csvfile, fieldnames=parsed_results.keys())
                        csv_writer.writerow(parsed_results)

                return True
            else:
                error_msg = result.stderr.strip()
                if "No contacts found for selection" in error_msg:
                    logging.warning(f"未找到 {filename} 中 P 链和 {protein_chain} 链的接触")
                else:
                    logging.error(f"分析 {filename} 失败: {error_msg}")
                return False

        except Exception as e:
            logging.error(f"分析 {pdb_file_path} 时发生异常: {e}")
            logging.error(traceback.format_exc())
            return False

    def process_pdb_files(self):
        """批量处理PDB文件"""
        total_files = 0
        successful_analyses = 0
        failed_files = []

        for pdb_file in os.listdir(self.input_folder):
            if not pdb_file.endswith('.pdb'):
                continue

            total_files += 1
            full_pdb_path = os.path.join(self.input_folder, pdb_file)

            # 获取蛋白质链
            protein_chains = self.get_protein_chains(full_pdb_path)

            # 尝试分析每个蛋白质链
            file_analyzed = False
            for chain in protein_chains:
                if self.run_prodigy_analysis(full_pdb_path, chain):
                    successful_analyses += 1
                    file_analyzed = True
                    break

            if not file_analyzed:
                failed_files.append(pdb_file)

        # 打印分析总结
        logging.info("\n=== 分析总结 ===")
        logging.info(f"总文件数: {total_files}")
        logging.info(f"成功分析: {successful_analyses}")
        logging.info(f"失败文件: {failed_files}")

def main():
    input_folder = config["step6_output_folder"]
    output_folder = os.path.join(config["output_folder"], config["step9_output_folder"])

    analyzer = PeptideProdigyAnalyzer(input_folder, output_folder)
    analyzer.process_pdb_files()

if __name__ == "__main__":
    main()



