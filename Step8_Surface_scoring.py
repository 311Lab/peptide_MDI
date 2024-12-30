import os
import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.pose import add_upper_terminus_type_to_pose_residue
import logging
from time import time
import csv
import yaml

# 读取配置文件
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)


def init_logging(log_file="interface_analysis.log"):
    """
    初始化日志功能。
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
    初始化 PyRosetta，允许忽略未知残基类型。
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
    准备PDB的pose对象，并为C端增加OXT原子。
    """
    try:
        pose = pyrosetta.pose_from_pdb(pdb_file)
        # 为C端增加OXT
        add_upper_terminus_type_to_pose_residue(pose, pose.total_residue())
        logging.info(f"OXT atom added to C-terminal of {pdb_file}")
        return pose
    except Exception as e:
        logging.error(f"Error preparing pose for {pdb_file}: {e}")
        raise


def get_chains_from_pdb(pdb_file, exclude_chain):
    """
    从PDB文件中提取所有链名，排除指定的配体链名（若存在）。
    """
    chains = set()
    try:
        with open(pdb_file, "r") as f:
            for line in f:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 22:
                    chain = line[21].strip()
                    # 排除指定链
                    if chain and chain != exclude_chain:
                        chains.add(chain)
        logging.info(f"Chains identified in {pdb_file}: {chains}")
    except Exception as e:
        logging.error(f"Error reading chains from {pdb_file}: {e}")
        raise
    return list(chains)


def analyze_interface(pdb_file, ligand_chain, distance_threshold):
    """
    分析单个PDB文件的界面特征。
    """
    try:
        logging.info(f"Analyzing interface for {pdb_file}...")

        # 准备Pose
        pose = prepare_pose(pdb_file)
        logging.info(f"Pose prepared successfully for {pdb_file}")

        # 在Pose中确认链信息
        pdbinfo = pose.pdb_info()
        pose_chains = set([pdbinfo.chain(i) for i in range(1, pose.total_residue() + 1)])
        # 确保ligand_chain存在于Pose中
        if ligand_chain not in pose_chains:
            raise ValueError(f"Ligand chain '{ligand_chain}' not found in {pdb_file}. Pose contains chains: {pose_chains}")

        # 从原始文件获取受体链信息（不含配体链）
        receptor_chains = get_chains_from_pdb(pdb_file, ligand_chain)
        if not receptor_chains:
            raise ValueError(f"No receptor chains found in {pdb_file}")

        logging.info(f"Receptor chains: {receptor_chains}, Ligand chain: {ligand_chain}")

        # 将受体链拼接为连续字符串
        chain_combination = f"{''.join(receptor_chains)}_{ligand_chain}"

        # 初始化InterfaceAnalyzerMover
        interface_analyzer = InterfaceAnalyzerMover(chain_combination, False)
        interface_analyzer.set_pack_input(False)
        interface_analyzer.set_pack_separated(False)
        interface_analyzer.set_compute_packstat(True)
        interface_analyzer.set_compute_interface_sc(True)
        interface_analyzer.set_compute_interface_delta_hbond_unsat(True)
        interface_analyzer.set_compute_separated_sasa(True)

        # 应用分析
        interface_analyzer.apply(pose)

        # 提取结果
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
    批量分析文件夹中的所有PDB文件并输出结果到CSV。
    """
    pdb_files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(".pdb")]
    total_files = len(pdb_files)
    logging.info(f"Total files to analyze: {total_files}")
    print(f"Total files to analyze: {total_files}")

    results = []
    start_time = time()

    for i, pdb_file in enumerate(pdb_files):
        logging.info(f"Starting analysis for {pdb_file} ({i + 1}/{total_files})...")
        print(f"Analyzing {pdb_file} ({i + 1}/{total_files})...")
        try:
            result = analyze_interface(pdb_file, ligand_chain, distance_threshold)
            results.append(result)
        except Exception as e:
            results.append({"description": pdb_file, "error": str(e)})

        elapsed_time = time() - start_time
        remaining_files = total_files - (i + 1)
        avg_time_per_file = elapsed_time / (i + 1)
        estimated_remaining_time = avg_time_per_file * remaining_files

        logging.info(f"Elapsed time: {elapsed_time:.2f}s. Estimated remaining time: {estimated_remaining_time:.2f}s.")
        print(f"Elapsed time: {elapsed_time:.2f}s. Estimated remaining time: {estimated_remaining_time:.2f}s.")

    # 写入CSV结果
    if results:
        fieldnames = results[0].keys()
        with open(output_csv, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        logging.info(f"Analysis complete. Results saved to {output_csv}")
        print(f"\nAnalysis complete. Results saved to {output_csv}")


if __name__ == "__main__":
    # 从配置文件中获取输入文件夹和输出文件夹
    input_folder = config["step6_output_folder"]
    ligand_chain = "P"  # 配体链名
    distance_threshold = 5.5  # 界面分析的距离阈值

    # 设置输出目录和CSV文件路径
    output_folder = os.path.join(config["output_folder"], config["step8_output_folder"])
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Output folder created: {output_folder}")

    # 定义CSV文件路径
    output_csv = os.path.join(output_folder, "interface_analysis_results.csv")

    # 初始化日志记录
    init_logging(log_file=os.path.join(output_folder, "interface_analysis.log"))

    try:
        # 初始化 PyRosetta
        init_pyrosetta()

        # 批量分析并保存结果到CSV
        batch_analyze_interfaces(input_folder, ligand_chain, distance_threshold, output_csv)
    except Exception as e:
        logging.critical(f"Critical error: {e}")
        print(f"Critical error: {e}")
