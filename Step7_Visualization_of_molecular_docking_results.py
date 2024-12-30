import os
import subprocess
import yaml

# 读取配置文件
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# 根文件夹路径
root_folder = config["step6_output_folder"]
output_root = os.path.join(config["output_folder"], config["step7_output_folder"])  # 结果存放目录

# 创建结果文件夹（如果不存在）
os.makedirs(output_root, exist_ok=True)

# 日志文件路径
log_file = os.path.join(output_root, "processing.log")


def process_pdb_file(pdb_file):
    """
    处理单个PDB文件，运行PLIP命令并记录结果。
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


def main():
    # 查找根文件夹中的 PDB 文件
    pdb_files = [f for f in os.listdir(root_folder) if f.endswith(".pdb")]

    if not pdb_files:
        print("No PDB files found in the root folder.")
        return

    # 逐个处理 PDB 文件
    with open(log_file, "w") as log:
        for pdb_file in pdb_files:
            result = process_pdb_file(pdb_file)
            log.write(result + "\n")
            print(result)

    print(f"Batch processing completed. Logs saved in {log_file}")


if __name__ == "__main__":
    # 检查路径是否存在
    if not os.path.exists(root_folder):
        raise FileNotFoundError(f"Root folder does not exist: {root_folder}")

    main()
