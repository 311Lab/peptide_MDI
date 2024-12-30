import os
import yaml
import subprocess
import glob
from multiprocessing import Pool
from functools import partial

# 读取 config.yaml 配置
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# 检查并创建 step5_output_folder
step5_output_path = os.path.join(config["output_folder"], config["step5_output_folder"])
if not os.path.exists(step5_output_path):
    os.makedirs(step5_output_path)
    print(f"Directory '{step5_output_path}' created.")

# 检查是否存在 pdbqt 文件
pdbqt_files = glob.glob(os.path.join(config["step4_output_folder"], "peptide_*.pdbqt"))
if not pdbqt_files:
    print("No pdbqt files found.")
    exit(1)

# 只选择前10个文件,这行代码到正式提交的时候要删除
pdbqt_files = pdbqt_files[:2]
print(f"Selected the first 10 pdbqt files: {pdbqt_files}")


# 创建存放日志文件的目录
log_dir = "log_file"

step5_log_dir = os.path.join(config["output_folder"], config["step5_output_folder"],log_dir)

if not os.path.exists(step5_log_dir):
    os.makedirs(step5_log_dir)
    print(f"Directory '{step5_log_dir}' created.")
else:
    print(f"Directory '{step5_log_dir}' already exists.")

# 定义执行对接任务的函数
def run_docking(f, config, step5_log_dir):
    base_name = os.path.splitext(os.path.basename(f))[0]
    print(f"Processing ligand {base_name}")
    ligand_dir = os.path.join(config["output_folder"], config["step5_output_folder"], base_name)
    os.makedirs(ligand_dir, exist_ok=True)

    output_file = os.path.join(ligand_dir, "out.pdbqt")
    log_file = os.path.join(step5_log_dir, f"{base_name}_out.log")

    config_file = os.path.abspath("./config.txt")
    command = [
               "./Autodock/Autodock_vina/vina",
               "--config", config_file,
               "--ligand", f,
               "--out", output_file]

    try:
        with open(log_file, "w") as log:
            subprocess.run(command, stdout=log, stderr=log, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred for ligand {base_name}: {e}")

# 使用多进程并行执行对接任务
run_docking_with_config = partial(run_docking, config=config, step5_log_dir=step5_log_dir)
with Pool(30) as pool:  # 设置并行的进程数
    pool.map(run_docking_with_config, pdbqt_files)

print("All ligands processed successfully.")







