from openbabel import openbabel
import os
import yaml
import logging
from concurrent.futures import ThreadPoolExecutor
import sys

# 配置日志，分离信息和错误日志
info_log_file = "process_info.log"
error_log_file = "process_error.log"

# 配置信息日志
info_handler = logging.FileHandler(info_log_file, mode='w')
info_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
info_handler.setLevel(logging.INFO)

# 配置错误日志
error_handler = logging.FileHandler(error_log_file, mode='w')
error_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
error_handler.setLevel(logging.ERROR)

# 设置日志记录器
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(info_handler)
logger.addHandler(error_handler)

# 将stderr重定向到error_log_file，避免写入Nextflow .err文件
sys.stderr = open(error_log_file, "a")

# 初始化 Open Babel
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("pdb", "pdbqt")
obConversion.SetOptions("K", openbabel.OBConversion.OUTOPTIONS)  # 忽略 kekulize 警告

obBuilder = openbabel.OBBuilder()

# 加载配置
config_path = "config.yaml"
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

input_folder = config["step3_output_folder"]
output_folder = os.path.join(config["output_folder"], config["step4_output_folder"])
os.makedirs(output_folder, exist_ok=True)

def ensure_3D_coordinates(mol, input_pdb):
    """
    确保分子具有3D坐标，如无则尝试生成。
    """
    if not mol.Has3D():
        logger.warning(f"{input_pdb} 没有3D坐标，尝试生成...")
        if not obBuilder.Build(mol):
            logger.error(f"{input_pdb} 生成3D坐标失败")
            return False
    return True

def convert_pdb_to_pdbqt(input_pdb, input_folder, output_folder):
    """
    将单个PDB文件转换为PDBQT文件。
    """
    mol = openbabel.OBMol()
    pdb_path = os.path.join(input_folder, input_pdb)
    logger.info(f"尝试读取文件: {pdb_path}")

    try:
        # 读取PDB文件
        if not obConversion.ReadFile(mol, pdb_path):
            logger.error(f"读取失败: {pdb_path}，请检查文件格式")
            return False

        # 添加氢原子并优化分子
        mol.AddHydrogens()
        if not obBuilder.Build(mol):
            logger.error(f"构建分子失败: {input_pdb}")
            return False

        # 确保3D坐标
        if not ensure_3D_coordinates(mol, input_pdb):
            return False

        # 写入PDBQT文件
        pdbqt_path = os.path.join(output_folder, f"{os.path.splitext(input_pdb)[0]}.pdbqt")
        if not obConversion.WriteFile(mol, pdbqt_path):
            logger.error(f"写入失败: {pdbqt_path}")
            return False

        logger.info(f"{input_pdb} 已成功转换为 {pdbqt_path}")
        return True

    except Exception as e:
        logger.error(f"处理文件 {input_pdb} 时发生错误: {e}")
        return False

def process_files(input_files, input_folder, output_folder):
    """
    并行处理文件转换。
    """
    success_count = 0
    failed_files = []

    def process_file(input_pdb):
        nonlocal success_count
        if convert_pdb_to_pdbqt(input_pdb, input_folder, output_folder):
            success_count += 1
        else:
            failed_files.append(input_pdb)

    with ThreadPoolExecutor() as executor:
        executor.map(process_file, input_files)

    return success_count, failed_files

# 获取PDB文件列表
input_files = [f for f in os.listdir(input_folder) if f.endswith(".pdb")]

# 执行转换
logger.info(f"开始处理 {len(input_files)} 个PDB文件...")
success_count, failed_files = process_files(input_files, input_folder, output_folder)

# 输出结果统计
logger.info(f"转换完成: 成功 {success_count} 个文件，失败 {len(failed_files)} 个文件")
if failed_files:
    logger.warning("失败文件列表:")
    for file in failed_files:
        logger.warning(f"  文件: {file}")
