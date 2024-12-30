import os
import yaml
import Bio.PDB
from PeptideBuilder import Geometry
import PeptideBuilder

# 配置输出目录和文件
config_path = "config.yaml"
try:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
except Exception as e:
    print(f"读取配置文件失败: {e}")
    exit(1)

# 读取输入文件路径
# 读取输入文件路径
input_file =  config["step2_output_file"]


# 配置输出文件夹
output_folder = os.path.join(config["output_folder"], config["step3_output_folder"])

# 如果输出文件夹不存在，则创建它
os.makedirs(output_folder, exist_ok=True)

OUT = Bio.PDB.PDBIO()

# 读取 Step 2 的输出文件 (已经去重的肽序列)
try:
    with open(input_file, "r") as file:
        lines = file.readlines()
except FileNotFoundError:
    print(f"文件未找到: {input_file}")
    exit(1)
except Exception as e:
    print(f"读取文件失败: {e}")
    exit(1)

# 逐个处理肽序列并生成 .pdb 文件
for line in lines:
    peptide_sequence = line.strip()  # 去除换行符和其他空白字符
    if peptide_sequence:
        print(f"Processing peptide sequence: {peptide_sequence}")

        # 使用 PeptideBuilder 生成肽的扩展结构
        try:
            structure = PeptideBuilder.make_extended_structure(peptide_sequence)
        except Exception as e:
            print(f"生成扩展结构时出错：{peptide_sequence}, 错误：{e}")
            continue

        # 将结构写入到指定的输出目录
        OUT.set_structure(structure)

        # 生成输出文件路径
        pdb_filename = f"peptide_{peptide_sequence}.pdb"
        pdb_filepath = os.path.join(output_folder, pdb_filename)

        # 保存 .pdb 文件
        try:
            OUT.save(pdb_filepath)
            print(f"{peptide_sequence} 已成功生成 .pdb 文件，保存路径：{pdb_filepath}")
        except Exception as e:
            print(f"保存 .pdb 文件失败：{pdb_filepath}, 错误：{e}")
            continue

print("Step3 转化成 .pdb 的程序已完成")
