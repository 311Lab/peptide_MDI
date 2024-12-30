import pandas as pd
import yaml
import os
import sys

# 检查命令行参数
if len(sys.argv) < 3:
    raise ValueError("Usage: python Step2_Peptide_Recognition.py <config.yaml> <step1_file.txt>")

config_path = sys.argv[1]
step1_file = sys.argv[2]

# 加载配置文件
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

# 定义输出路径
output_folder = config["output_folder"]
output_file = config["step2_output_file"]
res_path = os.path.join(output_folder, output_file)

# 确保输出目录存在
os.makedirs(output_folder, exist_ok=True)

# 读取 Step1 的输出文件，指定分隔符为制表符
try:
    df = pd.read_csv(step1_file, sep="\t")
    print(f"成功读取文件: {step1_file}")
except FileNotFoundError:
    raise FileNotFoundError(f"找不到文件: {step1_file}")
except Exception as e:
    raise RuntimeError(f"读取文件时发生错误: {e}")

# 去除 'Peptide sequence' 列中的重复值
if "Peptide sequence" not in df.columns:
    raise KeyError("列 'Peptide sequence' 不存在于输入文件中。")

unique_peptides = df["Peptide sequence"].drop_duplicates()

# 将去重后的序列保存到新的 TXT 文件
unique_peptides.to_csv(res_path, index=False, header=False, sep="\t")
print(f"去重后的序列已保存到 {res_path}")
