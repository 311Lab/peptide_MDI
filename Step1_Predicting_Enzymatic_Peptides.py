import os
import yaml
import json
import sys
from DeepDigest.the_main import DeepDigest

# 获取当前工作目录
current_work_dir = os.getcwd()
print(f"Current working directory: {current_work_dir}")

# 接收 JSON 目录路径
json_dir = sys.argv[2]  # 第二个命令行参数是 JSON 目录路径
print(f"JSON directory: {json_dir}")

# 加载 YAML 配置文件
config_path = sys.argv[1]  # 第一个命令行参数是 config.yaml
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

# 加载参数
data_path = os.path.join(current_work_dir, os.path.basename(config["data_path"]))

# 确保输出文件夹路径为当前工作目录下的 out_put
output_folder = os.path.join(current_work_dir, "out_put")
os.makedirs(output_folder, exist_ok=True)

# 更新结果文件路径
output_file = config["step1_output_file"]
res_path = os.path.join(output_folder, output_file)

# 调试输出
print(f"Input file path: {data_path}")
print(f"Output folder: {output_folder}")
print(f"Output file path: {res_path}")

# 动态加载 JSON 文件名
json_filename = config.get("protease", None) + ".json"  # 从配置中读取 JSON 文件名
json_file_path = os.path.join(json_dir, json_filename)

if not os.path.exists(json_file_path):
    raise FileNotFoundError(f"JSON file not found: {json_file_path}")

print(f"Using JSON file: {json_file_path}")

# 加载 JSON 文件内容
with open(json_file_path, "r") as json_file:
    json_data = json.load(json_file)

# 输出已加载的 JSON 数据
print(f"Loaded JSON data: {json_data}")

# 执行 DeepDigest
DeepDigest(
    data_path=data_path,
    res_path=res_path,
    regular=config["regular"],
    protease=config["protease"],
    missed_cleavages=config["missed_cleavages"],
    min_len=config["min_len"],
    max_len=config["max_len"]
)

print("Step1 completed successfully.")
