import os
from pymol import cmd
import yaml

# 读取配置文件
config_path = "config.yaml"
with open(config_path, "r") as f:
    config = yaml.safe_load(f)


def extract_first_conformer(ligand_file, output_file):
    """
    提取配体的第一个构象并保存为临时PDB文件。

    :param ligand_file: 原始配体PDBQT文件路径
    :param output_file: 提取后保存的PDB文件路径
    """
    with open(ligand_file, "r") as infile, open(output_file, "w") as outfile:
        write = False
        for line in infile:
            if line.startswith("MODEL"):
                # 开始写第一个模型
                if write:
                    break  # 第一个模型结束后停止
                write = True
            if write:
                outfile.write(line)


def merge_vina_results(receptor_file, ligand_files, output_file):
    """
    合并Vina的分子对接结果，将受体链设为A，配体链设为P，保存为单一PDB文件。

    :param receptor_file: 受体PDB文件路径
    :param ligand_files: 配体PDBQT文件路径列表
    :param output_file: 输出PDB文件路径
    """
    # 加载受体并设置链名为A
    cmd.load(receptor_file, "receptor")

    # 临时文件路径
    temp_files = []

    # 加载每个配体的第一个构象并设置链名为P
    for i, ligand_file in enumerate(ligand_files):
        temp_file = f"temp_ligand_{i + 1}.pdb"
        extract_first_conformer(ligand_file, temp_file)
        temp_files.append(temp_file)

        ligand_name = f"ligand_{i + 1}"
        cmd.load(temp_file, ligand_name)
        cmd.alter(ligand_name, "chain='P'")

    # 合并所有对象
    cmd.create("merged", "receptor or (chain P)")

    # 保存合并结果
    cmd.save(output_file, "merged")

    # 清除所有对象
    cmd.delete("all")

    # 删除临时文件
    for temp_file in temp_files:
        os.remove(temp_file)

    print(f"合并完成，结果保存至: {output_file}")


def batch_merge_vina_results(base_folder, receptor_file):
    """
    批量处理文件夹中的分子对接结果，合并后的结果保存到每个原始文件夹中，
    并将所有结果的副本保存到一个统一的result文件夹中。

    :param base_folder: 存放各个对接结果文件夹的主目录
    :param receptor_file: 受体PDB文件路径
    """
    # 创建统一的结果存放文件夹
    result_folder = os.path.join(config["output_folder"], config["step6_output_folder"])
    os.makedirs(result_folder, exist_ok=True)
    print(f"统一结果文件夹创建/存在于: {result_folder}")

    # 遍历 base_folder 中的子文件夹
    for folder in os.listdir(base_folder):
        folder_path = os.path.join(base_folder, folder)
        if os.path.isdir(folder_path):
            ligand_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".pdbqt")]
            if ligand_files:
                output_file = os.path.join(folder_path, f"{folder}.pdb")  # 保存到原始文件夹中
                merge_vina_results(receptor_file, ligand_files, output_file)

                # 将结果复制到 result_folder 中
                result_copy = os.path.join(result_folder, f"{folder}.pdb")
                with open(output_file, "r") as src, open(result_copy, "w") as dest:
                    dest.write(src.read())

                print(f"结果副本保存至: {result_copy}")


# 使用示例
batch_merge_vina_results(
    base_folder=config["step5_output_folder"],
    receptor_file="./protein.pdbqt"
)
