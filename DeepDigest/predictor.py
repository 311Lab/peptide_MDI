from keras.models import model_from_json
import numpy as np
import pandas as pd
import time
import re
import os

# 31-mer coding
def coding(mers):
    # amino acid dictionaries
    dic = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8,
           'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16,
           'V': 17, 'W': 18, 'Y': 19, 'Z': 20}
    return [[dic.get(aa, 20) for aa in mer] for mer in mers]  # Use 20 for unknown AA

# predicting
def predictor(protease, data, res_path, model_dir="./DeepDigest"):
    print("Loading model...")
    s3 = time.time()

    # 文件路径
    json_path = os.path.join(model_dir, f"{protease}.json")
    h5_path = os.path.join(model_dir, f"{protease}.h5")

    try:
        # 检查文件是否存在
        if not os.path.exists(json_path):
            raise FileNotFoundError(f"JSON file not found: {json_path}")
        if not os.path.exists(h5_path):
            raise FileNotFoundError(f"H5 file not found: {h5_path}")

        # 加载模型结构
        with open(json_path, 'r') as dd:
            model = dd.read()
        loaded_model = model_from_json(model)

        # 加载模型权重
        loaded_model.load_weights(h5_path)
        e3 = time.time()
        print("Time cost of loading model is %s seconds." % (e3 - s3))
    except FileNotFoundError as e:
        print(e)
        return

    # 预测部分
    print("Predicting! Please wait...")
    s4 = time.time()

    # 提取 31-mers
    mers = []
    for line in data:
        try:
            seqs = re.split('[\t,]', line[1])[1:]
            mers += [mer for mer in seqs if len(mer) > 1]
        except Exception as e:
            print(f"Error processing line: {line}")
            print(e)
            continue
    print("There are %s candidate cleavage sites." % len(mers))

    # predicting probabilities of all the sites
    x_mers = np.array(coding(mers))
    mers_pred = loaded_model.predict(x_mers)

    # calculating the peptide digestibilities
    print("Calculating peptide detectabilities.")
    results = []
    dig_ind = 0  # index of 31-mer
    for line in data:
        pro_id, seqs = line
        seqs_list = seqs.split('\t')
        if len(seqs_list) == 4:
            pep, left_mer, right_mer, missed_mers = seqs_list
            left_dig = 1 - mers_pred[dig_ind] if left_mer != '*' else 1
            dig_ind += (1 if left_mer != '*' else 0)
            right_dig = 1 - mers_pred[dig_ind] if right_mer != '*' else 1
            dig_ind += (1 if right_mer != '*' else 0)
            missed_probs = []
            missed_dig = 1
            for site in missed_mers.split(','):
                if len(site) == 31:
                    prob = mers_pred[dig_ind]
                    dig_ind += 1
                    missed_probs.append(str(1 - prob))
                    missed_dig *= prob
            dig_pro = left_dig * right_dig * missed_dig
            results.append([pro_id, pep, left_dig, right_dig, ','.join(missed_probs)])
        else:
            results.append([pro_id] + [''] * 4)

    e4 = time.time()
    print("Time cost of prediction is %s seconds." % (e4 - s4))

    # save results
    print("Saving results...")
    s5 = time.time()
    results_df = pd.DataFrame(results, columns=[
        "Protein id",
        "Peptide sequence",
        "Digestibility of the N-terminal site",
        "Digestibility of the C-terminal site",
        "Digestibility of the missed site(s)"
    ])
    results_df.to_csv(res_path, sep='\t', index=False)
    e5 = time.time()
    print("Time cost of saving results is %s seconds." % (e5 - s5))
