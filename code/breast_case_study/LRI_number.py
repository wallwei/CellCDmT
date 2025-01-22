import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches

# Modify the following parameters as needed.
# -------------------------
cancer = r'Breast_'  # File directory for cancer species
x_ytick = ['Breast cancer cells', 'Immune cells', 'Stromal cells', 'T cells', 'B cells', 'Myeloid cells'] # Melanoma->0...NK->6
cell_type = 5  # Modify the number of cell types here, such as melanoma ->6, colorectal cancer ->7
# -------------------------

sum_data = 0

for i in range(0, cell_type + 1):
    for j in range(0, cell_type + 1):
        list1 = pd.read_csv(cancer + "total_expr_list" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()  # 得分
        _range1 = np.max(list1[:, 1]) - np.min(list1[:, 1])
        aaa1 = (list1[:, 1] - np.min(list1[:, 1])) / _range1
        count = np.sum(aaa1 > 0.05)
        with open(cancer + "\\total_LRi_num.csv",mode="a") as f:
            f.write(str(i))
            f.write(str(j))
            f.write('____')
            f.write(str(count))
            f.write('\n')
        f.close()
for i in range(0, cell_type + 1):
    for j in range(0, cell_type + 1):
        list2 = pd.read_csv(cancer + "cell_expr_list" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()  # 得分
        _range2 = np.max(list2[:, 1]) - np.min(list2[:, 1])
        aaa2 = (list2[:, 1] - np.min(list2[:, 1])) / _range2
        count = np.sum(aaa2 > 0)
        with open(cancer + "\\cell_LRi_num.csv",mode="a") as f:
            f.write(str(i))
            f.write(str(j))
            f.write('____')
            f.write(str(count))
            f.write('\n')
        f.close()
for i in range(0, cell_type + 1):
    for j in range(0, cell_type + 1):
        list1 = pd.read_csv(cancer + "total_expr_list" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()
        _range1 = np.max(list1[:, 1]) - np.min(list1[:, 1])
        aaa1 = (list1[:, 1] - np.min(list1[:, 1])) / _range1

        list2 = pd.read_csv(cancer + "cell_expr_list" + str(i) + str(j) + ".csv", header=None, index_col=None).to_numpy()
        _range2 = np.max(list2[:, 1]) - np.min(list2[:, 1])
        aaa2 = (list2[:, 1] - np.min(list2[:, 1])) / _range2

        # 确保两个数组长度相同
        min_len = min(len(aaa1), len(aaa2))
        aaa1 = aaa1[:min_len]
        aaa2 = aaa2[:min_len]

        aaa = (aaa1 + aaa2) / 2
        gene = list1[:, 0]
        # 确保gene也与aaa的长度对齐
        gene = gene[:min_len]
        gene_val = np.column_stack((gene, aaa))
        df = pd.DataFrame(gene_val)
        output_path = cancer + "\Three"
        output_name = '{}{}.csv'
        output_file = os.path.join(output_path, output_name.format(i, j))
        df.to_csv(output_file, index=False, header=False)
        count = np.sum(aaa > 0)
        with open(cancer + "\\Three_LRi_num.csv", mode="a") as f:
            f.write(str(i))
            f.write(str(j))
            f.write('____')
            f.write(str(count))
            f.write('\n')
        f.close()
data = pd.read_csv(cancer + "\\Three_LRi_num.csv", header=None, index_col=None).to_numpy()

data2 = []
print("1", len(data2))

for j in range(len(data)):
    data11 = data[j][0]
    data1 = float(data11[6:])
    a1 = data11[:2]
    if a1[0] == a1[1]:
        data1 = 0.0
    data2.append(data1)

sum_data = pd.DataFrame(data2)

totol = sum_data.values
result = totol.reshape((cell_type + 1, cell_type + 1))
result = pd.DataFrame(result)

fig = plt.figure()
sns_plot = sns.heatmap(
    result,
    cmap='coolwarm',  # 设置渐变色，'coolwarm' 是从浅蓝到红色的颜色方案
    annot=True,       # 在单元格中显示数值
    fmt=".4g",        # 格式化显示的  数值
    xticklabels=x_ytick,  # 设置 x 轴的标签
    yticklabels=x_ytick,  # 设置 y 轴的标签
    linewidths=1      # 单元格间的线宽
)

plt.xticks(rotation=-45, size=15, ha='left')
plt.yticks(rotation=360, size=15)
plt.savefig(cancer + '\\Three_LRi_num.pdf', dpi=1080,bbox_inches = 'tight')   #
print("-----The number of LRIs mediating corresponding CCC----")
plt.show()