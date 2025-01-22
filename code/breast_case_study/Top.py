import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Modify the following parameters as needed.
# -----------------------------------------------
cancer = r'Breast_'  # File directory
cell_type = 5   # Modify the number of cell types here, such as melanoma ->6, colorectal cancer ->7
top = 20  # The number of top.
# Provide heat map X index names, from relevant cancer type 0 to 123456., and from 12345... to relevant cancer type 0 (The automatic generation of the Y index) (see example)
X_tick = pd.read_csv("X_ticklabels_top20.csv", header=None, index_col=None)
print(X_tick)

# -----------------------------------------------
cell_cell_num = cell_type * 2 + 1
data2 = []
top_three = pd.DataFrame()
for i in range(0, cell_type + 1):
    for j in range(0, cell_type + 1):
        if (i == 0 and j == 0) or (i == 0 and j != 0) or (i != 0 and j == 0):
            list1 = pd.read_csv(cancer + "total_expr_list" + str(i) + str(j) + ".csv",header=None, index_col=None)
            list2 = pd.read_csv(cancer + "cell_expr_list" + str(i) + str(j) + ".csv",header=None, index_col=None)
            _range = np.max(list1[1]) - np.min(list1[1])
            list11 = (list1[1] - np.min(list1[1])) / _range
            _range = np.max(list2[1]) - np.min(list2[1])
            list22 = (list2[1] - np.min(list2[1])) / _range
            list = (list11 + list22)/2
            a = [list1[0], list]
            train = pd.concat(a, axis=1, ignore_index=True)
            df1 = train.sort_values(by=1, ascending=False, ignore_index=True)
            three = df1[0].head(top)
            top_three = pd.concat([top_three, three], axis=0)
            top_three = top_three.drop_duplicates()
            top_three = top_three.sort_values(by=top_three.columns[0])
            top_three = top_three.reset_index(drop=True)
            df1.to_csv(cancer + "\Three\TOP\\" + str(i) + str(j) + ".csv", index = False,header = False)
for i in range(0, cell_type + 1):
    for j in range(0, cell_type + 1):
        if (i == 0 and j == 0) or (i == 0 and j != 0) or (i != 0 and j == 0):
            list = pd.read_csv(cancer + "\Three\TOP\\" + str(i) + str(j) + ".csv",header=None, index_col=None)
            for w in range(0, top_three.shape[0]):  #
                for x in range(0, list.shape[0]):  #
                    if top_three[0][w] == list[0][x]:
                        b1 = list[1][x]  #
                        b1 = b1.astype(float)
                        data2.append(b1)
                        break
# _range = np.max(data2) - np.min(data2)
# data2 = (data2 - np.min(data2)) / _range
sum_data = pd.DataFrame(data2)
totol = sum_data.values
result = totol.reshape((cell_cell_num, top_three.shape[0]))
result = result.T
result = pd.DataFrame(result)
result = result.astype(float)
ytick = top_three[0].tolist()
xtick = X_tick[0].tolist()
fig = plt.figure()
sns_plot = sns.heatmap(
    result,
    cmap='OrRd',  # 设置渐变色，'coolwarm' 是从浅蓝到红色的颜色方案
    xticklabels = xtick,
    yticklabels = ytick,
    linewidths=1      # 单元格间的线宽
)

plt.xticks(rotation=-45, size=10,ha='left')
plt.yticks(rotation=360, size=9)
plt.savefig(cancer + '\Three\TOP\\Top20.pdf', dpi=1080,bbox_inches = 'tight')   #
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "\Three\TOP\\Top_data_20.csv", index=True)
plt.show()