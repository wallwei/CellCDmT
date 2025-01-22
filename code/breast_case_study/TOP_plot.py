import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取已经计算好的数据文件
cancer = r'Breast_'  # 文件路径前缀
input_file = cancer + r'\Three\TOP\Top_data.csv'  # 数据文件路径
data = pd.read_csv(input_file, index_col=0)  # 读取数据，第一列作为索引

# 绘制热图
fig = plt.figure()

sns_plot = sns.heatmap(
    data,
    cmap='OrRd',  # 设置颜色渐变为浅蓝色到红色
    annot=False,      # 是否在单元格中显示数值，False 为不显示
    linewidths=1,     # 单元格之间的线宽
    xticklabels=True, # 显示 X 轴标签
    yticklabels=True  # 显示 Y 轴标签
)

# 设置 X 和 Y 轴的刻度样式
plt.xticks(rotation=-45, size=14, ha='left')  # X 轴标签倾斜
plt.yticks(rotation=360, size=11)             # Y 轴标签保持水平

# 保存热图到文件
output_file = cancer + r'\Three\TOP\Top_heatmap.pdf'  # 输出文件路径
plt.savefig(output_file, dpi=1080, bbox_inches='tight')  # 保存为高分辨率的 PDF 文件

# 显示热图
plt.show()
