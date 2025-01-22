import pandas as pd
import numpy as np
import seaborn as sns
from time import time
import networkx as nx
import matplotlib.pyplot as plt
import os
import psutil
itime = time()

# Modify the following parameters as needed.
folders_to_create = [
    r"breast\expression_thresholding",
    r"breast\expression_product",
    r"breast\cell_expression",
    r"breast\Three",
    r"breast\Three\TOP",
    r"breast"]
for folder in folders_to_create:
    os.makedirs(folder, exist_ok=True)
cancer = r'Breast_'  # File directory for cancer species
x_ytick = ['Breast cancer cells', 'Immune cells', 'Stromal cells', 'T cells', 'B cells', 'Myeloid cells'] # Melanoma->0...NK->6
cell_type = 5  # Modify the number of cell types here, such as melanoma ->8, colorectal cancer ->9

def min_max_scaling(data):
    return (data - data.min()) / (data.max() - data.min())


thrd = pd.read_csv(cancer + "thrdy.csv", header=None, index_col=None).to_numpy()
thrdy = min_max_scaling(thrd)
result1 = pd.DataFrame(thrdy)

pro = pd.read_csv(cancer + "multy.csv", header=None,index_col=None).to_numpy()
proy = min_max_scaling(pro)
result2 = pd.DataFrame(proy)

spec = pd.read_csv(cancer + "specy.csv", header=None,index_col=None).to_numpy()
specy = min_max_scaling(spec)
result3 = pd.DataFrame(specy)

total = pd.read_csv(cancer + "totaly.csv", header=None,index_col=None).to_numpy()
totaly = min_max_scaling(total)
result4 = pd.DataFrame(totaly)

reg_prod = pd.read_csv(cancer + "reg_prody.csv", header=None,index_col=None).to_numpy()
reg_prody = min_max_scaling(reg_prod)
result5 = pd.DataFrame(reg_prody)

cell = pd.read_csv(cancer + "cell_expy.csv", header=None,index_col=None).to_numpy()
celly = min_max_scaling(cell)
result6 = pd.DataFrame(celly)


# The three-point estimation method 1 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result1, result2)
result_med = np.median([result1, result2, result3], axis=0)
result_min = np.minimum(result1, result2)
result_matrix = np.maximum(result_max, result3) + result_med * 4 + np.minimum(result_min, result3)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score1.csv", index=True)
plt.savefig(cancer + 'score1-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score1.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score1-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 2 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result1, result2)
result_med = np.median([result1, result2, result4], axis=0)
result_min = np.minimum(result1, result2)
result_matrix = np.maximum(result_max, result4) + result_med * 4 + np.minimum(result_min, result4)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score2.csv", index=True)
plt.savefig(cancer + 'score2-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score2.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score2-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 3 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result2, result3)
result_med = np.median([result2, result3, result4], axis=0)
result_min = np.minimum(result2, result3)
result_matrix = np.maximum(result_max, result4) + result_med * 4 + np.minimum(result_min, result4)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score3.csv", index=True)
plt.savefig(cancer + 'score3-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score3.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score3-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 4 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result5, result3)
result_med = np.median([result5, result3, result4], axis=0)
result_min = np.minimum(result5, result3)
result_matrix = np.maximum(result_max, result4) + result_med * 4 + np.minimum(result_min, result4)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score4.csv", index=True)
plt.savefig(cancer + 'score4-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score4.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score4-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 5 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result1, result2)
result_med = np.median([result1, result2, result6], axis=0)
result_min = np.minimum(result1, result2)
result_matrix = np.maximum(result_max, result6) + result_med * 4 + np.minimum(result_min, result6)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score5.csv", index=True)
plt.savefig(cancer + 'score5-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score5.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score5-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 6 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result3, result4)
result_med = np.median([result3, result4, result6], axis=0)
result_min = np.minimum(result3, result4)
result_matrix = np.maximum(result_max, result6) + result_med * 4 + np.minimum(result_min, result6)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score6.csv", index=True)
plt.savefig(cancer + 'score6-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score6.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score6-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()


# The three-point estimation method 7 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result1, result4)
result_med = np.median([result1, result4, result6], axis=0)
result_min = np.minimum(result1, result4)
result_matrix = np.maximum(result_max, result6) + result_med * 4 + np.minimum(result_min, result6)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score7.csv", index=True)
plt.savefig(cancer + 'score7-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score7.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score7-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()


# The three-point estimation method 8 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result1, result3)
result_med = np.median([result1, result3, result6], axis=0)
result_min = np.minimum(result1, result3)
result_matrix = np.maximum(result_max, result6) + result_med * 4 + np.minimum(result_min, result6)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score8.csv", index=True)
plt.savefig(cancer + 'score8-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score8.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score8-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 9 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result3, result5)
result_med = np.median([result3, result5, result6], axis=0)
result_min = np.minimum(result3, result5)
result_matrix = np.maximum(result_max, result6) + result_med * 4 + np.minimum(result_min, result6)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score9.csv", index=True)
plt.savefig(cancer + 'score9-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score9.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score9-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 10 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result4, result5)
result_med = np.median([result4, result5, result6], axis=0)
result_min = np.minimum(result4, result5)
result_matrix = np.maximum(result_max, result6) + result_med * 4 + np.minimum(result_min, result6)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score10.csv", index=True)
plt.savefig(cancer + 'score10-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score10.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score10-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 11 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result1, result3)
result_med = np.median([result1, result3, result5], axis=0)
result_min = np.minimum(result1, result3)
result_matrix = np.maximum(result_max, result5) + result_med * 4 + np.minimum(result_min, result5)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score11.csv", index=True)
plt.savefig(cancer + 'score11-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score11.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score11-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 12 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result2, result3)
result_med = np.median([result2, result3, result6], axis=0)
result_min = np.minimum(result1, result3)
result_matrix = np.maximum(result_max, result6) + result_med * 4 + np.minimum(result_min, result6)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score12.csv", index=True)
plt.savefig(cancer + 'score12-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score12.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score12-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()


# The three-point estimation method 13 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result1, result5)
result_med = np.median([result1, result5, result6], axis=0)
result_min = np.minimum(result1, result5)
result_matrix = np.maximum(result_max, result6) + result_med * 4 + np.minimum(result_min, result6)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score13.csv", index=True)
plt.savefig(cancer + 'score13-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score13.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score13-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 14 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result1, result3)
result_med = np.median([result1, result3, result4], axis=0)
result_min = np.minimum(result1, result3)
result_matrix = np.maximum(result_max, result4) + result_med * 4 + np.minimum(result_min, result4)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score14.csv", index=True)
plt.savefig(cancer + 'score14-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score14.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score14-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 15 ------------------------------------------------------------------------------------------------------------------------------------

result_max = np.maximum(result2, result5)
result_med = np.median([result2, result5, result6], axis=0)
result_min = np.minimum(result2, result5)
result_matrix = np.maximum(result_max, result6) + result_med * 4 + np.minimum(result_min, result6)
result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score15.csv", index=True)
plt.savefig(cancer + 'score15-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score15.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score15-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()

# The three-point estimation method 16 ------------------------------------------------------------------------------------------------------------------------------------

result_max1 = np.maximum(result1, result2)
result_max2 = np.maximum(result_max1, result3)
result_max3 = np.maximum(result_max2, result4)
result_max4 = np.maximum(result_max3, result5)

result_med = np.median([result1, result2,result3, result4, result5, result6], axis=0)

result_min1 = np.minimum(result1, result2)
result_min2 = np.minimum(result_min1, result3)
result_min3 = np.minimum(result_min2, result4)
result_min4 = np.minimum(result_min3, result5)

result_matrix = np.maximum(result_max4, result6) + result_med * 4 + np.minimum(result_min4, result6)

result_matrix /= 6
result_matrix = pd.DataFrame(result_matrix)
# Generate heat map
fig = plt.figure()
sns_plot = sns.heatmap(result_matrix, cmap='Reds',
                       xticklabels=x_ytick,
                       yticklabels=x_ytick, linewidths=0.5  # , linecolor= 'black'
                       )
plt.xticks(rotation=-45, size=12, ha='left')
plt.yticks(rotation=360, size=12)
xticklabels = [label.get_text() for label in sns_plot.get_xticklabels()]
yticklabels = [label.get_text() for label in sns_plot.get_yticklabels()]
df = result_matrix.copy()
df.index = yticklabels
df.columns = xticklabels
df.index.name = 'cell_type'
df.to_csv(cancer + "score16.csv", index=True)
plt.savefig(cancer + 'score16-heat.pdf', dpi=1080,bbox_inches = 'tight')
print("-----Run Completed----")
plt.show()


# 读取数据
df = pd.read_csv(cancer + "score16.csv", index_col=0)
matrix = df.values.tolist()
a = df.shape[0]

# 定义节点颜色
#node_colors = ['hotpink', 'darkorange', 'b', 'y', 'm', 'c', 'gray', 'violet', 'r', 'khaki', 'darkred', 'saddlebrown'][:a]
node_colors = ['hotpink', 'orange', 'deepskyblue', 'mediumseagreen', 'thistle', 'indigo', 'pink', 'red']

# 创建有向图
G = nx.DiGraph()

# 添加节点
for i, name in enumerate(df.index):
    G.add_node(name, node_color=node_colors[i % len(node_colors)])

# 添加边，并设置边的颜色和权重
for i, row in enumerate(matrix):
    for j, weight in enumerate(row):
        if weight != 0:
            G.add_edge(df.index[i], df.columns[j], weight=weight, edge_color=node_colors[i % len(node_colors)])

# 使用 circular_layout 确保节点均匀分布
pos = nx.circular_layout(G)

# 获取边的宽度和颜色
edge_widths = [8 * d['weight'] for u, v, d in G.edges(data=True)]
edge_colors = [d['edge_color'] for u, v, d in G.edges(data=True)]

# 获取节点的颜色
node_colors = [d['node_color'] for n, d in G.nodes(data=True)]

# 设置节点大小，保持适中
node_size = 150  # 统一的节点大小，避免过大

G.remove_edges_from(nx.selfloop_edges(G))

# 创建画布
fig, ax = plt.subplots(figsize=(14, 14))  # 增加画布大小

# 绘制边
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.8, edge_color=edge_colors, arrowsize=20, connectionstyle='arc3,rad=0.1')

# 绘制节点
nx.draw_networkx_nodes(G, pos, node_shape='o', node_color=node_colors, node_size=node_size, alpha=0.8)

# 绘制标签，根据节点位置调整方向
for node, (x, y) in pos.items():
    # 根据节点位置调整偏移量
    if y > 0:  # 上方节点向上偏移
        label_y = y + 0.1
    else:  # 下方节点向下偏移
        label_y = y - 0.1

    if x > 0:  # 右方节点向右偏移
        label_x = x + 0.2
    else:  # 左方节点向左偏移
        label_x = x - 0.1

    plt.text(label_x, label_y, node, fontsize=20, ha='center', va='center')

# 设置画布的比例
ax.set_aspect('equal')

# 关闭坐标轴显示
plt.axis('off')

# 保存图像为 PDF
plt.savefig('breast\\score16-net.pdf', dpi=1080, bbox_inches='tight')

# 显示图形
plt.show()


print(' Time: {}秒, mem: {}MB'.format((time() - itime), psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024))
