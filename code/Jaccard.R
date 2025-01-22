# 加载必要的包
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

# 假设你已经读取数据并计算了jaccard_matrix
# 读取CSV文件
data <- read.csv("mela_upsetR.csv", header = TRUE)

# 计算Jaccard指数的函数
jaccard_index <- function(set1, set2) {
  intersect_len <- length(intersect(set1, set2))
  union_len <- length(union(set1, set2))
  return(intersect_len / union_len)
}

# 获取方法名
methods <- colnames(data)

# 初始化一个矩阵来存储Jaccard指数
jaccard_matrix <- matrix(0, nrow = length(methods), ncol = length(methods))
rownames(jaccard_matrix) <- methods
colnames(jaccard_matrix) <- methods

# 对每对方法计算Jaccard指数
for (i in 1:length(methods)) {
  for (j in i:length(methods)) {
    set1 <- data[[i]]
    set2 <- data[[j]]
    jaccard_matrix[i, j] <- jaccard_index(set1, set2)
    jaccard_matrix[j, i] <- jaccard_matrix[i, j]  # 因为Jaccard指数是对称的
  }
}

# 将Jaccard矩阵转换为数据框格式，方便ggplot绘制
jaccard_df <- melt(jaccard_matrix)

# 计算平均Jaccard指数
avg_jaccard <- rowMeans(jaccard_matrix)
avg_jaccard_df <- data.frame(Method = factor(names(avg_jaccard), levels = methods), AvgJaccard = avg_jaccard)

# 创建直立柱状图，与热图宽度一致，顺序一致
barplot_plot <- ggplot(avg_jaccard_df, aes(x = Method, y = AvgJaccard, fill = AvgJaccard)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = c(  "yellow", "red","lightblue", "skyblue", "blue","darkblue")) +  # 五种颜色
  theme_minimal() +
  labs(title = "Average Jaccard Index", x = "", y = "Average Jaccard") +
  theme(
    panel.background = element_blank(),      # 去掉背景
    panel.grid.major = element_blank(),      # 去掉主要网格线
    panel.grid.minor = element_blank(),      # 去掉次要网格线
    axis.text.x = element_blank(), # 去掉X轴标签
    axis.ticks.x = element_blank(), # 去掉X轴刻度线
    axis.title.x = element_blank(), # 去掉X轴标题
    axis.text.y = element_text(size = 14), #调整 Y 轴刻度标签的字体大小。
    axis.title.y = element_blank(),          # 去掉X轴标题
    legend.text = element_text(size = 12),   # 图例字体大小
    plot.title = element_text(size = 14),    # 图表标题字体大小
    plot.margin = margin(t = 40, r = 10, b = -80, l = 100) # 减少柱状图底部与热图顶部的间距
  )

# 创建热图，与柱状图X轴对齐
heatmap_plot <- ggplot(jaccard_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +  # 设置白色边框
  scale_fill_gradientn(colors = c("yellow", "red", "lightblue", "skyblue", "blue", "darkblue")) +  # 选择五种渐变色
  theme_minimal() + 
  labs(title = "Jaccard Similarity Heatmap", x = "", y = "") +
  coord_fixed(ratio = 1) +  # 确保每个单元格为正方形
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # 倾斜X轴标签，与柱状图对齐
    axis.text.y = element_text(size = 18),   # Y轴刻度字体大小
    axis.title.x = element_blank(),          # 去掉X轴标题
    axis.title.y = element_blank(),          # 去掉Y轴标题
    plot.title = element_blank(),            # 去掉热图顶部标题
    legend.text = element_text(size = 12),   # 图例字体大小
    plot.margin = margin(t = -140, r = 10, b = 20, l = 10), # 调整热图的边距，减少与柱状图的间距
    axis.ticks.x = element_blank(),  # 去掉X轴的刻度线
    axis.ticks.y = element_blank(),  # 去掉Y轴的刻度线
    axis.ticks = element_blank(),    # 确保去掉所有的刻度线
    legend.key = element_blank(),    # 去掉颜色条的背景
    panel.border = element_blank(),  # 去掉面板边框，防止产生外部刻度线
    panel.grid.major = element_blank(), # 去掉网格线
    panel.grid.minor = element_blank()  # 去掉次要网格线
  ) +
  guides(fill = guide_colorbar(ticks = FALSE))  # 去掉颜色条的刻度线


# 使用gridExtra包合并柱状图和热图
grid.arrange(barplot_plot, heatmap_plot, ncol = 1, heights = c(1, 15))
