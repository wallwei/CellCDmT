

library(UpSetR)

library(openxlsx)

library(RColorBrewer)

data <- read.csv("breast_upsetR.csv", header = TRUE, sep=',')
head(data)
#调整与美化后的集合图#

zero_indices <- which(data$freq == 0)

#总结#

upset(fromList(data),      
      
      nsets=length(data),#显示数据集的所有数据,nsets = 数值调整可视化数据集数量
      
      nintersects=30,#显示前多少个
      
      sets=c("CellCDmT","CellChat","NATMI","iTALK","CellPhoneDB"), # 指定集合或用keep.order = TRUE保持集合按输入的顺序排序
      
      number.angles = 0, #交互集合柱状图的柱标倾角
      
      point.size=3, #图中点的大小
      
      line.size=1, #图中连接线粗细
      
      mainbar.y.label="Overlapped LRIs", #y轴的标签
      
      main.bar.color = 'black', #y轴柱状图颜色
      
      matrix.color="black", #x轴点的颜色
      
      sets.x.label=expression(LRI[italic("sensitivity")]),   #x轴的标签
      
      sets.bar.color=brewer.pal(5,"Set1"),#x轴柱状图的颜色;Set1中只有9个颜色，Set3中有12个颜色，Paired中有12个颜色
      
      mb.ratio = c(0.7, 0.3), #bar plot和matrix plot图形高度的占比
      
      order.by = "freq", #y轴矩阵排序,如"freq"频率，"degree"程度
      
      decreasing = c(T,F), #以上排序是否降序c(FALSE,TRUE)
      
      text.scale=c(2,2,2,2,2,2), #6个参数intersection size title（y标题大小）,intersection size tick labels（y刻度标签大小）, set size title（set标题大小）, set size tick labels（set刻度标签大小）, set names（set 分类标签大小）, numbers above bars（柱数字大小）的设置
      
      shade.color="red", #图中阴影部分的颜色
      
      
      queries=list(list(query=intersects,params=list("CellCDmT","CellPhoneDB"),color="yellow",active=T),#设置自己想要展示的特定组的交集，通过queries参数进行设置，需要展示几个关注组合的颜色，就展示几个
                   
                   list(query=intersects,params=list("CellCDmT","iTALK"),color="blue",active=T),
                   
                   list(query=intersects,params=list("NATMI","iTALK"),color="green",active=T),
                   
                   list(query=intersects,params=list("CellCDmT","CellChat","NATMI","iTALK","CellPhoneDB"),color="orange",active=T),
                   
                   list(query=intersects,params=list("NATMI","CellCDmT"),color="lightgreen",active=T),
                   
                   list(query=intersects,params=list("CellCDmT","NATMI","iTALK"),color="orange",active=T) )
      
)

#----------------------------------------------------------------------------------------------------------------------------------------------------------#




