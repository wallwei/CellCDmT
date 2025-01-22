library(Seurat)

count = read.csv('melanoma1.csv',header=TRUE,row.names="cell",sep = ",",check.names=FALSE)
label = read.csv('melanoma2.csv',header=TRUE,row.names="cell",sep = ",",check.names=FALSE)

all(colnames(count) %in% rownames(label))

# 找出count的行名中不在label的列名中的值
missing_in_label <- colnames(count)[!colnames(count) %in% rownames(label)]
print(missing_in_label)



count = as.matrix(count)

which(is.na(count))
count[is.na(count)] <- 0




label = as.matrix(label)

which(is.na(label))

label[, 1] = label[,'labels']
colnames(label)[1] = 'subclass'
subclass = data.frame(as.factor(label[,'subclass']))
 # subclass = data.frame(as.matrix(label[,'subclass']))
colnames(subclass)[1] = 'subclass'

# 创建Seurat对象
testdata = Seurat::CreateSeuratObject(counts = count,
                                      meta.data = subclass)

testdata = Seurat::SetIdent(testdata,value = testdata@meta.data[["subclass"]])

# 标准化数据
testdata <- NormalizeData(testdata)


saveRDS(testdata,"melanoma.RDS")

