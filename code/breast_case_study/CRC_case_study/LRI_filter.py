import os
import warnings
warnings.filterwarnings("ignore")
from itertools import product as product
import pandas as pd
import numpy as np

nan = 0
path = '../out/human/pre+known.csv'
# 读取CSV文件，没有表头和索引列
df = pd.read_csv(path, header=None, index_col=None, sep=',')
print(df.shape)
# 使用 .iloc 来访问行和列
ligand = df.iloc[:, 0]
receptor = df.iloc[:, 1]

# 将Series转换为列表
l_gene = ligand.tolist()
r_gene = receptor.tolist()

# 将列表转换为NumPy数组，并垂直堆叠
LRI_gene = np.vstack((np.array(l_gene), np.array(r_gene))).T

dt = pd.read_csv('colorectal1.csv', index_col=0, header=None)
#a = np.array(dt.loc["Cell"])
dict = {}
dict["cell"] = np.array(dt.loc["cell"])
for i in range(1, dt.shape[0]):
    dict[dt.index[i]] = np.array(dt.loc[dt.index[i]], dtype=float)

savepath = 'HNSCC'

#0=colorectal cancer cells,1=B cells,2=T cells,3=Epithelial cells,4=Fibroblasts,5=Mast cells,6=Macrophages,7=Endothlial
malignant_index = np.where(dict['cell type'] == 0)[0]
B_index = np.where(dict['cell type'] == 1)[0]
T_index = np.where(dict['cell type'] == 2)[0]
Epithelial_index = np.where(dict['cell type'] == 3)[0]
Fibro_index = np.where(dict['cell type'] == 4)[0]
Mast_index = np.where(dict['cell type'] == 5)[0]
Macro_index = np.where(dict['cell type'] == 6)[0]
Endo_index = np.where(dict['cell type'] == 7)[0]

###
groupsize = []
groupsize.append(malignant_index.shape[0])
groupsize.append(Fibro_index.shape[0])
groupsize.append(B_index.shape[0])
groupsize.append(Epithelial_index.shape[0])
groupsize.append(Macro_index.shape[0])
groupsize.append(Endo_index.shape[0])
groupsize.append(T_index .shape[0])
groupsize.append(Mast_index.shape[0])
groupsize = pd.DataFrame(groupsize)
groupsize.to_csv('HNSCC/HNSCC.csv', header=None, index=None)


for i in range(8):
    for j in range(8):
        exec('mult_score{}{} = 0'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('thrd_score{}{} = 0'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('mult_list{}{} = []'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('mult_list_s{}{} = []'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('thrd_list{}{} = []'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('thrd_list_s{}{} = []'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('spec_score{}{} = 0'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('spec_list{}{} = []'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('spec_list_s{}{} = []'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('total_expr_score{}{} = 0'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('total_expr_list{}{} = []'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('total_expr_list_s{}{} = []'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('reg_prod_score{}{} = 0'.format(i, j))
        exec('reg_prod_list{}{} = []'.format(i, j))
        exec('reg_prod_list_s{}{} = []'.format(i, j))

# 细胞表达法的打分初始化
for i in range(8):
    for j in range(8):
        exec('cell_expr_score{}{} = 0'.format(i, j))
        exec('cell_expr_list{}{} = []'.format(i, j))
        exec('cell_expr_list_s{}{} = []'.format(i, j))

    # 定义细胞类型列表
cell_types = ['HNSCC cells', 'B cells', 'T cells',  'Epithelial cells','Fibroblasts', 'Mast cells','Macrophages', 'Endothelial cells']
cell_type_indices = {
    "HNSCC cells": malignant_index,
    "T cells": T_index,
    "B cells": B_index,
    "Epithelial cells": Epithelial_index,
    "Fibroblasts": Fibro_index,
    "Mast cells": Mast_index,
    "Macrophages": Macro_index,
    "Endothelial cells": Endo_index
}


g = 0
mu = np.mean([dict[gene].mean() for gene in dict if gene != "cell"])

for i in LRI_gene:
    if i[0] in dict and i[1] in dict:
        g= g+1
        print(g)
        malignant_l = 1 / malignant_index.shape[0] * sum(dict[i[0]][malignant_index])
        malignant_r = 1 / malignant_index.shape[0] * sum(dict[i[1]][malignant_index])
        Fibro_l = 1 / Fibro_index.shape[0] * sum(dict[i[0]][Fibro_index])
        Fibro_r = 1 / Fibro_index.shape[0] * sum(dict[i[1]][Fibro_index])
        B_l = 1 / B_index.shape[0] * sum(dict[i[0]][B_index])
        B_r = 1 / B_index.shape[0] * sum(dict[i[1]][B_index])
        Epithelial_l = 1 / Epithelial_index.shape[0] * sum(dict[i[0]][Epithelial_index])
        Epithelial_r = 1 / Epithelial_index.shape[0] * sum(dict[i[1]][Epithelial_index])
        Macro_l = 1 / Macro_index.shape[0] * sum(dict[i[0]][Macro_index])
        Macro_r = 1 / Macro_index.shape[0] * sum(dict[i[1]][Macro_index])
        Endo_l = 1 / Endo_index.shape[0] * sum(dict[i[0]][Endo_index])
        Endo_r = 1 / Endo_index.shape[0] * sum(dict[i[1]][Endo_index])
        T_l = 1 / T_index.shape[0] * sum(dict[i[0]][T_index])
        T_r = 1 / T_index.shape[0] * sum(dict[i[1]][T_index])
        Mast_l = 1 / Mast_index.shape[0] * sum(dict[i[0]][Mast_index])
        Mast_r = 1 / Mast_index.shape[0] * sum(dict[i[1]][Mast_index])
        l_list = [malignant_l, Fibro_l, B_l, Epithelial_l, Macro_l, Endo_l, T_l, Mast_l]
        r_list = [malignant_r, Fibro_r, B_r, Epithelial_r, Macro_r, Endo_r, T_r, Mast_r]

        a = b = 0
        for item in product(l_list, r_list):                    #product(A, B) 和 ((x,y) for x in A for y in B)一样
            # print("sigmoid:%f"%sigmoid(item[0]*item[1]))


            exec('mult_score{}{} += {}'.format(a, b, (item[0] * item[1])))
            exec('mult_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('mult_list_s{}{}.append({})'.format(a, b, (item[0] * item[1])))
            b += 1
            if b == 8:
                b = 0
                a += 1
            if a==8:
                break

        mean_l_malignant = np.mean(dict[i[0]][malignant_index])
        mean_l_Fibro = np.mean(dict[i[0]][Fibro_index])
        mean_l_B = np.mean(dict[i[0]][B_index])
        mean_l_Epithelial = np.mean(dict[i[0]][Epithelial_index])
        mean_l_Macro = np.mean(dict[i[0]][Macro_index])
        mean_l_Endo = np.mean(dict[i[0]][Endo_index])
        mean_l_T = np.mean(dict[i[0]][T_index])
        mean_l_Mast = np.mean(dict[i[0]][Mast_index])
        mean_l = np.mean((mean_l_malignant, mean_l_Fibro, mean_l_B, mean_l_Epithelial, mean_l_Macro, mean_l_Endo, mean_l_T, mean_l_Mast))
        std_l = np.std(dict[i[0]][np.concatenate((malignant_index, Fibro_index, B_index, Epithelial_index, Macro_index,
                                                  Endo_index, T_index, Mast_index))])
        sum_l = np.sum((mean_l_malignant, mean_l_Fibro, mean_l_B, mean_l_Epithelial, mean_l_Macro, mean_l_Endo, mean_l_T, mean_l_Mast))

        mean_r_malignant = np.mean(dict[i[1]][malignant_index])
        mean_r_Fibro = np.mean(dict[i[1]][Fibro_index])
        mean_r_B = np.mean(dict[i[1]][B_index])
        mean_r_Epithelial = np.mean(dict[i[1]][Epithelial_index])
        mean_r_Macro = np.mean(dict[i[1]][Macro_index])
        mean_r_Endo = np.mean(dict[i[1]][Endo_index])
        mean_r_T = np.mean(dict[i[1]][T_index])
        mean_r_Mast = np.mean(dict[i[0]][Mast_index])
        mean_r = np.mean((mean_r_malignant, mean_r_Fibro, mean_r_B, mean_r_Epithelial, mean_r_Macro, mean_r_Endo, mean_r_T, mean_r_Mast))
        std_r = np.std(dict[i[1]][np.concatenate((malignant_index, Fibro_index, B_index, Epithelial_index, Macro_index,
                                                  Endo_index, T_index, Mast_index))])
        sum_r = np.sum((mean_r_malignant, mean_r_Fibro, mean_r_B, mean_r_Epithelial, mean_r_Macro, mean_r_Endo, mean_r_T, mean_r_Mast))

        malignant_l = int(mean_l_malignant > mean_l + std_l)
        malignant_r = int(mean_r_malignant > mean_r + std_r)
        Fibro_l = int(mean_l_Fibro > mean_l + std_l)
        Fibro_r = int(mean_r_Fibro > mean_r + std_r)
        B_l = int(mean_l_B > mean_l + std_l)
        B_r = int(mean_r_B > mean_r + std_r)
        Epithelial_l = int(mean_l_Epithelial > mean_l + std_l)
        Epithelial_r = int(mean_r_Epithelial > mean_r + std_r)
        Macro_l = int(mean_l_Macro > mean_l + std_l)
        Macro_r = int(mean_r_Macro > mean_r + std_r)
        Endo_l = int(mean_l_Endo > mean_l + std_l)
        Endo_r = int(mean_r_Endo > mean_r + std_r)
        T_l = int(mean_l_T > mean_l + std_l)
        T_r = int(mean_r_T > mean_r + std_r)
        Mast_l = int(mean_l_Mast > mean_l + std_l)
        Mast_r = int(mean_r_Mast > mean_r + std_r)
        l_list = [malignant_l, Fibro_l, B_l, Epithelial_l, Macro_l, Endo_l, T_l, Mast_l]
        r_list = [malignant_r, Fibro_r, B_r, Epithelial_r, Macro_r, Endo_r, T_r, Mast_r]
        a = b = 0
        for item in product(l_list, r_list):
            exec('thrd_score{}{} += {}'.format(a, b, int(item[0] & item[1])))
            exec('thrd_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('thrd_list_s{}{}.append({})'.format(a, b, int(item[0] & item[1])))
            b += 1
            if b == 8:
                b = 0
                a += 1
            if a == 8:
                break
         #
        sp_l_malignant = mean_l_malignant/sum_l
        sp_l_Fibro = mean_l_Fibro/sum_l
        sp_l_B = mean_l_B/sum_l
        sp_l_Epithelial = mean_l_Epithelial/sum_l
        sp_l_Macro = mean_l_Macro/sum_l
        sp_l_Endo = mean_l_Endo/sum_l
        sp_l_T = mean_l_T/sum_l
        sp_l_Mast = mean_l_Mast / sum_l
        sp_l_list = [sp_l_malignant, sp_l_Fibro, sp_l_B,  sp_l_Epithelial, sp_l_Macro, sp_l_Endo, sp_l_T, sp_l_Mast]

        sp_r_malignant = mean_r_malignant / sum_r
        sp_r_Fibro = mean_r_Fibro / sum_r
        sp_r_B = mean_r_B / sum_r
        sp_r_Epithelial = mean_r_Epithelial / sum_r
        sp_r_Macro = mean_r_Macro / sum_r
        sp_r_Endo = mean_r_Endo / sum_r
        sp_r_T = mean_r_T / sum_r
        sp_r_Mast = mean_r_Mast / sum_r
        sp_r_list = [sp_r_malignant, sp_r_Fibro, sp_r_B,  sp_r_Epithelial, sp_r_Macro, sp_r_Endo, sp_r_T, sp_r_Mast]

        a = b = 0
        for item in product(sp_l_list, sp_r_list):                    #product(A, B) 和 ((x,y) for x in A for y in B)一样
            # print("sigmoid:%f"%sigmoid(item[0]*item[1]))

            exec('spec_score{}{} += {}'.format(a, b, (item[0] * item[1])))
            exec('spec_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('spec_list_s{}{}.append({})'.format(a, b, (item[0] * item[1])))
            b += 1
            if b == 8:
                b = 0
                a += 1
            if a == 8:
                break

        total_l_malignant = np.sum(dict[i[0]][malignant_index])
        total_l_T = np.sum(dict[i[0]][T_index])
        total_l_B = np.sum(dict[i[0]][B_index])
        total_l_Epithelial = np.sum(dict[i[0]][Epithelial_index])
        total_l_Fibroblasts = np.sum(dict[i[0]][Fibro_index])
        total_l_Mast = np.sum(dict[i[0]][Mast_index])
        total_l_Endothelial = np.sum(dict[i[0]][Endo_index])
        total_l_Macrophage = np.sum(dict[i[0]][Macro_index])
        total_l_list = [total_l_malignant, total_l_T, total_l_B, total_l_Epithelial, total_l_Fibroblasts, total_l_Mast, total_l_Endothelial, total_l_Macrophage]


        total_r_malignant = np.sum(dict[i[1]][malignant_index])
        total_r_T = np.sum(dict[i[1]][T_index])
        total_r_B = np.sum(dict[i[1]][B_index])
        total_r_Epithelial = np.sum(dict[i[1]][Epithelial_index])
        total_r_Fibroblasts = np.sum(dict[i[1]][Fibro_index])
        total_r_Mast = np.sum(dict[i[1]][Mast_index])
        total_r_Endothelial = np.sum(dict[i[1]][Endo_index])
        total_r_Macrophage = np.sum(dict[i[1]][Macro_index])
        total_r_list = [total_r_malignant, total_r_T, total_r_B, total_r_Epithelial, total_r_Fibroblasts, total_r_Mast, total_r_Endothelial, total_r_Macrophage]
        a = b = 0
        for item in product(total_l_list, total_r_list):
            exec('total_expr_score{}{} += {}'.format(a, b, item[0] * item[1]))
            exec('total_expr_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('total_expr_list_s{}{}.append({})'.format(a, b, item[0] * item[1]))
            b += 1
            if b == 8:
                b = 0
                a += 1
            if a == 8:
                break
        regpro_l_malignant = np.sum(dict[i[0]][malignant_index])
        regpro_l_T = np.sum(dict[i[0]][T_index])
        regpro_l_B = np.sum(dict[i[0]][B_index])
        regpro_l_Epithelial = np.sum(dict[i[0]][Epithelial_index])
        regpro_l_Fibroblasts = np.sum(dict[i[0]][Fibro_index])
        regpro_l_Mast = np.sum(dict[i[0]][Mast_index])
        regpro_l_Endothelial = np.sum(dict[i[0]][Endo_index])
        regpro_l_Macrophage = np.sum(dict[i[0]][Macro_index])
        regpro_l_list = [regpro_l_malignant, regpro_l_T, regpro_l_B, regpro_l_Epithelial, regpro_l_Fibroblasts, regpro_l_Mast, regpro_l_Endothelial, regpro_l_Macrophage]


        regpro_r_malignant = np.sum(dict[i[1]][malignant_index])
        regpro_r_T = np.sum(dict[i[1]][T_index])
        regpro_r_B = np.sum(dict[i[1]][B_index])
        regpro_r_Epithelial = np.sum(dict[i[1]][Epithelial_index])
        regpro_r_Fibroblasts = np.sum(dict[i[1]][Fibro_index])
        regpro_r_Mast = np.sum(dict[i[1]][Mast_index])
        regpro_r_Endothelial = np.sum(dict[i[1]][Endo_index])
        regpro_r_Macrophage = np.sum(dict[i[1]][Macro_index])
        regpro_r_list = [regpro_r_malignant, regpro_r_T, regpro_r_B, regpro_r_Epithelial, regpro_r_Fibroblasts, regpro_r_Mast, regpro_r_Endothelial, regpro_r_Macrophage]
        a = b = 0
        for item in product(regpro_l_list, regpro_r_list):
            score = np.sqrt(item[0] * item[1]) / (mu + np.sqrt(item[0] * item[1]))
            if score > 0.5:
                exec('reg_prod_score{}{} += {}'.format(a, b, score))
                exec('reg_prod_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
                exec('reg_prod_list_s{}{}.append({})'.format(a, b, score))
            b += 1
            if b == 6:
                b = 0
                a += 1
            if a == 6:
                break

        a = b = 0
        # 为每个配体-受体对初始化细胞类型表达列表
        cell_expr_values = []
        for cell_type in cell_types:
            c_i = dict[i[0]][cell_type_indices[cell_type]].sum()
            c_j = dict[i[1]][cell_type_indices[cell_type]].sum()
            cell_expr_values.append((c_i, c_j))

        # 计算细胞表达法的打分
        for item in product(cell_expr_values, repeat=2):
            c_i_k1, c_j_k1 = item[0]
            c_i_k2, c_j_k2 = item[1]
            n1 = len(cell_type_indices[cell_types[a]])
            n2 = len(cell_type_indices[cell_types[b]])
            score = (c_i_k1 / n1) * (c_j_k2 / n2)  # 注意：这里假设配体和受体在不同的细胞类型中

            if score > 0:  # 根据需要设置阈值
                exec('cell_expr_score{}{} += {}'.format(a, b, score))
                exec('cell_expr_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
                exec('cell_expr_list_s{}{}.append({})'.format(a, b, score))
            b += 1
            if b == 6:
                b = 0
                a += 1
            if a == 6:
                break


for i in range(8):
    for j in range(8):
        with open(savepath + "mult_score.txt", "a") as f:
            exec('f.write("mult_score{}{} = %f"%(mult_score{}{}))'.format(i, j, i, j))
            f.write('\n')
for i in range(8):
    for j in range(8):
        with open(savepath + "thrd_score.txt","a") as f:
            exec('f.write("thrd_score{}{} = %f"%(thrd_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(8):
    for j in range(8):
        exec('x = pd.DataFrame(mult_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(mult_list_s{}{},columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "mult_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(8):
    for j in range(8):
        exec('x = pd.DataFrame(thrd_list{}{})'.format(i,j))
        exec('y = pd.DataFrame(thrd_list_s{}{},columns=list("3"))'.format(i,j))
        x = x.join(y)
        exec('x.to_csv(savepath + "thrd_list{}{}.csv", header=None, index=None)'.format(i, j))
for i in range(8):
    for j in range(8):
        with open(savepath + "spec_score.txt", "a") as f:
            exec('f.write("spec_score{}{} = %f"%(spec_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(8):
    for j in range(8):
        exec('x = pd.DataFrame(spec_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(spec_list_s{}{},columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "spec_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(8):
    for j in range(8):
        with open(savepath + "total_expr_score.txt", "a") as f:
            exec('f.write("total_expr_score{}{} = %f"%(total_expr_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(8):
    for j in range(8):
        exec('x = pd.DataFrame(total_expr_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(total_expr_list_s{}{}, columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "total_expr_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(8):
    for j in range(8):
        with open(os.path.join(savepath + "reg_prod_score.txt"), "a") as f:
            exec('f.write("reg_prod_score{}{} = %f"%(reg_prod_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(8):
    for j in range(8):
        exec('x = pd.DataFrame(reg_prod_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(reg_prod_list_s{}{}, columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "reg_prod_list{}{}.csv", header=None, index=None)'.format(i, j))

# 将细胞表达法的打分结果写入文件
for i in range(8):
    for j in range(8):
        with open(savepath + "cell_expr_score.txt", "a") as f:
            exec('f.write("cell_expr_score{}{} = %f"%(cell_expr_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(8):
    for j in range(8):
        exec('x = pd.DataFrame(cell_expr_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(cell_expr_list_s{}{}, columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "cell_expr_list{}{}.csv".format(i, j), header=None, index=None)')

import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
mm = MinMaxScaler()
mult = []
thrd = []
spec = []
total_exp = []
reg_prod = []
cell_expr = []
for i in range(8):
    for j in range(8):
        exec('mult.append(mult_score{}{})'.format(i, j))
        exec('thrd.append(thrd_score{}{})'.format(i, j))
        exec('spec.append(spec_score{}{})'.format(i, j))
        exec('total_exp.append(total_expr_score{}{})'.format(i, j))
        exec('reg_prod.append(reg_prod_score{}{})'.format(i, j))
        exec('cell_expr.append(cell_expr_score{}{})'.format(i, j))
mult = np.array(mult).reshape(8, 8)
thrd = np.array(thrd).reshape(8, 8)
spec = np.array(spec).reshape(8, 8)
total_exp = np.array(total_exp).reshape(8, 8)
reg_prod = np.array(reg_prod).reshape(8, 8)
cell_exp = np.array(cell_expr).reshape(8, 8)

multy = pd.DataFrame(mult)
thrdy = pd.DataFrame(thrd)
specy = pd.DataFrame(spec)
totaly = pd.DataFrame(total_exp)
reg_prody = pd.DataFrame(reg_prod)
cell_expy = pd.DataFrame(cell_exp)

multy.to_csv(savepath + "multy.csv", header=None, index=None)
thrdy.to_csv(savepath + "thrdy.csv", header=None, index=None)
specy.to_csv(savepath + "specy.csv", header=None, index=None)
totaly.to_csv(savepath + "totaly.csv", header=None, index=None)
reg_prody.to_csv(savepath + "reg_prody.csv", header=None, index=None)
cell_expy.to_csv(savepath + "cell_expy.csv", header=None, index=None)