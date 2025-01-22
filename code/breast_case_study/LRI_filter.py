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

dt = pd.read_csv('Breast_1.csv', index_col=0, header=None)
#a = np.array(dt.loc["Cell"])
dict = {}
dict["cell"] = np.array(dt.loc["cell"])
for i in range(1, dt.shape[0]):
    dict[dt.index[i]] = np.array(dt.loc[dt.index[i]], dtype=float)

# def sigmoid(x):
#     return 1/(1+np.exp(-(x-6)))

savepath = 'Breast_'

# 0=cancer,1=Immune,2=Stromal,3=T,4=B,5=Myeloid
malignant_index = np.where(dict['cell type'] == 0)[0]
Immune = np.where(dict['cell type'] == 1)[0]
Stromal = np.where(dict['cell type'] == 2)[0]
T = np.where(dict['cell type'] == 3)[0]
B = np.where(dict['cell type'] == 4)[0]
Myeloid = np.where(dict['cell type'] == 5)[0]

###
groupsize = []
groupsize.append(malignant_index.shape[0])
groupsize.append(T.shape[0])
groupsize.append(B.shape[0])
groupsize.append(Myeloid.shape[0])
groupsize.append(Immune.shape[0])
groupsize.append(Stromal .shape[0])
groupsize = pd.DataFrame(groupsize)
groupsize.to_csv('Breast/Breast.csv', header=None, index=None)


for i in range(6):
    for j in range(6):
        exec('mult_score{}{} = 0'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('thrd_score{}{} = 0'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('mult_list{}{} = []'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('mult_list_s{}{} = []'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('thrd_list{}{} = []'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('thrd_list_s{}{} = []'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('spec_score{}{} = 0'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('spec_list{}{} = []'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('spec_list_s{}{} = []'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('total_expr_score{}{} = 0'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('total_expr_list{}{} = []'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('total_expr_list_s{}{} = []'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('reg_prod_score{}{} = 0'.format(i, j))
        exec('reg_prod_list{}{} = []'.format(i, j))
        exec('reg_prod_list_s{}{} = []'.format(i, j))

# 细胞表达法的打分初始化
for i in range(6):
    for j in range(6):
        exec('cell_expr_score{}{} = 0'.format(i, j))
        exec('cell_expr_list{}{} = []'.format(i, j))
        exec('cell_expr_list_s{}{} = []'.format(i, j))
# 定义细胞类型列表
cell_types = ["malignant", "Immune", "Stromal", "T", "B", "Myeloid"]
cell_type_indices = {
    "malignant": malignant_index,
    "T": T,
    "B": B,
    "Myeloid": Myeloid,
    "Immune": Immune,
    "Stromal": Stromal
}

g = 0
mu = np.mean([dict[gene].mean() for gene in dict if gene != "cell"])

for i in LRI_gene:
    if i[0] in dict and i[1] in dict:
        g= g+1
        print(g)
        malignant_l = 1 / malignant_index.shape[0] * sum(dict[i[0]][malignant_index])
        malignant_r = 1 / malignant_index.shape[0] * sum(dict[i[1]][malignant_index])
        T_l = 1 / T.shape[0] * sum(dict[i[0]][T])
        T_r = 1 / T.shape[0] * sum(dict[i[1]][T])
        B_l = 1 / B.shape[0] * sum(dict[i[0]][B])
        B_r = 1 / B.shape[0] * sum(dict[i[1]][B])
        Myeloid_l = 1 / Myeloid.shape[0] * sum(dict[i[0]][Myeloid])
        Myeloid_r = 1 / Myeloid.shape[0] * sum(dict[i[1]][Myeloid])
        Immune_l = 1 / Immune.shape[0] * sum(dict[i[0]][Immune])
        Immune_r = 1 / Immune.shape[0] * sum(dict[i[1]][Immune])
        Stromal_l = 1 / Stromal.shape[0] * sum(dict[i[0]][Stromal])
        Stromal_r = 1 / Stromal.shape[0] * sum(dict[i[1]][Stromal])
        l_list = [malignant_l, T_l, B_l, Myeloid_l, Immune_l, T_l, Stromal_l]
        r_list = [malignant_r, T_r, B_r, Myeloid_r, Immune_r, T_r, Stromal_r]

        a = b = 0
        for item in product(l_list, r_list):                    #product(A, B) 和 ((x,y) for x in A for y in B)一样
            # print("sigmoid:%f"%sigmoid(item[0]*item[1]))


            exec('mult_score{}{} += {}'.format(a, b, (item[0] * item[1])))
            exec('mult_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('mult_list_s{}{}.append({})'.format(a, b, (item[0] * item[1])))
            b += 1
            if b == 6:
                b = 0
                a += 1
            if a == 6:
                break


        mean_l_malignant = np.mean(dict[i[0]][malignant_index])
        mean_l_T = np.mean(dict[i[0]][T])
        mean_l_B = np.mean(dict[i[0]][B])
        mean_l_Myeloid = np.mean(dict[i[0]][Myeloid])
        mean_l_Immune = np.mean(dict[i[0]][Immune])
        mean_l_Stromal = np.mean(dict[i[0]][Stromal])
        mean_l = np.mean((mean_l_malignant, mean_l_T, mean_l_B, mean_l_Myeloid, mean_l_Immune, mean_l_Stromal))
        std_l = np.std(dict[i[0]][np.concatenate((malignant_index, T, B,  Myeloid, Immune, Stromal))])
        sum_l = np.sum((mean_l_malignant, mean_l_T, mean_l_B, mean_l_Myeloid, mean_l_Immune, mean_l_Stromal))

        mean_r_malignant = np.mean(dict[i[1]][malignant_index])
        mean_r_T = np.mean(dict[i[1]][T])
        mean_r_B = np.mean(dict[i[1]][B])
        mean_r_Myeloid = np.mean(dict[i[1]][Myeloid])
        mean_r_Immune = np.mean(dict[i[1]][Immune])
        mean_r_Stromal = np.mean(dict[i[1]][Stromal])
        mean_r = np.mean((mean_r_malignant, mean_r_T, mean_r_B, mean_r_Myeloid, mean_r_Immune,
                          mean_r_Stromal))
        std_r = np.std(dict[i[1]][np.concatenate((malignant_index, T, B,  Myeloid,
                                                  Immune, Stromal))])
        sum_r = np.sum((mean_r_malignant, mean_r_T, mean_r_B,  mean_r_Myeloid, mean_r_Immune, mean_r_Stromal))

        malignant_l = int(mean_l_malignant > mean_l + std_l)
        malignant_r = int(mean_r_malignant > mean_r + std_r)
        T_l = int(mean_l_T > mean_l + std_l)
        T_r = int(mean_r_T > mean_r + std_r)
        B_l = int(mean_l_B > mean_l + std_l)
        B_r = int(mean_r_B > mean_r + std_r)
        Myeloid_l = int(mean_l_Myeloid > mean_l + std_l)
        Myeloid_r = int(mean_r_Myeloid > mean_r + std_r)
        Immune_l = int(mean_l_Immune > mean_l + std_l)
        Immune_r = int(mean_r_Immune > mean_r + std_r)
        Stromal_l = int(mean_l_Stromal > mean_l + std_l)
        Stromal_r = int(mean_r_Stromal > mean_r + std_r)
        l_list = [malignant_l, T_l, B_l, Myeloid_l, Immune_l, Stromal_l, ]
        r_list = [malignant_r, T_r, B_r, Myeloid_r, Immune_r, Stromal_r, ]
        a = b = 0
        for item in product(l_list, r_list):
            exec('thrd_score{}{} += {}'.format(a, b, int(item[0] & item[1])))
            exec('thrd_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('thrd_list_s{}{}.append({})'.format(a, b, int(item[0] & item[1])))
            b += 1
            if b == 6:
                b = 0
                a += 1
            if a == 6:
                break
         #
        sp_l_malignant = mean_l_malignant/sum_l
        sp_l_T = mean_l_T/sum_l
        sp_l_B = mean_l_B/sum_l
        sp_l_Myeloid = mean_l_Myeloid/sum_l
        sp_l_Immune = mean_l_Immune/sum_l
        sp_l_Stromal = mean_l_Stromal / sum_l
        sp_l_list = [sp_l_malignant, sp_l_T, sp_l_B,  sp_l_Myeloid, sp_l_Immune, sp_l_Stromal]

        sp_r_malignant = mean_r_malignant / sum_r
        sp_r_T = mean_r_T / sum_r
        sp_r_B = mean_r_B / sum_r
        sp_r_Myeloid = mean_r_Myeloid / sum_r
        sp_r_Immune = mean_r_Immune / sum_r
        sp_r_Stromal = mean_r_Stromal / sum_r
        sp_r_list = [sp_r_malignant, sp_r_T, sp_r_B, sp_r_Myeloid, sp_r_Immune, sp_r_Stromal]

        a = b = 0
        for item in product(sp_l_list, sp_r_list):                    #product(A, B) 和 ((x,y) for x in A for y in B)一样
            # print("sigmoid:%f"%sigmoid(item[0]*item[1]))

            exec('spec_score{}{} += {}'.format(a, b, (item[0] * item[1])))
            exec('spec_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('spec_list_s{}{}.append({})'.format(a, b, (item[0] * item[1])))
            b += 1
            if b == 6:
                b = 0
                a += 1
            if a == 6:
                break

        total_l_malignant = np.sum(dict[i[0]][malignant_index])
        total_l_T = np.sum(dict[i[0]][T])
        total_l_B = np.sum(dict[i[0]][B])
        total_l_Myeloid = np.sum(dict[i[0]][Myeloid])
        total_l_Immune = np.sum(dict[i[0]][Immune])
        total_l_Stromal = np.sum(dict[i[0]][Stromal])
        total_l_list = [total_l_malignant, total_l_T, total_l_B, total_l_Myeloid, total_l_Immune, total_l_Stromal]


        total_r_malignant = np.sum(dict[i[1]][malignant_index])
        total_r_T = np.sum(dict[i[1]][T])
        total_r_B = np.sum(dict[i[1]][B])
        total_r_Myeloid = np.sum(dict[i[1]][Myeloid])
        total_r_Immune = np.sum(dict[i[1]][Immune])
        total_r_Stromal = np.sum(dict[i[1]][Stromal])
        total_r_list = [total_r_malignant, total_r_T, total_r_B, total_r_Myeloid, total_r_Immune, total_r_Stromal]

        a = b = 0
        for item in product(total_l_list, total_r_list):
            exec('total_expr_score{}{} += {}'.format(a, b, item[0] * item[1]))
            exec('total_expr_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('total_expr_list_s{}{}.append({})'.format(a, b, item[0] * item[1]))
            b += 1
            if b == 6:
                b = 0
                a += 1
            if a == 6:
                break

        # 计算每个细胞类型中配体和受体的表达量
        regpro_l_malignant = np.sum(dict[i[0]][malignant_index])
        regpro_l_T = np.sum(dict[i[0]][T])
        regpro_l_B = np.sum(dict[i[0]][B])
        regpro_l_Myeloid = np.sum(dict[i[0]][Myeloid])
        regpro_l_Immune = np.sum(dict[i[0]][Immune])
        regpro_l_Stromal = np.sum(dict[i[0]][Stromal])
        regpro_l_list = [regpro_l_malignant, regpro_l_T, regpro_l_B, regpro_l_Myeloid, regpro_l_Immune,
                         regpro_l_Stromal]

        regpro_r_malignant = np.sum(dict[i[1]][malignant_index])
        regpro_r_T = np.sum(dict[i[1]][T])
        regpro_r_B = np.sum(dict[i[1]][B])
        regpro_r_Myeloid = np.sum(dict[i[1]][Myeloid])
        regpro_r_Immune = np.sum(dict[i[1]][Immune])
        regpro_r_Stromal = np.sum(dict[i[1]][Stromal])


        regpro_r_list = [regpro_r_malignant, regpro_r_T, regpro_r_B, regpro_r_Myeloid, regpro_r_Immune, regpro_r_Stromal]

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


for i in range(6):
    for j in range(6):
        with open(savepath + "mult_score.txt", "a") as f:
            exec('f.write("mult_score{}{} = %f"%(mult_score{}{}))'.format(i, j, i, j))
            f.write('\n')
for i in range(6):
    for j in range(6):
        with open(savepath + "thrd_score.txt","a") as f:
            exec('f.write("thrd_score{}{} = %f"%(thrd_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(6):
    for j in range(6):
        exec('x = pd.DataFrame(mult_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(mult_list_s{}{},columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "mult_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(6):
    for j in range(6):
        exec('x = pd.DataFrame(thrd_list{}{})'.format(i,j))
        exec('y = pd.DataFrame(thrd_list_s{}{},columns=list("3"))'.format(i,j))
        x = x.join(y)
        exec('x.to_csv(savepath + "thrd_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(6):
    for j in range(6):
        with open(savepath + "spec_score.txt", "a") as f:
            exec('f.write("spec_score{}{} = %f"%(spec_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(6):
    for j in range(6):
        exec('x = pd.DataFrame(spec_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(spec_list_s{}{},columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "spec_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(6):
    for j in range(6):
        with open(savepath + "total_expr_score.txt", "a") as f:
            exec('f.write("total_expr_score{}{} = %f"%(total_expr_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(6):
    for j in range(6):
        exec('x = pd.DataFrame(total_expr_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(total_expr_list_s{}{}, columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "total_expr_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(6):
    for j in range(6):
        with open(savepath + "reg_prod_score.txt", "a") as f:
            exec('f.write("reg_prod_score{}{} = %f"%(reg_prod_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(6):
    for j in range(6):
        exec('x = pd.DataFrame(reg_prod_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(reg_prod_list_s{}{}, columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "reg_prod_list{}{}.csv", header=None, index=None)'.format(i, j))

# 将细胞表达法的打分结果写入文件
for i in range(6):
    for j in range(6):
        with open(savepath + "cell_expr_score.txt", "a") as f:
            exec('f.write("cell_expr_score{}{} = %f"%(cell_expr_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(6):
    for j in range(6):
        exec('x = pd.DataFrame(cell_expr_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(cell_expr_list_s{}{}, columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "cell_expr_list{}{}.csv".format(i, j), header=None, index=None)')

import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

def min_max_scaling(data):
    return (data - data.min()) / (data.max() - data.min())

mm = MinMaxScaler()
mult = []
thrd = []
spec = []
total_exp = []
reg_prod = []
cell_expr = []

for i in range(6):
    for j in range(6):
        exec('mult.append(mult_score{}{})'.format(i, j))
        exec('thrd.append(thrd_score{}{})'.format(i, j))
        exec('spec.append(spec_score{}{})'.format(i, j))
        exec('total_exp.append(total_expr_score{}{})'.format(i, j))
        exec('reg_prod.append(reg_prod_score{}{})'.format(i, j))
        exec('cell_expr.append(cell_expr_score{}{})'.format(i, j))
        
        
        
mult = np.array(mult).reshape(6, 6)
thrd = np.array(thrd).reshape(6, 6)
spec = np.array(spec).reshape(6, 6)
total_exp = np.array(total_exp).reshape(6, 6)
reg_prod = np.array(reg_prod).reshape(6, 6)
cell_exp = np.array(cell_expr).reshape(6, 6)


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
