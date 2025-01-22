import os
import warnings
warnings.filterwarnings("ignore")
from itertools import product as product
import pandas as pd
import numpy as np
import psutil
from time import time
itime = time()

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

dt = pd.read_csv('melanoma.csv', index_col=0, header=None)
#a = np.array(dt.loc["Cell"])
dict = {}
dict["Cell"] = np.array(dt.loc["Cell"])
for i in range(1, dt.shape[0]):
    dict[dt.index[i]] = np.array(dt.loc[dt.index[i]], dtype=float)

# def sigmoid(x):
#     return 1/(1+np.exp(-(x-6)))

savepath = 'GSE72056'

#0=m,1=T,2=B,3=myo,4=Macro,5=Endo,8=T
malignant_index = np.where(dict['malignant'] == 2)[0]
T = np.where(dict['non-cancer cell type'] == 1)[0]
B = np.where(dict['non-cancer cell type'] == 2)[0]
Macro = np.where(dict['non-cancer cell type'] == 3)[0]
Endo = np.where(dict['non-cancer cell type'] == 4)[0]
CAF = np.where(dict['non-cancer cell type'] == 5)[0]
NK = np.where(dict['non-cancer cell type'] == 6)[0]
###
groupsize = []
groupsize.append(malignant_index.shape[0])
groupsize.append(T.shape[0])
groupsize.append(B.shape[0])
groupsize.append(Macro.shape[0])
groupsize.append(Endo.shape[0])
groupsize.append(CAF.shape[0])
groupsize.append(NK.shape[0])
groupsize = pd.DataFrame(groupsize)
groupsize.to_csv('out/Melanoma.csv', header=None, index=None)


for i in range(7):
    for j in range(7):
        exec('mult_score{}{} = 0'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('thrd_score{}{} = 0'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('mult_list{}{} = []'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('mult_list_s{}{} = []'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('thrd_list{}{} = []'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('thrd_list_s{}{} = []'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('spec_score{}{} = 0'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('spec_list{}{} = []'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('spec_list_s{}{} = []'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('total_expr_score{}{} = 0'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('total_expr_list{}{} = []'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('total_expr_list_s{}{} = []'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('reg_prod_score{}{} = 0'.format(i, j))
        exec('reg_prod_list{}{} = []'.format(i, j))
        exec('reg_prod_list_s{}{} = []'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('cell_expr_score{}{} = 0'.format(i, j))
        exec('cell_expr_list{}{} = []'.format(i, j))
        exec('cell_expr_list_s{}{} = []'.format(i, j))

# 定义细胞类型列表
cell_types = ['Melanoma cancer cells', 'T cells', 'B cells', 'Macrophages', 'Endothelial cells', 'CAFs','NK cells']
cell_type_indices = {
    "Melanoma cancer cells": malignant_index,
    "T cells": T,
    "B cells": B,
    "Macrophages": Macro,
    "Endothelial cells": Endo,
    "CAFs": CAF,
    'NK cells': NK
}

g = 0

mu = np.mean([dict[gene].mean() for gene in dict if gene != "Cell"])

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
        Macro_l = 1 / Macro.shape[0] * sum(dict[i[0]][Macro])
        Macro_r = 1 / Macro.shape[0] * sum(dict[i[1]][Macro])
        Endo_l = 1 / Endo.shape[0] * sum(dict[i[0]][Endo])
        Endo_r = 1 / Endo.shape[0] * sum(dict[i[1]][Endo])
        CAF_l = 1 / CAF.shape[0] * sum(dict[i[0]][CAF])
        CAF_r = 1 / CAF.shape[0] * sum(dict[i[1]][CAF])
        NK_l = 1 / NK.shape[0] * sum(dict[i[0]][NK])
        NK_r = 1 / NK.shape[0] * sum(dict[i[1]][NK])
        l_list = [malignant_l, T_l, B_l, Macro_l, Endo_l, T_l, CAF_l, NK_l]
        r_list = [malignant_r, T_r, B_r, Macro_r, Endo_r, T_r, CAF_r, NK_r]

        a = b = 0
        for item in product(l_list, r_list):                    #product(A, B) 和 ((x,y) for x in A for y in B)一样
            # print("sigmoid:%f"%sigmoid(item[0]*item[1]))


            exec('mult_score{}{} += {}'.format(a, b, (item[0] * item[1])))
            exec('mult_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('mult_list_s{}{}.append({})'.format(a, b, (item[0] * item[1])))
            b += 1
            if b == 7:
                b = 0
                a += 1
            if a == 7:
                break


        mean_l_malignant = np.mean(dict[i[0]][malignant_index])
        mean_l_T = np.mean(dict[i[0]][T])
        mean_l_B = np.mean(dict[i[0]][B])
        mean_l_Macro = np.mean(dict[i[0]][Macro])
        mean_l_Endo = np.mean(dict[i[0]][Endo])
        mean_l_CAF = np.mean(dict[i[0]][CAF])
        mean_l_NK = np.mean(dict[i[0]][NK])
        mean_l = np.mean((mean_l_malignant, mean_l_T, mean_l_B, mean_l_Macro, mean_l_Endo, mean_l_CAF, mean_l_NK))
        std_l = np.std(dict[i[0]][np.concatenate((malignant_index, T, B,  Macro,
                                                  Endo, CAF, NK))])
        sum_l = np.sum((mean_l_malignant, mean_l_T, mean_l_B, mean_l_Macro, mean_l_Endo, mean_l_CAF, mean_l_NK))

        mean_r_malignant = np.mean(dict[i[1]][malignant_index])
        mean_r_T = np.mean(dict[i[1]][T])
        mean_r_B = np.mean(dict[i[1]][B])
        mean_r_Macro = np.mean(dict[i[1]][Macro])
        mean_r_Endo = np.mean(dict[i[1]][Endo])
        mean_r_CAF = np.mean(dict[i[0]][CAF])
        mean_r_NK = np.mean(dict[i[0]][NK])
        mean_r = np.mean((mean_r_malignant, mean_r_T, mean_r_B, mean_r_Macro, mean_r_Endo,
                          mean_r_CAF, mean_r_NK))
        std_r = np.std(dict[i[1]][np.concatenate((malignant_index, T, B,  Macro,
                                                  Endo, CAF, NK))])
        sum_r = np.sum((mean_r_malignant, mean_r_T, mean_r_B,  mean_r_Macro, mean_r_Endo, mean_r_CAF, mean_r_NK))

        malignant_l = int(mean_l_malignant > mean_l + std_l)
        malignant_r = int(mean_r_malignant > mean_r + std_r)
        T_l = int(mean_l_T > mean_l + std_l)
        T_r = int(mean_r_T > mean_r + std_r)
        B_l = int(mean_l_B > mean_l + std_l)
        B_r = int(mean_r_B > mean_r + std_r)
        Macro_l = int(mean_l_Macro > mean_l + std_l)
        Macro_r = int(mean_r_Macro > mean_r + std_r)
        Endo_l = int(mean_l_Endo > mean_l + std_l)
        Endo_r = int(mean_r_Endo > mean_r + std_r)
        CAF_l = int(mean_l_CAF > mean_l + std_l)
        CAF_r = int(mean_r_CAF > mean_r + std_r)
        NK_l = int(mean_l_NK > mean_l + std_l)
        NK_r = int(mean_r_NK > mean_r + std_r)
        l_list = [malignant_l, T_l, B_l, Macro_l, Endo_l, CAF_l, NK_l]
        r_list = [malignant_r, T_r, B_r, Macro_r, Endo_r, CAF_r, NK_r]
        a = b = 0
        for item in product(l_list, r_list):
            exec('thrd_score{}{} += {}'.format(a, b, int(item[0] & item[1])))
            exec('thrd_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('thrd_list_s{}{}.append({})'.format(a, b, int(item[0] & item[1])))
            b += 1
            if b == 7:
                b = 0
                a += 1
            if a == 7:
                break
         #
        sp_l_malignant = mean_l_malignant/sum_l
        sp_l_T = mean_l_T/sum_l
        sp_l_B = mean_l_B/sum_l
        sp_l_Macro = mean_l_Macro/sum_l
        sp_l_Endo = mean_l_Endo/sum_l
        sp_l_CAF = mean_l_CAF / sum_l
        sp_l_NK = mean_l_NK / sum_l
        sp_l_list = [sp_l_malignant, sp_l_T, sp_l_B,  sp_l_Macro, sp_l_Endo, sp_l_CAF, sp_l_NK]

        sp_r_malignant = mean_r_malignant / sum_r
        sp_r_T = mean_r_T / sum_r
        sp_r_B = mean_r_B / sum_r
        sp_r_Macro = mean_r_Macro / sum_r
        sp_r_Endo = mean_r_Endo / sum_r
        sp_r_CAF = mean_r_CAF / sum_r
        sp_r_NK = mean_r_NK / sum_r
        sp_r_list = [sp_r_malignant, sp_r_T, sp_r_B, sp_r_Macro, sp_r_Endo, sp_r_CAF,
                     sp_r_NK]

        a = b = 0
        for item in product(sp_l_list, sp_r_list):                    #product(A, B) 和 ((x,y) for x in A for y in B)一样
            # print("sigmoid:%f"%sigmoid(item[0]*item[1]))

            exec('spec_score{}{} += {}'.format(a, b, (item[0] * item[1])))
            exec('spec_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('spec_list_s{}{}.append({})'.format(a, b, (item[0] * item[1])))
            b += 1
            if b == 7:
                b = 0
                a += 1
            if a == 7:
                break
        total_l_malignant = np.sum(dict[i[0]][malignant_index])
        total_l_T = np.sum(dict[i[0]][T])
        total_l_B = np.sum(dict[i[0]][B])
        total_l_Endo = np.sum(dict[i[0]][Endo])
        total_l_Macro = np.sum(dict[i[0]][Macro])
        total_l_CAF = np.sum(dict[i[0]][CAF])
        total_l_NK = np.sum(dict[i[0]][NK])
        total_l_list = [total_l_malignant, total_l_T, total_l_B, total_l_Endo, total_l_Macro, total_l_CAF, total_l_NK]


        total_r_malignant = np.sum(dict[i[1]][malignant_index])
        total_r_T = np.sum(dict[i[1]][T])
        total_r_B = np.sum(dict[i[1]][B])
        total_r_Endo = np.sum(dict[i[1]][Endo])
        total_r_Macro = np.sum(dict[i[1]][Macro])
        total_r_CAF = np.sum(dict[i[1]][CAF])
        total_r_NK = np.sum(dict[i[1]][NK])
        total_r_list = [total_r_malignant, total_r_T, total_r_B, total_r_Endo, total_r_Macro, total_r_CAF, total_r_NK]

        a = b = 0
        for item in product(total_l_list, total_r_list):
            exec('total_expr_score{}{} += {}'.format(a, b, item[0] * item[1]))
            exec('total_expr_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('total_expr_list_s{}{}.append({})'.format(a, b, item[0] * item[1]))
            b += 1
            if b == 7:
                b = 0
                a += 1
            if a == 7:
                break

        regpro_l_malignant = np.sum(dict[i[0]][malignant_index])
        regpro_l_T = np.sum(dict[i[0]][T])
        regpro_l_B = np.sum(dict[i[0]][B])
        regpro_l_Endo = np.sum(dict[i[0]][Endo])
        regpro_l_Macro = np.sum(dict[i[0]][Macro])
        regpro_l_CAF = np.sum(dict[i[0]][CAF])
        regpro_l_NK = np.sum(dict[i[0]][NK])
        regpro_l_list = [regpro_l_malignant, regpro_l_T, regpro_l_B, regpro_l_Endo, regpro_l_Macro, regpro_l_CAF, regpro_l_NK]


        regpro_r_malignant = np.sum(dict[i[1]][malignant_index])
        regpro_r_T = np.sum(dict[i[1]][T])
        regpro_r_B = np.sum(dict[i[1]][B])
        regpro_r_Endo = np.sum(dict[i[1]][Endo])
        regpro_r_Macro = np.sum(dict[i[1]][Macro])
        regpro_r_CAF = np.sum(dict[i[1]][CAF])
        regpro_r_NK = np.sum(dict[i[1]][NK])
        regpro_r_list = [regpro_r_malignant, regpro_r_T, regpro_r_B, regpro_r_Endo, regpro_r_Macro, regpro_r_CAF, regpro_r_NK]

        a = b = 0
        for item in product(regpro_l_list, regpro_r_list):
            score = np.sqrt(item[0] * item[1]) / (mu + np.sqrt(item[0] * item[1]))
            if score > 0.5:
                exec('reg_prod_score{}{} += {}'.format(a, b, score))
                exec('reg_prod_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
                exec('reg_prod_list_s{}{}.append({})'.format(a, b, score))
            else:
                print("Score not greater than 0.5:", score)
            b += 1
            if b == 7:
                b = 0
                a += 1
            if a == 7:
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
            if b == 7:
                b = 0
                a += 1
            if a == 7:
                break

for i in range(7):
    for j in range(7):
        with open(savepath + "mult_score.txt", "a") as f:
            exec('f.write("mult_score{}{} = %f"%(mult_score{}{}))'.format(i, j, i, j))
            f.write('\n')
for i in range(7):
    for j in range(7):
        with open(savepath + "thrd_score.txt","a") as f:
            exec('f.write("thrd_score{}{} = %f"%(thrd_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(7):
    for j in range(7):
        exec('x = pd.DataFrame(mult_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(mult_list_s{}{},columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "mult_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(7):
    for j in range(7):
        exec('x = pd.DataFrame(thrd_list{}{})'.format(i,j))
        exec('y = pd.DataFrame(thrd_list_s{}{},columns=list("3"))'.format(i,j))
        x = x.join(y)
        exec('x.to_csv(savepath + "thrd_list{}{}.csv", header=None, index=None)'.format(i, j))
for i in range(7):
    for j in range(7):
        with open(savepath + "spec_score.txt", "a") as f:
            exec('f.write("spec_score{}{} = %f"%(spec_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(7):
    for j in range(7):
        exec('x = pd.DataFrame(spec_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(spec_list_s{}{},columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "spec_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(7):
    for j in range(7):
        with open(savepath + "total_expr_score.txt", "a") as f:
            exec('f.write("total_expr_score{}{} = %f"%(total_expr_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(7):
    for j in range(7):
        exec('x = pd.DataFrame(total_expr_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(total_expr_list_s{}{}, columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "total_expr_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(7):
    for j in range(7):
        with open(savepath + "reg_prod_score.txt", "a") as f:
            exec('f.write("reg_prod_score{}{} = %f"%(reg_prod_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(7):
    for j in range(7):
        exec('x = pd.DataFrame(reg_prod_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(reg_prod_list_s{}{}, columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "reg_prod_list{}{}.csv", header=None, index=None)'.format(i, j))

# 将细胞表达法的打分结果写入文件
for i in range(7):
    for j in range(7):
        with open(savepath + "cell_expr_score.txt", "a") as f:
            exec('f.write("cell_expr_score{}{} = %f"%(cell_expr_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(7):
    for j in range(7):
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

for i in range(7):
    for j in range(7):
        exec('mult.append(mult_score{}{})'.format(i, j))
        exec('thrd.append(thrd_score{}{})'.format(i, j))
        exec('spec.append(spec_score{}{})'.format(i, j))
        exec('total_exp.append(total_expr_score{}{})'.format(i, j))
        exec('reg_prod.append(reg_prod_score{}{})'.format(i, j))
        exec('cell_expr.append(cell_expr_score{}{})'.format(i, j))

mult = np.array(mult).reshape(7, 7)
thrd = np.array(thrd).reshape(7, 7)
spec = np.array(spec).reshape(7, 7)
total_exp = np.array(total_exp).reshape(7, 7)
reg_prod = np.array(reg_prod).reshape(7, 7)
cell_exp = np.array(cell_expr).reshape(7, 7)


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

# mult = mm.fit_transform(mult)
# thrd = mm.fit_transform(thrd)
# spec = mm.fit_transform(spec)
#
# da = []
# for i in range(9):
#     for j in range(9):
#         a = mult[i][j]
#         b = thrd[i][j]
#         c = spec[i][j]
#         if a > b:
#            if b > c:  # a>b，b>c
#                print("排列后：", a, b, c)
#                max = a
#                avg = b
#                min = c
#            elif a > c:  # a>b，b<c，a>c
#                print("排列后：", a, c, b)
#                max = a
#                avg = c
#                min = b
#            else:  # a>b，a<c
#                print("排列后：", c, a, b)
#                max = c
#                avg = a
#                min = b
#         elif a > c:  # a<b，a>c
#            print("排列后：", b, a, c)
#            max = b
#            avg = a
#            min = c
#         elif b > c:  # a<b，a<c，b>c
#            print("排列后：", b, c, a)
#            max = b
#            avg = c
#            min = a
#         else:  # a<b，a<c，b<c
#            print("排列后：", c, b, a)
#            max = c
#            avg = b
#            min = a
#         result = avg/2 + (max + min)/4
#         da.append(result)
# da = np.array(da).reshape(9, 9)
# da = pd.DataFrame(da)
# da.to_csv(savepath + "score_result.csv", header=None, index=None)

print(' Time: {}秒, mem: {}MB'.format((time() - itime), psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024))




