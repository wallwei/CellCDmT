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

dt = pd.read_csv('GSE103322.csv', index_col=0, header=None)
#a = np.array(dt.loc["Cell"])
dict = {}
dict["Cell"] = np.array(dt.loc["Cell"])
for i in range(1, dt.shape[0]):
    dict[dt.index[i]] = np.array(dt.loc[dt.index[i]], dtype=float)

def sigmoid(x):
    return 1/(1+np.exp(-(x-6)))

savepath = 'GSE103322'

#0=m,1=Fibro,2=B,3=myo,4=Macro,5=Endo,8=T
malignant_index = np.where(dict['classified  as cancer cell'] == 1)[0]
Fibro_index = np.where(dict['non-cancer cell type'] == 1)[0]
B_index = np.where(dict['non-cancer cell type'] == 2)[0]
myo_index = np.where(dict['non-cancer cell type'] == 3)[0]
Macro_index = np.where(dict['non-cancer cell type'] == 4)[0]
Endo_index = np.where(dict['non-cancer cell type'] == 5)[0]
T_index = np.where(dict['non-cancer cell type'] == 6)[0]
Den_index = np.where(dict['non-cancer cell type'] == 7)[0]
Mast_index = np.where(dict['non-cancer cell type'] == 8)[0]
###
groupsize = []
groupsize.append(malignant_index.shape[0])
groupsize.append(Fibro_index.shape[0])
groupsize.append(B_index.shape[0])
groupsize.append(myo_index.shape[0])
groupsize.append(Macro_index.shape[0])
groupsize.append(Endo_index.shape[0])
groupsize.append(T_index .shape[0])
groupsize.append(Den_index.shape[0])
groupsize.append(Mast_index.shape[0])
groupsize = pd.DataFrame(groupsize)
groupsize.to_csv('GSE103322/HNSCC.csv', header=None, index=None)


for i in range(9):
    for j in range(9):
        exec('mult_score{}{} = 0'.format(i, j))

for i in range(9):
    for j in range(9):
        exec('thrd_score{}{} = 0'.format(i, j))

for i in range(9):
    for j in range(9):
        exec('mult_list{}{} = []'.format(i, j))

for i in range(9):
    for j in range(9):
        exec('mult_list_s{}{} = []'.format(i, j))

for i in range(9):
    for j in range(9):
        exec('thrd_list{}{} = []'.format(i, j))

for i in range(9):
    for j in range(9):
        exec('thrd_list_s{}{} = []'.format(i, j))

for i in range(9):
    for j in range(9):
        exec('spec_score{}{} = 0'.format(i, j))

for i in range(9):
    for j in range(9):
        exec('spec_list{}{} = []'.format(i, j))

for i in range(9):
    for j in range(9):
        exec('spec_list_s{}{} = []'.format(i, j))
g = 0

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
        myo_l = 1 / myo_index.shape[0] * sum(dict[i[0]][myo_index])
        myo_r = 1 / myo_index.shape[0] * sum(dict[i[1]][myo_index])
        Macro_l = 1 / Macro_index.shape[0] * sum(dict[i[0]][Macro_index])
        Macro_r = 1 / Macro_index.shape[0] * sum(dict[i[1]][Macro_index])
        Endo_l = 1 / Endo_index.shape[0] * sum(dict[i[0]][Endo_index])
        Endo_r = 1 / Endo_index.shape[0] * sum(dict[i[1]][Endo_index])
        T_l = 1 / T_index.shape[0] * sum(dict[i[0]][T_index])
        T_r = 1 / T_index.shape[0] * sum(dict[i[1]][T_index])
        Den_l = 1 / Den_index.shape[0] * sum(dict[i[0]][Den_index])
        Den_r = 1 / Den_index.shape[0] * sum(dict[i[1]][Den_index])
        Mast_l = 1 / Mast_index.shape[0] * sum(dict[i[0]][Mast_index])
        Mast_r = 1 / Mast_index.shape[0] * sum(dict[i[1]][Mast_index])
        l_list = [malignant_l, Fibro_l, B_l, myo_l, Macro_l, Endo_l, T_l, Den_l, Mast_l]
        r_list = [malignant_r, Fibro_r, B_r, myo_r, Macro_r, Endo_r, T_r, Den_r, Mast_r]

        a = b = 0
        for item in product(l_list, r_list):                    #product(A, B) 和 ((x,y) for x in A for y in B)一样
            # print("sigmoid:%f"%sigmoid(item[0]*item[1]))


            exec('mult_score{}{} += {}'.format(a, b, (item[0] * item[1])))
            exec('mult_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('mult_list_s{}{}.append({})'.format(a, b, (item[0] * item[1])))
            b += 1
            if b == 9:
                b = 0
                a += 1

        mean_l_malignant = np.mean(dict[i[0]][malignant_index])
        mean_l_Fibro = np.mean(dict[i[0]][Fibro_index])
        mean_l_B = np.mean(dict[i[0]][B_index])
        mean_l_myo = np.mean(dict[i[0]][myo_index])
        mean_l_Macro = np.mean(dict[i[0]][Macro_index])
        mean_l_Endo = np.mean(dict[i[0]][Endo_index])
        mean_l_T = np.mean(dict[i[0]][T_index])
        mean_l_Den = np.mean(dict[i[0]][Den_index])
        mean_l_Mast = np.mean(dict[i[0]][Mast_index])
        mean_l = np.mean((mean_l_malignant, mean_l_Fibro, mean_l_B, mean_l_myo, mean_l_Macro, mean_l_Endo, mean_l_T,
                          mean_l_Den, mean_l_Mast))
        std_l = np.std(dict[i[0]][np.concatenate((malignant_index, Fibro_index, B_index, myo_index, Macro_index,
                                                  Endo_index, T_index, Den_index, Mast_index))])
        sum_l = np.sum((mean_l_malignant, mean_l_Fibro, mean_l_B, mean_l_myo, mean_l_Macro, mean_l_Endo, mean_l_T, mean_l_Den, mean_l_Mast))

        mean_r_malignant = np.mean(dict[i[1]][malignant_index])
        mean_r_Fibro = np.mean(dict[i[1]][Fibro_index])
        mean_r_B = np.mean(dict[i[1]][B_index])
        mean_r_myo = np.mean(dict[i[1]][myo_index])
        mean_r_Macro = np.mean(dict[i[1]][Macro_index])
        mean_r_Endo = np.mean(dict[i[1]][Endo_index])
        mean_r_T = np.mean(dict[i[1]][T_index])
        mean_r_Den = np.mean(dict[i[0]][Den_index])
        mean_r_Mast = np.mean(dict[i[0]][Mast_index])
        mean_r = np.mean((mean_r_malignant, mean_r_Fibro, mean_r_B, mean_r_myo, mean_r_Macro, mean_r_Endo, mean_r_T,
                          mean_r_Den, mean_r_Mast))
        std_r = np.std(dict[i[1]][np.concatenate((malignant_index, Fibro_index, B_index, myo_index, Macro_index,
                                                  Endo_index, T_index, Den_index, Mast_index))])
        sum_r = np.sum((mean_r_malignant, mean_r_Fibro, mean_r_B, mean_r_myo, mean_r_Macro, mean_r_Endo, mean_r_T, mean_r_Den, mean_r_Mast))

        malignant_l = int(mean_l_malignant > mean_l + std_l)
        malignant_r = int(mean_r_malignant > mean_r + std_r)
        Fibro_l = int(mean_l_Fibro > mean_l + std_l)
        Fibro_r = int(mean_r_Fibro > mean_r + std_r)
        B_l = int(mean_l_B > mean_l + std_l)
        B_r = int(mean_r_B > mean_r + std_r)
        myo_l = int(mean_l_myo > mean_l + std_l)
        myo_r = int(mean_r_myo > mean_r + std_r)
        Macro_l = int(mean_l_Macro > mean_l + std_l)
        Macro_r = int(mean_r_Macro > mean_r + std_r)
        Endo_l = int(mean_l_Endo > mean_l + std_l)
        Endo_r = int(mean_r_Endo > mean_r + std_r)
        T_l = int(mean_l_T > mean_l + std_l)
        T_r = int(mean_r_T > mean_r + std_r)
        Den_l = int(mean_l_Den > mean_l + std_l)
        Den_r = int(mean_r_Den > mean_r + std_r)
        Mast_l = int(mean_l_Mast > mean_l + std_l)
        Mast_r = int(mean_r_Mast > mean_r + std_r)
        l_list = [malignant_l, Fibro_l, B_l, myo_l, Macro_l, Endo_l, T_l, Den_l, Mast_l]
        r_list = [malignant_r, Fibro_r, B_r, myo_r, Macro_r, Endo_r, T_r, Den_r, Mast_r]
        a = b = 0
        for item in product(l_list, r_list):
            exec('thrd_score{}{} += {}'.format(a, b, int(item[0] & item[1])))
            exec('thrd_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('thrd_list_s{}{}.append({})'.format(a, b, int(item[0] & item[1])))
            b += 1
            if b == 9:
                b = 0
                a += 1
         #
        sp_l_malignant = mean_l_malignant/sum_l
        sp_l_Fibro = mean_l_Fibro/sum_l
        sp_l_B = mean_l_B/sum_l
        sp_l_myo = mean_l_myo/sum_l
        sp_l_Macro = mean_l_Macro/sum_l
        sp_l_Endo = mean_l_Endo/sum_l
        sp_l_T = mean_l_T/sum_l
        sp_l_Den = mean_l_Den / sum_l
        sp_l_Mast = mean_l_Mast / sum_l
        sp_l_list = [sp_l_malignant, sp_l_Fibro, sp_l_B,  sp_l_myo, sp_l_Macro, sp_l_Endo, sp_l_T, sp_l_Den, sp_l_Mast]

        sp_r_malignant = mean_r_malignant / sum_r
        sp_r_Fibro = mean_r_Fibro / sum_r
        sp_r_B = mean_r_B / sum_r
        sp_r_myo = mean_r_myo / sum_r
        sp_r_Macro = mean_r_Macro / sum_r
        sp_r_Endo = mean_r_Endo / sum_r
        sp_r_T = mean_r_T / sum_r
        sp_r_Den = mean_r_Den / sum_r
        sp_r_Mast = mean_r_Mast / sum_r
        sp_r_list = [sp_r_malignant, sp_r_Fibro, sp_r_B,  sp_r_myo, sp_r_Macro, sp_r_Endo, sp_r_T, sp_r_Den,
                     sp_r_Mast]

        a = b = 0
        for item in product(sp_l_list, sp_r_list):                    #product(A, B) 和 ((x,y) for x in A for y in B)一样
            # print("sigmoid:%f"%sigmoid(item[0]*item[1]))

            exec('spec_score{}{} += {}'.format(a, b, (item[0] * item[1])))
            exec('spec_list{}{}.append("{}" + "-" + "{}")'.format(a, b, i[0], i[1]))
            exec('spec_list_s{}{}.append({})'.format(a, b, (item[0] * item[1])))
            b += 1
            if b == 9:
                b = 0
                a += 1

for i in range(9):
    for j in range(9):
        with open(savepath + "mult_score.txt", "a") as f:
            exec('f.write("mult_score{}{} = %f"%(mult_score{}{}))'.format(i, j, i, j))
            f.write('\n')
for i in range(9):
    for j in range(9):
        with open(savepath + "thrd_score.txt","a") as f:
            exec('f.write("thrd_score{}{} = %f"%(thrd_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(9):
    for j in range(9):
        exec('x = pd.DataFrame(mult_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(mult_list_s{}{},columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "mult_list{}{}.csv", header=None, index=None)'.format(i, j))

for i in range(9):
    for j in range(9):
        exec('x = pd.DataFrame(thrd_list{}{})'.format(i,j))
        exec('y = pd.DataFrame(thrd_list_s{}{},columns=list("3"))'.format(i,j))
        x = x.join(y)
        exec('x.to_csv(savepath + "thrd_list{}{}.csv", header=None, index=None)'.format(i, j))
for i in range(9):
    for j in range(9):
        with open(savepath + "spec_score.txt", "a") as f:
            exec('f.write("spec_score{}{} = %f"%(spec_score{}{}))'.format(i, j, i, j))
            f.write('\n')

for i in range(9):
    for j in range(9):
        exec('x = pd.DataFrame(spec_list{}{})'.format(i, j))
        exec('y = pd.DataFrame(spec_list_s{}{},columns=list("3"))'.format(i, j))
        x = x.join(y)
        exec('x.to_csv(savepath + "spec_list{}{}.csv", header=None, index=None)'.format(i, j))


import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
mm = MinMaxScaler()
mult = []
thrd = []
spec = []
for i in range(9):
    for j in range(9):
        exec('mult.append(mult_score{}{})'.format(i, j))
        exec('thrd.append(thrd_score{}{})'.format(i, j))
        exec('spec.append(spec_score{}{})'.format(i, j))
mult = np.array(mult).reshape(9, 9)
thrd = np.array(thrd).reshape(9, 9)
spec = np.array(spec).reshape(9, 9)

multy = pd.DataFrame(mult)
thrdy = pd.DataFrame(thrd)
specy = pd.DataFrame(spec)

multy.to_csv(savepath + "multy.csv", header=None, index=None)
thrdy.to_csv(savepath + "thrdy.csv", header=None, index=None)
specy.to_csv(savepath + "specy.csv", header=None, index=None)

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




