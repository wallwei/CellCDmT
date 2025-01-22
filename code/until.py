import numpy as np
import pandas as pd
from pandas import DataFrame
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve,recall_score,precision_score
from sklearn.metrics import f1_score
from sklearn.model_selection import KFold

def perform(y_proba, y_true):  # Performance evaluation
    ACC = accuracy_score(y_true, y_proba.argmax(axis=1))
    fpr, tpr, threshold = roc_curve(y_true, y_proba[:, 1])
    pre, rec, _ = precision_recall_curve(y_true, y_proba[:, 1])
    RECALL = recall_score(y_true, y_proba.argmax(axis=1))
    AUC = auc(fpr, tpr)
    AUPR = auc(rec, pre)
    f1 = f1_score(y_true, y_proba.argmax(axis=1))
    Pre = precision_score(y_true, y_proba.argmax(axis=1))
    #plt
    curve_1 = np.vstack([fpr, tpr])
    curve_1 = pd.DataFrame(curve_1.T)
    curve_1.to_csv('out_new/1' + 'c' + "cv3" + '_au' + str(AUC) + '.csv', header=None, index=None)

    curve_2 = np.vstack([rec, pre])
    curve_2 = pd.DataFrame(curve_2.T)
    curve_2.to_csv('out_new/1' + 'c' + "cv3" + '_aupr' + str(AUPR) + '.csv', header=None, index=None)
    print("ACC=",ACC)
    print("AUC=", AUC)
    print("AUPR=", AUPR)
    print("RECALL=", RECALL)
    print("f1=", f1)
    print("Pre=", Pre)
    return ACC, AUC, AUPR, RECALL, f1, Pre


def coordinate_splicing(length,index_all,cv):
    temp= []
    for i in range(length): temp.append(i)
    index_x =[]
    index_y =[]
    for i in index_all:
        for j in temp:
                index_x.append(i)
                index_y.append(j)
    index_x = np.reshape(index_x, (-1, 1))
    index_y = np.reshape(index_y, (-1, 1))
    if(cv==1):
        index_end = np.hstack([index_x,index_y])
    else:
        index_end = np.hstack([index_y, index_x])
    return index_end

def Cv1_2_cross_validation_index(cv,Related):
    if cv==1:
        max = Related.shape[0]

    else:
        max = Related.shape[1]

    temp = []
    for i in range(max): temp.append(i)
    kf = KFold(n_splits=5,shuffle=True)
    train_index_all=[]
    test_index_all=[]
    for train_index, test_index in kf.split(temp):
        train_index_all.append(train_index)
        test_index_all.append(test_index)


    return train_index_all,test_index_all

# def CV3_Generate_data(Related,Line_feature,Row_feature,dim):
#     print("Generate sample coordinates...")
#     Associated_Index = np.argwhere(Related == 1)
#     Not_Associated_Index = np.argwhere(Related == 0)
#     np.random.shuffle(Not_Associated_Index)
#     Not_Associated_Index= Not_Associated_Index[: Associated_Index.shape[0]]
#     Related_index = np.vstack((Associated_Index, Not_Associated_Index))
#     np.random.shuffle(Related_index)
#     #Associated_Lable
#     Associated_Lable=[]
#     for x, y in Related_index: Associated_Lable.append(Related[x, y])
#     Associated_Lable=np.array(Associated_Lable).astype(np.int32)
#     if (len(Associated_Lable) == sum(Associated_Lable) * 2):
#         print("Positive and negative sample balance successfully")
#     else:
#         print("Failed to balance positive and negative samples")


    # Associated_feature = Generate_features(Related_index, Line_feature, Row_feature, dim)
    # return Associated_feature, Associated_Lable

def CV3_Generate_data(Related, Line_feature, Row_feature, dim):
    print("Generate sample coordinates...")
    Associated_Index = np.argwhere(Related == 1)
    Not_Associated_Index = np.argwhere(Related == 0)

    np.random.shuffle(Not_Associated_Index)
    Not_Associated_Index = Not_Associated_Index[:Associated_Index.shape[0]]
    Related_index = np.vstack((Associated_Index, Not_Associated_Index))
    np.random.shuffle(Related_index)

    # Associated_Label
    Associated_Label = np.array([Related[x, y] for x, y in Related_index]).astype(np.int32)

    if len(Associated_Label) == sum(Associated_Label) * 2:
        print("Positive and negative sample balance successfully")
    else:
        print("Failed to balance positive and negative samples")

    Associated_feature = Generate_features(Related_index, Line_feature, Row_feature, dim)
    return Associated_feature, Associated_Label

def Dim_reduction(feature,dim):
    from sklearn.decomposition import PCA
    PCA = PCA(n_components=dim)
    matrix = PCA.fit_transform(feature)
    return matrix

def Sample_balance(Related,index):
    Lable = []
    for x, y in index: Lable.append(Related[x,y])
    Associated_index = np.zeros((1, 2))
    No_Associated_index = np.zeros((1, 2))
    for i in range(len(Lable)):
        if Lable[i] == 1:
            Associated_index = np.vstack([Associated_index, index[i]])
        else:
            No_Associated_index = np.vstack([No_Associated_index,index[i]])

    Associated_index= Associated_index[1:]
    No_Associated_index=No_Associated_index[1:]
    # print("Associated_index.shape:",Associated_index.shape)
    np.random.shuffle(No_Associated_index)
    No_Associated_index=No_Associated_index[:Associated_index.shape[0]]
    # print("No_Associated_index.shape:", No_Associated_index.shape)
    index = np.vstack((Associated_index, No_Associated_index)).astype(np.int32)
    # print("index.shape:", index.shape)
    np.random.shuffle(index)

    Lable = []
    for x, y in index: Lable.append(Related[x,y])
    if(len(Lable)==sum(Lable)*2):print("Positive and negative sample balance successfully")
    else:print("Failed to balance positive and negative samples")

    return index

# def Generate_features(Related_index,Line_feature,Row_feature,dim):
#     # feature
#     feature_x = np.zeros((1, dim))
#     feature_y = np.zeros((1, dim))
#     Line_index, Row_index = Related_index[:, 0], Related_index[:, 1]
#     for j in Row_index: feature_x = np.vstack([feature_x, Row_feature[j]])
#     for i in Line_index: feature_y = np.vstack([feature_y, Line_feature[i]])
#
#     feature = np.hstack([feature_y[1:], feature_x[1:]]).astype(np.float64)
#     print("Generate feature vectors:",feature.shape)
#     return feature

def Generate_features(Related_index, Line_feature, Row_feature, dim):
    print("Generate features...")
    Line_index, Row_index = Related_index[:, 0], Related_index[:, 1]

    # 使用列表推导式收集特征向量
    feature_x = np.array([Row_feature[j] for j in Row_index])
    feature_y = np.array([Line_feature[i] for i in Line_index])

    # 水平堆叠特征向量
    feature = np.hstack((feature_y, feature_x)).astype(np.float64)
    print("Generate feature vectors:", feature.shape)
    return feature


def Case_Study_Generate_data(Related,Line_feature,Row_feature,dim):
    print("Generate sample coordinates...")
    Not_Associated_Index = np.argwhere(Related == 0)
    feature = Generate_features(Not_Associated_Index, Line_feature, Row_feature, dim)

    return feature,Not_Associated_Index

