import xgboost as xgb
from until import *
# path = os.path.abspath(os.path.join(os.getcwd(), "../"))
# Related = pd.read_csv(path+"/dataset/human/Related.csv", header=None).to_numpy()
# Line_feature= pd.read_csv(path+"/dataset/human/receptor_feature_CTD+pyfet.csv", header=None).to_numpy()
# Row_feature= pd.read_csv(path+"/dataset/human/ligand_feature_CTD+pyfet.csv", header=None).to_numpy()

dataname = 'mitab_4'
path= "../../dataset/"

Related = pd.read_csv(path+dataname+"/Related.csv", header=None).to_numpy()
Line_feature1= pd.read_csv(path+dataname+"/receptor_feature_CTD+pyfet.csv", header=None).to_numpy()
Line_feature2= pd.read_csv(path+dataname+"/receptor_PseudoAAC_feature.csv", header=None).to_numpy()
Row_feature1= pd.read_csv(path+dataname+"/ligand_feature_CTD+pyfet.csv", header=None).to_numpy()
Row_feature2 = pd.read_csv(path+dataname+"/ligand_PseudoAAC_feature.csv", header=None).to_numpy()


Line_feature=np.hstack([Line_feature1,Line_feature2])
Row_feature=np.hstack([Row_feature1,Row_feature2])

X_train, y_train = CV3_Generate_data(Related,Row_feature,Line_feature,377)

xgb_model=xgb.XGBClassifier()
xgbout=xgb_model.fit(X_train, y_train)
feature_importance=xgbout.feature_importances_
feature_number=-feature_importance
H1=np.argsort(feature_number)
mask=H1[:450]
X_train=X_train[:,mask]

print("After feature selection:",X_train.shape)

X_test,Associated_Index = Case_Study_Generate_data(Related,Row_feature,Line_feature,377)
X_test = X_test[:,mask]

cup = pd.DataFrame(Associated_Index)
cup.to_csv('../CD_case_study/out/mitab/pred_Associated_Index.csv', header=None, index=None)


from catboost import CatBoostClassifier
cat = CatBoostClassifier(boosting_type='Ordered',max_depth =8,n_estimators=3000, learning_rate = 0.01)
cat.fit(X_train, y_train)
cat_y_preds = cat.predict_proba(X_test)

from deepforest import CascadeForestClassifier
model = CascadeForestClassifier(n_estimators=3, n_trees=160)
model.fit(X_train, y_train)
cf_y_preds = model.predict_proba(X_test)

y_preds = (cat_y_preds*0.1+cf_y_preds*0.9)

cu = pd.DataFrame(y_preds)
cu.to_csv('out/mitab/Prediction_Associated_all.csv', header=None, index=None)
y_preds_out = y_preds[:, 1]
y_preds_out=np.array(y_preds_out)


Index = np.argwhere(y_preds_out > 0.5).reshape(1,-1)[0]
Associated = np.zeros((1, 2))

for i in Index:
    Associated = np.vstack([Associated, Associated_Index[i]])
Associated=Associated[1:]

curve = pd.DataFrame(Associated)
curve.to_csv('out/mitab/Associated_Index_out_pre_ture.csv', header=None, index=None)



# import yagmail
# import os
# print('运行完成')
# # 登录你的邮箱
# yag = yagmail.SMTP(user='451262502@qq.com', password='qdmqulccnkglbgig', host='smtp.qq.com')
# # 发送邮件
# yag.send(to=['867808507@qq.com'], subject='程序运行', contents=['--运行完毕'])








