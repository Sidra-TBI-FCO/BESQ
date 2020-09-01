# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost import plot_tree
import shap
from sklearn.model_selection import train_test_split, KFold
from lifelines.utils import concordance_index
from sklearn.model_selection import GridSearchCV
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import itertools
import operator
import collections
import graphviz
import random
import matplotlib as mpl

# +
#Load the dataset
survival_df = pd.read_csv("../Data/Feature_Matrix.csv",header='infer')
survival_status, icr_cluster_status, ethnicity_status,real_time,age = [],[],[],[],[]
for i in range(len(survival_df["Status"])):
    if survival_df["Status"][i]=="Dead":
        real_time.append((1.0*survival_df["Time"][i])/30)
        survival_status.append(0)
        #age.append(120)
    else:
        real_time.append(-(1.0*survival_df["Time"][i])/30)
        survival_status.append(1)
        #age.append(random.randint(30,50))

for i in range(len(survival_df["Ethnicity"])):
    if survival_df["Ethnicity"][i]=="Black":
        ethnicity_status.append(0)
    else:
        ethnicity_status.append(1)

for i in range(len(survival_df["ICR_Cluster"])):
    if survival_df["ICR_Cluster"][i]=="ICR_Low":
        icr_cluster_status.append(0)
    elif survival_df["ICR_Cluster"][i]=="ICR_Medium":
        icr_cluster_status.append(1)
    elif survival_df["ICR_Cluster"][i]=="ICR_High":
        icr_cluster_status.append(2)
survival_df["Status"]=survival_status
survival_df["Ethnicity"]=ethnicity_status
survival_df["ICR_Cluster"]=icr_cluster_status
survival_df["Time"]=real_time
#survival_df["Age"]=age

#Run the sensitivity analysis for AA population
survival_df2 = survival_df[survival_df["Ethnicity"]==0]
survival_df2.shape


# -
def c_statistic_harrell(pred, labels):
    total = 0
    matches = 0
    for i in range(len(labels)):
        for j in range(len(labels)):
            if labels[j] > 0 and abs(labels[i]) > labels[j]:
                total += 1
                if pred[j] > pred[i]:
                    matches += 1
    return matches/total


# +
#Build XGB model for cross-validation 
# use validation set to choose # of trees
params = {
    "eta": [0.001,0.01,0.1],
    "max_depth": [2,4],
    "colsample_bytree": [0.1,0.3,0.5],
    "min_child_weight": [3,5,7],
    "gamma": [0.0,0.5,1],
    "nthread":[4]
}
keys, values = zip(*params.items())
experiments = [dict(zip(keys, v)) for v in itertools.product(*values)]
per_param_performance=[]
    
for param in experiments:
    param["objective"]="survival:cox"
    param["booster"]="gbtree"
    param["mean_cv_c_index_many_runs"]=[]
        
    for i in range(1,16):
        #Shuffle the survival dataset
        if (i==9):
            continue
        temp_df = survival_df2.sample(frac=1, random_state=i).reset_index(drop=True)
        temp_survival_df2 = temp_df
        X_survival = temp_survival_df2.iloc[:,2:temp_survival_df2.shape[1]+1]
        y_survival = temp_survival_df2["Time"]
        survival_status = temp_survival_df2["Status"]
    
         
        kf = KFold(n_splits=5, random_state=7)
        c_index_list = []
        for train_index, test_index in kf.split(X_survival):
    
            #Get train, test: X, Y and survival info
            X_train, X_test, y_train, y_test = X_survival.iloc[train_index,:],X_survival.iloc[test_index,:],y_survival.iloc[train_index],y_survival.iloc[test_index]
            survival_status_train, survival_status_test = survival_status.iloc[train_index], survival_status.iloc[test_index]
            
            xgb_train = xgb.DMatrix(X_train, label=y_train)
            xgb_test = xgb.DMatrix(X_test, label=y_test)
            
            #Build optimal XGboost model
            model = xgb.train(param, xgb_train, num_boost_round=200, early_stopping_rounds=25, evals = [(xgb_test, "test")], verbose_eval=100)
            best_dict = model.attributes()
            best_num_round = int(best_dict["best_iteration"])
            
            #Test the predictions on the test set
            y_pred = model.predict(xgb_test, ntree_limit=best_num_round)
            c_index_value = concordance_index(y_test,y_pred,event_observed=survival_status_test)
            c_index_list.append(c_index_value)
        
        param["mean_cv_c_index_many_runs"].append(sum(c_index_list)/len(c_index_list))
    per_param_performance.append(param)
    
param_performance_df = pd.DataFrame(per_param_performance)
print(param_performance_df)
# -

many_performance_df = pd.DataFrame.from_records(param_performance_df["mean_cv_c_index_many_runs"].tolist())
max_value = np.max(np.array(many_performance_df.median(axis=1)))
max_id = np.argmax(np.array(many_performance_df.median(axis=1)))
print(max_id,max_value)

string_list = []
for i in range(many_performance_df.shape[0]):
    param_string = "P"+str(i)
    string_list.append(param_string)
tmp = many_performance_df
tmp.columns=["Exp1","Exp2","Exp3","Exp4","Exp5","Exp6","Exp7","Exp8","Exp9","Exp10","Exp11","Exp12","Exp13","Exp14"]
tmp.reset_index(drop=True)
tmp2 = tmp.copy()
tmp2['Params']=string_list
tmp2.set_index('Params')

tmp3 = tmp.T
tmp3.columns=string_list
fig = plt.figure(1, figsize=(20, 7))
mpl.rc('xtick', labelsize=10) 
mpl.rc('ytick', labelsize=16)
tmp3.boxplot(rot=90)
plt.xlabel("Parameter Setting",fontsize=20)
plt.ylabel("Mean 5-CV C-index Value",fontsize=20)
plt.ylim(0.4,0.7)
plt.axvline(x=max_id+1,color="red")
plt.title("Comparison of Mean 5-CV performance (15 randomizations) for grid of hyper-parameters for AA", fontsize=18)
plt.savefig("../Results/Mean_CV_Grid_AA.pdf", bbox_inches='tight')


fig = plt.figure(1, figsize=(10, 7))
rev_max_id = 'P'+str(max_id)
plt.boxplot(tmp3[rev_max_id].tolist())
plt.xlabel("Parameter Setting "+rev_max_id,fontsize=18)
plt.ylabel("Mean 5-CV C-index Value",fontsize=18)
plt.ylim(0.4,0.7)
plt.title("Sensitivity of optimal XGBoost model for AA")
plt.savefig("../Results/Optimal_XGB_AA.pdf", bbox_inches='tight')
