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
import shap
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import itertools
import operator
import collections
import graphviz

# +
#Load the dataset
survival_df = pd.read_csv("../Data/Feature_Matrix.csv",header='infer')
survival_status, icr_cluster_status, ethnicity_status,real_time = [],[],[],[]
for i in range(len(survival_df["Status"])):
    if survival_df["Status"][i]=="Dead":
        real_time.append(1.0*survival_df["Time"][i])
        survival_status.append(0)
    else:
        real_time.append(-1.0*survival_df["Time"][i])
        survival_status.append(1)

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

survival_df2 = survival_df[survival_df["Ethnicity"]==0]
X_black_survival = survival_df2.iloc[:,2:survival_df2.shape[1]+1]
y_black_survival = survival_df2["Time"]
print(X_black_survival.shape)

survival_df3 = survival_df[survival_df["Ethnicity"]==1]
X_white_survival = survival_df3.iloc[:,2:survival_df3.shape[1]+1]
y_white_survival = survival_df3["Time"]
print(X_white_survival.shape)

# +
#Load the saved optimal xgboost model
model_black = xgb.Booster()
model_black.load_model("../Models/Optimal_XGB_survival_black_model.m")
print(model_black)

model_white = xgb.Booster()
model_white.load_model("../Models/Optimal_XGB_survival_white_model.m")
print(model_white)

# +
#SHAP model to get importance for AA
shap_values = shap.TreeExplainer(model_black).shap_values(X_black_survival)

#Top pathways for AA people to determine survival
black_imp_df = pd.read_csv("../Results/black.csv",header='infer',sep=",")
flist = black_imp_df["Feature"].tolist()
#        "IPA:Myc_Mediated_Apoptosis_Signaling","IPA:EGF_Signaling","HM:Oxidative_phosphorylation",
#        "TPW:Immunogenic_Cell_Death_(ICD)","ICRscore","IPA:UVB_Induced_MAPK_Signaling","IPA:UVA_Induced_MAPK_Signaling"]

fids = [X_black_survival.columns.get_loc(c) for c in flist]
shap.summary_plot(shap_values[:,fids],X_black_survival.iloc[:,fids],sort=False)

# +
#SHAP model to get importance for black
shap_values = shap.TreeExplainer(model_white).shap_values(X_white_survival)

#Top pathways for white people to determine survival
white_imp_df = pd.read_csv("../Results/white.csv",header='infer',sep=",")
flist = white_imp_df["Features"][0:20]

#flist = ["IPA:Telomere_Extension_by_Telomerase","HM:PI3K_Akt_mTOR_signaling","LM:Proliferation",
#        "HM:Wnt_beta_catenin_signaling","TBI:Barrier_genes","IPA:AMPK_Signaling","IPA:PI3K_AKT_Signaling",
#        "HM:Angiogenesis","IPA:ErbB_Signaling","IPA:ERK5_Signaling","HM:G2M_checkpoint",
#        "HM:p53_pathway","HM:UV_response_down","IPA:UVC_Induced_MAPK_Signaling","IPA:HER_2_Signaling_in_Breast_Cancer",
#        "HM:Reactive_oxigen_species_pathway","IPA:VEGF_Signaling","IPA:Estrogen_Dependent_Breast_Cancer_Signaling",
#        "TBI:Phopholipase","IPA:Myc_Mediated_Apoptosis_Signaling","IPA:EGF_Signaling","HM:Oxidative_phosphorylation",
#        "TPW:Immunogenic_Cell_Death_(ICD)","ICRscore","IPA:UVB_Induced_MAPK_Signaling","IPA:UVA_Induced_MAPK_Signaling"]

fids = [X_white_survival.columns.get_loc(c) for c in flist]
shap.summary_plot(shap_values[:,fids],X_white_survival.iloc[:,fids],sort = False)
# -
shap.dependence_plot("IPA:EGF_Signaling", shap_values, X_white_survival, interaction_index="HM:DNA_repair")


