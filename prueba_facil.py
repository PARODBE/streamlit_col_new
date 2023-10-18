# from turtle import onclick
import streamlit as st
import numpy as np
import pandas as pd
from streamlit_ketcher import st_ketcher
from PIL import Image
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.neighbors import KNeighborsRegressor
import random
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
from xgboost import XGBRegressor
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score,cross_validate, KFold
from sklearn.model_selection import RepeatedKFold
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.metrics import make_scorer, mean_absolute_error
import time
import io
import os

# os.environ["R_HOME"] = 'C:/Program Files/R/R-4.3.1' 
# os.environ["PATH"] = 'C:/Program Files/R/R-4.3.1/bin/x64' + ";" + os.environ["PATH"] 

import rpy2.rinterface
from rpy2.robjects import pandas2ri, r
import rpy2.robjects as ro
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

pandas2ri.activate()

# os.environ['R_HOME']="C:\R"



st.image('juntas.png', use_column_width=False,width=300)


# molfile = st_ketcher(molecule_format="MOLFILE")
# st.markdown("molfile:")
# st.code(molfile)


def fingerprints_inputs(dataframe):

    X=np.array([AllChem.GetMorganFingerprintAsBitVect(mol,radius=2,nBits=2048,useFeatures=True) for mol in [Chem.MolFromSmiles(m) for m in list(dataframe.canonical_smiles)]])
    y=dataframe.pchembl_value.astype('float')
    return X,y

# Add a file uploader after Button 3 is clicked
file_uploaded = st.file_uploader("Select one file", type=["csv", "txt", "xlsx", "sdf"])

if file_uploaded is not None:  # Verifica si se ha cargado un archivo
    contenido = file_uploaded.read()
    
    if file_uploaded.type == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet":
        df_uploaded = pd.read_excel(io.BytesIO(contenido), index_col=0, engine='openpyxl')
    else:
        st.error("Unsupported file format. Please upload an Excel file.")
else:
    st.warning("Please upload a file before attempting to read it.")


# # Check and update the state of Button 1
# if st.button("Press to load transporter datasets"):
    
st.info("Loading transporter datasets...")

df_bcrp=PandasTools.LoadSDF(r"C:\Users\parod\OneDrive\Escritorio\Doctorado\etransafe_aop\Datasets\chembl_data_bcrp.sdf")
df_bcrp['inchi']=[Chem.MolToInchi(m) for m in df_bcrp.ROMol]
df_bcrp=df_bcrp[df_bcrp.standard_type=='IC50']
df_bcrp['pchembl_value']=df_bcrp['pchembl_value'].astype('float')
df_mrp2=PandasTools.LoadSDF(r"C:/Users/parod/OneDrive/Escritorio/Doctorado/etransafe_aop/Datasets/chembl_data_mrp2.sdf")
df_mrp2['inchi']=[Chem.MolToInchi(m) for m in df_mrp2.ROMol]
df_mrp2=df_mrp2[df_mrp2.standard_type=='IC50']
df_mrp2['pchembl_value']=df_mrp2['pchembl_value'].astype('float')
df_mrp3=PandasTools.LoadSDF(r"C:/Users/parod/OneDrive/Escritorio/Doctorado/etransafe_aop/Datasets/chembl_data_mrp3.sdf")
df_mrp3['inchi']=[Chem.MolToInchi(m) for m in df_mrp3.ROMol]
df_mrp3=df_mrp3[df_mrp3.standard_type=='IC50']
df_mrp3['pchembl_value']=df_mrp3['pchembl_value'].astype('float')
df_mrp4=PandasTools.LoadSDF(r"C:/Users/parod/OneDrive/Escritorio/Doctorado/etransafe_aop/Datasets/chembl_data_mrp4.sdf")
df_mrp4['inchi']=[Chem.MolToInchi(m) for m in df_mrp4.ROMol]
df_mrp4=df_mrp4[df_mrp4.standard_type=='IC50']
df_mrp4['pchembl_value']=df_mrp4['pchembl_value'].astype('float')
df_oat1=PandasTools.LoadSDF(r"C:/Users/parod/OneDrive/Escritorio/Doctorado/etransafe_aop/Datasets/chembl_data_OATP1b1.sdf")
df_oat1['inchi']=[Chem.MolToInchi(m) for m in df_oat1.ROMol]
df_oat1=df_oat1[df_oat1.standard_type=='IC50']
df_oat1['pchembl_value']=df_oat1['pchembl_value'].astype('float')
df_oat2=PandasTools.LoadSDF(r"C:/Users/parod/OneDrive/Escritorio/Doctorado/etransafe_aop/Datasets/chembl_data_OATP1b3.sdf")
df_oat2['inchi']=[Chem.MolToInchi(m) for m in df_oat2.ROMol]
df_oat2=df_oat2[df_oat2.standard_type=='IC50']
df_oat2['pchembl_value']=df_oat2['pchembl_value'].astype('float')
df_bsep=PandasTools.LoadSDF(r"C:/Users/parod/OneDrive/Escritorio/Doctorado/etransafe_aop/Datasets/chembl_data_bsep.sdf")
df_bsep['inchi']=[Chem.MolToInchi(m) for m in df_bsep.ROMol]
df_bsep=df_bsep[df_bsep.standard_type=='IC50']
df_bsep['pchembl_value']=df_bsep['pchembl_value'].astype('float')
df_pgp=PandasTools.LoadSDF(r"C:/Users/parod/OneDrive/Escritorio/Doctorado/etransafe_aop/Datasets/chembl_data_pgp.sdf")
df_pgp['inchi']=[Chem.MolToInchi(m) for m in df_pgp.ROMol]
df_pgp=df_pgp[df_pgp.standard_type=='IC50']
df_pgp['pchembl_value']=df_pgp['pchembl_value'].astype('float')

st.write('<p style="color:green; font-size:24px;">&#10003; Low-Level Models built successfully</p>', unsafe_allow_html=True)

st.info("Building Low-Level Models...")

X_bcrp,y_bcrp=fingerprints_inputs(df_bcrp)
X_mrp2,y_mrp2=fingerprints_inputs(df_mrp2)
X_mrp3,y_mrp3=fingerprints_inputs(df_mrp3)
X_mrp4,y_mrp4=fingerprints_inputs(df_mrp4)
X_oat1,y_oat1=fingerprints_inputs(df_oat1)
X_oat2,y_oat2=fingerprints_inputs(df_oat2)
X_bsep,y_bsep=fingerprints_inputs(df_bsep)
X_pgp,y_pgp=fingerprints_inputs(df_pgp)

random.seed(46)

model_bcrp=RandomForestRegressor(**{'criterion': 'squared_error', 'max_depth': None, 'min_samples_split': 2, 'n_estimators': 16},random_state=46).fit(X_bcrp,y_bcrp)

model_mrp2=RandomForestRegressor(**{'criterion': 'squared_error', 'max_depth': None, 'min_samples_split': 3, 'n_estimators': 16},random_state=46).fit(X_mrp2,y_mrp2)

model_mrp3=SVR(**{'C': 1, 'gamma': 0.001, 'kernel': 'linear'}).fit(X_mrp3,y_mrp3)

model_mrp4=RandomForestRegressor(**{'criterion': 'squared_error', 'max_depth': None, 'min_samples_split': 3, 'n_estimators': 80},random_state=46).fit(X_mrp4,y_mrp4)

model_oat1=RandomForestRegressor(**{'criterion': 'squared_error', 'max_depth': 2, 'min_samples_split': 5, 'n_estimators': 16},random_state=46).fit(X_oat1,y_oat1)

model_oat2=SVR(**{'C': 1, 'gamma': 0.0001, 'kernel': 'rbf'}).fit(X_oat2,y_oat2)

model_bsep=XGBRegressor(**{'colsample_bytree': 1, 'max_depth': 2, 'min_child_weight': 2, 'n_estimators': 16}).fit(X_bsep,y_bsep)

model_pgp=RandomForestRegressor(**{'criterion': 'squared_error', 'max_depth': 2, 'min_samples_split': 2, 'n_estimators': 16},random_state=46).fit(X_pgp,y_pgp)

st.write('<p style="color:green; font-size:24px;">&#10003; Low-Level Models built successfully</p>', unsafe_allow_html=True)


def fingerprints_inputs2(dataframe):
    X=np.array([AllChem.GetMorganFingerprintAsBitVect(mol,radius=2,nBits=2048,useFeatures=True) for mol in [Chem.MolFromSmiles(m) for m in list(dataframe.Smiles)]])
    y=dataframe.Activity.astype('int')
    return X,y

st.info("Calculating in vitro concentrations...")

X_test, y_test = fingerprints_inputs2(df_uploaded)

pred_bcrp = model_bcrp.predict(X_test)
pred_mrp2 = model_mrp2.predict(X_test)
pred_mrp3 = model_mrp3.predict(X_test)
pred_mrp4 = model_mrp4.predict(X_test)
pred_oat1 = model_oat1.predict(X_test)
pred_oat2 = model_oat2.predict(X_test)
pred_pgp = model_pgp.predict(X_test)
pred_bsep = model_bsep.predict(X_test)


df_in_vitro=pd.DataFrame(data={'BCRP':pred_bcrp,'MRP2':pred_mrp2,'MRP3':pred_mrp3,'MRP4':pred_mrp4,'OATP1B1':pred_oat1,'OATP1B3':pred_oat2,'BSEP':pred_bsep,
                                'PGP':pred_pgp,'Activity':y_test.values},index=df_uploaded.index)

df_nuevo_def3=pd.concat([df_in_vitro.iloc[:,:-1],df_uploaded[['FUB','CLint','Doses max','Activity']]],axis=1)

smi=df_uploaded.Smiles.tolist()

mwt=[rdkit.Chem.Descriptors.ExactMolWt(Chem.MolFromSmiles(m)) for m in smi]

df_nuevo_def3.insert(loc = 12,
        column = 'MW',
        value = mwt)

logp=[rdkit.Chem.Crippen.MolLogP(Chem.MolFromSmiles(m)) for m in smi]
df_nuevo_def3.insert(loc = 13,
        column = 'LOGP',
        value = logp)

df_nuevo_def3['CAS_fake']=[f'0-0-0-{i}' for i in range(10)]
df_nuevo_def3['DTXSID_fake']=[f'DTXSID-{i}' for i in range(10)]
df_nuevo_def3['fake_name']=[f'A-{i}' for i in range(10)]

# st.write (df_in_vitro.index)

st.write('<p style="color:green; font-size:24px;">&#10003; Transporter in vitro concentrations predicted successfully</p>', unsafe_allow_html=True)

st.info("Calculating in vivo doses with httk library...")
# Configura la conversión entre pandas y R dataframes
with (ro.default_converter + pandas2ri.converter).context():
    pandas_df_2_r = ro.conversion.py2rpy(df_nuevo_def3)

    r.assign("df_r", pandas_df_2_r)

    script_r = r('''
                    library(httk)
                    library(dplyr)
                    
                    my.new.data <- as.data.frame(df_r$fake_name,stringsAsFactors=FALSE)
                    my.new.data <- cbind(my.new.data,as.data.frame(df_r$CAS_fake,stringsAsFactors=FALSE))
                    my.new.data <- cbind(my.new.data,as.data.frame(c(df_r$DTXSID_fake),stringsAsFactors=FALSE))
                    my.new.data <- cbind(my.new.data,as.data.frame(c(df_r$MW)))
                    my.new.data <- cbind(my.new.data,as.data.frame(c(df_r$LOGP)))
                    my.new.data <- cbind(my.new.data,as.data.frame(c(df_r$FUB)))
                    my.new.data <- cbind(my.new.data,as.data.frame(c(df_r$CLint)))
                
                    colnames(my.new.data) <- c("Name","CASRN","DTXSID","MW","LogP","Fup","CLint")
                    chem.physical_and_invitro.data <- add_chemtable(my.new.data,
                    current.table=
                    chem.physical_and_invitro.data,
                    data.list=list(Compound="Name",CAS="CASRN",DTXSID="DTXSID",MW="MW",logP="LogP",
                    Funbound.plasma="Fup",Clint="CLint"),overwrite=TRUE,species="Human",reference="httk|chembl|medscape|NIH")
                    
                    set.seed(42)
                    idx_list <- list(10^(-df_r$BCRP)*10^6, df_r$DTXSID_fake)
                    idx_vect <- c(1:length(idx_list[[1]]))
                    answer_bcrp<-data.frame()
                    equiv_dose<-for(i in idx_vect){
                        x <- idx_list[[1]][i]
                        j <- idx_list[[2]][i]
                        output = calc_mc_oral_equiv(conc = x,dtxsid=j,which.quantile = c(0.9),model="1compartment")
                        answer_bcrp<-rbind(answer_bcrp,output)
                    }
                    colnames(answer_bcrp) <- c("0.9_oral_dose_bcrp")
                    answer_bcrp$DTXSID<-df_r$DTXSID_fake

                    set.seed(42)
                    idx_list <- list(10^(-df_r$MRP2)*10^6, df_r$DTXSID_fake)
                    idx_vect <- c(1:length(idx_list[[1]]))
                    answer_mrp2<-data.frame()
                    equiv_dose<-for(i in idx_vect){
                        x <- idx_list[[1]][i]
                        j <- idx_list[[2]][i]
                        output = calc_mc_oral_equiv(conc = x,dtxsid=j,which.quantile = c(0.9),model="1compartment")
                        answer_mrp2<-rbind(answer_mrp2,output)
                    }
                    colnames(answer_mrp2) <- c("0.9_oral_dose_mrp2")
                    answer_mrp2$DTXSID<-df_r$DTXSID_fake
                    
                    set.seed(42)
                    idx_list <- list(10^(-df_r$MRP3)*10^6, df_r$DTXSID_fake)
                    idx_vect <- c(1:length(idx_list[[1]]))
                    answer_mrp3<-data.frame()
                    equiv_dose<-for(i in idx_vect){
                        x <- idx_list[[1]][i]
                        j <- idx_list[[2]][i]
                        output = calc_mc_oral_equiv(conc = x,dtxsid=j,which.quantile = c(0.9),model="1compartment")
                        answer_mrp3<-rbind(answer_mrp3,output)
                    }
                    colnames(answer_mrp3) <- c("0.9_oral_dose_mrp3")
                    answer_mrp3$DTXSID<-df_r$DTXSID_fake
                    
                    set.seed(42)
                    idx_list <- list(10^(-df_r$MRP4)*10^6, df_r$DTXSID_fake)
                    idx_vect <- c(1:length(idx_list[[1]]))
                answer_mrp4<-data.frame()
                    equiv_dose<-for(i in idx_vect){
                        x <- idx_list[[1]][i]
                        j <- idx_list[[2]][i]
                        output = calc_mc_oral_equiv(conc = x,dtxsid=j,which.quantile = c(0.9),model="1compartment")
                        answer_mrp4<-rbind(answer_mrp4,output)
                    }
                    colnames(answer_mrp4) <- c("0.9_oral_dose_mrp4")
                    answer_mrp4$DTXSID<-df_r$DTXSID_fake

                    set.seed(42)
                    idx_list <- list(10^(-df_r$OATP1B1)*10^6, df_r$DTXSID_fake)
                    idx_vect <- c(1:length(idx_list[[1]]))
                    answer_oat1<-data.frame()
                    equiv_dose<-for(i in idx_vect){
                        x <- idx_list[[1]][i]
                        j <- idx_list[[2]][i]
                        output = calc_mc_oral_equiv(conc = x,dtxsid=j,which.quantile = c(0.9),model="1compartment")
                        answer_oat1<-rbind(answer_oat1,output)
                    }
                    colnames(answer_oat1) <- c("0.9_oral_dose_oat1")
                    answer_oat1$DTXSID<-df_r$DTXSID_fake

                    set.seed(42)
                    idx_list <- list(10^(-df_r$OATP1B3)*10^6, df_r$DTXSID_fake)
                    idx_vect <- c(1:length(idx_list[[1]]))
                    answer_oat2<-data.frame()
                    equiv_dose<-for(i in idx_vect){
                        x <- idx_list[[1]][i]
                        j <- idx_list[[2]][i]
                        output = calc_mc_oral_equiv(conc = x,dtxsid=j,which.quantile = c(0.9),model="1compartment")
                        answer_oat2<-rbind(answer_oat2,output)
                    }
                    colnames(answer_oat2) <- c("0.9_oral_dose_oat2")
                    answer_oat2$DTXSID<-df_r$DTXSID_fake


                    set.seed(42)
                    idx_list <- list(10^(-df_r$BSEP)*10^6, df_r$DTXSID_fake)
                    idx_vect <- c(1:length(idx_list[[1]]))
                    answer_bsep<-data.frame()
                    equiv_dose<-for(i in idx_vect){
                        x <- idx_list[[1]][i]
                        j <- idx_list[[2]][i]
                        output = calc_mc_oral_equiv(conc = x,dtxsid=j,which.quantile = c(0.9),model="1compartment")
                        answer_bsep<-rbind(answer_bsep,output)
                    }
                    colnames(answer_bsep) <- c("0.9_oral_dose_bsep")
                    answer_bsep$DTXSID<-df_r$DTXSID_fake

                    set.seed(42)
                    idx_list <- list(10^(-df_r$PGP)*10^6, df_r$DTXSID_fake)
                    idx_vect <- c(1:length(idx_list[[1]]))
                    answer_pgp<-data.frame()
                    equiv_dose<-for(i in idx_vect){
                        x <- idx_list[[1]][i]
                        j <- idx_list[[2]][i]
                        output = calc_mc_oral_equiv(conc = x,dtxsid=j,which.quantile = c(0.9),model="1compartment")
                        answer_pgp<-rbind(answer_pgp,output)
                    }
                    colnames(answer_pgp) <- c("0.9_oral_dose_pgp")
                    answer_pgp$DTXSID<-df_r$DTXSID_fake

                    dataframe_together_pred<-Reduce(function(x, y) left_join(x, y, by = "DTXSID"), list(answer_bcrp,answer_mrp2,answer_mrp3,answer_mrp4,answer_bsep,answer_oat1,answer_oat2,answer_pgp))
                    
                    ''')
    
df_qivive = pd.DataFrame(script_r).\
            set_index(df_nuevo_def3.index).\
            rename(columns={'DTXSID':'DTXSID_fake'})

df_combi = pd.concat([df_uploaded[['name','ID','Smiles','Doses max']],df_qivive.drop(columns='DTXSID_fake'),df_uploaded[['Activity']]],axis=1)

st.write('<p style="color:green; font-size:24px;">&#10003; Transporter in vivo doses extrapolated successfully</p>', unsafe_allow_html=True)

class LogicalOrEstimatorpk(BaseEstimator, TransformerMixin):
    
    def __init__(self, k=1):
        self.k = k

    def fit(self, X, y):
        return self
    
    def predict(self, X):
        # Apply the logical OR rule to the test data
        pred = np.where(X['Doses max'].values>self.k*X.iloc[:,4:-1].max(1).values,1,0)
        return pred
    def get_params(self, deep=True):
        return {"k": self.k}
    
    def set_params(self, **parameters):
        for parameter, value in parameters.items():
            setattr(self, parameter, value)
        return self
    
ensembl_modelpk=LogicalOrEstimatorpk(k=4.696969696969697).fit(df_combi,df_combi.Activity)

predictions = None
real_class = None
mols = None

if 'df_uploaded' in st.session_state:

    predictions=ensembl_modelpk.predict(df_combi)
    real_class=df_uploaded['Activity'].values

    mols=[Chem.MolFromSmiles(m) for m in df_combi.Smiles.values]

classes = {'Non-cholestatic':0, 'Cholestatic':1}
rclasses = {0:'Non-cholestatic', 1:'Cholestatic'}

def mol2fp(mol):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    arr = np.zeros((0,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def make_data(sdf):
    mols = [m for m in Chem.SDMolSupplier(sdf) if m !=None]
    X = np.array([mol2fp(m) for m in mols])
    Y = [classes[m.GetProp('SOL_classification')] for m in mols]
    Y = np.array(Y)
    return (X, Y)

option = st.selectbox('Please select index of test molecules',[i for i in range(len(mols))])
    
if mols is not None and option is not None:

    st.write('you selected:', Chem.MolToSmiles(mols[option]))
    fp = [mol2fp(mols[option])]

    img = Draw.MolToImage(mols[option])
    bio = BytesIO()
    img.save(bio, format='png')
    st.image(img)


    st.write("Real class: ", rclasses[real_class[0]])
    st.write("Predicted class: ", rclasses[predictions[0]])

    
    # with (ro.default_converter + pandas2ri.converter).context():
    #     r_from_pd_df = ro.conversion.get_conversion().py2rpy(dataframe_together_pred)    
    # # pandas_dataframe = pandas2ri.ri2py(dataframe_together_pred)

    # # with (ro.default_converter + pandas2ri.converter).context():
    # #     r_from_pd_df = ro.conversion.get_conversion().py2rpy(dataframe_together_pred)

    

    #     # Execute the R script
    #     # result = robjects.r(colnames)

    # st.write(r_from_pd_df.head())
    

# # Check and update the state of Button 4 (Calculate with R)
# if st.session_state.button3 and not st.session_state.button4:
#     if st.button("Calculating in vivo doses"):
#         st.session_state.button4 = not st.session_state.button4

# if st.session_state.button4:
#     st.info("Calculating in vivo doses with httk library...")
#     # Configura la conversión entre pandas y R dataframes
#     pandas2ri.activate()

#     # Inicializa R
#     rpy2.rinterface.initr()

#     # Convierte df_nuevo_def3 a un objeto de datos R
#     r_dataframe = pandas2ri.py2ri(df_nuevo_def3)

#     # Envía el objeto de datos R a R
#     r.assign("df_r", r_dataframe)

#     script_r=r('''colnames(df_r)''')

#     st.write(script_r)
