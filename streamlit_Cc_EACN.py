import streamlit as st
import pickle
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
import sys
import os
import pandas as pd
import numpy as np





def EACN_calculation():
    sys.path.append(os.path.abspath('EACN'))
    from function import calculate_EACN
    #Import all necessary files
    with open('EACN/trained_pipe.pickle', 'rb') as inp:
        pipe = pickle.load(inp)
    with open('EACN/imputer.pickle', 'rb') as inp:
        imputer = pickle.load(inp) 
    with open('EACN/app_domain.pickle', 'rb') as inp:
        app_domain = pickle.load(inp) 
    with open('EACN/features_to_drop_Padel.pickle', 'rb') as inp:
        features_to_drop = pickle.load(inp)
    st.title('EACN calculation')
      
    st.write('In the field below, enter SMILES of your molecule')
    smi = st.text_input('Enter SMILES of your molecule, if the molecule is valid, it will be drawn')
    st.write(smi)
    if smi:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            img = Draw.MolToImage(mol)
            st.image(img, smi)
            st.write('Would you like to use SGIR predictor? :red[Torch-moleule package must be installed, as welll as CUDA should be installed and configured]')
            gpu_available = st.button('Yes', type = 'primary')
            gpu_non_available = st.button('No', type = 'secondary')
            if gpu_available:
                from torch_molecule import SGIRMolecularPredictor
                sgir = SGIRMolecularPredictor()
                sgir.load_from_local('EACN/sgir_final.pt')
                st.write(':green[SGIR model loaded sucessfully]')
                svr_result, sgir_result, flag_AD = calculate_EACN(smi, app_domain = app_domain,  imputer = imputer, features_to_drop = features_to_drop, svr_model = pipe, use_sgir = True, sgir_model = sgir)
                if flag_AD:
                    st.write(":red[You're outside the applicability domain!")
                st.write('EACN of the molecule according to SVR is **{:.3g}**, according to SGIR is **{:.3g}**'.format(svr_result , sgir_result))
            if gpu_non_available:
                svr_result, flag_AD = calculate_EACN(smi, app_domain = app_domain,  imputer = imputer, features_to_drop = features_to_drop, svr_model = pipe)
                if flag_AD:
                    st.write(":red[You're outside the applicability domain!]")
                st.write('EACN of the molecule according to SVR is **{:.3g}**'.format(svr_result))
        
        
        
    

def Cc_calculation():
    st.title('Cc calculation')
    
    sys.path.append(os.path.abspath('Critical curvature'))
    from functions import embed_optimize, calculate_cc, calculate_padel_descriptors, calculate_rdkit_descriptors, calculate_leverage
    st.write('In the field below, enter SMILES of your molecule')
    smi = st.text_input('Enter SMILES of your molecule, if the molecule is valid, it will be drawn')
    if smi:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            img = Draw.MolToImage(mol)
            st.image(img, smi)
    st.write('Surfactant class selection')
    type = st.radio('Select type of your surfactant', ['Cationic', 'Anionic', 'Nonionic', 'Extended'], index = None)
    st.write(type)
    sys.path.append(os.path.abspath('Critical curvature'))
    if type == 'Cationic':
        with open('Critical curvature/Appl_domains_data/Cationic.pickle', 'rb') as inp:
            AData = pickle.load(inp)
            warning = AData['warning']
            invxtx = AData['invxtx']
        with open('Critical curvature/models_for_prediction_cationic.pickle', 'rb') as inp:
            cationic_model = pickle.load(inp)[0]
        with open('Critical curvature/Features_PaDEL/featured_to_drop_padel.pickle', 'rb') as inp:
            features_to_drop_padel = pickle.load(inp)
            features_to_drop_cationic_padel = features_to_drop_padel['PaDEL_cationic']
        with open('Critical curvature/Features_RDKit/features_to_drop_rdkit.pickle', 'rb') as inp:    
            features_to_drop_cationic_rdkit = pickle.load(inp)['RDKit_cationic']
        with open('Critical curvature/Features_United/features_to_drop_united.pickle', 'rb') as inp:
            features_to_drop_cationic = pickle.load(inp)['United_cationic']
        with open('Critical curvature/Features_PaDEL/imputers_for_padel.pickle', 'rb') as inp:
            imputers  = pickle.load(inp)
            imputer_cationic = imputers['PaDEL_cationic']
            
        array_of_cationics = []
    #feats = calculate_padel_descriptors(smi = smi, imputer = imputer_cationic, features_to_drop = features_to_drop_cationic, extended = False)
        for i in range(5):
            feats_padel = calculate_padel_descriptors(smi = smi, imputer = imputer_cationic, features_to_drop = features_to_drop_cationic_padel, extended = False)
            feats_rdkit = calculate_rdkit_descriptors(smi = smi,  features_to_drop = features_to_drop_cationic_rdkit, if_non_ionic = False)
            feats_cationic = pd.concat([feats_padel, feats_rdkit], axis = 1)
            feats_cationic.drop(columns = features_to_drop_cationic, inplace=True)
            array_of_cationics.append(feats_cationic)
        average_descriptors = sum(array_of_cationics)/len(array_of_cationics)
        leverage = calculate_leverage(average_descriptors, invxtx)
        if leverage > warning:
            st.write(":red[You're outside the applicability domain! Leverage is {:.3g}, whereas warning level is {:.3g}]".format(leverage, warning))
        st.write('Cc acc. to Ridge-based model is **{:.4g}** '.format(calculate_cc(cationic_model, feats_cationic),  ))

    if type == 'Anionic':
        with open('Critical curvature/Appl_domains_data/Anionic.pickle', 'rb') as inp:
            AData = pickle.load(inp)
            warning = AData['warning']
            invxtx = AData['invxtx']
        with open('Critical curvature/models_for_prediction_anionic.pickle', 'rb') as inp:
            models = pickle.load(inp)
            anionic_model_svr = models[0]
            anionic_model_ridge = models[1]
        with open('Critical curvature/Features_PaDEL/featured_to_drop_padel.pickle', 'rb') as inp:
            features_to_drop_padel = pickle.load(inp)
            features_to_drop_anionic = features_to_drop_padel['PaDEL_anionic']
        with open('Critical curvature/Features_PaDEL/featured_to_drop_padel.pickle', 'rb') as inp:
            features_to_drop_padel = pickle.load(inp)
            features_to_drop_anionic = features_to_drop_padel['PaDEL_anionic']
        with open('Critical curvature/Features_PaDEL/imputers_for_padel.pickle', 'rb') as inp:
            imputers  = pickle.load(inp)
            imputer_anionic = imputers['PaDEL_anionic']
            
        feats = calculate_padel_descriptors(smi = smi, imputer = imputer_anionic, features_to_drop = features_to_drop_anionic, extended = False)
        leverage = calculate_leverage(feats, invxtx)
        if leverage > warning:
            st.write(":red[You're outside the applicability domain! Leverage is {:.3g}, whereas warning level is {:.3g}]".format(leverage, warning))
        st.write('Cc for  SVR-based model is {:.3g}, Cc for  Ridge-based model is {:.3g} '.format(calculate_cc(anionic_model_svr, feats), calculate_cc(anionic_model_ridge, feats) ))

    if type == 'Nonionic':
        with open('Critical curvature/Appl_domains_data/Nonionic.pickle', 'rb') as inp:
            AData = pickle.load(inp)
            warning = AData['warning']
            invxtx = AData['invxtx']
        with open('Critical curvature/models_for_prediction_non_ionic.pickle', 'rb') as inp:
            models = pickle.load(inp)
            non_ionic_model_svr = models[0]
            non_ionic_model_ridge = models[1]
        with open('Critical curvature/Features_RDKit/features_to_drop_rdkit.pickle', 'rb') as inp:
            features_to_drop_nonionic = pickle.load(inp)['RDKIt_non_ionic']
        array_svr = []
        array_linear = []
        array_of_nonionics = []
        for i in range(5):
            feats = calculate_rdkit_descriptors(smi = smi, features_to_drop = features_to_drop_nonionic)
            array_svr.append(calculate_cc(non_ionic_model_svr, feats))
            array_linear.append(calculate_cc(non_ionic_model_ridge, feats))
            array_of_nonionics.append(feats)
        average_descriptors = sum(array_of_nonionics)/len(array_of_nonionics)
        leverage = calculate_leverage(average_descriptors, invxtx)
        if leverage > warning:
            st.write(":red[You're outside the applicability domain! Leverage is {:.3g}, whereas warning level is {:.3g}]".format(leverage, warning))
        st.write('Cc for SVR-based model is {:.3g}, for Ridge-based model is {:.3g} '.format(sum(array_svr)/10, sum(array_linear)/10))
        st.write('Standard deviation for SVR is {:.3f}, for linear regression is {:.3f} '.format(np.std(array_svr), np.std(array_linear)))

    if type == 'Extended':
        with open('Critical curvature/Appl_domains_data/Extended.pickle', 'rb') as inp:
            AData = pickle.load(inp)
            warning = AData['warning']
            invxtx = AData['invxtx']
        with open('Critical curvature/models_for_prediction_extended.pickle', 'rb') as inp:
            extended_model = pickle.load(inp)[0]
        with open('Critical curvature/Features_PaDEL/featured_to_drop_padel.pickle', 'rb') as inp:
            features_to_drop_padel = pickle.load(inp)
            features_to_drop_extended = features_to_drop_padel['PaDEL_extended']
        with open('Critical curvature/Features_PaDEL/imputers_for_padel.pickle', 'rb') as inp:
            imputers  = pickle.load(inp)
            imputer_extended = imputers['PaDEL_extended']
        
        feats = calculate_padel_descriptors(smi = smi, imputer = imputer_extended, features_to_drop = features_to_drop_extended, extended = True)
        leverage = calculate_leverage(feats, invxtx)
        if leverage > warning:
            st.write(":red[You're outside the applicability domain! Leverage is {:.3g}, whereas warning level is {:.3g}]".format(leverage, warning))
        st.write('Cc for SVR-based model {:.3g}'.format(calculate_cc(extended_model, feats)))
        
            
pg = st.navigation([EACN_calculation, Cc_calculation])

pg.run()
