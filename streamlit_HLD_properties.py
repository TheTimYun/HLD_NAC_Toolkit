import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem import Draw
st.set_page_config(layout='wide')


from functions_for_surfactants import draw_type_1, draw_type_2, draw_type_3, draw_tubes, HLD_calculation, len_calculation, xi_calculation





        
def HLD_tubes():
    #First app - to calculate HLD and draw tubes of Winsor Type I, II, III emulsions
    st.title('HLD calculation and drawing tubes')
    st.write('Calculation of HLD for the range of tubes and plotting their phase state schematically')
    #Selection of surfactant type, which defines coefficients for HLD calculation
    surfactant_selection = st.selectbox('Type of surfactants', options = ('Ionic', 'Nonionic'))
    #HLD property that is varied in some range
    prop = st.selectbox('Property to be varied', options = (None, 'Cc', 'EACN', 'Temperature', 'Salinity'), placeholder = None)
    if prop:
        #when property is selected, we enter its range and fixed values of others - 4 different scenarios depending on the selected property
        st.write(prop)
        if prop == 'EACN':
            col1, col2, col3 = st.columns(3)
            with col1:
                Cc = st.number_input('Enter Cc', value = 0.00)
                st.write('Cc is {}'.format(Cc))
            with col2:
                Temp = st.number_input('Enter temperature', value = 25.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a range of EACN", -15.00, 15.00, (0.00, 10.00))
            EACN_range = np.linspace(values[0], values[1], 10)
            #calculates HLD and draw tubes
            HLDs = HLD_calculation(EACN = EACN_range, Cc = Cc, Sal = Sal, Temp = Temp, Type = surfactant_selection)
            fig = draw_tubes(HLDs, EACN_range, 'EACN')
            #shows figure in streamlit
            st.pyplot(fig)
        #do the same stages if other properties are selected
        elif prop == 'Cc':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Temp = st.number_input('Enter temperature', value = 25.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a range of Cc", -15.00, 15.00, (-5.00, 5.00))
            Cc_range = np.linspace(values[0], values[1], 10)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc_range, Sal = Sal, Temp = Temp, Type = surfactant_selection)
            fig = draw_tubes(HLDs, Cc_range, 'Cc')
            st.pyplot(fig)
        elif prop == 'Temperature':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Cc = st.number_input('Enter Cc', value = 0.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a temperature range", 5.0, 95.0, (25.0, 50.0))
            temp_range = np.linspace(values[0], values[1], 10)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc, Sal = Sal, Temp = temp_range, Type = surfactant_selection)
            fig = draw_tubes(HLDs, temp_range, 'T')
            st.pyplot(fig)
        elif prop == 'Salinity':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Cc = st.number_input('Enter Cc', value = 0.00)
            with col3:
                Temp = st.number_input('Enter temperature', value = 25.00)
            values = st.slider("Select a salinity range, g NaCl/100 mL", 0.00, 100.00, (0.00, 5.00))
            Sal_range = np.linspace(values[0], values[1], 10)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc, Sal = Sal_range, Temp = Temp, Type = surfactant_selection)
            fig = draw_tubes(HLDs, Sal_range, 'S')
            st.pyplot(fig)

def Calculation_of_xi():
    #calculates xi parameter for HLD-NAC concept
    st.title('Calculation of xi parameter in HLD-NAC equation')
    #Enter surafactant properties
    col1, col2 = st.columns(2)
    with col1:
       #Select the input type - with SMILES (automatic MWt and tail length calculation) or manual input
        type_of_input = st.selectbox('How do you want to input tail length of surfactant', options = (None, 'SMILES', 'Manual'))
        if type_of_input == 'Manual':
            L = st.slider('Select length of carbon chain in Å', 0, 100, 20)
            st.write('Tail length is {}Å'.format(L))
        #Input surfactant molecular structure (SMILES) with subsequent calculation of tail length and MWt 
        elif type_of_input == 'SMILES':
            smi = st.text_input('Enter SMILES of you surfactant - only surface-active part, without counterion!')
            #transforms SMILES into RDKit mol object
            mol = Chem.MolFromSmiles(smi)
            #if molecule is valid, then calculation proceeds
            if mol:
                L = len_calculation(mol)
                st.write('Tail length is {}Å'.format(L))
                img = Draw.MolToImage(mol)
                st.image(img, smi)
        #Input of interfacial area of surfactant 
        with col2:
            interfacial_area = st.slider('Select interfacial area in Å2', 0, 200, 45)
            st.write('Interfacial area is {} '.format(interfacial_area))
    #Input of oil properties - moleular wight and density
    col1, col2 = st.columns(2)
    with col1:
        #Type of input - with SMILES (automatic MWt calculation) or manual MWT value
        type_of_input = st.selectbox('How do you want to input molecular weight of oil', options = (None, 'SMILES', 'Manual'))
        if type_of_input == 'Manual':
            MWt = st.slider('Select oil molecular weight', 0, 400, 95)
            st.write('Molecular weight of oil is {}'.format(MWt))
        elif type_of_input == 'SMILES':
            smi_oil = st.text_input('Enter SMILES of your oil')
            #SMILES is transformed into molecule rdkit object
            mol = Chem.MolFromSmiles(smi_oil)
            #If molecule is valid, that molecular weight is calculated
            if mol:
                MWt = MolWt(mol)
                st.write('Moleular weight is {}'.format(MWt))
                img = Draw.MolToImage(mol)
                st.image(img, smi_oil)
        with col2:
            density = st.slider('Select density, g/cc', 0.0, 1.2, 0.8)
            st.write('Density is{} '.format(density))
        
        #Final calculation of xi using special function
        surfactant_type = st.selectbox('Type of surfactants', options = (None, 'Ionic', 'Nonionic'))
        if surfactant_type:
            xi =  xi_calculation(L, interfacial_area, density, MWt, surfactant_type)
            st.write('**ξ is {}**'.format(xi))
        
       
            
    

def Phase_diagram():
    #This app plots the fishtail diagram - concentration vs different HLD properties (temperature, salinity, HLD, Cc) vs ME type
    st.title('Phase diagram plotting')
    st.write('Calculation of phase diagram on the base of HLD-NAC framework: concentration vs external parameters')
    #Three blocks - surfactant properties, HLD properties, system properties
    st.write('**Input of surfactant intrinsic properties**')
    #Select the input type - with SMILES (automatic MWt and tail length calculation) or manual input
    type_of_input = st.selectbox('How do you want to input properties of your surfactants (tail length, molecular weight?', options = (None, 'SMILES', 'Manual'))
    if type_of_input == 'Manual':
        col1, col2 = st.columns(2)
        with col1:
            L = st.number_input('Enter tail length in Å', value = 0.00)
            st.write('Tail length is {}Å'.format(L))
        with col2:
            MWt = st.number_input('Enter molecular weight', value = 0.00)
            st.write('Molecular weight is {} '.format(MWt))
    #Input surfactant molecular structure (SMILES) with subsequent calculation of tail length and MWt 
    elif type_of_input == 'SMILES':
        smi = st.text_input('Enter SMILES of you surfactant - only surface-active part, without counterion!')
        #SMILES is transformed into molecule object
        mol = Chem.MolFromSmiles(smi)
        if mol:
            #if molecule is valid, MWt and tail length are calculated
            L = len_calculation(mol)
            MWt = MolWt(mol)
            col1, col2 = st.columns(2)
            with col1:
                st.write('Tail length is {}'.format(L))
            with col2:
                st.write('Molecular weight is {} '.format(MWt))
            img = Draw.MolToImage(mol)
            st.image(img, smi)

    #Second block - input of HLD properties
    st.write('**Input of HLD properties**')
    #Selection of surfactant type, which defines HLD calculation
    surfactant_selection = st.selectbox('Type of surfactants', options = ('Ionic', 'Nonionic'))
    #HLD property that is varied in some range (for fishtail plot of nonionics usually - temperature)
    prop = st.selectbox('Property to be varied', options = (None, 'Cc', 'EACN', 'Temperature', 'Salinity'), placeholder = None)
    if prop:
        st.write(prop)
        #when property is selected, we enter its range and fixed values of others - 4 different scenarios depending on the selected property
        if prop == 'EACN':
            col1, col2, col3 = st.columns(3)
            with col1:
                Cc = st.number_input('Enter Cc', value = 0.00)
                st.write('Cc is {}'.format(Cc))
            with col2:
                Temp = st.number_input('Enter temperature', value = 25.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a range of EACN", -15.00, 15.00, (0.00, 10.00))
            EACN_range = np.linspace(values[0], values[1], 100)
            #calculates HLD range
            HLDs = HLD_calculation(EACN = EACN_range, Cc = Cc, Sal = Sal, Temp = Temp, Type = surfactant_selection)
        #do the same stages if other properties are selected
        elif prop == 'Cc':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Temp = st.number_input('Enter temperature', value = 25.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a range of Cc", -15.00, 15.00, (-5.00, 5.00))
            Cc_range = np.linspace(values[0], values[1], 100)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc_range, Sal = Sal, Temp = Temp, Type = surfactant_selection)
        elif prop == 'Temperature':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Cc = st.number_input('Enter Cc', value = 0.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a temperature range", 5.0, 95.0, (25.0, 50.0))
            temp_range = np.linspace(values[0], values[1], 100)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc, Sal = Sal, Temp = temp_range, Type = surfactant_selection)
        elif prop == 'Salinity':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Cc = st.number_input('Enter Cc', value = 0.00)
            with col3:
                Temp = st.number_input('Enter temperature', value = 25.00)
            values = st.slider("Select a salinity range, g NaCl/100 mL", 0.00, 100.00, (0.00, 5.00))
            Sal_range = np.linspace(values[0], values[1], 100)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc, Sal = Sal_range, Temp = Temp, Type = surfactant_selection)
        #Third block - input of system parameters
        st.write('**Input of system parameters**')
        col1, col2, col3 = st.columns(3)
        with col1:
            #Select volume fraction of water with subsequent calculation of Vo
            Vw = st.slider("Enter the volume fraction of water ", 0.01, 0.99, 0.5, step = 0.01)
            Vo = 1 - Vw
        with col2:
            #Enter weight concentration and recalculate to mol/ml
            wt_concentrations_range = st.slider("Select a concentration in % wt",0.1 , 20.1, (1.1, 10.1), step = 0.1)
            wt_concentrations = np.linspace(wt_concentrations_range[0], wt_concentrations_range[1], 100)
            if MWt != 0:
                C_mol_range = ((1000 * wt_concentrations / 100)/MWt)/1000
        with col3:
                #Enter interfacial area of surfactant on O/W boundary
            interfacial_area = st.number_input('Enter interfacial area of surfactant on O/W boundary in Å2', value = 45)
            st.write('Interfacial area is {} '.format(interfacial_area))
        #xi calculation - manual input (from app) or automatic calculation
        calculation_type = st.selectbox('How do you want to calculate xi', options = (None, 'Automatically', 'Manually'), placeholder = None)
        if calculation_type == 'Automatically':
            #Only oil parameters are needed for successful calculation, as surfactant parameters were given earlier
            col1, col2 = st.columns(2)
            with col1:
                density = st.number_input('Enter density of your oil', value = 0.8)
                st.write('Density is {}'.format(density))
            with col2:
                MWt_oil = st.number_input('Enter molecular weight of your oil', value = 100)
                st.write('Molecular weight of the oil is {} '.format(MWt_oil))
            #xi calculation with special function
            xi = xi_calculation(L, interfacial_area, density, MWt_oil, surfactant_selection)
            st.write('**ξ** is {}'.format(xi))
        elif calculation_type == 'Manually':
            #Manual input if user wants to do it
            xi = st.number_input('Enter xi', value = 100)
            st.write('**ξ** is {}'.format(xi))
        #if xi value is sucessfully calculated and all other inout parameters are defined, then other calculations take place
        if calculation_type:
                        #Calculation of As - interfacial area of surfactant molecules in water phase, Rw_red, Ro_ref 
            As = {C: Vw*C*6E23*interfacial_area/1E24 for C in C_mol_range}
            Rw_ref_range = {C:3*Vw/As[C] for C in C_mol_range}
            Ro_ref_range = {C:3*Vo/As[C] for C in C_mol_range}
            from itertools import product
            #Calculation of 1/Rw - 1/Ro = H
            Hs = -(HLDs/L)
            #two zero arrays for subsequent filling
            Ros = np.zeros((len(Hs), len(C_mol_range)), dtype = 'float')
            Rws = np.zeros((len(Hs), len(C_mol_range)), dtype = 'float')
            #Calculation of Ro, Rw and filling arrays
            for (i, C), (j, H) in product(enumerate(C_mol_range), enumerate(Hs)):
                if H > 0:
                    Rw = Rw_ref_range[C]
                    Ro = 1/(H+1/Rw)
                    Ros[i,j] = Ro
                    Rws[i,j] = Rw
                else:
                    Ro = Ro_ref_range[C]
                    Rw = 1/(np.abs(H)+1/Ro)
                    Ros[i,j] = Ro
                    Rws[i,j] = Rw
            #NAC calculation and determination of ME type (Rw > Ro - I, Ro > Rw - II, 1/Ro + 1/Rw < 1/xi - III). Emulsion type is coded into NACs_bin
            NACs = 0.5*((1/Ros)+(1/Rws))
            NACs_bin = Ros > Rws
            NACs_bins  = (Ros > Rws).astype(int)
            mask = (NACs<1/xi)
            NACs_bins[mask] = 2
            #Meshgrid for plotting
            if prop == 'EACN':
                HH, CC = np.meshgrid(EACN_range, wt_concentrations)
            elif prop == 'Cc':
                HH, CC = np.meshgrid(Cc_range, wt_concentrations)
            elif prop == 'Salinity':
                HH, CC = np.meshgrid(Sal_range, wt_concentrations)
            elif prop == 'Temperature':
                HH, CC = np.meshgrid(temp_range, wt_concentrations)
            #Plotting phase diagram
            fig, ax = plt.subplots()
            NACs_bins_fl = NACs_bins.flatten()
            color = {0:'blue', 1:'red', 2:'green'}
            labels = {0:'Type I', 1:'Type II', 2:'Type III'}
            
            
            for val in (0,1,2):
                mask = (NACs_bins_fl==val)
                scatter = ax.scatter(x = HH.flatten()[mask], y = CC.flatten()[mask], c = color[val], label = labels[val])
            
            ax.set_xlabel(prop)
            ax.set_ylabel('Concentration - weight')
            ax.legend()
            st.pyplot(fig)
            
        
def Volumes_of_phases():
    #This app calculates volume of all phases in the given system under equilibrium
    st.title('Phases volumes and boundaries')
    st.write('Calculation of oil, water and ME boundaries in tubes')
    st.write('**Input of surfactant intrinsic properties**')
    #Select the input type - with SMILES (automatic MWt and tail length calculation) or manual input
    type_of_input = st.selectbox('How do you want to input properties of your surfactants (tail length, molecular weight?', options = (None, 'SMILES', 'Manual'))
    if type_of_input == 'Manual':
        col1, col2 = st.columns(2)
        with col1:
            L = st.number_input('Enter tail length in Å', value = 0)
            st.write('Tail length is {}Å'.format(L))
        with col2:
            MWt = st.number_input('Enter molecular weight', value = 0)
            st.write('Molecular weight is {} '.format(MWt))
    #Input surfactant molecular structure (SMILES) with subsequent calculation of tail length and MWt
    elif type_of_input == 'SMILES':
        #SMILES is transformed into molecule object
        smi = st.text_input('Enter SMILES of you surfactant - only surface-active part, without counterion!')
        mol = Chem.MolFromSmiles(smi)
        #if molecule is valid, MWt and tail length are calculated
        if mol:
            L = len_calculation(mol)
            MWt = MolWt(mol)
            col1, col2 = st.columns(2)
            with col1:
                st.write('Tail length is {}'.format(L))
            with col2:
                st.write('Molecular weight is {} '.format(MWt))
            img = Draw.MolToImage(mol)
            st.image(img, smi)
    #Second block - input of HLD properties
    #HLD property that is varied in some range 
    st.write('**Input of HLD properties**')
    #Selection of surfactant type, which defines HLD calculation
    surfactant_selection = st.selectbox('Type of surfactants', options = ('Ionic', 'Nonionic'))
    prop = st.selectbox('Property to be varied', options = (None, 'Cc', 'EACN', 'Temperature', 'Salinity'), placeholder = None)
    if prop:
        st.write(prop)
        #when property is selected, we enter its range and fixed values of others - 4 different scenarios depending on the selected property
        if prop == 'EACN':
            col1, col2, col3 = st.columns(3)
            with col1:
                Cc = st.number_input('Enter Cc', value = 0.00)
                st.write('Cc is {}'.format(Cc))
            with col2:
                Temp = st.number_input('Enter temperature', value = 25.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a range of EACN", -15.00, 15.00, (0.00, 10.00))
            EACN_range = np.linspace(values[0], values[1], 20)
            #calculates HLD range
            HLDs = HLD_calculation(EACN = EACN_range, Cc = Cc, Sal = Sal, Temp = Temp, Type = surfactant_selection)
        #do the same stages if other properties are selected
        elif prop == 'Cc':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Temp = st.number_input('Enter temperature', value = 25.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a range of Cc", -15.00, 15.00, (-5.00, 5.00))
            Cc_range = np.linspace(values[0], values[1], 20)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc_range, Sal = Sal, Temp = Temp, Type = surfactant_selection)
        elif prop == 'Temperature':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Cc = st.number_input('Enter Cc', value = 0.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a temperature range", 5.00, 95.00, (25.00, 50.00))
            temp_range = np.linspace(values[0], values[1], 20)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc, Sal = Sal, Temp = temp_range, Type = surfactant_selection)
        elif prop == 'Salinity':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Cc = st.number_input('Enter Cc', value = 0.00)
            with col3:
                Temp = st.number_input('Enter temperature', value = 25.00)
            values = st.slider("Select a salinity range, g NaCl/100 mL", 0.00, 100.00, (0.00, 5.00))
            Sal_range = np.linspace(values[0], values[1], 20)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc, Sal = Sal_range, Temp = Temp, Type = surfactant_selection)
        #Third block - input of system parameters
        st.write('**Input of system paramters**')
        col1, col2, col3 = st.columns(3)
        with col1:
            #Enter volume of water and oil in the system
            Vw = st.slider("Enter the water fraction in the system", 0.01, 0.99, 0.5, step = 0.01)
            Vo = 1 - Vw
        with col2:
             #Enter weight concentration and recalculate to mol/ml
            wt_concentration = st.slider("Select a concentration of surfactant in water in % wt",0.1 , 20.1, step = 0.1)
            if MWt != 0:
                mol_concentration = ((1000 * wt_concentration / 100)/MWt)/1000
            Vs = (wt_concentration/100)*Vw
        with col3:
            interfacial_area = st.number_input('Enter interfacial area in Å2', value = 45)
            st.write('Interfacial area is {} '.format(interfacial_area)) 
        #Calculation of As - interfacial area of surfactant molecules in water phase, Rw_red, Ro_ref 
        As =  Vw*mol_concentration*6E23*interfacial_area/1E24  
        Rw_ref = 3*Vw/As
        Ro_ref = 3*Vo/As
        
        #xi calculation
        calculation_type = st.selectbox('How do you want to calculate xi', options = (None, 'Automatically', 'Manually'), placeholder = None)
        if calculation_type == 'Automatically':
            #Only oil parameters are needed for successful calculation, as surfactant parameters were given earlier
            col1, col2 = st.columns(2)
            with col1:
                density = st.number_input('Enter density of your oil', value = 0.8)
                st.write('Density is {}'.format(density))
            with col2:
                MWt_oil = st.number_input('Enter molecular weight of your oil', value = 100)
                st.write('Molecular weight of the oil is {} '.format(MWt_oil))
            #calculates xi with the spoecial function
            xi = xi_calculation(L, interfacial_area, density, MWt_oil, surfactant_selection)
            st.write('**ξ** is {}'.format(xi))
        elif calculation_type == 'Manually':
            xi = st.number_input('Enter xi', value = 100)
            st.write('**ξ** is {}'.format(xi))
        
        
        #if xi value is calculated (i.e. calculation type is chosen) and other properties are defined, then calculations of volume take place
        if calculation_type:
            #calculation of NACs and emulsion types
            Hs = -(HLDs/L)
            #Creating zero arrays for radiia of water and oil droplets
            Ros = np.zeros((len(Hs)), dtype = 'float')
            Rws = np.zeros((len(Hs)), dtype = 'float')
            #Calculation of Ro, Rw and filling arrays
            for i,H in enumerate(Hs):
                if H > 0:
                    Rw = Rw_ref
                    Ro = 1/(H+1/Rw)
                else:
                    Ro = Ro_ref
                    Rw = 1/(np.abs(H)+1/Ro)
                Ros[i] = Ro
                Rws[i] = Rw
            
            #NAC calculation and determination of ME type (Rw > Ro - I, Ro > Rw - II, 1/Ro + 1/Rw < 1/xi - III). Emulsion type is coded into NACs_bin
            NACs = 0.5*((1/Ros)+(1/Rws))
            NACs_bin = (Ros > Rws).astype(int)
            NACs_bin[NACs<1/xi] = 2
            mask = NACs_bin==2
            #Calculation of oil, water droplets in MEs of different Types
            Rome  = Ros
            Rwme = Rws
            Rwme[mask] = 1/((HLDs[mask]/(2*L))+1/xi)
            Rome[mask] = 1/(1/Rwme[mask] - (HLDs[mask]/L))
            #Calculation of volumes of oil and water in ME
            Vome = Rome*As/3
            Vwme = Rwme*As/3
            Vtot = Vome + Vwme
            
            # Calculation of volume fractions of oil and water in ME
            phio = np.zeros(20)
            phiw = np.zeros(20)
            
            for i in range(20):
                if HLDs[i] < 0:
                    phio[i] = Vome[i]/Vtot[i]
                    phiw[i] = 1 - phio[i]
                else:
                    phiw[i] = Vwme[i]/Vtot[i]
                    phio[i] = 1 - phiw[i]
            
            #calculation of volume of oil, water and ME in system
            Vo_range = np.zeros(20)
            Vw_range = np.zeros(20)
            Vme_range = np.zeros(20)
            
            for i in range(20):
                if NACs_bin[i] == 0:
                    Vo_range[i] = Vo - Vome[i]
                    Vw_range[i] = Vw + Vs + Vome[i] 
                elif NACs_bin[i] == 1:
                    Vo_range[i] = Vo + Vs + Vwme[i]
                    Vw_range[i] = Vw - Vwme[i]
                elif (NACs_bin[i] == 2):
                    Vo_range[i] = Vo - Vome[i]
                    Vw_range[i] = Vw - Vwme[i]
                    Vme_range[i] = Vs + Vome[i] + Vwme[i]
            
            #Calculation of voume fractions of oil, water and ME in system 
            Vt = Vo_range+Vw_range+Vme_range
            fr_o = Vo_range/Vt 
            fr_w = Vw_range/Vt
            fr_me = Vme_range/Vt  
            
            #Boundaries calculatuion
            boundary_1 = np.zeros(20)
            boundary_2 = np.zeros(20)
            for i in range(20):
                if fr_me[i] > 0:
                    boundary_1[i] = fr_w[i] + fr_me[i]
                    boundary_2[i] = fr_w[i]
                else:
                    if fr_w[i] > fr_o[i]:
                        boundary_1[i] = fr_w[i]
                        boundary_2[i]  = 0
                    else:
                        boundary_2[i] = fr_w[i]
                        boundary_1[i] = 1
            #Plotting
            from functions_for_surfactants import draw_tubes_volumes
            if prop == "EACN":
                varied_value_range = EACN_range
            elif prop=='Cc':
                 varied_value_range = Cc_range
            elif prop=='Temperature':
                varied_value_range =  temp_range
            elif prop=='Salinity':
                varied_value_range = Sal_range
            fig = draw_tubes_volumes(boundary_1, boundary_2, fr_o, fr_w, fr_me,  varied_value_range, prop, HLDs)
            
            st.pyplot(fig, width =  1600)
            col1, col2, col3 = st.columns(3)
            with col1:
                fig, ax = plt.subplots(1,1, figsize = (12,12))
                ax.plot(varied_value_range, Vome, label = 'Volume of oil in ME')
                ax.plot(varied_value_range, Vwme, label = 'Volume of water in ME')
                ax.set_xlabel(prop)
                ax.set_ylabel('Volume, mL')
                ax.legend()
                ax.grid(visible = True)
                st.pyplot(fig)
                st.write('Volume of oil and ater in ME')
            with col2:
                fig, ax = plt.subplots(1,1, figsize = (12,12))
                ax.plot(varied_value_range, Vo_range, label = 'Volume of oil')
                ax.plot(varied_value_range, Vw_range, label = 'Volume of water ')
                ax.plot(varied_value_range, Vme_range, label = 'Volume of ME ')
                ax.set_xlabel(prop)
                ax.set_ylabel('Volume, mL')
                ax.legend()
                ax.grid(visible = True)
                st.pyplot(fig)
                st.write('Volumes of oil, water and ME')
            with col3:
                fig, ax = plt.subplots(1,1, figsize = (12,12))
                ax.plot(varied_value_range, boundary_1, label = 'Upper boundary')
                ax.plot(varied_value_range, boundary_2, label  ='Lower boundary')
                ax.set_xlabel(prop)
                ax.set_ylabel('Volume, mL')
                ax.legend()
                ax.grid(visible = True)
                st.pyplot(fig)
                st.write('Upper and lower boundaries')
                                     

def Salts_additives_calculator():
    st.title('Calculation of equivalent water salinity')
    st.write('Calculates equivalent salinity in terms of g Nacl / 100 mL of the aqueous phase')
    #Selection of surfactant type
    surfactant_type = st.selectbox('Type of your surfactant', options = [None, 'Cationic', 'Anionic', 'Nonionic'])
    #For ionic surfactant
    if (surfactant_type == 'Cationic')|(surfactant_type=='Anionic'):
        sum_eq = 0
        import pickle
        #Loading data
        with open('Ionic_salts_data.pickle', 'rb') as inp:
            cat_dict, an_dict = pickle.load(inp)
        #Each iteration of this function creates a row with input windows for cation, anion and concentration of salt
        def render_row(i):
            col1, col2, col3, col4 = st.columns(4)
            #Salt selection
            with col1:
                cation = st.selectbox('Cation', options = list(cat_dict.keys()), index = 2, key = 'A{}'.format(i))
            with col2:
                anion = st.selectbox('Anion', options =  list(an_dict.keys()), index = 8,key =  'B{}'.format(i))
            #Getting MWt of salt and relative affinity of cation and anion
            MW = cat_dict[cation]['MW, g/mol']*an_dict[anion]['n'] + an_dict[anion]['MW, g/mol']*cat_dict[cation]['p']
            cat_aff = cat_dict[cation]['Rel.Affin'] 
            an_aff = an_dict[anion]['Rel.Affin'] 
            #If surfactant is anionic, only cationic affect it and vice versa. Proper relative affinity is chosen
            if surfactant_type == 'Cationic':
                aff = (1/(an_dict[anion]['n']*cat_dict[cation]['p']))/an_aff
            elif surfactant_type == 'Anionic':
                aff = aff = (1/(an_dict[anion]['n']*cat_dict[cation]['p']))/cat_aff
            #Salinity input
            with col3:
                Salinity = st.number_input('Enter salt concentration in g/100 mL', value = 0.00, key = 'C{}'.format(i) )
            #Calculation of mole concentration and equivalent salinity
            M = Salinity/MW
            S_eq = M/aff*58.5
            with col4:
                st.metric('Equivalent salinity is', value = round(S_eq,2))
            return S_eq
        #Five iterations of the row and final result
        for i in range(5):
            sum_eq += render_row(i)
        st.write('Summed equivalent salinity is **{:.2f}** g NaCl/ 100 mL'.format(sum_eq))
    #For nonionic surfactant
    elif surfactant_type == 'Nonionic':
        sum_eq = 0
        import pickle
        #Loading data
        with open('Nonionic_salts_data.pickle', 'rb') as inp:
            cat_dict, an_dict = pickle.load(inp)
        #Each iteration of this function creates a row with input windows for cation, anion and concentration of salt
        def render_row(i):
            col1, col2, col3, col4 = st.columns(4)
            #Salt selection
            with col1:
                cation = st.selectbox('Cation', options = list(cat_dict.keys()), index = 2, key = 'A{}'.format(i))
            with col2:
                anion = st.selectbox('Anion', options =  list(an_dict.keys()), index = 8,key =  'B{}'.format(i))
            #Getting MWt of salt and relative affinity of cation and anion
            MW = cat_dict[cation]['MW, g/mol']*an_dict[anion]['n'] + an_dict[anion]['MW, g/mol']*cat_dict[cation]['p']
            cat_bf = cat_dict[cation]['B']
            an_bf = an_dict[anion]['B']
            #Dmod is an equivalent of relative affinity
            BMod = (cat_bf * cat_dict[cation]['p'] + an_bf * an_dict[anion]['n'])/cat_dict[cation]['p']
            #Salinity input
            with col3:
                Salinity = st.number_input('Enter salt concentration in g/100 mL', value = 0.00, key = 'C{}'.format(i) )
            #Calculation of mole concentration and equivalent salinity
            b = 96*BMod/MW
            S_eq = (Salinity*b)/0.13
            with col4:
                st.metric('Equivalent salinity is', value = round(S_eq,2))
            return S_eq
        #Five iterations of the row and final result
        for i in range(5):
            sum_eq += render_row(i)
        st.write('Summed equivalent salinity is **{:.2f}** g NaCl/ 100 mL'.format(sum_eq))
            
def IFT_and_viscosity_calculation():
    st.title('IFT and viscosities')
    st.write('Calculation of IFT and microemulsion viscosity')
    #Three blocks - surfactant properties, HLD properties, system properties
    st.write('**Input of surfactant intrinsic properties**')
    #Select the input type - with SMILES (automatic MWt and tail length calculation) or manual input
    type_of_input = st.selectbox('How do you want to input properties of your surfactants (tail length, molecular weight?', options = (None, 'SMILES', 'Manual'))
    if type_of_input == 'Manual':
        col1, col2 = st.columns(2)
        with col1:
            L = st.number_input('Enter tail length in Å', value = 0)
            st.write('Tail length is {}Å'.format(L))
        with col2:
            MWt = st.number_input('Enter molecular weight', value = 0)
            st.write('Molecular weight is {} '.format(MWt))
    #Input surfactant molecular structure (SMILES) with subsequent calculation of tail length and MWt
    elif type_of_input == 'SMILES':
        smi = st.text_input('Enter SMILES of you surfactant - only surface-active part, without counterion!')
        #SMILES is transformed into molecule object
        mol = Chem.MolFromSmiles(smi)
        #if molecule is valid, MWt and tail length are calculated
        if mol:
            L = len_calculation(mol)
            MWt = MolWt(mol)
            col1, col2 = st.columns(2)
            with col1:
                st.write('Tail length is {}'.format(L))
            with col2:
                st.write('Molecular weight is {} '.format(MWt))
            img = Draw.MolToImage(mol)
            st.image(img, smi)
    #Second block - input of HLD properties
    st.write('**Input of HLD properties**')
    #Selection of surfactant type, which defines HLD calculation
    surfactant_selection = st.selectbox('Type of surfactants', options = ('Ionic', 'Nonionic'))
    #Surfactant type and "is extended" defines the coefficient before kbT for IFT calculation 
    is_extended = st.checkbox('Is surfactant extended')
    if is_extended:
        coeff = 7
    elif (surfactant_selection == 'Ionic') & (not is_extended) :
        coeff = 1
    elif (surfactant_selection == 'Nonionic') & (not is_extended):
        coeff = 4
    #HLD property that is varied in some range and hthe app draws tubes against its value
    prop = st.selectbox('Property to be varied', options = (None, 'Cc', 'EACN', 'Temperature', 'Salinity'), placeholder = None)
    if prop:
        st.write(prop)
        if prop == 'EACN':
            col1, col2, col3 = st.columns(3)
            with col1:
                Cc = st.number_input('Enter Cc', value = 0.00)
                st.write('Cc is {}'.format(Cc))
            with col2:
                Temp = st.number_input('Enter temperature', value = 25.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a range of EACN", -15.00, 15.00, (0.00, 10.00))
            EACN_range = np.linspace(values[0], values[1], 20)
            #calculates HLD range
            HLDs = HLD_calculation(EACN = EACN_range, Cc = Cc, Sal = Sal, Temp = Temp, Type = surfactant_selection)
        #do the same stages if other properties are selected
        elif prop == 'Cc':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Temp = st.number_input('Enter temperature', value = 25.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a range of Cc", -15.00, 15.00, (-5.00, 5.00))
            Cc_range = np.linspace(values[0], values[1], 20)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc_range, Sal = Sal, Temp = Temp, Type = surfactant_selection)
        elif prop == 'Temperature':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Cc = st.number_input('Enter Cc', value = 0.00)
            with col3:
                Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
            values = st.slider("Select a temperature range", 5.00, 95.00, (25.00, 50.00))
            temp_range = np.linspace(values[0], values[1], 20)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc, Sal = Sal, Temp = temp_range, Type = surfactant_selection)
        elif prop == 'Salinity':
            col1, col2, col3 = st.columns(3)
            with col1:
                EACN = st.number_input('Enter EACN', value = 0.00)
                st.write('EACN is {}'.format(EACN))
            with col2:
                Cc = st.number_input('Enter Cc', value = 0.00)
            with col3:
                Temp = st.number_input('Enter temperature', value = 25.00)
            values = st.slider("Select a salinity range, g NaCl/100 mL", 0.00, 100.00, (0.00, 5.00))
            Sal_range = np.linspace(values[0], values[1], 20)
            HLDs = HLD_calculation(EACN = EACN, Cc = Cc, Sal = Sal_range, Temp = Temp, Type = surfactant_selection)
        #Third block - input of system parameters
        st.write('**Input of system paramters**')
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            Vw = st.slider("Enter the volume % of water", 0.01, 0.99, 0.5, step = 0.01)
            Vo = 1 - Vw
        with col2:
             #Enter weight concentration and recalculate to mol/ml
            wt_concentration = st.slider("Select a concentration of surfactant in water in % wt",0.1 , 20.1, step = 0.1)
            if MWt != 0:
                mol_concentration = ((1000 * wt_concentration / 100)/MWt)/1000
            Vs = (wt_concentration/100)*Vw
        with col3:
            interfacial_area = st.number_input('Enter interfacial area in Å2', value = 45)
            st.write('Interfacial area is {} '.format(interfacial_area)) 
        with col4:
            mu_water = st.number_input('Enter water viscosity (cP)', value = 1)
            mu_oil = st.number_input('Enter oil viscosity (cP)', value = 1)
        As =  Vw*mol_concentration*6E23*interfacial_area/1E24  
        Rw_ref = 3*Vw/As
        Ro_ref = 3*Vo/As
        
        
        #xi calculation - manual input (from app) or automatic calculation
        calculation_type = st.selectbox('How do you want to calculate xi', options = (None, 'Automatically', 'Manually'), placeholder = None)
        #Only oil parameters are needed for successful automatic calculation, as surfactant parameters were typed in earlier
        if calculation_type == 'Automatically':
            col1, col2 = st.columns(2)
            with col1:
                density = st.number_input('Enter density of your oil', value = 0.8)
                st.write('Density is {}'.format(density))
            with col2:
                MWt_oil = st.number_input('Enter molecular weight of your oil', value = 100)
                st.write('Molecular weight of the oil is {} '.format(MWt_oil))
            xi = xi_calculation(L, interfacial_area, density, MWt_oil, surfactant_selection)
            st.write('**ξ** is {}'.format(xi))
        elif calculation_type == 'Manually':
            xi = st.number_input('Enter xi', value = 100)
            st.write('**ξ** is {}'.format(xi))

        #if xi value is sucessfully calculated, then other calculations take place        
        if calculation_type:
            #calculation of NACs and emulsion types
            Hs = -(HLDs/L)
            Ros = np.zeros((len(Hs)), dtype = 'float')
            Rws = np.zeros((len(Hs)), dtype = 'float')
            #Calculation of Ro, Rw and filling arrays
            for i,H in enumerate(Hs):
                if H > 0:
                    Rw = Rw_ref
                    Ro = 1/(H+1/Rw)
                else:
                    Ro = Ro_ref
                    Rw = 1/(np.abs(H)+1/Ro)
                Ros[i] = Ro
                Rws[i] = Rw
            #NAC calculation and determination of ME type (Rw > Ro - I, Ro > Rw - II, 1/Ro + 1/Rw < 1/xi - III). Emulsion type is coded into NACs_bin
            NACs = 0.5*((1/Ros)+(1/Rws))
            NACs_bin = (Ros > Rws).astype(int)
            NACs_bin[NACs<1/xi] = 2
            mask = NACs_bin==2
            #Calculation of dispersed phase radii in emulsion
            Rome  = Ros
            Rwme = Rws
            Rwme[mask] = 1/((HLDs[mask]/(2*L))+1/xi)
            Rome[mask] = 1/(1/Rwme[mask] - (HLDs[mask]/L))
            #Calculatuion of oil and water volume in ME
            Vome = Rome*As/3
            Vwme = Rwme*As/3
            Vtot = Vome + Vwme
            
            # Calculation of volume fractions of oil and water in ME
            phio = np.zeros(20)
            phiw = np.zeros(20)
            
            for i in range(20):
                if HLDs[i] < 0:
                    phio[i] = Vome[i]/Vtot[i]
                    phiw[i] = 1 - phio[i]
                else:
                    phiw[i] = Vwme[i]/Vtot[i]
                    phio[i] = 1 - phiw[i]
            # IFT calculation
            #Temperature from Kelvin to Celsius, two states: tempreature is a range or varied HLD parameter (set earlier)
            if prop == "Temperature":
                T = temp_range+273
            else:
                T = Temp + 273
            #Main calculation
            kb = 1.38E-23
            kbT = kb*T
            IFTome = coeff*kbT/(4*3.14*((Ros*1E-10)**2))
            IFTwme = coeff*kbT/(4*3.14*((Rws*1E-10)**2))

            #Viscosity calculation
            #Masking to get emulsion type I (0) (0) and type II (1) and arrays to store viscosities
            mask_I = NACs_bin == 0
            mask_II = NACs_bin == 1
            viscosities_I = np.zeros(sum(mask_I))
            viscosities_II = np.zeros(sum(mask_II))
            #SymPy to solve equations and get Ld, Rd - diameters of rigid rods (droplets) of emulsion - see Keiran, Acosta, see Ind. Engg. Chem. Res, 2010, 49, 3424-3432. 
            from sympy import symbols, solve
            Ld, Rd = symbols('Ld Rd')
            #Solving equations, relating Ha, NAC, Ld, Rd and calculating viscosities for Type I MEs
            for i, (Ha, NAC) in enumerate(zip(Hs[mask_I], NACs[mask_I])):
                system = [2/(Ld+2*Rd) + (Ld/Rd)*(1/(2*Ld+4*Rd)) - Ha/2, (2*Ld + 4 * Rd)/(3*Ld*Rd + (4 * Rd**2)) - NAC]
                solution = solve(system, Ld, Rd)
                for tup in solution:
                    #solution may have negative value
                    if (tup[0]>0) & (tup[1] > 0):
                        Ld, Rd = tup
                mu = mu_water*(1+(4*phio[mask_I][i]*(Ld**2)/(3.14*1*((2*Rd)**2))))
                viscosities_I[i] = mu
            #the same stage for Type II MEs
            for i, (Ha, NAC) in enumerate(zip(Hs[mask_II], NACs[mask_II])):
                system = [2/(Ld+2*Rd) + (Ld/Rd)*(1/(2*Ld+4*Rd)) - Ha/2, (2*Ld + 4 * Rd)/(3*Ld*Rd + (4 * Rd**2)) - NAC]
                solution = solve(system, Ld, Rd)
                for tup in solution:
                    if (tup[0]>0) & (tup[1] > 0):
                        Ld, Rd = tup
                mu = mu_oil*(1+(4*phio[mask_II][i]*(Ld**2)/(3.14*1*((2*Rd)**2))))
                viscosities_II[i] = mu

            
            #Plotting
            if prop == "EACN":
                varied_value_range = EACN_range
            elif prop=='Cc':
                 varied_value_range = Cc_range
            elif prop=='Temperature':
                varied_value_range =  temp_range
            elif prop=='Salinity':
                varied_value_range = Sal_range
            #Plotting IFT and viscosity vs varied HLD parameter (viscosity - only for Type II and I MEs)    
            col1, col2  = st.columns(2)
            with col1:
                fig, ax = plt.subplots(1,1, figsize = (12,12))
                ax.plot(varied_value_range, IFTome*1000, label = 'IFT oil/ME')
                ax.plot(varied_value_range, IFTwme*1000, label = 'IFT water/ME')
                plt.plot(varied_value_range, 500*(IFTome+IFTwme), label  = 'IFT oil/water')
                ax.set_xlabel(prop)
                ax.set_ylabel('IFT, mN/m')
                ax.set_yscale('log')
                ax.grid(visible = True)
                ax.legend()
                st.pyplot(fig)
                st.write('IFT between phases')
            with col2:
                fig, ax = plt.subplots(1,1, figsize = (12,12))
                ax.plot(varied_value_range[mask_I], viscosities_I, label = 'Type I ME viscosity')
                ax.plot(varied_value_range[mask_II], viscosities_II, label = 'Type II ME viscosity')
                ax.set_ylabel('Viscosity, cP')
                ax.set_xlabel(prop)
                ax.legend()
                ax.grid(visible = True)
                st.pyplot(fig)
                st.write('Viscosity of microemulsions')
            
def Ternary_phase_diagram():
    #Optimizer will be needed to calculate end value of phi_o_s and phi_w_s
    from scipy.optimize import minimize_scalar
    import ternary
    st.title('Plotting ternary phase diagrams')
    st.write('Calculation and plotting of O-W-S ternary phase diagrams')
    #Select the input type - with SMILES (automatic MWt and tail length calculation) or manual input
    type_of_input = st.selectbox('How do you want to input properties of your surfactants (tail length, molecular weight?', options = (None, 'SMILES', 'Manual'))       
    #Session state code block is used to prevent immediate plotting with unexisting parameters tight after the app is launched 
    if "L" not in st.session_state:
        st.session_state.L = None
    if type_of_input == 'Manual':
        col1, col2 = st.columns(2)
        with col1:
            L = st.number_input('Enter tail length in Å', value = 0)
            st.session_state.L = L
            st.write('Tail length is {}Å'.format(L))
        with col2:
            MWt = st.number_input('Enter molecular weight', value = 0)
            #Vs_as is parameter, defined by Acosta 
            Vs_as = MWt*10/(6.022*100)
            st.write('Molecular weight is {} '.format(MWt))
    #Input surfactant molecular structure (SMILES) with subsequent calculation of tail length and MWt
    elif type_of_input == 'SMILES':
        smi = st.text_input('Enter SMILES of you surfactant - only surface-active part, without counterion!')
        #SMILES is transformed into molecule object
        mol = Chem.MolFromSmiles(smi)
        #if molecule is valid, MWt and tail length are calculated
        if mol:
            L = len_calculation(mol)
            st.session_state.L = L
            MWt = MolWt(mol)
            col1, col2 = st.columns(2)
            with col1:
                st.write('Tail length is {}'.format(L))
            with col2:
                st.write('Molecular weight is {} '.format(MWt))
            img = Draw.MolToImage(mol)
            st.image(img, smi)
    #Second block - input of HLD properties. Nothing is varied, you just enter properties and calculate HLD.
    #Defines type of the surfactant, which determines HLD equation
    surfactant_selection = st.selectbox('Type of surfactants', options = ('Ionic', 'Nonionic')) 
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        EACN = st.number_input('Enter EACN', value = 0.00)
    with col2:
        Temp = st.number_input('Enter temperature', value = 25.00)
    with col3:
         Sal = st.number_input('Enter salinity, g NaCl/100 mL', value = 1.00)
    with col4:
        Cc = st.number_input('Enter Cc value of your surfactant', value = -1.0)
    #HLD value calculation
    HLD = HLD_calculation(EACN = EACN, Cc = Cc, Sal = Sal, Temp = Temp, Type = surfactant_selection)
    
    #xi calculation - manual input (from app) or automatic calculation
    calculation_type = st.selectbox('How do you want to calculate xi', options = (None, 'Automatically', 'Manually'), placeholder = None)
    #Session state is used to prevent preliminary calculation
    if "xi" not in st.session_state:
        st.session_state.xi = None
    if calculation_type == 'Automatically':
        col1, col2 = st.columns(2)
        #Only oil parameters are needed for successful calculation, as surfactant parameters were typed in earlier
        with col1:
            density = st.number_input('Enter density of your oil', value = 0.8)
            st.write('Density is {}'.format(density))
        with col2:
            MWt_oil = st.number_input('Enter molecular weight of your oil', value = 100)
            st.write('Molecular weight of the oil is {} '.format(MWt_oil))
        xi = xi_calculation(L, interfacial_area, density, MWt_oil, surfactant_selection)
        st.session_state.xi = xi
        st.write('**ξ** is {}'.format(xi))
    elif calculation_type == 'Manually':
        xi = st.number_input('Enter xi', value = 100)
        st.session_state.xi = xi
        st.write('**ξ** is {}'.format(xi))
    col1, col2 = st.columns(2)
    #fso and fsw - fraction of surfactant in oil and water respectively
    with col1:
        fso = st.number_input('Enter surfactant fraction in the oil', value = 0.2)
    with col2:
        fsw = st.number_input('Enter surfactant fraction in the water', value = 0.2)

    #Calculation of ternrary diagram   
    if (st.session_state.L != None) & (st.session_state.xi != None):
        #HLD that defines transition betwee emulsion types
        HLDi_iii =  -2*L/xi
        HLDiii_ii =  2*L/xi
        if HLD <  HLDi_iii:
            Ro_d = -L/HLD
            Rw_d = 1000
        else:
            if HLD > HLDiii_ii:
                Ro_d = 1000
                Rw_d = L/HLD
            else:
                Ro_d = 2/(2/xi - HLD/L)
                Rw_d = 2/(2/xi + HLD/L)   
        Vs_as = MWt*10/(6.022*100)
    
        #Calculates R of dispersed phase if surfactant continuum is taken into account
        def R_calculation(phi_s_fl, fs, R_d):
            Rws = 3 * (1-phi_s_fl+fs*phi_s_fl)*Vs_as/phi_s_fl
            Rw = 3 * (1-phi_s_fl)*Vs_as/phi_s_fl
            Ros = 1/(1/Rws-1/Rw+1/R_d)
            return (Rws, Rw, Ros)
        #objective function, which is used for minimization
        def obj_fn(phi_s_fl, fs, R_d):
            return 1000/np.abs(R_calculation(phi_s_fl, fs, R_d)[2])
        # ME + excess water phase boundary line - 1) find the final value of surfactant water ratio, where Ros -> infinite 2) in the whole range calclulate Ros 3)calculate surfactant/oil raio 3) calculate fractions of surfactant, water and oil
        res = minimize_scalar(obj_fn, args = (fsw, Ro_d), bounds  = (0,1), method='bounded')
        phi_s_aq = np.linspace(0.01, res.x, 30)           
        Ros = R_calculation(phi_s_aq, fsw, Ro_d)[2]
        phi_s_o = 1/(Ros/(3*Vs_as)+1-fso)
        phi_o = (1/phi_s_o-1)/((1/phi_s_o-1)+(1/phi_s_aq-1)+1)
        phi_s = 1/(1+(1/phi_s_o-1)+(1/phi_s_aq-1))
        phi_w = 1 - phi_o - phi_s      
        ME_oil_exc_points = []
        for w, o, s in zip(phi_w, phi_o, phi_s):
            ME_oil_exc_points.append((o, s, w))
        #ME + excess oil phase bundary line - the same procedure as above
        res = minimize_scalar(obj_fn, args = (fso, Rw_d), bounds  = (0,1), method='bounded')
        phi_s_o = np.linspace(0.01, res.x, 30)
        Rws = R_calculation(phi_s_o, fso, Rw_d)[2]
        phi_s_aq = 1/(Rws/(3*Vs_as)+1-fsw)
        phi_o = (1/phi_s_o-1)/((1/phi_s_o-1)+(1/phi_s_aq-1)+1)
        phi_s = 1/(1+(1/phi_s_o-1)+(1/phi_s_aq-1))
        phi_w = 1 - phi_o - phi_s
        ME_water_exc_points = []
        for w, o, s in zip(phi_w, phi_o, phi_s):
            ME_water_exc_points.append((o, s, w))
     #plot everything on the ternary diagram
        fig, ax = plt.subplots(1, 1, figsize = (6,6))
        fig1, tax = ternary.figure(scale = 1, ax = ax)
        
        tax.boundary(linewidth=1.5)
        tax.gridlines(color="black", multiple=0.1)
        tax.gridlines(color="blue", multiple=0.1, linewidth=0.1)
        
        #tax.left_axis_label("Left label $\\alpha^2$", fontsize=12)
        #tax.right_axis_label("Right label $\\beta^2$", fontsize=12)
        tax.top_corner_label('Surfactant', fontsize = 12)
        tax.right_corner_label('Oil', fontsize = 12)
        tax.left_corner_label('Water', fontsize = 12)
        
        tax.bottom_axis_label("Oil, % volume", fontsize = 12)
        tax.plot(ME_oil_exc_points, label = 'ME + escess oil')
        tax.ticks(multiple = 0.2, tick_formats="%.1f", offset = 0.01)
        tax.plot(ME_water_exc_points, label = 'ME + excess water')
        tax.legend()
        tax.clear_matplotlib_ticks()
        tax.get_axes().axis('off')
        st.pyplot(fig, width = 500)



pg = st.navigation([Salts_additives_calculator, HLD_tubes, Calculation_of_xi, Phase_diagram, Volumes_of_phases, IFT_and_viscosity_calculation, Ternary_phase_diagram])

pg.run()
