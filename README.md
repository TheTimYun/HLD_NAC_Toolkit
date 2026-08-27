# HLD-NAC toolkit

This toolkit was created to ease the formulation of surfactant-containing products using the HLD-NAC theory. The essence of the theory and supporting materials can be found in ["Practical Surfactants Science"](https://www.stevenabbott.co.uk/practical-surfactants/) by Steven Abbott. 

There are two main parts of the toolkit:
1. Prediction of Cc parameter of surfactant and EACN parameter of oil from their molecular structure using shallow machine learning methods;
2. Calculation of different surfactant-containing system properties with HLD-NAC

## How to use

The best way to use these apps for Python newcomers is:

1. Download and install [Anaconda Navigator](https://www.anaconda.com/download);
2. Using GUI interface of the Anaconda Navigator ("Environments" tab), create new environment ("Create" button) and call it any way you like, for example "HLD_prediction" or "Toolkit_env". Choose Python version 3.11. Newer versions will work likely, too, however, they have not been checked. 
3. Go again to "Home" tab, choose "JupyterLab" app and click "Install";
4. Launch "JupyterLab" app, it will open in your browser;
5. Select necessary folder with streamlit files in the left tab of the JupyterLab app. Launch "Terminal" below the "Other" title;
6. Type in necessary installation command (starting with ``` pip install ```) - see below and press "Enter", packages will be downloaded and installed;
7. Check, that you're working in the right directory. If so, you can run streamlit file you want with one of the following commands in Terminal. 

```
nohup streamlit run streamlit_Cc_EACN.py

nohup streamlit run streamlit_HLD_properties.py
```

And you can enjoy working with HLD-NAC toolkit! If you're an experienced user, you can choose your own way. Of course, you can modify all files to meet your needs (see below).



## Prediction of HLD parameters (streamlit_Cc_EACN.py)

This Streamlit file was created for easy calculation of simple HLD properties of surfactants and oils using Machine Learning methods. Characteristic Curvature (Cc) for surfactants and Equivalent Alkane Carbon Number (EACN) for oils are calculated using Streamlit instrument. If you want to understand, what models are used, how they are obtained and trained, you can go to the corresponding folders **EACN** and **Critical curvature** and dive deeply into the process. Workflow for predicting HLD for surfactants is described in paper [Predicting surfactants characteristic curvature for HLD framework using shallow machine learning methods](https://www.sciencedirect.com/science/article/abs/pii/S0009250926002411) .


First, you need to install required packages in new Anaconda environment. 

```
pip install datamol==0.12.5 molfeat==0.10.1 padelpy molfeat-padel numpy==1.26.4 scikit-learn streamlit
```
To make calculation of PaDEL descriptors possible, Java JRE 8+ should be installed and added to the system PATH.

If you want to use deep learning with Graph Neural Networks for EACN prediction, you can install torch-molecule with the following command

```
pip install torch-molecule
```

Please ensure, that CUDA is installed and configured for use in your system. However,  prediction with DL method dows not add much to the accuracy. 

After launch of Streamlit app, you can choose the property that you want to predict (Cc or EACN) on the left tab. After selection of neccesary property, you will be prompted to enter molecular structure of your molecule (oil or surfactant) in [SMILES](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System) format. After typing the SMILES in, select surfactant type (if Cc is calculated) or if you want to use a SGIR predictor (if EACN is calculated). The calculation will proceed then and take some time. Results of Cc and EACN will be shown below. If your moleule lies outside the applicability domain, a red warning will be shown 

**Please note**, that SMILES of ionic surfactants should be typed in without counterion, only surface-active ion.

## Calculation of system properties (streamlit_HLD_properties.py)

This file contains numerous apps, that failitate the formulation process. Several of them  apps resemble the ones made by [Steven Abbott](https://www.stevenabbott.co.uk/practical-surfactants/) and several of them follow  Excel spreadsheetds, provided by [E. J. Acosta](https://www.researchgate.net/profile/Edgar-Acosta-4/research) All these apps they can be launched offline and use simple cheinformatics methods to skip annoying stages of manual molecular weight and tail length calculation. Required packages are


```
pip install numpy==1.26.4  streamlit rdkit sympy python-ternary 
```

Following apps are contained here: 

1. **CC mixture calulator** - to calculate Cc of surfactants mixture, indicidual Cc value, ratios in mixture and molecular weights/SMILES are required
2. **Salts additives calculator**. It re-calculates water solution with several different salts into units of g NaCl/100 mL. You have to select type of your surfactants and then choose salt (as cation and anion) and its concentration. Equivalent salinity will be shown on the left side, summed equivalent salinity will be shown below;
3. **HLD tubes** - an equivalent of [this app](https://www.stevenabbott.co.uk/practical-surfactants/HLD-Tubes.php). You choose  HLD parameter (Cc, EACN, salinity or temperature), that you want to vary, select its min and max values and enter other parameters as well as type of surfactant. You will get a simple picture, that shows phase behvaiour of your system, depdnding on the varied value;
4. **$\xi$** calculation - equivalent of [this app](https://www.stevenabbott.co.uk/practical-surfactants/xi.php). Calculates $\xi$ parameter from HLD-NAC equations;
5. **Phase diagram** - equivalent of [this app](https://www.stevenabbott.co.uk/practical-surfactants/fishtail.php). First, you should enter properties of the surfactant (manually or swith SMILES) and select its type, then enter HLD properties (again with one varied value and with other being fixed), then enter properties of the system - oil-water ratio, range of surfactant **weight** concentration (min and max value) and $\xi$. $\xi$ can be entered manually (from previous app) or calculated directly here, using molecular weight and density of oil. If all stages are done, you will get a phase diagram, that relates concentration, varied HLD parameter and type of microemulsion;
6. **Volumes of phases** - equivalent of [this app](https://www.stevenabbott.co.uk/practical-surfactants/phase-volumes.php). Again, you must enter intrinsic properties of surfactant, HLD parameters, oil-water ratio and **weight** concetration of surfactant, $\xi$. Then you will get several plots, that shows relative volumes of phases and boundaries positions of microemulsion in the system;
7. **IFT and viscosity calculation** - requires the same input, as the previous apps. Calculates the IFT between microemulsions and excess phase and viscosity of Type I and Type II microemulsions, using [this file](https://www.researchgate.net/publication/369655812_Formulation_Engineering_with_the_Hydrophilic-Lipophilic-Difference_HLD_and_Net-Average_Curvature_NAC_HLD-NAC_Tutorial?_sg%5B0%5D=Y5p97dX2clIeXcf1I69AInRqi3P5LI0w4D44wfmxJfKCAAd03GvWqfKmnF-nf6AgV6m6IhqeqOi5W_LBE1HcMJCI9p-grBjWCJ2x0Z1w.OjfuzotsmK3wXyYMpMyOKFCe80U048pOqlQgRX8iLAq-mKoAv78Hd0ut-pXWOJESkqT3ZRYsUhUz860ocmzGSw&_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6ImhvbWUiLCJwYWdlIjoicHJvZmlsZSIsInBvc2l0aW9uIjoicGFnZUNvbnRlbnQifX0) and [this paper](https://pubs.acs.org/doi/abs/10.1021/ie9013106).
8. **Ternary phase diagrams** - builds ternary diagram of surfactant-oil-water system. After you choose surfactant intrinsic properties, HLD parameters (non of them is varied here!), $\xi$, you choose surfactant fraction in oil and water. These parameters are somewhat tricky and right now they should be considered as adjustable parameters of the system - and then you get  ternary phase diagam in S-O-W coordinates with two curves (microemulsion and excess oil, microemulsion and excess water). 

Two **very important rules** must be followed for the correct work of apps:
1. SMILES of ionic surfactants should be typed in without counterion, only surface-active ion;
2. All actions must be executed consequtively. Please don't start with HLD parameter, skipping entering molecular weight. Otherwise, an error will occur. In this case, just step back, enter all the values consecutively, and highly likely, the error will disappear. 

If you want to dive deeply into coding of every app, you are welcome to the Jupyter Notebook files in HLD-NAC folder. These notebooks perform the same actions, however, they cannot operate in a live mode, as Streamlit app does. 

If you have any questions, suggestions of complaints, you can always write me to timyun96@gmail.com 
