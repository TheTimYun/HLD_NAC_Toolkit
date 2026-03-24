# Calculation of EACN and Cc for oils and surfactants

This Streamlit file was created for easy calculation of simple HLD properties of surfactants and oils using Machine Learning methods. Characteristic Curvature (Cc) for surfactants and Equivalent Alkane Carbon Number (EACN) are calculated using Streamlit instrument.

First, you need to install required packages in new Anaconda environment (I actually used Anaconda Navigator with pre-installed packages)

```
pip install datamol==0.12.5 molfeat==0.10.1 padelpy molfeat-padel numpy==12.6.4 scikit-learn streamlit
```

After installation, you can run streamlit file

```
nohup streamlit run streamlit_Cc_EACN.py
```

If you want to use deep learning with Graph Neural Networks for EACN prediction, you can install torch-molecule with the following command

```
pip install torch-molecule
```

Please ensure, that CUDA is installed and configured for use in your system. However, accuracy of EACN prediction with DL method is on the similar leverl with shallow methods.

If you want to get, what models are used, how they are obtained and trained, you can go to the corresponding folders **EACN** and **Critical curvature** and dive deeply into the process. Workflow for preidcting HLD for surfactants is described in paper [Predicting surfactants characteristic curvature for HLD framework using shallow machine learning methods](https://www.sciencedirect.com/science/article/abs/pii/S0009250926002411) in Chemical Engineering Science.

If you have any questions, you can always contact me on timyun96@gmail.com 
