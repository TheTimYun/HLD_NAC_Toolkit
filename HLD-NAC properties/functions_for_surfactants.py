import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from rdkit import Chem
import rdkit


def draw_type_3(x, y):
    """
    Drawing type III emulsion tube. x, y - lower left corner
    
    """ 
    return (patches.Rectangle((x,y), width = 4, height = 4, facecolor = 'cyan', edgecolor = 'black'), 
    patches.Rectangle((x,y+4), width = 4, height = 2, facecolor = 'white', edgecolor = 'black'),
    patches.Rectangle((x,y+4+2), width = 4, height = 4, facecolor = 'yellow', edgecolor = 'black'))

def draw_type_2(x, y):
    """
    Drawing type II emulsion tube. x, y - lower left corner
    """ 
    return (patches.Rectangle((x,y), width = 4, height = 6, facecolor = 'cyan', edgecolor = 'black'), 
    patches.Rectangle((x,y+6), width = 4, height = 4, facecolor = 'white', edgecolor = 'black'))


def draw_type_1(x, y):
    """
    Drawing type I emulsion tube. x, y - lower left corner
    """ 
    return (patches.Rectangle((x,y), width = 4, height = 6, facecolor = 'white', edgecolor = 'black'), 
    patches.Rectangle((x,y+6), width = 4, height = 4, facecolor = 'yellow', edgecolor = 'black'))


def draw_tubes(HLDs, variables, variable_name):
    """
    Drawing a row of tubes. Variable - the range of HLD parameter, which is varied (np.array), variable_name - it's name.
    """ 
    fig, ax = plt.subplots(1,1, figsize = (12,5))
    
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 20)
    x = 1
    y = 3
    for HLD, variable in zip(HLDs, variables):
        if (HLD > -0.5) & (HLD < 0.5):
            for patch in draw_type_3(x,y):
                ax.add_patch(patch)
            ax.text(x = x+0.2, y = 16, s = 'Type III')
            ax.text(x = x+0.2, y = 14, s = 'HLD = {:.2g}'.format(HLD))
            ax.text(x = x+0.2, y = 1, s = '{} = {:.2g}'.format(variable_name, variable))
        elif (HLD >= 0.5):
            for patch in draw_type_2(x,y):
                ax.add_patch(patch)
            ax.text(x = x+0.2, y = 16, s = 'Type II')
            ax.text(x = x+0.2, y = 14, s = 'HLD = {:.2g}'.format(HLD))
            ax.text(x = x+0.2, y = 1, s = '{} = {:.2g}'.format(variable_name, variable))
        elif (HLD <= 0.5):
            for patch in draw_type_1(x,y):
                ax.add_patch(patch)
            ax.text(x = x+0.2, y = 16, s = 'Type I')
            ax.text(x = x+0.2, y = 14, s = 'HLD = {:.2g}'.format(HLD))
            ax.text(x = x+0.2, y = 1, s = '{} = {:.2g}'.format(variable_name, variable))
        x += 6
    
    #ax.set_axis_off()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    return fig

def draw_tubes_volumes(boundary_1, boundary_2, fr_o, fr_w, fr_me,  range, prop, HLDs):
    """
    Drawing a row of tubes with relative voleume of phases
    Boundary 1, boundary 2 - np.arrays with boundaries of ME
    fr_o, fr_w, fr_me - np.arrays with volume fractions of eacch component 
    range - np.array with the range of varied HLD property
    prop - name of HLD property
    HLDs - np.array with HLD values
    """ 
    fig, ax = plt.subplots(1,1,figsize = (40,15))
    ax.set_xlim(2, 120)
    ax.set_ylim(-.2, 1.60)
    
    x, y = 4, 0
    for  bound_1, bound_2, o, w, me, var, HLD  in zip(boundary_1, boundary_2, fr_o, fr_w, fr_me, range, HLDs):
        if bound_2 == 0:
            rec_1 = patches.Rectangle(xy = (x,y), width = 2, height = bound_1, facecolor = 'blue', edgecolor = 'black')
            rec_2 = patches.Rectangle(xy = (x,y+bound_1), width = 2, height = 1 - bound_1, facecolor = 'yellow', edgecolor = 'black')
            ax.add_patch(rec_1)
            ax.add_patch(rec_2)
            ax.text(x, 1.3, 'HLD=\n{:.2f}'.format(HLD), fontsize = 20)
            ax.text(x, 1.1, '{:.2f}o\n{:.2f}w'.format(o, w), fontsize = 20)
            ax.text(x, -0.1, '{:.2f}'.format(var), fontsize = 20)
            x += 5
            
        if (bound_2 != 0) & (bound_1 != 1):
            rec_1  = patches.Rectangle(xy = (x,y), width = 2, height = bound_2, facecolor = 'cyan', edgecolor = 'black')
            rec_2 = patches.Rectangle(xy = (x,y+bound_2), width = 2, height = bound_1 -bound_2, facecolor = 'white', edgecolor = 'black')
            rec_3 = patches.Rectangle(xy = (x,y+bound_1), width = 2, height = 1 - bound_1, facecolor = 'yellow', edgecolor = 'black')
            ax.add_patch(rec_1)
            ax.add_patch(rec_2)
            ax.add_patch(rec_3)
            ax.text(x, 1.3, 'HLD=\n{:.2f}'.format(HLD), fontsize = 20)
            ax.text(x, 1.1, '{:.2f}o\n{:.2f}w\n{:.2f}me'.format(o, w, me), fontsize = 20)
            ax.text(x, -0.1, '{:.2f}'.format(var), fontsize = 20)
            x += 5
        if (bound_2 != 0) & (bound_1 == 1):
            rec_1 = patches.Rectangle(xy = (x,y), width = 2, height = bound_2, facecolor = 'cyan', edgecolor = 'black')
            rec_2 = patches.Rectangle(xy = (x,y+bound_2), width = 2, height = 1 - bound_2, facecolor = 'orange', edgecolor = 'black')
            ax.add_patch(rec_1)
            ax.add_patch(rec_2)
            ax.text(x, 1.3, 'HLD=\n{:.2f}'.format(HLD), fontsize = 20)
            ax.text(x, 1.1, '{:.2f}o\n{:.2f}w'.format(o, w), fontsize = 20)
            ax.text(x, -0.1, '{:.2f}'.format(var), fontsize = 20)
            x += 5
        
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.text(50, -0.3, prop, fontsize = 28)
        
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Volume fraction')
    return fig

def HLD_calculation(EACN, Cc, Sal, Temp, Type, Cor):
    """
    Calculation of HLD for the selected surfactant
    """
    if Type == 'Ionic':
        HLD = Cc - 0.17 * EACN - 0.01 * (Temp - 25) + np.log(Sal + Cor)
    elif Type == 'Nonionic':
        HLD = Cc - 0.17 * EACN + 0.06 * (Temp - 25) + 0.13 * Sal
    elif Type == 'Zwitterionic':
        HLD = Cc - 0.17 * EACN - 0.00 * (Temp - 25) + 0.13 * Sal
    return HLD

def len_calculation(mol):
    """
    Calculation of tail length of surfactant. Requires "mol" RDKit object
    """
    
    #Structural fragments, described as SMARTS
    p1 = '[CX4H3][CX4H2]'
    p2 = '[CX4H3][CX4H1]([#6])[#6]'
    secondary = '[CX4H2]'
    tertiary_2 = '[CH3][CH1;X4]([CH3])'
    tertiary = '[CH1;X4]'
    benzene = 'c1ccccc1'
    EO = 'O[CH2][CH2]O'
    PO = 'OC(C)CO'
    #getting the numbers of each fragment
    n_tertiary = len(mol.GetSubstructMatches(Chem.MolFromSmarts(tertiary)))
    n_tertiary_2 = len(mol.GetSubstructMatches(Chem.MolFromSmarts(tertiary_2)))
    n_PO = len(mol.GetSubstructMatches(Chem.MolFromSmarts(PO)))
    n_EO = len(mol.GetSubstructMatches(Chem.MolFromSmarts(EO)))
    n_p1 = len(mol.GetSubstructMatches(Chem.MolFromSmarts(p1)))
    n_p2 = len(mol.GetSubstructMatches(Chem.MolFromSmarts(p2)))
    n_secondary = len(mol.GetSubstructMatches(Chem.MolFromSmarts(secondary)))
    n_benzene = len(mol.GetSubstructMatches(Chem.MolFromSmarts(benzene)))
    n_tertiary_non_2 = n_tertiary-n_tertiary_2-n_PO
    #Calculation of tail length, according to Abbott
    number_of_carbons = (n_p1 + n_p2 +n_secondary)+6.9*n_tertiary_2 + 0.71 * n_tertiary_non_2 +3.5 * n_benzene
    Lc = 1.5 + 1.265 * number_of_carbons
    L = 1.3 * Lc
    return L

def xi_calculation(L, interfacial_area, density, MWt, surfactant_type):
    """
    Calculates xi value for ionic and nonionics
    L - tail length - surfactant property (Angstroms)
    Interfacial area - surfactant property (squared Angstroms)
    Density - oil property, g/cc
    MWt - oil property, g/mol
    """
    if surfactant_type == 'Ionic':
        OF = 13 * (L/interfacial_area)/((MWt/density)**(1/3))
        xi = 40 * np.exp(1.2*OF)
        return xi
    elif surfactant_type == 'Nonionic':
        OF = 13 * (L/interfacial_area)/((MWt/density)**(1/3))
        xi = 0.4 * np.exp(6.6*OF)
        return xi

    