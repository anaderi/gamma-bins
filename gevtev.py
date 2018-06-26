# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 18:34:49 2018

@author: Admin
"""
import numpy as np
import pandas as pd
from astropy.io import fits

_epsilon = 1e-2
_path_gev = 'data/gll_psc_v16.fit'
_path_tev = 'data/gammacat.fits.gz'
_names_gev = [
    'CLASS1',    
    'RAJ2000',
    'DEJ2000', 
    'GLON', 
    'GLAT', 
    'Variability_Index', 
    'Flux1000', 
    'Flux10000_100000', 
    'Flux1000_3000', 
    'Flux100_300', 
    'Flux3000_10000', 
    'Flux300_1000',  
    'Flux30_100',
    ] 
_names_tev = [
    'classes', 
    'glat', 
    'glon', 
    'morph_pa', 
    'pos_ra',
    'pos_dec',
    'sed_dnde', 
    'sed_dnde_err', 
    'sed_e_ref', 
    'spec_dnde_1TeV', 
    'spec_dnde_1TeV_err', 
    'spec_eflux_1TeV_10TeV', 
    'spec_eflux_1TeV_10TeV_err', 
    'spec_flux_1TeV', 
    'spec_flux_1TeV_crab', 
    'spec_flux_1TeV_crab_err'
    ]
_names_common = [
    'glat',
    'glon',
    'morph_pa',
    'pos_ra',
    'pos_dec',
    'sed_dnde',
    'sed_dnde_err',
    'sed_e_ref', 
    'spec_dnde_1TeV', 
    'spec_dnde_1TeV_err', 
    'spec_eflux_1TeV_10TeV', 
    'spec_eflux_1TeV_10TeV_err', 
    'spec_flux_1TeV', 
    'spec_flux_1TeV_crab', 
    'spec_flux_1TeV_crab_err', 
    'ASSOC_TEV', 
    'Variability_Index', 
    'Flux1000', 
    'Flux10000_100000', 
    'Flux1000_3000', 
    'Flux100_300', 
    'Flux3000_10000', 
    'Flux300_1000', 
    'Flux30_100', 
    'Flux1000', 
    'Flux10000_100000', 
    'Flux1000_3000', 
    'Flux100_300', 
    'Flux3000_10000', 
    'Flux300_1000', 
    'Flux30_100', 
    'CLASS1'
    ]
_gevToTev = {'BLL': 'blazar', 
            'FRSQ': 'frsq', 
            'HMB': 'bin' , 
            'BIN': 'bin', 
            'GAL': 'galaxy', 
            'PSR': 'psr', 
            'PWN': 'pwn', 
            'SNR': 'snr', 
            '': 'unid'}
_tevToGev = {v:k for k, v in _gevToTev.items()}

_d_coord = {'GLON' : 'glon', 
            'GLAT' : 'glat', 
            'RAJ2000' : 'pos_ra', 
            'DEJ2000' : 'pos_dec'}


def cat_gev_tev(path_gev, path_tev):    
    hdul_tev = fits.open(path_tev)
    cat_tev = hdul_tev[1].data
    hdul_gev = fits.open(path_gev)
    cat_gev = hdul_gev[1].data
    return cat_gev, cat_tev

def common(cat_gev, cat_tev, epsilon):
    """This function looks for the same objects in GeV and TeV catalogs
    
    Return: dictionary with corresponding values
    
    feature1(list) - features of GeV 
    feature2(list) - features of TeV 
    epsilon(double) - distance accepted as equivalence
    """
    d = {}

    class_gev = cat_gev['CLASS1']
    class_tev = cat_tev['classes']
    glat_gev = cat_gev['GLAT']
    glat_tev = cat_tev['glat']
    glon_gev = cat_gev['GLON']
    glon_tev = cat_tev['glon']
    
    for i in range(len(glat_gev)):
        for j in range(len(glat_tev)):
            classGeV = class_gev[i]
            start = -1
            try: 
                start = (class_tev[j].find(_gevToTev[classGeV]))
            except KeyError:
                continue
            if (start != -1):
                if ((np.abs(glat_gev[i] - glat_tev[j])/np.abs(glat_gev[i]) < epsilon) and (np.abs(glon_gev[i] - glon_tev[j])//np.abs(glat_tev[j]) < epsilon)) :
                    d.update({j : i})
    return d

def create_common_data(cat_gev, cat_tev, D, namefinal):
    """The fonction adds objects found both in GeV and TeV.
    
    D(dictionary) - a dictionary with repeated objects from TeV and GeV
    namefinal(list) - names of required features to fill in data
    """
    data = []
    k = 0
    for j in D.keys():
        data.append([]) 
        for i in range(15): 
            data[k].append(cat_tev[namefinal[i]][j]) 
        for i in range(15, len(namefinal)): 
            data[k].append(cat_gev[namefinal[i]][D[j]])
        k = k+1  
    df_common = pd.DataFrame(data = data, columns = namefinal)
    return df_common 


def create_only_tev_data(cat_tev, D, name_tev):
    """The fonction adds objects found only in TeV.
    
    D(dictionary) - a dictionary with repeated objects from TeV and GeV
    name_tev(list) - names of TeV columns
    """
    k = 0
    data = []
    
    for j in range(len(cat_tev)):
        if not(j in D.keys()):
            data.append([])
            for i in range(len(name_tev)):
                data[k].append(cat_tev[name_tev[i]][j])
            k = k+1
    df_only_tev = pd.DataFrame(data = data, columns = name_tev)
    df_only_tev = df_only_tev.rename(columns = {'classes' : 'CLASS1'})    
    return df_only_tev


def create_only_gev_data(cat_gev, D, names_gev):
    """The fonction adds objects found only in GeV.
    
    D(dictionary) - a dictionary with repeated objects from TeV and GeV
    name_gev(list) - names of GeV columns
    """
    k = 0
    data = []
    
    for j in range(len(cat_gev)):
        if not(j in D.values()):
            data.append([])
            for i in range(len(names_gev)):
                data[k].append(cat_gev[names_gev[i]][j])
            k = k+1            
    df_only_gev = pd.DataFrame(data = data, columns = names_gev)  
    df_only_gev = df_only_gev.rename(columns = _d_coord)        
    return df_only_gev

def gev_tev_data():
    cat_gev, cat_tev = cat_gev_tev(_path_gev, _path_tev)
    D = common(cat_gev, cat_tev, _epsilon)
    common_data = create_common_data(cat_gev, cat_tev, D, _names_common)
    only_tev_data = create_only_tev_data(cat_tev, D, _names_tev)
    only_gev_data = create_only_gev_data(cat_gev, D, _names_gev)
    return common_data, only_tev_data, only_gev_data

#common_data, only_tev_data, only_gev_data = gev_tev_data()
#print(common_data.head())
#print(common_data.info)
#print(only_tev_data.head())
#print(only_gev_data.head())