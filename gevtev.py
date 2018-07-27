# -*- coding: utf-8 -*-
"""
Created on Sun Jun 17 18:34:49 2018

@author: Admin
"""
import numpy as np
import pandas as pd
from astropy.io import fits

_epsilon = 1.2e-1
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
    'ASSOC_TEV',
    '2FGL_Name',
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
    'spec_flux_1TeV_crab_err',
    'tevcat_name',
    'gamma_names',
    'other_names',
    'fermi_names',
    'common_name',
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

_interesting_types = []

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

def create_common(cat_gev, cat_tev, epsilon):
    """
    Returns 2 vectors for gev and tev respectively.
    C_associations_gev coordinate is equal 
    
    -1, if no similar objects in TEV catalog found
    i, where i is a corresponding index of a similar object from TEV
        to the object from GEV with an index equal to number 
        of the coordinate.
    
    Arguments:
    cat_gev -- rec array with GEV data
    cat_tev -- rec array with TEV data
    epsilon - threshold for similarity
    
    Returns:
    C_associations_gev - numpy array (n,)
    C_associations_tev - numpy array (m, )
    
    n - number of examples in GEV
    m - number of examples in TEV
    """
    
    glat_gev = cat_gev['GLAT']
    glat_tev = cat_tev['glat']
    glon_gev = cat_gev['GLON']
    glon_tev = cat_tev['glon']
    
    glat_dif_matrix = np.dot(np.vstack((glat_gev, -np.ones_like(glat_gev))).T,
                      np.vstack((np.ones_like(glat_tev), glat_tev)))
    glon_dif_matrix = np.dot(np.vstack((glon_gev, -np.ones_like(glon_gev))).T,
                      np.vstack((np.ones_like(glon_tev), glon_tev))) 
    pairs_matrix = np.logical_and(np.abs(glat_dif_matrix) < epsilon,
                                  np.abs(glon_dif_matrix) < epsilon)
    return pairs_matrix

def create_pandas_frames(cat_gev, cat_tev):   
    """
    Creates pandas dataframes with the same values as in cat_gev 
    and cat_tev, adding to columns names "gev_" and "tev_" 
    respectively.
    
    Arguments:
    cat_gev -- rec array with GEV data
    cat_tev -- rec array with TEV data
    
    Returns:
    data_gev -- pandas DataFrame with GEV data
    data_tev -- pandas DataFrame with TEV data
    """
    #for i in range(len(cat_tev.dtype.names)):
    #    print(cat_tev.dtype.names[i])
    #    print(cat_tev.tolist()[0][i])
    data_gev = pd.DataFrame.from_records(cat_gev.tolist(), columns=cat_gev.dtype.names)
    gev_match_names = {}
    for i in data_gev.columns:
        gev_match_names.update({i : "gev_" + i})
    data_gev = data_gev.rename(columns = gev_match_names)

    data_tev = pd.DataFrame.from_records(cat_tev.tolist(), columns=cat_tev.dtype.names)
    tev_match_names = {}
    for i in data_tev.columns:
        tev_match_names.update({i : "tev_" + i})
    data_tev = data_tev.rename(columns = tev_match_names)
    return data_gev, data_tev

def create_common_data(data_gev, data_tev, pairs_matrix):
    """
    The fonction adds objects found both in GeV and TeV.
    
    Arguments:
    data_gev -- pandas DataFrame with GEV data
    data_tev -- pandas DataFrame with TEV data
    pairs matrix - (n, m) associations
    
    n - number of examples in GEV
    m - number of examples in TEV
    
    Returns:
    pd_common_gevtev - pandas DataFrame with all chosen columns 
    from GEV and TEV
    """
    vector_association = np.where(np.sum(pairs_matrix, axis=0) > 0)[0]
    pd_common_gevtev = pd.DataFrame()
    for i in vector_association:
        data_gev_join = (pairs_matrix[:, i] > 0)*(i + 1) - 1
        data_gev["join"] = data_gev_join
        pd_common_gevtev0 = pd.merge(data_gev, data_tev, right_index=True, left_on='join', how='inner')
        if (len(pd_common_gevtev)):
            pd_common_gevtev = pd_common_gevtev.append(pd_common_gevtev0)
            del pd_common_gevtev["join"]
        else:
            pd_common_gevtev = pd_common_gevtev0.copy()
    array_non_duplicate = ['tev_glon', 'gev_GLAT', 'gev_GLON', 'tev_glat', 'gev_CLASS1', 'tev_classes']
    pd_common_gevtev = pd_common_gevtev.drop_duplicates(array_non_duplicate)
    pd_common_gevtev = pd_common_gevtev.reset_index()
    #df_common = pd.DataFrame(data = data, columns = namefinal)
    return pd_common_gevtev

def create_only_tev_data(data_tev):
    """The fonction adds objects found only in TeV.
    """
    #data_tev['join'] = C_associations_tev
    #data_only_tev = data_tev[data_tev['join'] < 0]
    #del data_only_tev['join']
    data_only_tev = data_tev
    data_only_tev = data_only_tev.reset_index()
    return data_only_tev

def create_only_gev_data(data_gev):
    """The fonction adds objects found only in GeV.
    """
    #data_gev['join'] = C_associations_gev
    #data_only_gev = data_gev[data_gev['join'] < 0]
    #del data_only_gev['join']
    data_only_gev = data_gev
    data_only_gev = data_only_gev.reset_index()
    return data_only_gev


def compare_gev_tev_data(epsilon):
    """
    The fonction returns common objects for GEV and TEV,
    only GEV and only TEV objects.
    
    Returns:
    common_data - pandas DataFrame with common for GEV and TEV objects 
    only_tev_data - pandas DataFrame of only TEV objects 
    only_gev_data - pandas DataFrame of only GEV objects 
    """
    cat_gev, cat_tev = cat_gev_tev(_path_gev, _path_tev)
    pairs_matrix = create_common(cat_gev, cat_tev, epsilon)
    data_gev, data_tev = create_pandas_frames(cat_gev, cat_tev)
    common_data = create_common_data(data_gev, data_tev, pairs_matrix)
    only_tev_data = create_only_tev_data(data_tev)
    only_gev_data = create_only_gev_data(data_gev)
    return common_data, only_tev_data, only_gev_data


#print(only_tev_data.head())
#print(only_gev_data.head())