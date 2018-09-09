# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 13:13:07 2018

@author: Admin
"""

def get_epsilon():
    return 0.12

def list_xmm_spectra_columns():
    xmm_spectrum_columns = [
            "SC_EP_1_FLUX",
            "SC_EP_2_FLUX",
            "SC_EP_3_FLUX",
            "SC_EP_4_FLUX",
            "SC_EP_5_FLUX",
    ]
    xmm_spectrum_columns = list(map(lambda x: "xmm_" + x, xmm_spectrum_columns))
    return xmm_spectrum_columns


def list_gev_spectrum_columns():
    gev_spectrum_columns = [
            'gev_nuFnu10000_100000',
            'gev_nuFnu1000_3000',
            'gev_nuFnu100_300',
            'gev_nuFnu3000_10000',
            'gev_nuFnu300_1000',
            'gev_nuFnu30_100',
            ]
    return gev_spectrum_columns


def list_s_spectrum_columns():
    s_spectrum_columns = [
            's_FLUX_U', 
            's_FLUX_B', 
            's_FLUX_V', 
            's_FLUX_R', 
            's_FLUX_I',
            's_FLUX_G', 
            's_FLUX_J', 
            's_FLUX_H',
            's_FLUX_K', 
            's_FLUX_u',
            's_FLUX_g', 
            's_FLUX_r',
            's_FLUX_i',    
            's_FLUX_z',
    ]
    return s_spectrum_columns


def list_tev_spectrum_columns():
    tev_spectrum_columns =  [
            'tev_0.3TeV', 
            'tev_1TeV', 
            'tev_3TeV',  
            'tev_10TeV',  
            'tev_30TeV' 
    ]
    return tev_spectrum_columns

def list_s_other_columns():
    s_other_columns = ['s_MAIN_ID',
                       's_class', 
                       's_nan_spectra'  
                 ]
    return s_other_columns

def list_gev_other_columns():
    gev_other_columns = ['gev_1FGL_Name',
                         'gev_CLASS1',
                         ]
    return gev_other_columns

def list_tev_other_columns():
    tev_other_columns = ['tev_fermi_names', 
                         'tev_classes',
                         ]
    return tev_other_columns

def list_xmm_other_columns():
    xmm_other_columns = ['xmm_IAUNAME',
                         'xmm_WEBPAGE_URL',
                         ]
    return xmm_other_columns
