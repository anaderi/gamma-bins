# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 13:13:07 2018

@author: Admin
"""
import numpy as np
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord

def get_epsilon():
    return 0.01

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

def change_units_in_radec(glat, glon): 
    ra = np.zeros((len(glat), 1))
    dec = np.zeros((len(glat), 1))
    for i in range(len(glat)):
        if (str(glat[i]) != ''):
            c_icrs = SkyCoord(l = glat[i]*u.degree, b = glon[i]*u.degree,  frame='galactic')
            c_new = c_icrs.transform_to('fk5')
            ra[i][0] = c_new.ra.deg
            dec[i][0] = c_new.dec.deg
        else:
            ra[i][0] = None
            dec[i][0] = None
    return ra, dec

def angsep(ra1, dec1, ra2, dec2):
    SEP = np.cos(dec1 * np.pi / 180) * np.cos(dec2 * np.pi / 180) * np.cos((ra1 - ra2) * np.pi / 180)
    SEP += np.sin(dec1 * np.pi / 180) * np.sin(dec2 * np.pi / 180) #returns values between 0 and pi radians
    SEP = np.arccos(SEP)
    return SEP * 180. / np.pi

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
                         'gev_GLON', 
                         'gev_GLAT', 
                         ]
    return gev_other_columns

def list_tev_other_columns():
    tev_other_columns = ['tev_fermi_names', 
                         'tev_classes',
                         'tev_glat', 
                         'tev_glon', 
                         ]
    return tev_other_columns

def list_xmm_other_columns():
    xmm_other_columns = ['xmm_IAUNAME',
                         'xmm_WEBPAGE_URL',
                         ]
    return xmm_other_columns
