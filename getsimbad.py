import astropy
import astropy.units as u
import pandas as pd
import numpy as np
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from definitions import *

def change_spectra_columns_units(s_data):
    s_spectrum_columns =  list_s_spectrum_columns()
    if len(s_spectrum_columns) > 0:
        s_data[s_spectrum_columns] = 10 ** (-s_data[s_spectrum_columns])
    return s_data

def get_simbad_data() :
	#queries for required objects
    qry_Be = ("sptypes=Be")
    qry_O = ("sptypes=O")
    qry_B = ("sptypes=B")
    qry_plsr = ("sptypes=plsr")
    qry_B0Ve = ("otypes=HMXB")
    qry_BO = ("sptypes=B0")


    #adding spectral information
    Simbad.add_votable_fields('flux(U)',
                              'flux(B)',
                              'flux(V)',
                              'flux(R)',
                              'flux(I)',
                              'flux(G)',
                              'flux(J)',
                              'flux(H)',
                              'flux(K)',
                              'flux(u)',
                              'flux(g)',
                              'flux(r)',
                              'flux(i)',
                              'flux(z)')
	
    astropy.utils.data.conf.remote_timeout = 1000
    #tables with these objects
    table_Be = Simbad.query_criteria(qry_Be)
    table_plsr = Simbad.query_criteria(qry_plsr)
    table_O = Simbad.query_criteria(qry_O)
    table_B = Simbad.query_criteria(qry_B)
    table_B0Ve = Simbad.query_criteria(qry_B0Ve)
    table_BO = Simbad.query_criteria(qry_BO)
    print(len(table_BO))
    print(len(table_O))

	#create pandas tables and add column "class"
    simbad_plsr = table_plsr.to_pandas()
    simbad_plsr["class"] = "plsr"
    simbad_Be = table_Be.to_pandas()
    simbad_Be["class"] = "Be"
    simbad_B = table_B.to_pandas()
    simbad_B["class"] = "B"
    simbad_O = table_O.to_pandas()
    simbad_O["class"] = "O" 
    simbad_B0Ve = table_B0Ve.to_pandas()
    simbad_B0Ve["class"] = "B0Ve" 
    simbad_BO = table_BO.to_pandas()
    simbad_BO["class"] = "B0Ve" 
    for i in simbad_O.columns:
        print(i)

	#concatenating tables into one
    simbad_tables = [simbad_plsr, 
                     simbad_Be, 
                     simbad_B, 
                     simbad_O, 
                     simbad_BO, 
                     simbad_B0Ve]
    s_data = pd.concat(simbad_tables)

	#adding columns containing glat and glon(takes 5 minutes)
    pos_ra_simbad = s_data['RA'].values
    pos_dec_simbad = s_data['DEC'].values
    glat_simbad = np.zeros(len(pos_ra_simbad))
    glon_simbad = np.zeros(len(pos_ra_simbad))
    for i in range(len(pos_ra_simbad)):
       if i % 100 == 0:
            print(str(i) + '/' + str(len(pos_ra_simbad)))
       if (str(pos_ra_simbad[i]) != ''):
            c_icrs = SkyCoord(pos_ra_simbad[i], pos_dec_simbad[i], unit=(u.hourangle, u.deg), frame='icrs')
            glat_simbad[i] = c_icrs.galactic.b.deg
            glon_simbad[i] = c_icrs.galactic.l.deg
       else:
            glat_simbad[i] = None
            glon_simbad[i] = None
    s_data['glat'] = glat_simbad
    s_data['glon'] = glon_simbad


	#renaming columns : adding "s_" to the names
    s_match_names = {}
    for i in s_data.columns:
        s_match_names.update({i : "s_" + i})
    s_data = s_data.rename(columns = s_match_names)
    s_data = change_spectra_columns_units(s_data)
    
    s_data.to_csv("data/simbad.txt", sep='\t', encoding='utf-8')
    return s_data