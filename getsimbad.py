import astropy
import astropy.units as u
import pandas as pd
import numpy as np
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse_single_table, parse

def simbad_data() :
	#queries for required objects
    qry_Be = ("sptypes=Be")
    qry_O = ("sptypes=O")
    qry_B = ("sptypes=B")
    qry_plsr = ("sptypes=plsr")

	#tables with these objects
    table_Be = Simbad.query_criteria(qry_Be)
    table_plsr = Simbad.query_criteria(qry_plsr)
    table_O = Simbad.query_criteria(qry_O)
    table_B = Simbad.query_criteria(qry_B)

	#create pandas tables and add column "class"
    simbad_plsr = table_plsr.to_pandas()
    simbad_plsr["class"] = "plsr"
    simbad_Be = table_Be.to_pandas()
    simbad_Be["class"] = "Be"
    simbad_B = table_B.to_pandas()
    simbad_B["class"] = "B"
    simbad_O = table_O.to_pandas()
    simbad_O["class"] = "O" 

	#concatenating tables into one
    simbad_tables = [simbad_plsr, simbad_Be, simbad_B, simbad_O]
    s_data = pd.concat(simbad_tables)

	#adding columns containing glat and glon(takes 5 minutes)
    pos_ra_simbad = s_data['RA'].values
    pos_dec_simbad = s_data['DEC'].values
    glat_simbad = np.zeros(len(pos_ra_simbad))
    glon_simbad = np.zeros(len(pos_ra_simbad))
    for i in range(len(pos_ra_simbad)):
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
    
    s_data.to_csv("data/symbad.txt", sep='\t', encoding='utf-8')

s_data = pd.read_csv("data/symbad.txt", sep='\t', encoding='utf-8')
print(s_data.columns)