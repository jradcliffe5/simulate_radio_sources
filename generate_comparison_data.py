import numpy as np
import pandas as pd
import scipy
import os
from VLBI_analysis_functions import *


### Inputs ###
run_BANE = False
get_blobcat_params = True

if run_BANE == True:
	for i in os.listdir('./'):
		if i.endswith('casa.fits'):
			os.system('BANE --grid 50 50 %s' % i)

if get_blobcat_params == True:
	full_cat = pd.read_csv('target_uvsub_2_CASA_uv_delt_SN100.0_input_model.csv')
	full_cat = full_cat[['x','y','x_deg','y_deg']]
	full_cat.to_csv('int.csv')
	for i in os.listdir('./'):
		if (i.startswith('catalogue_BLOBCAT_')):
			SN = i.split('catalogue_BLOBCAT_target_uvsub_2_CASA_uv_delt_SN')[1].split('.image.csv')[0]
			postfix = i.split('catalogue_BLOBCAT_')[1].split('.csv')[0]

			master_cat = pd.read_csv('target_uvsub_2_CASA_uv_delt_SN%s_input_model.csv'%SN)
			blobcat_cat = pd.read_csv('catalogue_BLOBCAT_target_uvsub_2_CASA_uv_delt_SN%s.image.csv'%SN)
			df = match_catalogues(cat1=master_cat,cat2=blobcat_cat,\
							 RA1='x_deg',Dec1='y_deg',\
							 RA2='RA_p_%s'%postfix,Dec2='Dec_p_%s' %postfix,\
							 unit1=('deg','deg'),unit2=('deg','deg'),\
							 name1='master_cat',name2='%s'%SN,\
							 columns1=master_cat.columns.values,columns2=blobcat_cat.columns.values,\
							 distance=1*u.arcsecond,keep_cat1=True,keep_cat2=True)
			temp_cat = pd.DataFrame({'S_p_mod_SN%s'%SN:master_cat['mode_flux'], 'S_p_SN%s'%SN:df['S_p_%s'%postfix], 'S_i_SN%s'%SN:df['S_int_%s'%postfix]})
			full_cat = pd.concat([full_cat,temp_cat],axis=1)
			#full_cat['S_p_mod_SN%s'%SN] = master_cat['mode_flux']
			#full_cat['S_p_SN%s'%SN] = df['S_p_%s'%postfix]
			#full_cat['S_i_SN%s'%SN] = df['S_int_%s'%postfix]
			del master_cat, blobcat_cat
	full_cat.to_csv('test.csv')

#catalogue_BLOBCAT_target_uvsub_2_CASA_uv_delt_SN10.0.image.csv
