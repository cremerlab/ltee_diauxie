#%%
import numpy as np 
import pandas as pd 
import cmdstanpy 
import arviz as az

# Load the data and restrict to what we care about
calib_data = pd.read_csv('../../data/calibration/2021-04-05_hplc_calibration/processed/2021-04-05_NC_DM_calibration_relative_areas.csv')

# restrict the  data to the relevant quantites
calib_data = calib_data[(calib_data['buffer_base']=='DM') & 
                        (calib_data['compound'].isin(['glucose', 'acetate']))
                        ][['carbon_conc_mM', 'rel_area_phosphate', 'compound']]

