#%%
import numpy as np 
import pandas as pd 
import cmdstanpy
import arviz as az

# Load the datasets 
calib_data = pd.read_csv('../../data/calibration/2021-04-05_hplc_calibration/processed/2021-04-05_NC_DM_calibration_relative_areas.csv')
glucose_data = pd.read_csv('../../data/metabolite_turnover/2021-04-04_REL606_glucose_turnover/processed/2021-04-04_REL606_glucose_turnover_relative_areas.csv')
acetate_data = pd.read_csv('../../data/metabolite_turnover/2021-04-27_REL606_acetate_turnover/processed/2021-04-27_REL606_acetate_turnover_relative_areas.csv')

# restrict the  data to the relevant quantites
calib_data = calib_data[(calib_data['buffer_base']=='DM') & 
                        (calib_data['compound'].isin(['glucose', 'acetate']))
                        ][['carbon_conc_mM', 'rel_area_phosphate', 'compound']]
glucose_data = glucose_data[glucose_data['compound'].isin(['glucose', 'acetate'])
                            ][['replicate', 'od_600nm', 'compound', 'rel_area_phosphate',
                               'date', 'carbon_source']]
acetate_data = acetate_data[acetate_data['compound']=='acetate'
                            ][['replicate', 'od_600nm', 'compound', 'rel_area_phosphate',
                               'date', 'carbon_source']]

# Merge the turnover measurements and save
turnover = pd.concat([glucose_data, acetate_data], sort=False)
turnover.to_csv('../../data/collated_turnover_measurements.csv', index=False)
# %%
# Load and compile the inferrential model
model = cmdstanpy.CmdStanModel(stan_file='../stan/hierarchical_yield_coefficient.stan')

# %%

# Define the percentiles to compute
percs = [(2.5, 97.5), (12.5, 87.5), (25, 75), (37.5, 62.5), (47.5, 52.5)]
perc_labels = [95, 75, 50, 25, 5]

# Group by strain, carbon source, and then compound
param_samples, conc_samples = [], []
param_summary, conc_summary = pd.DataFrame([]), pd.DataFrame([]) 
yield_fits, calib_fits = [], []
for g, d in turnover.groupby(['carbon_source', 'compound']):
    # Get the correct calibration data
    calib = calib_data[calib_data['compound']==g[1]]
    
    # Define the data dictionary
    data_dict = {'J':d['replicate'].max(),
                 'N_yield': len(d),
                 'N_calib': len(calib),
                 'idx': d['replicate'].values.astype(int),
                 'calib_conc': calib['carbon_conc_mM'].values.astype(float),
                 'calib_rel_areas':calib['rel_area_phosphate'].values.astype(float),
                 'optical_density':d['od_600nm'].values.astype(float),
                 'yield_rel_areas':d['rel_area_phosphate'].values.astype(float)}

    # Sample the inferrential model
    samps = model.sample(data=data_dict)
    samps = az.from_cmdstanpy(samps)
    samps = samps.posterior.to_dataframe().reset_index()

    # Tidy low-level parameters
    params = samps[['yield_inter_dim_0', 'yield_inter', 'yield_slope']]
    params.drop_duplicates(inplace=True)
    params = params.melt('yield_inter_dim_0', var_name='parameter')
    params.rename(columns={'yield_inter_dim_0':'level'}, inplace=True)
    params['level'] = [f'replicate {i+1}' for i in params['level'].values]

    # Tidy the hyper parameters
    hyperparams = samps[['yield_inter_mu',  'yield_slope_mu']]
    hyperparams.drop_duplicates(inplace=True)
    hyperparams['level'] = 'hyperparameter'
    hyperparams = hyperparams.melt('level', var_name='parameter')
    hyperparams.loc[hyperparams['parameter']=='yield_slope_mu', 'parameter'] = 'yield_slope'
    hyperparams.loc[hyperparams['parameter']=='yield_inter_mu', 'parameter'] = 'yield_inter'

    # Tidy the calibration samples
    calib_params = samps[['calib_slope', 'calib_inter', 'calib_sigma']]
    calib_params.drop_duplicates(inplace=True)
    calib_params['level'] = 'calibration'
    calib_params = calib_params.melt('level', var_name='parameter')

    # Save the parameter samples
    params = pd.concat([params, hyperparams, calib_params], sort=False)
    params['carbon_source'] = g[0]
    params['compound_turnover'] = g[1]
    params['strain'] = 'REL606'
    param_samples.append(params)

    # Tidy the concentration samples
    yield_concs = samps[['yield_concs_dim_0', 'yield_concs']]
    yield_concs.drop_duplicates(inplace=True)

    # Map the dimension to the od 
    yield_concs['od_600nm'] = [data_dict['optical_density'][i] for i in yield_concs['yield_concs_dim_0'].values]
    yield_concs['replicate'] = [data_dict['idx'][i] for i in yield_concs['yield_concs_dim_0'].values]
    yield_concs['strain'] = 'REL606'
    yield_concs['carbon_source'] = g[0]
    yield_concs['compound_turnover'] = g[1]
    yield_concs.rename(columns={'yield_concs':'compound_conc_mM'}, inplace=True)
    yield_concs.drop(columns=['yield_concs_dim_0'], inplace=True)
    conc_samples.append(yield_concs)

    # Compute the summary statistics for yield
    for _g, _d in yield_concs.groupby(['replicate', 'od_600nm']):
        mean_val = _d['compound_conc_mM'].mean()
        median_val = _d['compound_conc_mM'].median()
        ci_95_upper, ci_95_lower = np.percentile(_d['compound_conc_mM'], (97.5, 2.5))
        ci_75_upper, ci_75_lower = np.percentile(_d['compound_conc_mM'], (87.5, 12.5))
        conc_summary = conc_summary.append({
                            'strain':'REL606',
                            'carbon_source':g[0],
                            'compound_turnover': g[1],
                            'replicate':_g[0],
                            'od_600nm':_g[1],
                            'mean_val_mM':mean_val,
                            'median_val_mM':median_val,
                            'ci_95th_upper_mM': ci_95_upper,
                            'ci_95th_lower_mM': ci_95_lower,
                            'ci_75th_upper_mM': ci_75_upper,
                            'ci_75th_lower_mM':ci_75_lower
                            }, ignore_index=True)

    # Compute the summary statistics for the parameters
    for _g, _d in params.groupby(['level', 'parameter']):
        mean_val = _d['value'].mean()
        median_val = _d['value'].median()
        ci_95_upper, ci_95_lower = np.percentile(_d['value'], (97.5, 2.5))
        ci_75_upper, ci_75_lower = np.percentile(_d['value'], (87.5, 12.5))
        param_summary = param_summary.append({
                            'strain':'REL606',
                            'carbon_source':g[0],
                            'compound_turnover': g[1],
                            'level':_g[0],
                            'parameter':_g[1],
                            'mean_val':mean_val,
                            'median_val':median_val,
                            'ci_95th_upper': ci_95_upper,
                            'ci_95th_lower': ci_95_lower,
                            'ci_75th_upper': ci_75_upper,
                            'ci_75th_lower':ci_75_lower
                            }, ignore_index=True)

    
    # Compute the dependent variable ranges
    conc_range = np.linspace(0, 30, 100)
    
    # Isolate the calibration parameters
    calib_slope = params[params['parameter']=='calib_slope']['value'].values
    calib_inter = params[params['parameter']=='calib_inter']['value'].values

    # Compute the calibration ranges
    for perc, lab in zip(percs, perc_labels):
        ci = np.zeros((2, len(conc_range))) 
        for i, c in enumerate(conc_range):
            fit = calib_inter + calib_slope * c
            ci[:, i] = np.percentile(fit, perc)
        _df = pd.DataFrame([])
        _df['conc_range_mM'] = conc_range
        _df['rel_area_lower'] = ci[0, :]
        _df['rel_area_upper'] = ci[1, :]
        _df['percentile'] = lab
        _df['strain'] = 'REL606'
        _df['carbon_source'] = g[0]
        _df['compound_turnover'] = g[1]
        calib_fits.append(_df)

    # Compute the yield fit range
    min_od = 0.9 * d['od_600nm'].min() 
    max_od = 1.1 * d['od_600nm'].max()
    od_range = np.linspace(min_od, max_od, 100)
    for _g, _d in params.groupby(['level']):
        if _g != 'calibration':
            # Get the relevant parameters
            yield_slope = _d[_d['parameter']=='yield_slope']['value'].values
            yield_inter = _d[_d['parameter']=='yield_inter']['value'].values
            for perc, lab in zip(percs, perc_labels):
                ci = np.zeros((2, len(conc_range))) 
                for i, od in enumerate(od_range):
                    fit =  yield_inter + yield_slope * od 
                    ci[:, i] = np.percentile(fit, perc)
                _df = pd.DataFrame([])
                _df['od_600nm'] = od_range
                _df['conc_lower_mM'] = ci[0, :]
                _df['conc_upper_mM'] = ci[1, :]
                _df['percentile'] = lab
                _df['level'] = _g
                _df['strain'] = 'REL606'
                _df['carbon_source'] = g[0]
                _df['compound_turnover'] = g[1]
                yield_fits.append(_df)

# Concateate everything 
param_samples = pd.concat(param_samples, sort=False)
conc_samples = pd.concat(conc_samples, sort=False)
calib_fits = pd.concat(calib_fits, sort=False)
yield_fits = pd.concat(yield_fits, sort=False)

# Save results to disk
param_samples.to_csv('../../data/collated_turnover_parameter_MCMC_samples.csv', 
                     index=False)
conc_samples.to_csv('../../data/collated_turnover_concentration_MCMC_samples.csv', 
                    index=False)
param_summary.to_csv('../../data/collated_turnover_parameter_MCMC_summary.csv',
                    index=False) 
conc_summary.to_csv('../../data/collated_turnover_concentration_MCMC_summary.csv',
                    index=False) 
calib_fits.to_csv('../../data/collated_turnover_calibration_MCMC_percentiles.csv',
                    index=False)
yield_fits.to_csv('../../data/collated_turnover_yield_MCMC_percentiles.csv',
                    index=False)

# %%
