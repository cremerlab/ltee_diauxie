#%%
import numpy as np 
import pandas as pd 
import arviz as az 
import cmdstanpy

# Load the data 
data = pd.read_csv('../../data/collated_growth_measurements.csv')

# Load the inferrential model 
model = cmdstanpy.CmdStanModel(stan_file='../stan/hierarchical_growth_rate.stan')

# Define the percentiles to compute for the fits. 
percs = [(2.5, 97.5), (12.5, 87.5), (25, 75), (37.5, 62.5), (47.5, 52.5)]
perc_labels = [95, 75, 50, 25, 5]
# %%

# Iterate through each carbon source and strain and perform the inference
samples, summary, fits = [], [], []
for g, d in data.groupby(['strain', 'carbon_source', 'medium_base']):
    data_dict = {'J':d['replicate'].max(),
                 'N':len(d),
                 'idx':d['replicate'].values.astype(int),
                 'time':d['elapsed_time_hr'].values.astype(float),
                 'OD':d['od_600nm'].values.astype(float)}
    _samples = model.sample(data=data_dict)
    _samples = az.from_cmdstanpy(_samples)
    _samples = _samples.posterior.to_dataframe().reset_index()

    # Get the low level parameters
    params = _samples[['lam_dim_0', 'OD_init', 'lam']]  
    params.rename(columns={'lam_dim_0':'level'}, inplace=True)
    params['level'] = [f'replicate {i+1}' for i in params['level'].values]
    params = params.melt('level', var_name='parameter')
    
    # hyper_params
    hyper_params = _samples[['OD_init_mu', 'lam_mu', 'OD_init_sigma', 'lam_sigma', 'sigma']]
    hyper_params['level'] = 'hyperparameter'
    hyper_params = hyper_params.melt('level', var_name='parameter')
    hyper_params.loc[hyper_params['parameter']=='lam_mu', 'parameter']  = 'lam'
    hyper_params.loc[hyper_params['parameter']=='OD_init_mu', 'parameter']  = 'OD_init'

    # Concatenate the longform
    params = pd.concat([params, hyper_params], sort=False)
    
    # Add identifiers and store
    params['strain'] = g[0]
    params['carbon_source'] = g[1]
    params['medium_base'] = g[2]
    samples.append(params)

    # Compute the parameter summaries
    summaries = pd.DataFrame([])
    for _g, _d in params.groupby(['level', 'parameter']):
        mean = _d['value'].mean()
        median = _d['value'].median()
        ci_95_up, ci_95_low = np.percentile(_d['value'].values, (97.5, 2.5))
        summaries = summaries.append({
                    'level': _g[0],
                    'parameter': _g[1],
                    'mean_val': mean ,
                    'median_val': median ,
                    'ci_95th_upper':ci_95_up,
                    'ci_95th_lower':ci_95_low,
                    'strain':g[0],
                    'carbon_source':g[1],
                    'medium_base':g[2]
                    },
                    ignore_index=True)
    summary.append(summaries)

    # Compute the fits
    min_time = d['elapsed_time_hr'].min()
    max_time = d['elapsed_time_hr'].max()
    if min_time == 0:
        min_time = 0
    else:
        min_time *= 0.9
    max_time *= 1.1 
    time_range = np.linspace(min_time, max_time, 200)
    for _g, _d in params.groupby(['level']):
        lam = _d[_d['parameter']=='lam']['value'].values
        od_init = _d[_d['parameter']=='OD_init']['value'].values
        for perc, lab in zip(percs, perc_labels):
            ci = np.zeros((2, len(time_range)))
            for j, t in enumerate(time_range): 
                fit = od_init*np.exp(lam * t)
                ci[:, j] = np.percentile(fit, perc)
            perc_df = pd.DataFrame([]) 
            perc_df['elapsed_time_hr'] = time_range
            perc_df['level'] = _g
            perc_df['percentile'] = lab
            perc_df['od_600nm_lower'] = ci[0, :]
            perc_df['od_600nm_upper'] = ci[1, :]
            perc_df['strain'] = g[0]
            perc_df['carbon_source'] = g[1]
            perc_df['medium_base'] = g[2]
            fits.append(perc_df)

# Concatenate and save the data
samples = pd.concat(samples, sort=False)
summary = pd.concat(summary, sort=False)
fits = pd.concat(fits, sort=False)
samples.to_csv('../../data/collated_growth_MCMC_parameter_samples.csv', index=False)
summary.to_csv('../../data/collated_growth_MCMC_parameter_summary.csv', index=False)
fits.to_csv('../../data/collated_growth_MCMC_percentiles.csv', index=False)

# %%


