#%% 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import altair as alt
from altair_saver import save
import diaux.viz
import cremerlab.growth
colors, palette = diaux.viz.matplotlib_style()
_ = diaux.viz.altair_style()

# Define the data path and load raw data
DATA_PATH = '../../../data/growth_curves/2021-04-27_REL606_acetate_growth'
raw_data = pd.read_csv(f'{DATA_PATH}/raw/2021-04-27_REL606_acetate_growth.csv')
data, params = cremerlab.growth.infer_growth_rate(raw_data, 
                                                  od_bounds=[0.09, 0.4],
                                                  groupby=['replicate'],
                                                  viz=False)
_data = []
for g, d in data.groupby(['replicate']):
        d.sort_values('elapsed_time_hr', inplace=True)
        d['rel_od'] = d['od_600nm'].values / d['od_600nm'].values[0]
        _data.append(d)
_data = pd.concat(_data, sort=False)
#%% 
# compute fits
time_range = np.linspace(0, 4, 200)
fit_dfs = []
for g, d in params.groupby(['replicate']):
        lam = d[d['parameter']=='growth_rate']['map_val'].values[0]
        inter = d[d['parameter']=='od_init']['map_val'].values[0]
        fit = inter * np.exp(time_range * lam)
        _df = pd.DataFrame([])
        _df['elapsed_time_hr'] = time_range
        _df['od_600nm'] = fit
        _df['rel_od'] = _df['od_600nm'].values / fit[0]
        _df['replicate'] = g
        fit_dfs.append(_df)
fit_df = pd.concat(fit_dfs, sort=False)

#%%
# Generate a plot with all the replicates overlaid with the fits
points = alt.Chart(_data).mark_point(size=80, opacity=0.5).encode(
                x=alt.X('elapsed_time_hr:Q', title='elapsed time [hr]'),
                y=alt.Y('od_600nm:Q', title='optical density',
                        scale=alt.Scale(type='log')),
                color=alt.Color('replicate:N', title='biological replicate'))

fits = alt.Chart(fit_df).mark_line(size=1).encode(
                x=alt.X('elapsed_time_hr:Q', title='elapsed time [hr]'),
                y=alt.Y('od_600nm:Q', title='optical density',
                        scale=alt.Scale(type='log')),
                color=alt.Color('replicate:N', title='biological replicate'))

agged = params.groupby(['parameter'])['map_val'].agg(('mean', 'sem')).reset_index()
mean_lam, sem_lam = agged[agged['parameter']=='growth_rate'][['mean', 'sem']].values[0]

out = (fits + points).properties(title=f'REL606 acetate growth, λ = {mean_lam:0.3f} ± {sem_lam:0.3f} per hr.')
save(out, './output/2021-04-27_REL606_acetate_growth.pdf')
out

#%%
# Transform the boolean sampled column to annotated time points
dfs = []
for g, d in data.groupby(['replicate']):
        d['time_idx'] = 0
        sample_idx = []
        iter = 0
        for s in d['sampled'].values:
                if s == True:
                        iter += 1
                        sample_idx.append(iter)
                else:
                        sample_idx.append(0)
        d['time_idx'] = sample_idx
        dfs.append(d)

data = pd.concat(dfs, sort=False)
data.drop(columns=['sampled'], inplace=True)

# Save the pruned data and growth rate parameters
data.to_csv(f'{DATA_PATH}/processed/2021-04-27_REL606_acetate_growth.csv')
params.to_csv(f'{DATA_PATH}/processed/2021-04-27_REL606_acetate_growth_parameters.csv')
#%%
