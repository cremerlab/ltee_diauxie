#%%
import numpy as np 
import pandas as pd 
import altair as alt 
from altair_saver import save 
import diaux.viz
colors, palette = diaux.viz.altair_style()
alt.data_transformers.disable_max_rows()
# %% 
#  Load the collated data
data = pd.read_csv('../../data/collated_turnover_measurements.csv')
calib = pd.read_csv('../../data/calibration/2021-04-05_hplc_calibration/processed/2021-04-05_NC_DM_calibration_relative_areas.csv')
calib = calib[(calib['buffer_base']=='DM') & 
                        (calib['compound'].isin(['glucose', 'acetate']))
                        ][['carbon_conc_mM', 'rel_area_phosphate', 'compound']]
concs = pd.read_csv('../../data/collated_turnover_concentration_MCMC_summary.csv')
calib_percs = pd.read_csv('../../data/collated_turnover_calibration_MCMC_percentiles.csv')
turn_percs = pd.read_csv('../../data/collated_turnover_yield_MCMC_percentiles.csv')
params = pd.read_csv('../../data/collated_turnover_parameter_MCMC_samples.csv')

#%%
calib_data_base = alt.Chart(calib)
calib_points = calib_data_base.mark_point(size=80).encode(
                x=alt.X('carbon_conc_mM:Q', title='concentration [mM]'),
                y=alt.Y('rel_area_phosphate:Q', title='relative peak area of phosphate'),
                color=alt.Color('compound:N', title='compound')
)


# calib_points

#%%
_percs = calib_percs[(calib_percs['carbon_source']=='acetate') & 
                     (calib_percs['compound_turnover']=='acetate')]
perc_plot = alt.Chart(_percs).mark_area(color=colors['primary_black']).encode(
                x=alt.X('conc_range_mM:Q', title='concentration [mM]'),
                y=alt.Y('rel_area_lower:Q', title = 'relative peak area of phosphate'),
                y2='rel_area_upper:Q',
                opacity=alt.Opacity('percentile:O', sort='descending')
)
(perc_plot + calib_points).interactive()
# %%
g = ['glucose', 'glucose']
ac_data = concs[(concs['carbon_source']==g[0]) & 
                (concs['compound_turnover']==g[1])]
_turn = turn_percs[(turn_percs['carbon_source']==g[0]) & 
                    (turn_percs['compound_turnover']==g[1]) &
                   (turn_percs['level']=='hyperparameter')]
_params = params[(params['carbon_source']==g[0]) & 
                 (params['compound_turnover']==g[1])]
yield_slope = _params[(_params['parameter']=='yield_slope')]
yield_inter = _params[(_params['parameter']=='yield_inter')]

# Set up the yield plot
conc_base = alt.Chart(ac_data)
points = conc_base.mark_point(size=80, color=colors['primary_black']).encode(
            x=alt.X('od_600nm:Q', title='optical density [a.u.]',
                     scale=alt.Scale(zero=False)),
            y=alt.Y('mean_val_mM:Q', title=f'{g[1]} concentration [mM]',
                     scale=alt.Scale(zero=False)),
            shape=alt.Shape('replicate:N', title='replicate'),
)
errors = conc_base.mark_errorbar(color=colors['primary_black']).encode(
        x=alt.X('od_600nm:Q', title='optical density [a.u.]'),
        y=alt.Y('ci_95th_lower_mM:Q', title=f'{g[1]} concentration [mM]'),
        y2='ci_95th_upper_mM:Q',
)

percs = alt.Chart(_turn).mark_area(opacity=0.45).encode(
        x=alt.X('od_600nm:Q', title='optical density [a.u.]'),
        y=alt.Y('conc_lower_mM:Q', title=f'{g[1]} concentration [mM]'),
        y2='conc_upper_mM:Q',
        color=alt.Color('percentile:N', sort='descending',
                        scale=alt.Scale(scheme='blues')))

# Plot the stripplots
slope_reps =  alt.Chart(yield_slope).mark_tick(opacity=0.1).encode(
                    x=alt.X('value:Q', title=f' {g[1]} yield coefficient [mM / OD]'),
                    y=alt.Y('level:N', title=None)
).properties(height=100, width=200)

# (percs + errors + points).resolve_scale(color='independent')
slope_reps
# %%
