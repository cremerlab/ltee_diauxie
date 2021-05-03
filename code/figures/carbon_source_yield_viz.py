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
summary = pd.read_csv('../../data/collated_turnover_parameter_MCMC_summary.csv')

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
for g, _ in data.groupby(['carbon_source', 'compound']):

        # Get the parameter summary
        _summary = summary[(summary['carbon_source']==g[0]) & 
                           (summary['compound_turnover'] == g[1]) &
                           (summary['level']=='hyperparameter') & 
                           (summary['parameter']=='yield_slope')]

        # Define the plot title
        mean_val = np.abs(_summary['mean_val'].values[0])
        ci_95_low, ci_95_up = np.sort(np.abs(np.array(
                                        [_summary['ci_95th_upper'].values[0], 
                                        _summary['ci_95th_lower'].values[0]])))
        ci_95_up = np.abs(ci_95_up) - mean_val
        ci_95_low = mean_val - np.abs(ci_95_low)
        title = f'{g[0]} growth, {g[1]} yield = {mean_val:0.0f} [+{ci_95_up:0.0f}, -{ci_95_low:0.0f}] mM / OD'
        _data = concs[(concs['carbon_source']==g[0]) & 
                        (concs['compound_turnover']==g[1])]
        _turn = turn_percs[(turn_percs['carbon_source']==g[0]) & 
                            (turn_percs['compound_turnover']==g[1]) &
                           (turn_percs['level']=='hyperparameter')]
        _params = params[(params['carbon_source']==g[0]) & 
                         (params['compound_turnover']==g[1])]
        yield_slope = _params[(_params['parameter']=='yield_slope')]
        yield_inter = _params[(_params['parameter']=='yield_inter')]

        # Set up the yield plot
        conc_base = alt.Chart(_data, height=250)
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
        slope_reps =  alt.Chart(yield_slope
                               ).mark_tick(
                                  opacity=0.4, 
                                  thickness=0.01
                               ).encode(
                                  x=alt.X('value:Q', title=f' {g[1]} yield coefficient [mM / OD]',
                                    scale=alt.Scale(zero=False)),
                                  y=alt.Y('level:N', title=None)
                                ).properties(
                                  height=100, 
                                  width=200)
        int_reps =  alt.Chart(yield_inter
                             ).mark_tick(
                                opacity=0.4, 
                                thickness=0.01
                            ).encode(
                                x=alt.X('value:Q', title=f' {g[1]} intercept [mM]',
                                     scale=alt.Scale(zero=False)),
                                y=alt.Y('level:N', title=None)
        ).properties(height=100, width=200)


        yield_plot = (percs + errors + points).resolve_scale(
                                        color='independent'
                                        ).properties(
                                                title=title
                                        )
        rep_plot = slope_reps & int_reps
        out = rep_plot | yield_plot
        save(out, f'../../figures/inferred_REL606_{g[0]}_growth_{g[1]}_turnover.pdf')


# (percs + errors + points).resolve_scale(color='independent')

# %%
