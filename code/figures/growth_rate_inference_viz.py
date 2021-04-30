#%%
import numpy as np 
import pandas as pd 
import altair as alt
from altair_saver import save
import diaux.viz
alt.data_transformers.disable_max_rows()
colors, palette = diaux.viz.altair_style()

# Load the experimental data, samping information, and fits
data = pd.read_csv('../../data/collated_growth_measurements.csv')
samples = pd.read_csv('../../data/collated_growth_MCMC_parameter_samples.csv')
fits = pd.read_csv('../../data/collated_growth_MCMC_percentiles.csv')
summary = pd.read_csv('../../data/collated_growth_MCMC_parameter_summary.csv')
summary = summary[(summary['level']=='hyperparameter') &
                  (summary['parameter']=='lam')]

# %%
for g, d in data.groupby(['strain', 'carbon_source']):
    _samps = samples[(samples['strain']==g[0]) & 
                     (samples['carbon_source']==g[1])]
    _summary = summary[(summary['strain']==g[0]) & 
                        (summary['carbon_source']==g[1])]
    _fits = fits[(fits['strain']==g[0]) & (fits['carbon_source']==g[1])]
    avg = _summary['mean_val'].values[0]
    ci_upper = _summary['ci_95th_upper'].values[0] - avg
    ci_lower = avg - _summary['ci_95th_lower'].values[0]
    # Define the plot title

    title = f'growth of {g[0]} in {g[1]}; λ = {avg:0.3f} [+{ci_upper:0.3f}, -{ci_lower:0.3f}] per hr'

    # Define the chart bases
    lam_samp_base = alt.Chart(_samps[_samps['parameter']=='lam'])
    od_init_samp_base = alt.Chart(_samps[_samps['parameter']=='OD_init'])
    fit_base = alt.Chart(_fits[_fits['level']=='hyperparameter'])
    data_base = alt.Chart(d)

    # Set up the tick plots
    lam_ticks = lam_samp_base.mark_tick(
                        strokeWidth=0.005,
                        opacity=0.1,
                        ).encode(
                        x=alt.X('value:Q', title='growth rate [per hr]',
                                scale=alt.Scale(zero=False)),
                        y=alt.Y('level:N', title=None),
                        color=alt.Color('level:N', sort='descending', legend=None) 
                        ).properties(height=100, width=200, title='MCMC growth rate samples')
    od_init_ticks = od_init_samp_base.mark_tick(
                        strokeWidth=0.005,
                        opacity=0.1
                    ).encode(
                        x=alt.X('value:Q', title='optical density [a.u.]',
                                scale=alt.Scale(zero=False)),
                        y=alt.Y('level:N', title=None),
                        color=alt.Color('level:N', sort='descending', legend=None)
                    ).properties(height=100, width=200, title='MCMC initial optical density samples')

    points = data_base.mark_point(size=80, opacity=0.75).encode(
                x=alt.X('elapsed_time_hr:Q', title='elapsed time [hr]'),
                y=alt.Y('od_600nm:Q', title='optical density [a.u.]',
                    scale=alt.Scale(type='log')),
                color=alt.Color('replicate:N', title='replicate')
                ).properties(width=350, height=250)
    fit = fit_base.mark_area(color=colors['primary_black']).encode(
                x=alt.X('elapsed_time_hr:Q', title='elapsed time [hr]'),
                y=alt.Y('od_600nm_lower:Q', title='optical density [a.u.]'),
                y2='od_600nm_upper:Q',
                opacity=alt.Opacity('percentile:O', title='percentile',
                                    sort='descending')
        )

    comparison = (fit + points).properties(title=title)
    layer = (lam_ticks & od_init_ticks) | comparison
    layer.properties(title=title)
    layer.resolve_scale(color='independent')
    save(layer, f'../../figures/inferred_{g[0]}_{g[1]}_growth.pdf')

# %%
