#%%
import numpy as np 
import pandas as pd 
import diaux.viz 
import altair as alt 
from altair_saver import save
colors, palette = diaux.viz.altair_style()

# Load the parameter samples for the yields
samples = pd.read_csv('../../data/collated_turnover_parameter_MCMC_samples.csv')

# Restrict only to the yield coefficients of the hyper parameters
samples = samples[(samples['level']=='hyperparameter') & 
                  (samples['parameter']=='yield_slope')]

# Compute the absolute value of each 
samples['value'] = np.abs(samples['value'])

# Tidy
gluc_yield = samples[(samples['carbon_source']=='glucose') & 
                              (samples['compound_turnover']=='glucose')]['value']
ace_yield = samples[(samples['carbon_source']=='acetate') & 
                              (samples['compound_turnover']=='acetate')]['value']
ace_secr = samples[(samples['carbon_source']=='glucose') & 
                              (samples['compound_turnover']=='acetate')]['value']

# Convert to biomass yield 
# from Meyers et al 2013, OD 1 of 1 L approx 0.45 g DW
gluc_mass = 180 * 1E-3 #  per mM
ace_mass = 60 * 1E-3 # g per mM
OD_CONV = 0.45 * 0.55 # g protein per L at OD 1 

# Compute the yield produce/consumed in mass of substrate per unit biommass
gluc_yield = gluc_yield * gluc_mass / OD_CONV
ace_yield = ace_yield * ace_mass / OD_CONV
ace_secr = ace_secr * ace_mass /OD_CONV


# Compute the secreted acetate mass per mass of consumed glucose
rand_ace_secr = np.random.choice(ace_secr, replace=True, size=len(ace_secr))
rand_gluc_yield = np.random.choice(gluc_yield, replace=True, size=len(ace_secr))
ace_secr_gluc =  rand_ace_secr / rand_gluc_yield 
gluc_yield_df = pd.DataFrame(1/gluc_yield, columns=['value'])
ace_yield_df = pd.DataFrame(1/ace_yield, columns=['value'])
ace_secr_df = pd.DataFrame(ace_secr_gluc, columns=['value'])
gluc_yield_df
    

#%%
# Set up baic growth parameters
VOL = 0.01 # in L
GLUC_0 = 0.025 # in g / L
m_glucose = GLUC_0 * VOL

biomass_prod_gluc = m_glucose / gluc_yield # Biomass produced by consuming glucose
ace_prod_gluc = m_glucose * ace_secr_gluc # Mass of acetate secreted

# Do subsampling to get the biomass produced from the secreted acetate
rand_ace_yield = np.random.choice(ace_yield, replace=True, size=len(biomass_prod_gluc))
rand_ace_prod = np.random.choice(ace_prod_gluc, replace=True, size=len(biomass_prod_gluc))
biomass_prod_ace = rand_ace_prod / rand_ace_yield # in units of biomass

# Compute the extremes of fitness from subsampling
fitness = (biomass_prod_ace + biomass_prod_gluc) / (biomass_prod_gluc)
fit_pred = pd.DataFrame([])
fit_pred['fitness'] = fitness
fit_pred['label'] = r'prediction (95% CI)'


# Load the Lanski and Travisano fitnes data 
fit_data = pd.read_csv('../../data/literature/Lenski_Travisano_1994/Lenski_Travisano_1994_relative_fitness.csv')
fit_data['label'] = 'Lenski & Travisano, 1994'
# %%
# Set up the yield histograms
gluc_yield_strip = alt.Chart(gluc_yield_df, height=50, width=300).mark_tick(clip=True).encode(
                    x=alt.X('value:Q', title='biomass produced / mass glucose',
                            scale=alt.Scale(domain=[0.1, 1])),
).properties(title=f'glucose yield, Ω glucose ≈ {np.mean(1/gluc_yield):0.2f}')

ace_yield_strip = alt.Chart(ace_yield_df, height=50, width=300).mark_tick(clip=True).encode(
                    x=alt.X('value:Q', title='biomass produced / mass acetate',
                        scale=alt.Scale(domain=[0.1, 1]))
).properties(title=f'acetate yield, Ω acetate ≈ {np.mean(1/ace_yield):0.2f} ', height=50)

ace_secr_strip = alt.Chart(ace_secr_df, height=50, width=300).mark_tick(clip=True).encode(
                    x=alt.X('value:Q', title='mass acetate secreted / mass glucose consumed',
                       scale=alt.Scale(domain=[0.1, 1])),
).properties(title=f'acetate secretion, ω acetate ≈ {np.mean(ace_secr_gluc):0.2f}')

strips = gluc_yield_strip & ace_yield_strip & ace_secr_strip
strips
#%%
# Set up the data plots
data_base = alt.Chart(fit_data)
points = data_base.mark_point(size=80, color=colors['primary_blue']).encode(
            x=alt.X('generation:Q', title='generation'),
            y=alt.Y('mean(fitness):Q', title='relative fitness',
                    scale=alt.Scale(zero=False)),
            color=alt.Color('label:N', title=None)
)
errs = data_base.mark_errorbar(extent='stdev').encode(
            x=alt.X('generation:Q', title='generation'),
            y=alt.Y('fitness:Q', title='relative fitness'),
            color=alt.Color('label:N', title=None)
)

fit_window = alt.Chart(fit_pred).mark_errorband(extent='ci').encode(
                y=alt.Y('fitness:Q', title='relative fitness'),
                color=alt.Color('label:N', title=None)
)

data_plot = (points + errs + fit_window).properties(title='fitness measurements and predicted diauxic fitnesss gain')

save(strips | data_plot, '../../figures/yields_and_fitness_gain.pdf')

# %%
