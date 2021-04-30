#%%
import pandas as pd 

folders = ['2021-04-04_REL606_glucose_growth',
           '2021-04-27_REL606_acetate_growth']
dfs = []
for i, f in enumerate(folders):
    data = pd.read_csv(f'../../../data/growth_curves/{f}/processed/{f}.csv') 
    dfs.append(data)
data = pd.concat(dfs, sort=False)
if 'time_idx' in data.keys():
    data.drop(columns='time_idx', inplace=True)
data.to_csv('../../../data/collated_growth_measurements.csv', index=False)




# %%
