#%%
import numpy as np 
import pandas as pd 
import cremerlab.hplc
import diaux.hplc
import glob

# Conver the files
DATA_PATH='../../../data/metabolite_turnover/2021-04-27_REL606_acetate_turnover/'
raw_files = sorted(glob.glob(f'{DATA_PATH}/raw/*.txt'))
cremerlab.hplc.convert(raw_files, peak_table=True)

# Load the peak table and assemble as necessary
files = sorted(glob.glob(f'{DATA_PATH}/raw/converted/*peak_table.csv'))
#%%
peaks = []
for i, f in enumerate(files):
    # Parse the file name
    _, _, _, rep, _, time_idx, _, _ = f.split('/')[-1].split('_') 
    rep = int(rep)
    time_idx = int(time_idx)
    
    # Load peaks and assign the replicate and time point identifier
    peak_df = pd.read_csv(f, comment='#')
    peak_df['replicate'] = rep
    peak_df['time_idx'] = time_idx
    peaks.append(peak_df)

# Concatenate 
peaks = pd.concat(peaks, sort=False)

#%%
# Assign the peaks
peaks = diaux.hplc.assign_peaks(peaks)
peaks.to_csv(f'{DATA_PATH}/processed/2021-04-27_REL606_acetate_peak_table.csv',
             index=False)

# %%
