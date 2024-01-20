## Merge the individual hyperparameter performance files generated using the
## xgb_hyper_param_tuning_parallelized_new.py script. 

import glob
import os
import pandas as pd
import click

@click.command()
@click.option('--input_dir', type=str, required=True, 
help='Name of the input directory that contains output subdirectories for FM, VM and FVM from running the xgb_hyper_param_tuning_parallelized_new.py script. Should be the same as the value for outdir used for the xgb_hyper_param_tuning_parallelized_new.py script.')
@click.option('--outfile', type=str, required=True, help='Name of the file to which the consolidated output/performance will be written into.')

# Empty dataframe to store performance for all data
df_perf_all = pd.DataFrame()
if not input_dir.endswith('/'):
    input_dir += '/'

for model_i in ['FM', 'VM', 'FVM']:
    inputdir_model_i = input_dir + model_i
    subdirs = glob.glob(inputdir_model_i + '/seed_*')
    count = 1
    for subdir_i in subdirs: # Each train/test iteration
        seed_i = subdir_i.split('/')[-1].split('_')[1]
        print (f"Running for directory {count}...", end="\r", flush=True)
        # Read all the files containing performance metrics
        # for each hyperparameter combination (total 3^6+2^2 = 2916 combinations)
        for file_i in glob.glob(subdir_i + '/' +  'comb*'):
            comb_id = int(file_i.split('/')[-1].split('.')[0].split('_')[1])
            df_perf_comb_i = pd.read_csv(file_i)
            # Include a column for seed
            df_perf_comb_i['Seed'] = int(seed_i)
            # Include a column for model
            df_perf_comb_i['Model'] = model_i
            # Include a column to include the hyper-param combination ID
            df_perf_comb_i['Hyperparam_comb_id'] = comb_id
            # Append to df_perf_all
            df_perf_all = pd.concat([df_perf_all, df_perf_comb_i])
        count += 1    
df_perf_all.reset_index(drop=True, inplace=True)
df_perf_all.to_csv(outfile, index=False)
 








