##########################################################################
# Author: Sambit K. Mishra, 
#         Cancer Genomics Research Laboratory
#         Frederick National Laboratory
#         Division of Cancer Epidemiology and Genetics
#         National Cancer Institute
# Created: Aug. 25, 2022
# 
# The input to the script is the csv file generated by 
# scripts/split_multi_allelic_var.py.
# This script will parse through a set of VCF files and predict
# which single nucleotide variants (SNVs) are true. Depending on
# the model type selected by the user, the prediction will be
# restricted only to the SNVs having specific low/intermediate VAFs.
# E.g., for VM models there will be no specific filter imposed for low VAF
# but, a filter for high VAF (i.e., VAF < 0.6) will be imposed.
##########################################################################

import pandas as pd
import os
import click
from vcf_parser import VCFParser
import glob
import numpy as np
from pathlib import Path
from Bio import SeqIO
import sys
import joblib
import warnings
warnings.filterwarnings("ignore") # Will suppress any unnecessary warnings

def filter_snvs(snv_file, model_type):
    """
    Will filter SNVs by their FDP, FAO and VAF depending on the model type.
    """
    df = pd.read_csv(snv_file)
    # First filter SNVs by FAO (>= 1) and FDP (>= 100)
    df = df.loc[(df['FAO'] >= 1) & (df['FDP'] >= 100)]
    df.reset_index(drop=True, inplace=True)
    # The cutoffs used here are the cutoffs that rendered optimal performance
    # for each model.
    if model_type == 'FM':
        vaf_upper = 0.6
        vaf_lower = 0.01
    elif model_type == 'VM':
        vaf_upper = 0.6
        vaf_lower = None
    else:
        vaf_upper = 0.01
        vaf_lower = 0.5
    if model_type == 'VM':
        df = df.loc[df['AF'] < vaf_upper]
    else:
        df = df.loc[(df['AF'] > vaf_lower) &
                    (df['AF'] < vaf_upper)]
    df.reset_index(drop=True, inplace=True)
    return df

def add_mutation_context(df_data, seq_file):
    """
    For each variant, add the tri-nucleotide context - nucleotide at the 5'end and
    the nucleotide at the 3' end.
    
    The reference sequence must be unpadded.
    """
    ref_seq = str(SeqIO.read(seq_file, "fasta").seq)
    nucl_numeric_map = {'A': 1, 
                       'G': 2,
                       'T': 3,
                       'C': 4}
    ref_len = len(ref_seq)
    tri_nucl_context = []
    nucl_5_prime = []
    nucl_3_prime = []
    mutation_list = []
    for idx_i, row_i in df_data.iterrows():
        pos, ref, alt = row_i['POS'], row_i['REF'], row_i['ALT']
        pos = int(pos)
        # Remember - position is the actual serial position, not the index.
        if pos == 1:
            nt_5 = ref_seq[-1]
            nt_3 = ref_seq[1]
        elif pos == ref_len:
            nt_5 = ref_seq[ref_len-2]
            nt_3 = ref_seq[0]
        else:
            nt_5 = ref_seq[pos-2] # index of the 5' nucleotide
            nt_3 = ref_seq[pos] # index of the 3' nucleotide
        nucl_5_prime.append(nucl_numeric_map[nt_5]) # Convert the nucleotide alphabets into numerics
        nucl_3_prime.append(nucl_numeric_map[nt_3]) # Convert the nucleotide alphabets into numerics
        tri_nucl_context.append(nt_5 + '.'  + nt_3)
        mutation_list.append(ref + '>' + alt)
    df_data['TRI_NUCLEOTIDE_CONTEXT'] = tri_nucl_context
    df_data['Mutation'] = mutation_list
    df_data['5_PRIME_NUCLEOTIDE_CONTEXT'] = nucl_5_prime
    df_data['3_PRIME_NUCLEOTIDE_CONTEXT'] = nucl_3_prime

    return df_data   

def create_feature_dataframe(df_snvs, hpv_ref_file, model_type):
    """
    Add the remaining features (e.g, 3' context, 5' context) to the dataframe
    and extract the model-specific features for prediction.
    """
    # Add the mutation context
    df_snvs = add_mutation_context(df_snvs, hpv_ref_file)
    
    # Filter features depending on model type
    # feat_cat_dict = {'Moderate': ['FSAF','FSAR','FSRF','FSRR','FWDB','FXX','GQ','MLLD','QUAL','REFB',
    #                             'REVB','SAF','SAR','SRF','SRR','SSSB','STB','VARB', '5_PRIME_NUCLEOTIDE_CONTEXT',
    #                             '3_PRIME_NUCLEOTIDE_CONTEXT'],
    #                 'Strict': ['FSAF','FSAR','FSRF','FSRR','FWDB','FXX','MLLD','QUAL','REFB','REVB','SSSB','VARB', 
    #                        '5_PRIME_NUCLEOTIDE_CONTEXT', '3_PRIME_NUCLEOTIDE_CONTEXT'],
    #                 'Exhaustive': ['AO','DP','FAO','FDP','FRO','FSAF','FSAR','FSRF','FSRR','FWDB',
    #                         'FXX','GQ','HRUN','LEN','MLLD','QD','QUAL','RBI','REFB','REVB',
    #                         'RO','SAF','SAR','SRF','SRR','SSSB','STB','STBP','VARB', '5_PRIME_NUCLEOTIDE_CONTEXT',
    #                         '3_PRIME_NUCLEOTIDE_CONTEXT']
    #                 }
    # if model_type == 'FM':
    #     feature_cols = feat_cat_dict['Exhaustive']
    # elif model_type == 'VM':
    #     feature_cols = feat_cat_dict['Moderate']
    # else:
    #     feature_cols = feat_cat_dict['Moderate']
    # df_snvs = df_snvs[feature_cols]
    return df_snvs 

def get_ensemble_prediction(y_scores_all_models):
    """
    For a given testing point, get the median score across all the models.
    If the median score is >= 0.5, then label 1 else label 0. 
    """
    median_scores = list(y_scores_all_models.median(axis=1))
    median_labels = list(map(lambda x: 0 if (round(x,2) < 0.5) else 1, median_scores))
    return median_labels
    

@click.command()
@click.option('--snv_file', required=True, help='The csv file generated by the split_multi_allelic_var.py script.\
For VM and FVM models, remember to filter the vcf files with VCFgenie before running the split_multi_allelic_var.py script.')
@click.option('--model_type', type=click.Choice(['VM', 'FM', 'FVM']), required = True, default='VM',
              show_default=True, help='The model type to be used for making the prediction')
@click.option('--hpv_ref_file', required = True, help='The reference genome file for the HPV type. File must contain the unpadded \
genome sequence in fasta format (i.e., genome sequence should not include a duplication of bases in the start region.)')
@click.option('--outfile', required=True, help='Name of the output file that will contain the prediction results')

def predict_true_snvs(snv_file, model_type, outfile, hpv_ref_file):
    if not os.path.isfile(hpv_ref_file):
        sys.exit(f"Error! HPV reference genome file {hpv_ref_file} not found!")
    print ("Extracting SNVs...", end="", flush=True)
    df_snvs = filter_snvs(snv_file, model_type)
    print ("Done!")
    print ("Calculating features ...", end="", flush=True)
    df_snvs = create_feature_dataframe(df_snvs, hpv_ref_file, model_type)
    print ("Done!")
    # Make predictions using the ensemble of trained XGB models
    print ("Running predictions...", end="", flush=True)
    model_dir = '../' + model_type + '_models/'
    all_models = glob.glob(model_dir + "xgb*")
    df_y_pred_all_models = pd.DataFrame()
    df_y_scores_all_models = pd.DataFrame()
    # Filter features depending on model type
    feat_cat_dict = {'Moderate': ['FSAF','FSAR','FSRF','FSRR','FWDB','FXX','GQ','MLLD','QUAL','REFB',
                                'REVB','SAF','SAR','SRF','SRR','SSSB','STB','VARB', '5_PRIME_NUCLEOTIDE_CONTEXT',
                                '3_PRIME_NUCLEOTIDE_CONTEXT'],
                    'Strict': ['FSAF','FSAR','FSRF','FSRR','FWDB','FXX','MLLD','QUAL','REFB','REVB','SSSB','VARB', 
                           '5_PRIME_NUCLEOTIDE_CONTEXT', '3_PRIME_NUCLEOTIDE_CONTEXT'],
                    'Exhaustive': ['AO','DP','FAO','FDP','FRO','FSAF','FSAR','FSRF','FSRR','FWDB',
                            'FXX','GQ','HRUN','LEN','MLLD','QD','QUAL','RBI','REFB','REVB',
                            'RO','SAF','SAR','SRF','SRR','SSSB','STB','STBP','VARB', '5_PRIME_NUCLEOTIDE_CONTEXT',
                            '3_PRIME_NUCLEOTIDE_CONTEXT']
                    }
    if model_type == 'FM':
        feature_cols = feat_cat_dict['Exhaustive']
    elif model_type == 'VM':
        feature_cols = feat_cat_dict['Moderate']
    else:
        feature_cols = feat_cat_dict['Moderate']
    # Only retain the required features for prediction purposes
    df_snv_features = df_snvs[feature_cols]
    for model_i in all_models:
        # Load model i and predict whether true or false SNV
        xgb_mdl_i = joblib.load(model_i)
        model_num = model_i.split('/')[-1].split('_')[-1]        
        y_pred = xgb_mdl_i.predict(df_snv_features)
        y_pred_scores = xgb_mdl_i.predict_proba(df_snv_features)[:, 1]
        df_y_pred_all_models[f"Model_{str(model_num)}"] = y_pred
        df_y_scores_all_models[f"Model_{model_i}"] = y_pred_scores
    # Get the consensus predictions from all models
    y_pred_consensus = get_ensemble_prediction(df_y_scores_all_models)
    pred_label = [True if pred_i == 1 else False for pred_i in y_pred_consensus]
    df_snvs['Predicted_SNV_Type'] = pred_label
    df_snvs.to_csv(outfile, index=False)
    print ("Done!")
    print (f"Wrote predictions to {outfile}")


if __name__ == '__main__':
    predict_true_snvs()    