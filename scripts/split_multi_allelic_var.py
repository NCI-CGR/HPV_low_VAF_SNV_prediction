################################################################################
# Author: Sambit K. Mishra, 
#         Cancer Genomics Research Laboratory
#         Frederick National Laboratory
#         Division of Cancer Epidemiology and Genetics
#         National Cancer Institute
# Created: Aug. 25, 2022
# 
# This script will parse through a given folder containing VCF files
# using the parse vcf repo at https://github.com/moonso/vcf_parser.
# Will create a single VCF file with all the multi-allelic variants split.
# 
# Important: The script assumes the python code from the above repo to be already installed
# in the current environment. If the code is not installed, then please install it first
# and then run the script.
################################################################################

import pandas as pd
import os
import click
from vcf_parser import VCFParser
import glob
import numpy as np

def split_info(info_str, info_cols):
    """
    Split the info string and return only values for fields indicated by info_cols
    """
    info_list_sel = []
    info_arr = info_str.split(';')
    info_dict = {}
    for elem_i in info_arr:
        if '=' not in elem_i:
            # print(elem_i)
            continue
        k,v = elem_i.split('=')
        info_dict[k] = v
    for col_i in info_cols:
        if col_i in info_dict.keys():
            info_list_sel.append(info_dict[col_i])
        else:
            info_list_sel.append('')
    return info_list_sel                

def split_format(format_header_str, format_val_str, format_cols):
    """
    Split the format_header and format_val string and return 
    only values for fields indicated by format_cols
    """
    format_dict = {}
    format_list_sel = []
    header_arr = format_header_str.split(':')
    val_arr = format_val_str.split(':')
    for h, v in zip(header_arr, val_arr):
        format_dict[h] = v
    for col_i in format_cols:    
        if col_i in format_dict.keys():
            format_list_sel.append(format_dict[col_i])
        else:
            format_list_sel.append('')    
    return format_list_sel      


def split_multiple_alleles(vcf_file_i, num_metadata_lines):
    """
    For a given vcf file split multiple variant alleles for a position
    into separate rows.
    """
    # print (vcf_file_i)
    # First get the name of the sample from the vcf file
    df_temp = pd.read_csv(vcf_file_i, delimiter='\t', skiprows=num_metadata_lines) # For vcf genie generated vcf files
    # print (df_temp)
    sample_name = df_temp.columns[9]
    # print (f"sample name = {sample_name}")

    # Parse the vcf file into a list of dictionary items using the VCF parser.
    # Each line in the vcf file is dictionary.
    vcf_parser = VCFParser(infile=vcf_file_i, split_variants=True, check_info=True)
    # print (vcf_parser)
    # Declare a pandas dataframe that will store the vcf information
    df_vcf_1 = pd.DataFrame() # INFO fields and FORMAT fileds will be a string of text
    df_vcf_2 = pd.DataFrame() # INFO and format fields will be further broken down
                                 # into individual columns
    col_names_1 = ['sample', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
        'INFO', 'FORMAT', 'SAMPLE_COL']
    info_cols = ['AF','AO','DP','FAO','FDP','FR','FRO','FSAF','FSAR','FSRF','FSRR',
    'FWDB','FXX','HRUN','LEN','MLLD','QD','RBI', 'REFB','REVB','RO','SAF','SAR','SRF',
    'OPOS', 'OREF', 'OALT', 'OID', 'OMAPALT', 'SRR','SSEN','SSEP','SSSB','STB',
    'STBP','TYPE','VARB']
    format_cols = ['GT','GQ']            
    col_names_2 = ['sample', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER'] + info_cols + format_cols
    for vcf_i in vcf_parser:
        # print (vcf_i)
        # print (f"sample name = {sample_name}")
        # vcf_i is a dictionary where the keys are the different columns 
        # Create a dataframe for each record
        data_list_1 = [sample_name, vcf_i['CHROM'], vcf_i['POS'], vcf_i['ID'],
            vcf_i['REF'], vcf_i['ALT'], vcf_i['QUAL'], vcf_i['FILTER'],
            vcf_i['INFO'], vcf_i['FORMAT'], vcf_i[sample_name]]
        # print (data_list_1)
        # print (f"Length of data_list_1 = {len(data_list_1)}")    
        df_1 = pd.DataFrame(data=[data_list_1], columns=col_names_1)
        df_vcf_1 = pd.concat([df_vcf_1, df_1])
        # print (df_vcf_1.head())
        
        # Create another dataframe that splits the INFO and FORMAT strings
        info_values = split_info(vcf_i['INFO'], info_cols)
        format_values = split_format(vcf_i['FORMAT'], vcf_i[sample_name], format_cols)
        data_list_2 = [sample_name, vcf_i['CHROM'], vcf_i['POS'], 
            vcf_i['ID'], vcf_i['REF'], vcf_i['ALT'], vcf_i['QUAL'], vcf_i['FILTER']] + info_values + format_values # Append info_values and format_values
        df_2 = pd.DataFrame(data=[data_list_2], columns=col_names_2)
        df_vcf_2 = pd.concat([df_vcf_2,df_2])
    return df_vcf_1, df_vcf_2        


def concat_vcfs(vcf_dir, num_metadata_lines):
    """
    Get concatenated tables of vcf records after splitting the multiple alleles
    """
    all_vcf_files = glob.glob(vcf_dir + '/*.vcf')
    df_1_all = pd.DataFrame()
    df_2_all = pd.DataFrame()
    count = 1
    for vcf_file_i in all_vcf_files:
        print (f"Processing file {vcf_file_i} ({count}/{len(all_vcf_files)})...", end='', flush=True)
        df_1, df_2 = split_multiple_alleles(vcf_file_i, num_metadata_lines)
        if df_1_all.empty:
            df_1_all = df_1
            df_2_all = df_2
        else:
            df_1_all = pd.concat([df_1_all,df_1])
            df_1_all.reset_index(drop=True)    
            df_2_all = pd.concat([df_2_all,df_2])
            df_2_all.reset_index(drop=True)
        print ("Done!")
        count += 1    
    outfile = 'vcfs_multallele_info_format_split.csv'
    # Drop the columns from df_2_all that are not split properly 
    # between the multiple alleles.
    drop_cols_df2 = ['OPOS', 'OREF', 'OALT', 'OID', 'OMAPALT']
    df_2_all.drop(columns=drop_cols_df2, inplace=True)
    df_2_all.to_csv(outfile, index=False)
    print (f"Wrote parsed vcf data to {outfile}")             

@click.command()
@click.option('--vcf_dir', required=True, help='The directory containing the vcf files (with the headers and meta data)')
@click.option('--num_metadata_lines', required=True, type=int, default=70, help='The number of meta data lines \
in the vcf files (excluding the header line). Meta data lines begin with "##". Default is 70 for the raw vcfs \
generated by our pipeline. However, for the filtered VCF files generated by VCFgenie it is 86.')

def main(vcf_dir, num_metadata_lines):
    if not vcf_dir.endswith('/'):
        vcf_dir += '/'
    # print(f"vcf_dir = {vcf_dir}, outfile={outfile}")
    concat_vcfs(vcf_dir, num_metadata_lines)

if __name__ == '__main__':
    main()