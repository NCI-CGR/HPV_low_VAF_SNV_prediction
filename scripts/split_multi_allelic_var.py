################################################################################
# Author: Sambit K. Mishra, 
#         Cancer Genomics Research Laboratory
#         Frederick National Laboratory
#         Division of Cancer Epidemiology and Genetics
#         National Cancer Institute
# Created: Aug. 18, 2022
# 
# This script will parse through a given folder containing VCF files
# using the parse vcf repo at https://github.com/moonso/vcf_parser.
# Will create a single VCF file with all the multi-allelic variants split. 
# In addition, it will also split the multinucleotide variants into separate 
# rows. E.g., AC > AT at position 10 is a C > T change at position 11.
# 
# The script reports the total number of SNVs which can be provided
# as input to VCFgenie for the multiple testing correction. 
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
from pathlib import Path

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
    For a given vcf file split the multiple alleles
    """
    df_temp = pd.read_csv(vcf_file_i, delimiter='\t', skiprows=num_metadata_lines)
    sample_name = df_temp.columns[9]
    # Parse the vcf file into a list of dictionary items using the VCF parser.
    # Each line in the vcf file is dictionary.
    vcf_parser = VCFParser(infile=vcf_file_i, split_variants=True, check_info=True)
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
        # vcf_i is a dictionary where the keys are the different columns 
        # Create a dataframe for each record
        data_list_1 = [sample_name, vcf_i['CHROM'], vcf_i['POS'], vcf_i['ID'],
            vcf_i['REF'], vcf_i['ALT'], vcf_i['QUAL'], vcf_i['FILTER'],
            vcf_i['INFO'], vcf_i['FORMAT'], vcf_i[sample_name]]
        df_1 = pd.DataFrame(data=[data_list_1], columns=col_names_1)
        df_vcf_1 = pd.concat([df_vcf_1, df_1])
        
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
    # Drop the columns from df_2_all that are not split properly 
    # between the multiple alleles.
    drop_cols_df2 = ['OPOS', 'OREF', 'OALT', 'OID', 'OMAPALT']
    df_2_all.drop(columns=drop_cols_df2, inplace=True)
    return df_2_all

def split_mnp(df_data):
    """
    Split the Multi-nucletide polymorphisms into separate rows: e.g., REF = AC, ALT = AT
    """
    df_mnv = pd.DataFrame()
    df_all_new = pd.DataFrame()
    for ids_i, rec_i in df_data.iterrows():
        if ((len(str(rec_i['REF'])) > 1 and len(str(rec_i['ALT'])) > 1) and 
            (len(rec_i['REF']) == len(rec_i['ALT'])) and
            (hamming_distance(str(rec_i['REF']),str(rec_i['ALT'])) == 1)):
            df_mnv = pd.concat([df_mnv, rec_i.to_frame().transpose()])
            ref_i = str(rec_i['REF'])
            alt_i = str(rec_i['ALT'])
            diff_ind, chr1, chr2 = string_diff(ref_i,alt_i)
            new_rec = rec_i
            new_rec['REF'] = chr1
            new_rec['ALT'] = chr2
            new_rec['POS'] = int(rec_i['POS']) + diff_ind
            df_all_new = pd.concat([df_all_new, new_rec.to_frame().transpose()])
        else:
            df_all_new = pd.concat([df_all_new, rec_i.to_frame().transpose()])
    df_all_new.reset_index(drop=True, inplace=True)    
    return df_all_new, df_mnv

def hamming_distance(str1, str2):
    """
    Get the hamming distance between 2 strings
    """
    return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(str1, str2))))

def string_diff(str1, str2):
    """
    Get the position where str1 and str2 differ by a single nucleotide
    """
    if hamming_distance(str1,str2) > 1:
        print (f"Error! Hamming distance between {str1} and {str2} > 1")
        return
    if len(str1) != len(str2):
        print (f"Error! Length of {str1} does not equal length of {str2}")
        return
    else:
        for ind_i, chr1, chr2 in zip(range(0,len(str1)), list(str1), list(str2)):
            if chr1 != chr2:
                return ind_i, chr1, chr2             

def get_snp_variants(df_vcf):
    """
    From a dataframe of vcfs, retain only monoallelic SNPs
    """
    df_vcf_snps = df_vcf.loc[(df_vcf['ALT'].apply(len) == 1) & (df_vcf['REF'].apply(len) == 1)]
    df_vcf_snps.reset_index(drop=True, inplace=True)
    return df_vcf_snps


@click.command()
@click.option('--vcf_dir', required=True, help='The directory containing the vcf files (with the headers)')
@click.option('--output_dir', required=True, help='The directory to which the parsed files will be written')
@click.option('--num_metadata_lines', required=True, type=int, default=70, help='The number of meta data lines \
in the vcf files (excluding the header line). Meta data lines begin with "##". Default is 70 for the raw vcfs \
generated by our pipeline. However, for the filtered VCF files generated by VCFgenie it is 86.')


def main(vcf_dir, output_dir, num_metadata_lines):
    if not vcf_dir.endswith('/'):
        vcf_dir += '/'
    if not output_dir.endswith('/'):
        output_dir += '/'
    # Create output directory if not already present
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # print(f"vcf_dir = {vcf_dir}, outfile={outfile}")
    
    # Split multi-allelic SNPs
    df_2_all_split1 = concat_vcfs(vcf_dir, num_metadata_lines)
    df_2_all_split1.reset_index(drop=True, inplace=True)
    
    # Split SNPs with a multi-nucleotide context (AT > AC)
    print ("Splitting MNPs...", end="", flush=True)
    df_2_all_split2, _ = split_mnp(df_2_all_split1)
    df_2_all_split2.reset_index(drop=True, inplace=True)
    print ("Done!")
    
    # Retain only SNPs
    print ("Extracting only SNVs...", end='', flush=True)
    df_snps = get_snp_variants(df_2_all_split2)
    print ("Done!")

    # Write the SNPs to a .csv file
    df_snps.to_csv(output_dir + 'snvs.csv')

    print (f"Total SNVs = {len(df_snps)}. Wrote SNVs to {output_dir +  'snvs.csv'}")

if __name__ == '__main__':
    main()