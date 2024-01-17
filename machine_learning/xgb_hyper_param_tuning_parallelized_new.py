### Train and test models using all the hyperparameter combinations.
### Includes testing for 50 iterations using
### 50 seeds. 
### The training/testing split is performed at the sample-level: 25 samples for trainig
### and 6 for testing.

### It takes a lot of time for this script to run - I have therefore parallelized items
### to run on multiple cores. The user will need to run separate instances of this script
### for each of the FM, VM and FVM prediction strategies.

import pandas as pd
import seaborn as sns
import math
import os
import matplotlib.pyplot as plt
import glob
import random
from scipy import stats
from pathlib import Path
import warnings
from itertools import product
import joblib
from joblib import Parallel, delayed
import click
import psutil

warnings.filterwarnings("ignore") # Will suppress any unnecessary warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# ML libraries
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import cross_validate
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_auc_score, roc_curve, auc, classification_report
from sklearn.metrics import f1_score, precision_score, recall_score, matthews_corrcoef
from sklearn.metrics import mean_squared_error, cohen_kappa_score, make_scorer
from sklearn.metrics import confusion_matrix, accuracy_score, average_precision_score, ConfusionMatrixDisplay


def remove_low_coverage_samples(df_snv_w_cov, min_cov_cutoff, coverage_column):
    """
    Updated to exclude minimum depth requirement.
    Remove iSNVs that do not meet quality control criteria: FAO > 0
        df_snv_w_cov: The dataframe with all SNVs and their features
        min_cov_cutoff: The minimum value for coverage
        coverage_column: The coverage column (FDP or FAO for Ion Torrent data)
    """
    df_snv_high_cov = df_snv_w_cov.copy()
    df_snv_high_cov = df_snv_high_cov.loc[df_snv_high_cov[coverage_column] >= min_cov_cutoff]
    df_snv_high_cov.reset_index(drop=True, inplace=True)
    return df_snv_high_cov
    
def create_balanced_datasets(X_train, y_train, true_label=1, false_label=0):
    """
    Create balanced datasets for training. The test data will remain imbalanced.
        seed_value : The seed for random state selection
        t_size : Size of the test data. Note that this is 30% of the entire data 
                (true label + false label), not 30% of each label. You may notice that 
                the test data for the individual labels (true label or false label) are not exactly 
                30% of the entire data. 
        true_label : The numeric label for the true variants
        false_label: The numeric label for the false variants
    """
    train_samples = {} # Dictionary in which we will store the training samples
    train_index_list = [] # Will keep track of the indices of the sampled rows in the training data.
                          # Depending on the value set for alpha, the number of unique indices
                          # should equal to the size of the training data.
    X_train.reset_index(inplace=True, drop=True)
    
    # If the training data is not balanced, then undersample the over-represented class
    # and create multiple training datasets.
    num_true = y_train.count(true_label) # Number of true labels
    num_false = y_train.count(false_label) # Number of false labels
    
    # Identify the over-represented class 
    if num_true > num_false:
        over_rep_class = true_label
        over_rep_count = num_true
        under_rep_class = false_label
        under_rep_count = num_false
    elif num_false > num_true:
        over_rep_class = false_label
        over_rep_count = num_false
        under_rep_class = true_label
        under_rep_count = num_true
    else:
        # Both classes are balanced
        train_samples['1'] = {}
        train_samples['1']['X_train'] = X_train
        train_samples['1']['y_train'] = y_train
        return train_samples, train_index_list
    
    # # Commenting out this debug info statement for now. Uncomment when needed.
    #print ("\tTrain-test split resulted in imbalanced training data", flush=True)
    #print (f"\tOver-represented class = {over_rep_class}, count = {over_rep_count}", flush=True)
    #print (f"\tUnder-represented class = {under_rep_class}, count = {under_rep_count}", flush=True)
    #print ("\tWill proceed by creating balanced training data.", flush=True)
    
    alpha = 5 # A sampling weight constant that determines the number of times
            # the over-represented class will be sampled.
    num_sample_iter = int(round(over_rep_count/under_rep_count)*alpha)
    sample_size = under_rep_count
    seed_list = list(range(1,num_sample_iter+1)) # The seeds we will use for each iteration of sampling.
                                                 # This is to make the results reproducible.
    # Identify the indices of the over-represented and under-represented class
    # and get the respective data in X_train.
    over_rep_indices = list(np.where(np.array(y_train) == over_rep_class)[0])
    X_train_over_rep = X_train.iloc[over_rep_indices]
    under_rep_indices = list(np.where(np.array(y_train) == under_rep_class)[0])
    X_train_under_rep = X_train.iloc[under_rep_indices]
   
    # Perform sampling
    for seed_i in seed_list:
#         print (f"Sampling iteration {seed_i}...", flush=True, end='')
#         print (sample_size)
#         print (len(X_train_over_rep))
#         print (len(X_train_under_rep))
        X_train_sample_i_over_rep= X_train_over_rep.sample(n=sample_size, replace=False, random_state=seed_i)
        y_train_sample_i_over_rep = [over_rep_class] * sample_size
        # Consolidate the training feature data for the under represented class
        # and over-represented class into a single data frame
        #X_train_sample_i = X_train_sample_i_over_rep.append(X_train_under_rep)
        X_train_sample_i = pd.concat([X_train_sample_i_over_rep, X_train_under_rep])
        
        ind_i_list = X_train_sample_i.index.tolist()
    #     print (ind_i_list)
    #     break
        if len(train_index_list) == 0:
            train_index_list = ind_i_list
        else:    
            train_index_list.extend(ind_i_list)
        
        # Turn off reset index for de-bugging
        X_train_sample_i.reset_index(inplace=True, drop=True)
        # Consolidate the training labels for the under represented
        # and over represented classes
        y_train_sample_i = y_train_sample_i_over_rep  + [under_rep_class] * sample_size
        
        # Shuffle the rows. Otherwise the Top N rows will be over-rep class and 
        # bottom N rows will be under-rep class.
        X_train_sample_i_shuffled = X_train_sample_i.sample(frac=1, random_state=seed_i)
        shuffled_indices = X_train_sample_i_shuffled.index.to_list()
        #
        y_train_sample_i_shuffled = [y_train_sample_i[ind_i] for ind_i in shuffled_indices]
        train_samples[seed_i] = {}
        train_samples[seed_i]['X_train'] = X_train_sample_i_shuffled
        train_samples[seed_i]['y_train'] = y_train_sample_i_shuffled
    return train_samples, train_index_list
 
def get_ensemble_prediction(y_scores_all_models):
    """
    For a given testing point, get the median score across all the models.
    If the median score is >= 0.5, then label 1 else label 0. 
    """
    median_scores = list(y_scores_all_models.median(axis=1))
    median_labels = list(map(lambda x: 0 if (round(x,2) < 0.5) else 1, median_scores))
    return median_labels


def get_ensemble_prediction_scores(y_scores_all_models):
    """
    Return the median scores across all the models for a given test 
    data point.
    """
    median_scores = list(y_scores_all_models.median(axis=1))
    median_scores = [round(score_i,2) for score_i in median_scores]
    return median_scores


# Instead of doing a 5-fold CV, we will perform training and testing on the given set of samples.
# That is, instead of further breaking the training set into training and validation, we will use the
# entire training set for training on the given hyperparameter combination and testing on the given testing data.
def train_test_comb_i (comb_id, comb_i, df_snvs_training_comb_i, df_snvs_testing_comb_i, true_var, false_var, af_lower, af_upper, feature_cat, target_cols, outdir):
    print (f"comb id = {comb_id}, combination = {comb_i}", flush=True)
    eta_i, n_trees_i, gamma_i, max_depth_i, subsample_i, colsample_bytree_i, reg_lambda_i, reg_alpha_i = comb_i
    # Skip combination if output file already exists
    if os.path.isfile(outdir + "comb_" + str(comb_id) + '.csv'):
        print ("Skipping as performance metric file already exists")
        return
    
    df_perf_comb_i = pd.DataFrame()
    feat_cat_dict = {'Moderate': ['FSAF','FSAR','FSRF','FSRR','FWDB','FXX','GQ','MLLD','QUAL','REFB',
                    'REVB','SAF','SAR','SRF','SRR','SSSB','STB','VARB', '5_PRIME_NUCLEOTIDE_CONTEXT', '3_PRIME_NUCLEOTIDE_CONTEXT'],
                    'Strict': ['FSAF','FSAR','FSRF','FSRR','FWDB','FXX','MLLD','QUAL','REFB','REVB','SSSB','VARB', 
                               '5_PRIME_NUCLEOTIDE_CONTEXT', '3_PRIME_NUCLEOTIDE_CONTEXT'],
                    'Exhaustive': ['AO','DP','FAO','FDP','FRO','FSAF','FSAR','FSRF','FSRR','FWDB',
                    'FXX','GQ','HRUN','LEN','MLLD','QD','QUAL','RBI','REFB','REVB',
                    'RO','SAF','SAR','SRF','SRR','SSSB','STB','STBP','VARB', '5_PRIME_NUCLEOTIDE_CONTEXT',
                    '3_PRIME_NUCLEOTIDE_CONTEXT']
                    }
    feature_cols = feat_cat_dict[feature_cat]
    X_train = df_snvs_training_comb_i[feature_cols]
    y_train = df_snvs_training_comb_i[target_cols].tolist()
    
    X_testing = df_snvs_testing_comb_i[feature_cols]
    y_testing = df_snvs_testing_comb_i[target_cols].tolist()
                                                      
    # Balance the training set
    train_sample, train_index_list = create_balanced_datasets(X_train,y_train)

    # Train XGB model
    xgb_model = XGBClassifier(use_label_encoder=False,
                          booster='gbtree', # boosting algorithm to use, default gbtree, othera: gblinear, dart
                          n_estimators=n_trees_i, # number of trees, default = 100
                          eta=eta_i, # this is learning rate, default = 0.3
                          max_depth=max_depth_i, # maximum depth of the tree, default = 6
                          gamma = gamma_i, # used for pruning, if gain < gamma the branch will be pruned, default = 0
                          reg_lambda = reg_lambda_i, # regularization parameter, defautl = 1
                          eval_metric = 'logloss',
                          colsample_bytree = colsample_bytree_i,
                          alpha = reg_alpha_i,
                          subsample = subsample_i    
                             )

    fpr_tpr = [] # Store all the fpr, tpr values for each model 
    auc_scores = []
    df_y_pred_all_models = pd.DataFrame()
    df_y_scores_all_models = pd.DataFrame()
    df_all_cv_scores = pd.DataFrame()
    df_feature_importances_all_models = pd.DataFrame()
    for model_i in list(train_sample.keys()):
        X_train_model_i = train_sample[model_i]['X_train']
        y_train_model_i = train_sample[model_i]['y_train']
#         print (f"Running CV and predictions for model {model_i}...", end='', flush=True)
        X_train_model_i = X_train_model_i.apply(pd.to_numeric)
        xgb_model.fit(X_train_model_i,y_train_model_i)

        # Predict using model_i on the testing data
        y_pred = xgb_model.predict(X_testing)
        df_y_pred_all_models[f"Model_{model_i}"] = y_pred
        y_testing_pred_scores = xgb_model.predict_proba(X_testing)[:, 1]
        df_y_scores_all_models[f"Model_{model_i}"] = y_testing_pred_scores
        false_positive, true_positive, _ = roc_curve(y_testing, y_testing_pred_scores)
        fpr_tpr.append([false_positive,true_positive])
        # Calculate ROC only if both true and false labels are present.
        # Otherwise, the function roc_auc_score returns error.
        if len(set(y_testing)) == 2:
            auc_test = roc_auc_score(y_testing, y_pred)
        else:
            auc_test = math.nan
        auc_scores.append(auc_test)
    # Calculate metrics only if y_testing has both true and false iSNVs labels.
    # Otherwise, the AUC cannot be calculated as the function roc_auc_score returns error.
    if len(set(y_testing)) == 2:
        # Get the consensus predictions from all models
        y_pred_consensus = get_ensemble_prediction(df_y_scores_all_models)

        # Calculate metrics on the consensus prediction
        auc_test = roc_auc_score(y_testing, y_pred_consensus)

        # Calculate the F1 score
        F1_score_test = f1_score(y_testing, y_pred_consensus)

        # Calculate the mean-squared error
        mse_test = mean_squared_error(y_testing, y_pred_consensus)

        # Calculate the accuracy
        accuracy_test = accuracy_score(y_testing, y_pred_consensus)

        # Calculate the MCC
        mcc_test = matthews_corrcoef(y_testing, y_pred_consensus)

        # Print the performance metrics on test data
        # Uncomment to see the metrics for each combination
        #print ("AUC =", round(auc_test, 3), ", F1 score =", F1_score_test, ", Mean-squared error =", mse_test, ", Accuracy =", accuracy_test, ", Matthews correlation coefficient =", mcc_test)

        # Get the confusion matrix
        cm = confusion_matrix(y_testing, y_pred_consensus)

        # Get the indices of false positives
        false_positive_ind = list(set(np.where(np.array(y_testing) == 0)[0]).intersection(set(np.where(np.array(y_pred_consensus) == 1)[0])))
        # Get the indices of true positives
        true_positive_ind = list(set(np.where(np.array(y_testing) == 1)[0]).intersection(set(np.where(np.array(y_pred_consensus) == 1)[0])))
        # Get the indices of true negatives
        true_negative_ind = list(set(np.where(np.array(y_testing) == 0)[0]).intersection(set(np.where(np.array(y_pred_consensus) == 0)[0])))
        # Get the indices of false negatives
        false_negative_ind = list(set(np.where(np.array(y_testing) == 1)[0]).intersection(set(np.where(np.array(y_pred_consensus) == 0)[0])))

        fp = len(false_positive_ind)
        tp = len(true_positive_ind)
        tn = len(true_negative_ind)
        fn = len(false_negative_ind)
    else:
        # Get the indices of false positives
        false_positive_ind = list(set(np.where(np.array(y_testing) == 0)[0]).intersection(set(np.where(np.array(y_pred_consensus) == 1)[0])))
        # Get the indices of true positives
        true_positive_ind = list(set(np.where(np.array(y_testing) == 1)[0]).intersection(set(np.where(np.array(y_pred_consensus) == 1)[0])))
        # Get the indices of true negatives
        true_negative_ind = list(set(np.where(np.array(y_testing) == 0)[0]).intersection(set(np.where(np.array(y_pred_consensus) == 0)[0])))
        # Get the indices of false negatives
        false_negative_ind = list(set(np.where(np.array(y_testing) == 1)[0]).intersection(set(np.where(np.array(y_pred_consensus) == 0)[0])))

        fp = len(false_positive_ind)
        tp = len(true_positive_ind)
        tn = len(true_negative_ind)
        fn = len(false_negative_ind)
        auc_test = math.nan
        mcc_test = math.nan
        F1_score_test = math.nan
        accuracy_test = math.nan
        mse_test = math.nan

    print (f"Hyper-param comb = {comb_id}, AUC = {round(auc_test , 3)}, F1 score = {F1_score_test}, MSE = {mse_test}, Accuracy = {accuracy_test}, MCC = {mcc_test}", flush=True)
    cols = ['ETA', 'Num_trees', 'Gamma', 'Max_depth', 'Sub_sample', 'Colsample_bytree', 'Reg_lambda', 'Reg_alpha', 
                'True_Var_Def', 'False_Var_Def', 'VAF_Lower_Limit', 'VAF_Upper_Limit',  'Feature_Category',
                'Num_True_Var_Testing', 'Num_False_Var_Testing', 'TP', 'FP', 'TN', 'FN', 'AUC', 'MCC', 
                'F1_score', 'Accuracy', 'MSE']

    df_perf_comb_i = pd.DataFrame(data=[[eta_i, n_trees_i, gamma_i, max_depth_i, subsample_i, colsample_bytree_i, 
                                     reg_lambda_i, reg_alpha_i, true_var, false_var, af_lower, af_upper, feature_cat, 
                                     tp+fn, tn+fp, tp, fp, tn, fn, auc_test, mcc_test, F1_score_test,
                                     accuracy_test, mse_test]], columns=cols)
    
    df_perf_comb_i.to_csv(outdir + "comb_" + str(comb_id) + '.csv', index=False)


@click.command()
@click.option('--model_type', type=click.Choice(['VM', 'FM', 'FVM']), required = True, default='FM',
              show_default=True, help='The model type for which hyper-parameters need to be tuned')
@click.option('--num_jobs', type=int, required=True, help='Total tuning jobs to be run in parallel. Each job will correspond to a single hyperparameter combination')
@click.option('--outdir', type=str, required=True, help='The name of the output directory to which results will be written')
              

def run_hyper_param_tuning(model_type, num_jobs, outdir):
    if model_type == 'FM':
        # The input vcf .csv file
        snv_sbs_file = '../data/SNV_data_wo_vcf_genie.csv'
    else:
        snv_sbs_file = '../data/SNV_data_with_vcf_genie.csv'
  
    if model_type == 'FM':
        af_lower = 0.01
        af_upper = 0.5
        feature_cat = 'Strict'
        true_var = 1
        false_var = 0.33
    elif model_type == 'VM':
        af_upper = 0.5
        feature_cat = 'Moderate'
        true_var = 1
        false_var = 0.33
    else:
        af_lower = 0.01
        af_upper = 0.5
        feature_cat = 'Strict'
        true_var = 1
        false_var = 0.33

    eta_list = [0.1, 0.2, 0.3]
    n_trees_list = [50, 100, 150]
    gamma_list = [0, 1, 2]
    max_depth_list = [3, 6, 10]
    subsample_list = [0, 0.5, 1]
    colsample_bytree_list =[0, 0.5, 1]
    reg_lambda_list = [0, 1, 2]
    reg_alpha_list = [0, 1, 2]

    hyper_param_comb = list(product(eta_list, n_trees_list, gamma_list, max_depth_list, 
                                    subsample_list, colsample_bytree_list, reg_lambda_list, reg_alpha_list))
    print (f"Total combinations = {len(hyper_param_comb)}")

    df_snv_sbs = pd.read_csv(snv_sbs_file)

    # Create a column to include the sample id (excluding the well id)
    df_snv_sbs['sample_id'] = df_snv_sbs['sample'].apply(lambda x: x.split('_')[0])



    # Remove variants with FAO = 0
    fao_cutoff = 1
    df_snv_sbs_fao_filtered = remove_low_coverage_samples(df_snv_sbs, fao_cutoff, 'FAO')

    random.seed(10)
    
    # Using the above seed we will select 50 random integers between 1 and 500.
    # Using this set seed above will help us generate the same set of random integers
    # and make the results reproducible.    
    int_list = list(range(1,500))
    seed_list = random.sample(int_list, 50)
        
    # Each seed will correspond to one iteration of train/test split and 
    # 5-fold CV on the selected training set.
    for seed_i in seed_list:
        # We will consider the SNVs in each replicate set only once and perform
        # training and testing for each replicate set independently.
        df_perf_out_xgb = pd.DataFrame()
        
        ## Create output directory
        if not outdir.endswith('/'):
            outdir += '/'
        # Create a sub-directory to include results from this seed    
        outdir_seed_i = outdir + 'seed_' + str(seed_i) + '/'    
        Path(outdir_seed_i).mkdir(parents=True, exist_ok=True)
        
        # Get all the sample IDs
        all_sample_ids = df_snv_sbs_fao_filtered['sample_id'].drop_duplicates().tolist()

        # Set the random sampling seed to seed_i
        random.seed(seed_i)

        # We will pick 6 testing samples and make sure SNVs from these samples
        # are not included in the testing set.
        testing_set_samples = random.sample(all_sample_ids, k=6)
        print (f"testing samples = {testing_set_samples}")
        training_set_samples = [sample_id for sample_id in all_sample_ids if sample_id not in testing_set_samples]

        # Get the SNVs for corresponding to the training and testing data
        df_snvs_testing = df_snv_sbs_fao_filtered.loc[df_snv_sbs_fao_filtered['sample_id'].isin(testing_set_samples)]
        df_snvs_training = df_snv_sbs_fao_filtered.loc[~df_snv_sbs_fao_filtered['sample_id'].isin(testing_set_samples)]
        df_snvs_testing.reset_index(drop=True, inplace=True)
        df_snvs_training.reset_index(drop=True, inplace=True)

        # Prepare training and testing data for the given parameter-feature combination
        df_snvs_training_comb_i = df_snvs_training.copy()
        if model_type != 'VM':
            df_snvs_training_comb_i = df_snvs_training_comb_i.loc[(df_snvs_training_comb_i['AF'] > af_lower) & 
                                                       (df_snvs_training_comb_i['AF'] < af_upper)]
        else:
            df_snvs_training_comb_i = df_snvs_training_comb_i.loc[(df_snvs_training_comb_i['AF'] < af_upper)]    
        
        df_snvs_training_comb_i = df_snvs_training_comb_i.loc[(df_snvs_training_comb_i['PERCENT_OVERLAP'] == true_var) | (df_snvs_training_comb_i['PERCENT_OVERLAP'] == false_var)]
        df_snvs_training_comb_i.reset_index(drop=True, inplace=True)
        overlap_cat = list(df_snvs_training_comb_i['PERCENT_OVERLAP'].apply(lambda x: 1 if (x == true_var) else 0 ))

        # Create a binary response column to identify a variant as true variant (1) or false variant (0)
        df_snvs_training_comb_i['OVERLAP_CATEGORY'] = overlap_cat
        target_cols = 'OVERLAP_CATEGORY' # Response column 

        df_snvs_testing_comb_i = df_snvs_testing.copy()
        if model_type != 'VM':
            df_snvs_testing_comb_i = df_snvs_testing_comb_i.loc[(df_snvs_testing_comb_i['AF'] > af_lower) & 
                                                       (df_snvs_testing_comb_i['AF'] < af_upper)]
        else:
            df_snvs_testing_comb_i = df_snvs_testing_comb_i.loc[(df_snvs_testing_comb_i['AF'] < af_upper)]
        df_snvs_testing_comb_i = df_snvs_testing_comb_i.loc[(df_snvs_testing_comb_i['PERCENT_OVERLAP'] == true_var) | (df_snvs_testing_comb_i['PERCENT_OVERLAP'] == false_var)]
        df_snvs_testing_comb_i.reset_index(drop=True, inplace=True)
        overlap_cat = list(df_snvs_testing_comb_i['PERCENT_OVERLAP'].apply(lambda x: 1 if (x == true_var) else 0 ))

        # Create a binary response column to identify a variant as true variant (1) or false variant (0)
        df_snvs_testing_comb_i['OVERLAP_CATEGORY'] = overlap_cat
        target_cols = 'OVERLAP_CATEGORY' # Response column 
        
        # If performance file for testing data is already present,
        # then skip this iteration.
        if os.path.isfile(outdir_seed_i + model_type + '_' + str(seed_i) + '_hyperparameter_tuning_perf.csv'):
            continue

        # Iterate over the hyperparameter combinations and evaluate performance for the given set data
        # We will parallelize this step.
        if model_type == 'VM':
            af_lower = None
        Parallel(n_jobs=num_jobs)(delayed(train_test_comb_i) (index_i, comb_i, df_snvs_training_comb_i, df_snvs_testing_comb_i, true_var, false_var, af_lower, af_upper, feature_cat, target_cols, outdir_seed_i) for index_i, comb_i in enumerate(hyper_param_comb))
        
        # Combine metrics from individual files into one
        for comb_num in list(range(0,len(hyper_param_comb))):
            file_i = outdir_seed_i + "comb_" + str(comb_num) + '.csv'
            df_i = pd.read_csv(file_i)
            df_perf_out_xgb = pd.concat([df_perf_out_xgb, df_i])
        
        # The output file that will include the performance metrics
        # over all the combinations.
        perf_out_file_xgb = outdir_seed_i + model_type + '_' + str(seed_i) + '_hyperparameter_tuning_perf.csv'    
        df_perf_out_xgb.to_csv(perf_out_file_xgb, index=False)


if __name__ == '__main__':
    run_hyper_param_tuning()
    
