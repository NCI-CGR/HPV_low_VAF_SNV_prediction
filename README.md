
# Predicting True Low-VAF SNVs in the Human Papillomavirus 


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#repository-contents">Repository Contents</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
## About The Project
Identifying true single nucleotide variants (SNVs) existing at low variant allele fractions (VAFs) is challenging, particularly because it is tricky to tease apart SNVs that are sequencing artifacts from genuine low-VAF SNVs. In this work, we present a novel machine learning technique to identify SNVs occurring at low-intermediate VAFs by using the human papillomavirus (HPV) as a case study.

We sequenced the HPV whole genome from 31 HPV18 positive samples in <b> triplicates (i.e., 31 samples x 3 reps per samples = 93 sequenced samples in total) </b> and used an in-house variant calling pipeline to detect variants in each sequenced sample. The variants identified were with respect to the reference HPV18 genome (NCBI Ref. ID AY262282). We labelled each variant with its <b> replicate frequency (the number of replicates of the sample in which the variant was detected) </b> and used the replicate frequency to identify true SNVs. Based on our underlying hypothesis, true SNVs would be detected in all (3/3) replicates, while false SNVs in 1/3 replicates of a sample. Consistent with our goal of detecting low-VAF variants, we restricted our study only to variants with VAF 1% - 60%. We then extracted sequencing features for each variant from the VCF file, supplemented them with features related to the 5' (5-prime) and 3' (3-prime) nucleotide contexts of the variant and trained Xtreme Gradient Boosting (XGBoost) binary classification models to predict true and false SNVs. We further complemented our supervised XGBoost models with an unsupervised approach [VCFgenie](https://github.com/chasewnelson/VCFgenie.git) and observed improved prediction performance. 

## Installation

The code in this repository is written in Python (v3.8) You will need access to a Linux terminal (even a Windows 10 having Windows Subsystem for Linux enabled and Ubuntu installed) to be able to execute the code. The specific instructions below work best in a Ubuntu 18.02/Ubuntu20.04 platform. Follow the steps below to install the required packages and run the scripts.

1. Clone the repo
```sh
  git clone https://github.com/NCI-CGR/HPV_low_VAF_SNV_prediction.git
```
2. Install miniconda
```sh
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  chmod +x Miniconda3-latest-Linux-x86_64.sh
  ./Miniconda3-latest-Linux-x86_64.sh
```
&ensp;&ensp;Follow the on-screen instructions to complete installation of miniconda and to initialize conda.

3. Create a conda environment

&ensp;&ensp;Create a conda environment that will include all the required packages necessary for running the repository code. We will use the <b>hpv_triplicate_env.yaml</b> file for this purpose.

```sh
  conda env create -f psp_gnm_env.yaml
```

4. Activate the conda environment
```sh
  conda activate hpv_triplicate_env
```
5. Install the vcf_parser package from https://github.com/moonso/vcf_parser

6. Make sure you are good to go
```sh
  cd HPV_low_VAF_SNV_prediction # go to the repo directory that you just cloned in Step 1
  python scripts/split_multi_allelic_var.py --help
```
&ensp;&ensp;The following help menu should be displayed
```sh
Usage: split_multi_allelic_var.py [OPTIONS]

Options:
  --vcf_dir TEXT                The directory containing the vcf files (with
                                the headers)  [required]
  --output_dir TEXT             The directory to which the parsed files will
                                be written  [required]
  --num_metadata_lines INTEGER  The number of meta data lines in the vcf files
                                (excluding the header line). Meta data lines
                                begin with "##". Default is 70 for the raw
                                vcfs generated by our pipeline. However, for
                                the filtered VCF files generated by VCFgenie
                                it is 86.  [required]
  --help                        Show this message and exit.
```


## Repository Contents
The contents of the repository are described below.
- <b> data/ </b>- The data directory includes the raw and the processed variant data that is required for training and testing the machine learning models. Within this directory, the <b> vcf_files/raw_vcf </b> and <b> vcf_files/vcfgenie_vcf </b> subdirectories include the raw VCFs generated by our in-house variant calling pipeline and the filtered VCFs generated upon post-filtering of the raw VCFs by [VCFgenie](https://github.com/chasewnelson/VCFgenie.git), respectively. The <b> SNV_data_wo_vcf_genie.csv </b> contains all the raw SNVs from our in-house pipeline and their features. The <b> SNV_data_with_vcf_genie.csv </b> includes SNVs extracted from the processed VCF files generated by VCFgenie. The PERCENT_OVERLAP indicates the replicate frequency of each SNV (1 => 3/3, 0.67 => 2/3, 0.33 => 1/3) and the AF is the VAF for a SNV. The HPV18_not_padded.fasta file includes the reference genome sequence for HPV18.

- <b> scripts/ </b> - The split_multi_allelic_var.py script is used to parse through a directory having VCF files, split multi-allelic positions into separate rows, only extract SNVs and their sequencing features and write the SNVs and their features into a csv file. <b> Note: to run this script, you must first download and install the vcf_parser package from https://github.com/moonso/vcf_parser </b>. The data/SNV_data_wo_vcf_genie.csv and data/SNV_data_with_vcf_genie.csv files were generated using this script. Checkout the command `python scripts/split_multi_allelic_var.py --help` to explore how to use this script.

- <b> machine_learning/ </b> - This directory includes two pieces of code related to machine learning. The <b> xgb_model_training_and_performance.ipynb </b> is an interactive jupyter notebook intended to guide the user through step-by-step instructions on how to reproduce the study results. It is intended to demonstrate to the user how training and testing were performed for the 3 scenarios: <b> a. Machine learning with VAF filters (FM models) </b>, <b> b. Machine learning with VCFgenie without any low-VAF filters (VM models) </b>, and <b> c. Machine learning with VCFgenie with low-VAF filters (FVM models). </b> The <b> predict_true_low_vaf_snvs.py </b> script is a more independent script that can be used for predicting true/false SNVs using the pre-trained models (doesn't require prior model training). More on this script is elaborated in the <a href="#usage">Usage</a> section of this document.

- <b> FM_models/ </b> - Includes trained XGB models for the FM scenario. Only includes models trained for the optimum parameter/feature combination.
- <b> VM_models/ </b> - Includes trained XGB models for the VM scenario. Only includes models trained for the optimum parameter/feature combination.
- <b> FVM_models/ </b> - Includes trained XGB models for the FVM scenario. Only includes models trained for the optimum parameter/feature combination.

- <b> results/ </b> - Includes performance files for models trained in each of the 3 scenarios - FM, VM and FVM. Files are generated upon running the jupyter notebook (machine_learning/xgb_model_training_and_performance.ipynb). For each scenario, the performance metrics are reported for models trained for each parameter/feature combination. Note that for each combination, performance metrics are reported for 10 testing datasets generated by random sampling using 10 different seeds. When calculating overall performance for a combination, we consider the median performance across all the 10 testing datasets.

<!-- USAGE EXAMPLES -->
## Usage

Now we will go through some specific examples of how to use this repository.

### Example 1: Reproducing study results
 To reproduce the results included in our study, simply run the jupyter notebook (machine_learning/xgb_model_training_and_performance.ipynb). For each scenario (FM, VM and FVM), the notebook will train and test XGBoost models on all possible parameter/feature combinations and write the performance metrics into results directory. <i> To reduce the execution time, you can skip running the sections that correspond to saving models for the best performing parameter/feature combination in each scenario (i.e., the section titled <i>"Save the models for the best performing scenario - in this case the VM models" </i> and all the code following this section) </i>. These sections have already been executed and the models saved in their respective directories.

Following is an example of the performance output. Below is a snapshot of the performance results for the FM model as observed in the <b> performance_FM.csv </b> file. Specifically, I am showing the performance for the optimal parameter/feature combination for FM by applying the following filters: <i> True_Var_Def = 1 (i.e., 3/3 replicates), False_Var_Def = 0.33 (i.e., 1/3 replicates), VAF_Lower_Limit = 0.01, VAF_Upper_Limit = 0.6, Feature_Category = Exhaustive. </i>

As you will observe there are 10 rows after applying the above filters, each row reports performance on a separate testing dataset obtained by randomly sampling with a different seed. The overall performance for this combination is obtained by calculating the median of each metric column.

 | Seed | True_Var_Def | False_Var_Def | VAF_Lower_Limit | VAF_Upper_Limit | Feature_Category | Num_True_Var_Testing | Num_False_Var_Testing | TP  | FP  | TN  | FN | AUC      | MCC   | F1_score | Accuracy | MSE   |
| ---- | ------------ | ------------- | --------------- | --------------- | ---------------- | -------------------- | --------------------- | --- | --- | --- | -- | -------- | ----- | -------- | -------- | ----- |
| 10   | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 159                  | 1055                  | 127 | 228 | 827 | 32 | 0.791314 | 0.432 | 0.494    | 0.786    | 0.214 |
| 20   | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 137                  | 1077                  | 112 | 251 | 826 | 25 | 0.792232 | 0.404 | 0.448    | 0.773    | 0.227 |
| 3247 | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 119                  | 1095                  | 95  | 224 | 871 | 24 | 0.796877 | 0.401 | 0.434    | 0.796    | 0.204 |
| 24   | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 126                  | 1088                  | 94  | 249 | 839 | 32 | 0.758586 | 0.35  | 0.401    | 0.769    | 0.231 |
| 4501 | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 144                  | 1070                  | 109 | 231 | 839 | 35 | 0.770528 | 0.39  | 0.45     | 0.781    | 0.219 |
| 9879 | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 136                  | 1078                  | 113 | 240 | 838 | 23 | 0.804124 | 0.422 | 0.462    | 0.783    | 0.217 |
| 878  | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 122                  | 1092                  | 100 | 235 | 857 | 22 | 0.802235 | 0.407 | 0.438    | 0.788    | 0.212 |
| 76   | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 133                  | 1081                  | 102 | 203 | 878 | 31 | 0.789564 | 0.417 | 0.466    | 0.807    | 0.193 |
| 187  | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 135                  | 1079                  | 102 | 246 | 833 | 33 | 0.763783 | 0.367 | 0.422    | 0.77     | 0.23  |
| 299  | [1]          | [0.33]        | 0.01            | 0.6             | Exhaustive       | 136                  | 1078                  | 98  | 210 | 868 | 38 | 0.762892 | 0.381 | 0.441    | 0.796    | 0.204 |



### Example 2: Predicting true SNVs using the raw VCF files
In this example, we will learn how to predict true/false SNVs using the raw VCF files generated by our in-house script. We will use the raw vcf files provided in the data/vcf_files/raw_vcf directory for this purpose.

As the first step, you will need to process the raw vcf files and convert them into a csv file using <i> scripts/split_multi_allelic_var.py </i> script. Make sure your current directory is <i> HPV_low_VAF_SNV_prediction </i> and you have activated the conda environment as described in Step 4 of the <a href="#installation">Installation</a> section. Then use the following commands.

```sh
mkdir example2
cd example2
python ../scripts/split_multi_allelic_var.py --vcf_dir ../data/vcf_files/raw_vcf/ --output_dir raw_vcf_parsed --num_metadata_lines 70
```

You should see the following output (I have truncated the output as it was long!)
```sh
Processing file ../data/vcf_files/raw_vcf/sample10_D4.tvc.vcf (1/93)...Done!
Processing file ../data/vcf_files/raw_vcf/sample10_D5.tvc.vcf (2/93)...Done!
Processing file ../data/vcf_files/raw_vcf/sample10_D6.tvc.vcf (3/93)...Done!
Processing file ../data/vcf_files/raw_vcf/sample11_H1.tvc.vcf (4/93)...Done!
Processing file ../data/vcf_files/raw_vcf/sample11_H2.tvc.vcf (5/93)...Done!
...
...
...
Processing file ../data/vcf_files/raw_vcf/sample9_H7.tvc.vcf (91/93)...Done!
Processing file ../data/vcf_files/raw_vcf/sample9_H8.tvc.vcf (92/93)...Done!
Processing file ../data/vcf_files/raw_vcf/sample9_H9.tvc.vcf (93/93)...Done!
Splitting MNPs...Done!
Extracting only SNVs...Done!
Total SNVs = 9225. Wrote SNVs to raw_vcf_parsed/snvs.csv
```

The total number of SNVs identified is 9,225 and the SNVs and their features have been written to the file <i>raw_vcf_parsed/snvs.csv</i>

Each row in raw_vcf_parsed/snvs.csv is a SNV. We will now consider that the raw_vcf_parsed/snvs.csv file is our hypothetical independent test dataset and our goal is to predict true and false low-VAF SNVs on this dataset. We will achieve this goal by using the <i>**machine_learning/predict_true_low_vaf_snvs.py**</i> script. Checkout the usage of this script as follows.

```sh
python ../machine_learning/predict_true_low_vaf_snvs.py --help

Usage: predict_true_low_vaf_snvs.py [OPTIONS]

Options:
  --snv_file TEXT           The csv file generated by the
                            split_multi_allelic_var.py script.For VM and FVM
                            models, remember to filter the vcf files with
                            VCFgenie before running the
                            split_multi_allelic_var.py script.  [required]
  --model_type [VM|FM|FVM]  The model type to be used for making the
                            prediction  [default: VM; required]
  --hpv_ref_file TEXT       The reference genome file for the HPV type. File
                            must contain the unpadded genome sequence in fasta
                            format (i.e., genome sequence should not include a
                            duplication of bases in the start region.)
                            [required]
  --outfile TEXT            Name of the output file that will contain the
                            prediction results  [required]
  --help                    Show this message and exit.
```

We will make predictions using the FM model. Run the following command and you will see the output that follows the command.

```sh
python ../machine_learning/predict_true_low_vaf_snvs.py --snv_file raw_vcf_parsed/snvs.csv --model_type FM --hpv_ref_file ../data/HPV18_not_padded.fasta --outfile FM_predictions.csv

Extracting SNVs...Done!
Calculating features ...Done!
Running predictions...Done!
Wrote predictions to FM_predictions.csv
```
The prediction output is written to <i> FM_predictions.csv </i>. As you will observe, the script only made predictions for 4,590/9,225 SNVs. That's because the prediction is only performed for SNVs with VAF 1% - 60% (the optimal VAF range for FM models), having a minimum FDP of 100 and FAO of 1. The <b>Predicted_SNV_Type</b> column is a binary column that includes the prediction : TRUE => Predicted true SNVs, FALSE => Predicted false SNVs. The <b> True_SNV_Probability </b> column provides a probability value for the SNV to be true - higher probability means higher likelihood of the SNV to be a true SNV. SNVs with probabilities < 0.5 are labelled as false SNVs.

### Example 3: Predicting true SNVs using processed VCF files generated by VCFgenie
In this example, we will learn how to predict true/false SNVs using the filtered VCF files generated by VCFgenie. We will use the raw vcf files provided in the data/vcf_files/raw_vcf directory for this purpose. 

The first step is to run VCFgenie on the raw vcf files. VCFgenie requires the user to provide a p value threshold as a runtime argument. The p value in this scenario is the Bonferroni-corrected p value and can be calculated by dividing 0.05 by the total number of SNVs (in our case it's 9,225 as calculated by the split_multi_allelic_var.py script). The corrected p value would therefore be <b>0.000005420054201</b>. We will also need the platform error rate - for Ion Torrent it's 0.00012936. We will now use the following commands to install VCFgenie, execute it and generate the filtered VCF files. Before running these commands make sure that you are in the <i> HPV_low_VAF_SNV_prediction </i> directory.

```sh
mkdir example3
cd example3
git clone https://github.com/chasewnelson/VCFgenie.git
# Copy the raw vcf files to the current dir
cp -rf ../data/vcf_files/raw_vcf/ .
cd raw_vcf
python ../VCFgenie/VCFgenie.py --error_rate=0.00012936 --p_cutoff=0.000005420054201 --AC_key=FAO --AC_key_new=FAO --AF_key=AF --AF_key_new=AF --DP_key=FDP --overwrite_INFO --VCF_files sample*.vcf > vcfgenie.out
```

You will notice that once VCFgenie has finished running, a directory *VCFgenie_output* is created that contains the filtered/processed vcf files. The log output from VCFgenie can be found in the file vcfgenie.out.

Next, we will compile a csv file from the vcf files that will only contain SNVs. As in [Example 2](#example-2-predicting-true-snvs-using-the-raw-vcf-files), we will use the <i> split_multi_allelic_var.py </i> script to perform this operation.

```sh
# Go back to the parent directory
cd ..
python ../scripts/split_multi_allelic_var.py --vcf_dir ./raw_vcf/VCFgenie_output --output_dir VCFgenie_SNVs --num_metadata_lines 86
```

Note that we changed the `num_metadata_lines` to 86 compared to example 2 as the VCFgenie-processed VCF files have 16 additional metadata lines.
You should see the following output (I have truncated the output as it was long!)
```sh
Processing file VCFgenie_output/sample10_D4.tvc_filtered.vcf (1/93)...Done!
Processing file VCFgenie_output/sample10_D5.tvc_filtered.vcf (2/93)...Done!
Processing file VCFgenie_output/sample10_D6.tvc_filtered.vcf (3/93)...Done!
Processing file VCFgenie_output/sample11_H1.tvc_filtered.vcf (4/93)...Done!
Processing file VCFgenie_output/sample11_H2.tvc_filtered.vcf (5/93)...Done!
Processing file VCFgenie_output/sample11_H3.tvc_filtered.vcf (6/93)...Done!
...
...
...
Processing file VCFgenie_output/sample8_G7.tvc_filtered.vcf (88/93)...Done!
Processing file VCFgenie_output/sample8_G8.tvc_filtered.vcf (89/93)...Done!
Processing file VCFgenie_output/sample8_G9.tvc_filtered.vcf (90/93)...Done!
Processing file VCFgenie_output/sample9_H7.tvc_filtered.vcf (91/93)...Done!
Processing file VCFgenie_output/sample9_H8.tvc_filtered.vcf (92/93)...Done!
Processing file VCFgenie_output/sample9_H9.tvc_filtered.vcf (93/93)...Done!
Splitting MNPs...Done!
Extracting only SNVs...Done!
Total SNVs = 9225. Wrote SNVs to VCFgenie_SNVs/snvs.csv
```

Finally, we will predict true SNVs with the VM models and observe the output below.
```sh
python ../machine_learning/predict_true_low_vaf_snvs.py --snv_file VCFgenie_SNVs/snvs.csv --model_type VM --hpv_ref_file ../data/HPV18_not_padded.fasta --outfile VM_predictions.csv

Extracting SNVs...Done!
Calculating features ...Done!
Running predictions...Done!
Wrote predictions to VM_predictions.csv
```


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Sambit Mishra - sambit.mishra@nih.gov

Project Link: [https://github.com/NCI-CGR/HPV_low_VAF_SNV_prediction](https://github.com/NCI-CGR/HPV_low_VAF_SNV_prediction.git)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
The work is an outcome of a collaborative effort between the researchers at the Cancer Genomics Research Laboratory ([CGR](https://dceg.cancer.gov/about/organization/cgr)), [Frederick National Laboratory](https://frederick.cancer.gov/) and the principal investigators and postdoctoral fellows at the Division of Cancer Epidemiology and Genetics ([DCEG](https://dceg.cancer.gov/)). I would like to acknowledge the scientific efforts and research guidance received from Dr. Lisa Mirabello, Dr. Meredith Yeager, Dr. Laurie Burdette, Dr. Michael Dean, Dr. Bin Zhu, Dr. Chase Nelson and Dr. Maisa Pinheiro without which this work wouldn't be complete.
