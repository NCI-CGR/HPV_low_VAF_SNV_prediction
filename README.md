
# Predicting True Low-VAF SNVs in the Human Papillomavirus 


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
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


<!-- GETTING STARTED -->
## Getting Started

The code in this repository required for training and testing the machine learning models is written in Python (v3.8) You will need access to a Linux terminal (even a Windows 10 having Windows Subsystem for Linux enabled and Ubuntu installed) to be able to execute the code. The specific instructions below work best in a Linux (Ubuntu 18.02/Ubuntu20.04) platform. Follow the steps below to install the required packages and run the scripts.


### Installation

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

5. Make sure you are good to go
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



<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_




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
