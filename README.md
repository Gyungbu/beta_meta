# Beta-Meta: a meta-analysis application considering heterogeneity among genome-wide association studies

![Python](https://img.shields.io/badge/Python-v3.9.7-blue.svg?style=flat&logo=python)&nbsp;
![Pandas](https://img.shields.io/badge/pandas-v1.4.2-blue.svg?style=flat&logo=pandas)&nbsp;
![Numpy](https://img.shields.io/badge/NumPy-v1.22.4-blue.svg?style=flat&logo=numpy)&nbsp;
![Scipy](https://img.shields.io/badge/SciPy-v1.8.1-blue.svg?style=flat&logo=scipy)&nbsp;
![scikit-learn](https://img.shields.io/badge/scikit--learn-v1.1.2-blue.svg?style=flat&logo=scikit-learn)&nbsp;
![GitHub](https://img.shields.io/badge/GitHub-grey.svg?style=flat&logo=github)&nbsp;
![Exe](https://img.shields.io/badge/exe-grey.svg?style=flat&logo=exe)&nbsp;

![beta_meta_icon](https://user-images.githubusercontent.com/106565330/177929679-992af204-532d-4a59-8330-27c37ed96208.png)

Beta-Meta is a meta-analysis application considering heterogeneity among GWAS studies. It provides a step-by-step meta-analysis of GWAS in the following order: heterogeneity test, two different calculations of an effect size and a p-value based on heterogeneity, and the Benjamini-Hochberg (BH) p-value adjustment. It uses the pandas library to handle dataframes from an excel input file and only requires the single file to conduct a meta-analysis.

## Installation

Beta-Meta comes in two versions: `script` version and `exe` version.

You can install the Beta-Meta with following command.
	
	git clone https://github.com/Gyungbu/beta_meta.git

The list of required packages for `script` is shown in the `requirements.txt` file. When you run the command below, the required packages will be downloaded. (version : `python 3.9.7`)
	
	conda create -n env_beta
	conda activate env_beta
	conda install pip  
	conda install python=3.9.7
	pip install -r ./beta_meta/requirements.txt

# Beta_Meta_LD : (Optional) Suggest Correlated SNP with LD=1
## How to use

### 1. Prepare Input data
Place the `input_SNPs.txt` file of your input data in the `input` folder (`./beta_meta/script/beta_meta_LD_script/input/` or `./beta_meta/exe/beta_meta_LD_exe/input/` for script and exe version respectively).

Caveats: 

1. Write down rs_num to search for correlated SNPs.
2. Separate rs_num by enter.
3. You have to install the R program. 

### 2. Run Beta_Meta_LD
To run Beta_Meta_LD,

- For script version:
    
    Run the command below:

		python ./beta_meta/script/beta_meta_LD_script/haploR.py
    
- For exe version:
    
    When you double-click the `./beta_meta/exe/beta_meta_LD_exe/haploR.exe` file, the application will be executed.

When Beta_Meta_LD is executed as above, the file `correlated_with_{rs_num}.txt` will be created or modified in the `output` folder (`./beta_meta/script/beta_meta_LD_script/output/` or `./beta_meta/exe/beta_meta_LD_exe/output/` for script and exe version respectively).

# Beta_Meta : Meta-analysis 
## How to use

### 1. Prepare Input data
Place the excel file of your input data in the `input` folder (`./beta_meta/script/beta_meta_script/input/` or `./beta_meta/exe/beta_meta_exe/input/` for script and exe version respectively).

Caveats: 

1. All values of (`PHENOTYPE`, `SNP`, `EFFECT_ALLELE`, `NON_EFFECT_ALLELE`, `P_VAL`) should be written in the input file.
2. Either (`BETA`, `BETA_SE`) or (`OR`, `OR_95%CI_LOWER`, `OR_95%CI_UPPER`) values must be also written.
3. As Beta-Meta calculates SNP-phenotype associations separately, it is acceptable to include as many phenotypes as desired in the single input file. The ones (in `PHENOTYPE` and `SNP` columns) to be integrated must have exactly the same spelling (including spaces) as Beta-Meta matches exact phenotype and SNP.
4. Beta-Meta deals with strand flipping and provides the direction of effect size relative to the same allele. When the effect and the non-effect allele are inverted between the individual studies, this can also be resolved automatically by changing the sign of the normalized effect.
5. When only one study for a certain SNP-phenotype association is provided, ‘No Meta’ will be shown in the `I_SQUARE`, `Q_HET` columns in the output file.
6. If above not possible, Beta-Meta will remove them from file and will not conduct a meta-analysis for them.

### 2. Run Beta_Meta
To run Beta_Meta,

- For script version:
    
    Run the command below:

		python ./beta_meta/script/beta_meta_script/beta_meta.py
    
- For exe version:
    
    When you double-click the `./beta_meta/exe/beta_meta_exe/beta_meta.exe` file, the application will be executed.
    

When Beta-Meta is executed as above, the file `meta_output.xlsx` & `meta_forestplot.png` will be created or modified in the `output` folder (`./beta_meta/script/beta_meta_script/output/` or `./beta_meta/exe/beta_meta_exe/output/` for script and exe version respectively).

## Data Description

### <Input file - Column Description>

* `PHENOTYPE` (str) : Phenotype 	
* `SNP` (str) : SNP 
* `EFFECT_ALLELE` (str) : Effect allele of SNP	
* `NON_EFFECT_ALLELE` (str) : Non-effect allele of SNP	
* `BETA` (float) : Effect size of SNP from an individual study
* `BETA_SE` (float) : Standard error of the beta from an individual study	
* `OR` (float) : Odds ratio for an association from an individual study	
* `OR_95%CI_LOWER` (float) : Lower bound of 95% confidence interval of odds ratio	
* `OR_95%CI_UPPER` (float) : Upper bound of 95% confidence interval of odds ratio 	
* `P_VAL` (float) : P-value of the effect size

### <Output file - Column Description>
	
* `PHENOTYPE` : Phenotype
* `SNP` : SNP 
* `EFFECT_ALLELE` : Effect allele from our beta-meta analysis with the lowest BH adjusted p-value for each SNP-phenotype association
* `NON_EFFECT_ALLELE` : Non-effect allele from our beta-meta analysis with the lowest BH adjusted p-value for each phenotype-SNP association
* `BETA` : Weighted average of the effect sizes
* `BETA_SE` : Standard error of the integrated beta
* `P_VAL` : Meta p-value
* `BH_P_VAL` : BH adjusted meta p-value
* `I_SQUARE` : Higgin’s heterogeneity metric
* `Q_HET` : Cochran’s Q statistic
