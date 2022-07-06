# Beta-Meta: a meta-analysis application considering heterogeneity among genome-wide association studies

![Python](https://img.shields.io/badge/-Python-05122A?style=flat&logo=python)&nbsp;
![Pandas](https://img.shields.io/badge/-Pandas-05122A?style=flat&logo=pandas)&nbsp;
![Numpy](https://img.shields.io/badge/-Numpy-05122A?style=flat&logo=numpy)&nbsp;
![Scipy](https://img.shields.io/badge/-Scipy-05122A?style=flat&logo=scipy)&nbsp;
![SVG](https://img.shields.io/badge/-SVG-05122A?style=flat&logo=svg)&nbsp;
![GitHub](https://img.shields.io/badge/-GitHub-05122A?style=flat&logo=github)&nbsp;

![beta_meta_icon](https://user-images.githubusercontent.com/106565330/177439764-bd1d00c4-a553-4ec0-b75a-3222bef155c6.png)

Beta-Meta is a meta-analysis application considering heterogeneity among GWAS studies. It provides a step-by-step meta-analysis of GWAS in the following order: heterogeneity test, two different calculations of an effect size and a p-value based on heterogeneity, and the Benjamini-Hochberg (BH) p-value adjustment. It uses the pandas library to handle dataframes from an excel input file and only requires the single file to conduct a meta-analysis.

## Download

Beta-Meta comes in two versions: `beta_meta_script` for Linux and `beta_meta_exe` for Windows.

It is not necessary to download all of the folders above to run Beta-Meta (*but it is fine if you do*). 

For Windows users, you may download `beta_meta_exe` folder only. 

For Linux users, you may download `beta_meta_script` folder and `requirements.txt` only.

The list of required packages for `beta_meta_script` is shown in the `requirements.txt` file. When you run the command below, the required packages will be downloaded.

	$ pip install -r requirements.txt

## How to use

### 1. Prepare Input data
Place the excel file of your input data in the `input` folder.

Caveats: 

1. All values of (`PHENOTYPE`, `SNP`, `EFFECT_ALLELE`, `NON_EFFECT_ALLELE`, `P_VAL`) should be written in the input file.
2. Either (`BETA`, `BETA_SE`) or (`OR`, `OR_95%CI_LOWER`, `OR_95%CI_UPPER`) values must be also written.
3. As Beta-Meta calculates SNP-phenotype associations separately, it is acceptable to include as many phenotypes as desired in the single input file. The ones (in `PHENOTYPE` and `SNP` columns) to be integrated must have exactly the same spelling (including spaces) as Beta-Meta matches exact phenotype and SNP.
4. Beta-Meta deals with strand flipping and provides the direction of effect size relative to the same allele. When the effect and the non-effect allele are inverted between the individual studies, this can also be resolved automatically by changing the sign of the normalized effect.
5. When only one study for a certain SNP-phenotype association is provided, ‘No Meta’ will be shown in the `I_SQUARE`, `Q_HET` columns in the output file.
6. If above not possible, Beta-Meta will remove them from file and will not conduct a meta-analysis for them.

### 2. Run Beta-Meta
To run Beta-Meta,

- For Linux:
    
    Run the command below:

		$ python beta_meta.py
    
- For Windows:
    
    When you double-click the `BetaMeta.exe` file, the application will be executed.
    

When Beta-Meta is executed as above, the file `meta_output.xlsx` will be created or modified in the `output` folder.

## Data Description

### <Input file - Column Description>

* `PHENOTYPE` : Phenotype 	
* `SNP` : SNP 
* `EFFECT_ALLELE` : Effect allele of SNP	
* `NON_EFFECT_ALLELE` : Non-effect allele of SNP	
* `BETA` : Effect size of SNP from an individual study
* `BETA_SE` :Standard error of the beta from an individual study	
* `OR` : Odds ratio for an association from an individual study	
* `OR_95%CI_LOWER` : Lower bound of 95% confidence interval of odds ratio	
* `OR_95%CI_UPPER` : Upper bound of 95% confidence interval of odds ratio 	
* `P_VAL` : P-value of the effect size

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
