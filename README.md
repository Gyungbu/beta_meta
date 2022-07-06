# Beta-Meta: a meta-analysis application considering heterogeneity among genome-wide association studies

<img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=Python&logoColor=white">

![KakaoTalk_20220608_101041824_01](https://user-images.githubusercontent.com/106565330/177439137-cde237e1-abbf-4ff2-902e-6033e27f2621.png)

Beta-Meta is a meta-analysis application considering heterogeneity among GWAS studies. It provides a step-by-step meta-analysis of GWAS in the following order: heterogeneity test, two different calculations of an effect size and a p-value based on heterogeneity, and the Benjamini-Hochberg (BH) p-value adjustment. It uses the pandas library to deal with dataframes from excel input files.

Beta-Meta consists of two versions: python script version and exe application.

Please just download the folder below that corresponds to the version you want.
	
	1. python script version : 'beta_meta_script'
	
	2. python exe application version : 'beta_meta_exe'

The list of required packages for `python script version` is shown in the `requirements.txt` file. When you run the command below, the package will be downloaded.

	$ pip install -r requirements.txt

## How to use

When Beta-Meta is executed as follows, the file `meta_output.xlsx` is created or modified in the `output` folder by the files in the `input` floder.

1. python script version

When you run the command below, the code will be executed.

	$ python beta_meta.py

2. python exe application version

When you click the `BetaMeta.exe` file, the application will be executed. 


## Data Description

### <Input file - Precautions>
1. All values (`PHENOTYPE`, `SNP`, `EFFECT_ALLELE`, `NON_EFFECT_ALLELE`, `P_VAL`) should be written.
2. (`BETA`, `BETA_SE`) or (`OR`, `OR_95%CI_LOWER`, `OR_95%CI_UPPER`) values must be written in the input file.
3. For the same phenotype and SNP to be meta-analyzed, the values of the `PHENOTYPE` and `SNP` columns must have the same spacing and spelling.


### <Input file - Column Description>

* `PHENOTYPE` : Phenotype of the individual study.	
* `SNP` : SNP of the individual study.
* `EFFECT_ALLELE` : Effect allele of the individual study.	
* `NON_EFFECT_ALLELE` : Non-effect allele of the individual study.	
* `BETA` : Effect size of the individual study.	
* `BETA_SE` : Standard error of the Beta of the individual study.	
* `OR` : Odds ratio of the individual study.	
* `OR_95%CI_LOWER` : Lower bound of 95% Confidence interval of Odds ratio	
* `OR_95%CI_UPPER` : Upper bound of 95% Confidence interval of Odds ratio  	
* `P_VAL` : P-value of the individual study.

### <Output file - Column Description>
	
* `PHENOTYPE` : Phenotype of the integrated study.
* `SNP` : SNP of the integrated study.
* `EFFECT_ALLELE` : Effect allele of the study with the lowest p value for a specific phenotype, SNP
* `NON_EFFECT_ALLELE` : Non-effect allele of the study with the lowest p value for a specific phenotype, SNP
* `BETA` : Weighted average of the Effect sizes for specific phenotypes, SNPs.
* `BETA_SE` : Standard error of the integrated beta for specific phenotypes, SNPs.
* `P_VAL` : Meta P-value for specific phenotypes, SNPs..
* `BH_P_VAL` : BH adjusted Meta p-value for specific phenotypes, SNPs..
* `I_SQUARE` : Higgin’s heterogeneity metric for specific phenotypes, SNPs.
* `Q_HET` : Cochran’s Q statistic for specific phenotypes, SNPs.
