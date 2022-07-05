# Beta-Meta: a meta-analysis application considering heterogeneity among genome-wide association studies

Beta-Meta is a meta-analysis application considering heterogeneity among GWAS studies. It uses the pandas library to deal with dataframes from excel input files.

Beta-Meta consists of two versions: python script version and exe application.

Please just download the folder below that corresponds to the version you want.
	
	1. python script version : 'beta_meta_script'
	
	2. python exe application version : 'beta_meta_exe'

The list of required packages for 'python script version' is shown in the 'requirements.txt' file. When you run the command below, the package will be downloaded.

	$ pip install -r 'The Absolute path of beta_meta.py file'

## How to use

When Beta-Meta is executed in the following way, the 'meta_output.xlsx' file in the 'output' folder is created or modified.

1. python script version

When you run the command below, the code will be executed.

	$ find -name beta_meta.py | xargs python

	or

	$ python 'The Absolute path of beta_meta.py file'

2. python exe application version

When you click the 'BetaMeta.exe' file, the application will be executed. 


## Data Description

### <Input file - Precautions>
1. All values ('PHENOTYPE', 'SNP', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'P_VAL') should be written.
2. ('BETA', 'BETA_SE') or ('OR', 'OR_95%CI_LOWER', 'OR_95%CI_UPPER') values must be written in the input file.
3. For the same phenotype and SNP to be meta-analyzed, the values of the 'PHENOTYPE' and 'SNP' columns must have the same spacing and spelling.


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
	
	
<span style="color: #ff0000;">빨강</span>,
