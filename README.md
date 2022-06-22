# beta_meta
Meta-analysis with Heterogeneity Test

The whole beta_meta folder should be downloaded.

The list of required packages was shown in the requiremets.txt file.

When you run the command below, the code will be executed.

$ find -name beta_meta.py | xargs python

or

$ python 'Absolute path of beta_meta.py'

### <Input file - Precautions>
1. All values ('PHENOTYPE', 'SNP', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'P_VAL') should be written.
2. ('BETA', 'BETA_SE') or ('OR', 'OR_95%CI_LOWER', 'OR_95%CI_UPPER') values must be written in the input file.
3. For the same phenotype and SNP to be meta-analyzed, the values of the 'PHENOTYPE' and 'SNP' columns must have the same spacing and spelling.


### <Input file - Column Description>
1. PHENOTYPE : Phenotype of the individual study.
2. SNP : SNP of the individual study.
3. EFFECT_ALLELE : Effect allele of the individual study.
4. NON_EFFECT_ALLELE : Non-effect allele of the individual study.
5. BETA : Effect size of the individual study.
6. BETA_SE : Standard error of the Beta of the individual study.
7. OR : Odds ratio of the individual study.
8. OR_95%CI_LOWER : Lower bound of 95% Confidence interval of Odds ratio
9. OR_95%CI_UPPER : Upper bound of 95% Confidence interval of Odds ratio  
10. P_VAL : P-value of the individual study.

### <Output file - Column Description>
1. PHENOTYPE : Phenotype of the integrated study.
2. SNP : SNP of the integrated study.
3. EFFECT_ALLELE : Effect allele of the study with the lowest p value for a specific phenotype, SNP
4. NON_EFFECT_ALLELE : Non-effect allele of the study with the lowest p value for a specific phenotype, SNP
5. BETA : Weighted average of the Effect sizes for specific phenotypes, SNPs.
6. BETA_SE : Standard error of the integrated beta for specific phenotypes, SNPs.
7. P_VAL : Meta P-value for specific phenotypes, SNPs..
8. BH_P_VAL : BH adjusted Meta p-value for specific phenotypes, SNPs..
9. I_SQUARE : Higgin’s heterogeneity metric for specific phenotypes, SNPs.
10. Q_HET : Cochran’s Q statistic for specific phenotypes, SNPs.
