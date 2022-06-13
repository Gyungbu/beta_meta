# genome_meta_beta
Meta-analysis with Heterogeneity Test

The whole meta_data folder should be downloaded.

The list of required packages was displayed in the requiremets.txt file.

When you run the command below, the code will be executed.

$ find -name meta_beta.py | xargs python

or

$ python 'Absolute path of meta_beta.py'

### <Input file - Precautions>
1. All values ('PHENOTYPE', 'SNP', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'P_VAL') should be written.
2. ('BETA', 'BETA_SE') or ('OR', 'OR_95%CI_LOWER', 'OR_95%CI_UPPER') values must be written in the input file.
3. For meta-analysis, all values in the input file must have the same spacing and spelling.
4. In the case of SNPs with the same phenotype, the 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE' must be the same. (순서는 상관 없음)
5. When the Effect allele and the Non-effect allele were reversed, the sign of the normalized effect β_i was calculated to change.

### <Input file - Column Description>
1. BETA : Effect size of the individual study.
2. BETA_SE : Standard error of the Beta of the individual study.
3. OR : Odds ratio of the individual study.
4. OR_95%CI_LOWER : Lower bound of 95% Confidence interval of Odds ratio
5. OR_95%CI_UPPER : Upper bound of 95% Confidence interval of Odds ratio  
6. P_VAL : P-value of the individual study.

### <Output file - Column Description>
1. BETA : Weighted average of the Effect size.
2. BETA_SE : Standard error of the Beta.
3. P_VAL : Meta P-value.
4. BH_P_VAL : BH adjusted Meta p-value.
5. I_SQUARE : Higgin’s heterogeneity metric
6. Q_HET : Cochran’s Q statistic
