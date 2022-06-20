import os
import pandas as pd
import numpy as np
from scipy.stats import norm
    
# Input Folder 
# file_list : List of File name in the Input Folder

path_meta_data_dir = os.path.dirname(os.path.abspath(__file__)) + "/input/"
file_list = os.listdir(path_meta_data_dir)

# Column Extraction / Concat the Dataframes / Deduplication
# li_column : List of column name in the Input File
# df_meta : Data frame of Input Files to be Meta-Analyzed

li_column = ['PHENOTYPE', 'SNP', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'BETA', 'BETA_SE', 'OR', 'OR_95%CI_LOWER', 'OR_95%CI_UPPER','P_VAL']
df_meta = pd.DataFrame(columns = li_column)

for i in range(len(file_list)):
  path_meta_data = path_meta_data_dir + file_list[i]
  df_meta_data = pd.read_excel(path_meta_data)
  df_meta_data = df_meta_data.loc[:, li_column]
  df_meta_data = df_meta_data.replace({'EFFECT_ALLELE':{1:'A', 2:'C', 3:'G', 4:'T'},'NON_EFFECT_ALLELE':{1:'A', 2:'C', 3:'G', 4:'T'}})

  df_meta = pd.concat([df_meta, df_meta_data])

df_meta = df_meta.drop_duplicates(li_column)

# Calculation - BETA & SE
# li_BETA : List of Beta corresponding to df_meta
# li_BETA_SE : List of Beta Standard Error corresponding to df_meta

li_BETA = []
li_BETA_SE = []

for idx, row in df_meta.iterrows():
  if (np.isnan(row['OR']) == False) & (np.isnan(row['OR_95%CI_LOWER']) == False) & (np.isnan(row['OR_95%CI_UPPER']) == False):
    Beta = np.log(row['OR'])
    Beta_SE = (np.log(row['OR_95%CI_UPPER']) - np.log(row['OR_95%CI_LOWER']))/3.92 
    
    li_BETA.append(Beta)
    li_BETA_SE.append(Beta_SE)
    
  else:
    li_BETA.append(row['BETA'])
    li_BETA_SE.append(row['BETA_SE'])
    
df_meta['BETA'] = li_BETA
df_meta['BETA_SE'] = li_BETA_SE

# Create the list [PHENOTYPE, SNP]
# li_PHENOTYPE_SNP : List of Phenotype & SNP 

li_PHENOTYPE_SNP = []

for idx, row in df_meta.iterrows(): 
  if [row['PHENOTYPE'], row['SNP']] not in li_PHENOTYPE_SNP:
    li_PHENOTYPE_SNP.append([row['PHENOTYPE'], row['SNP']])

# Effect Direction Correction
# li_EFFECT_ALLELE : List of Effect Allele corresponding to li_PHENOTYPE_SNP
# li_NON_EFFECT_ALLELE : List of Non-Effect Allele corresponding to li_PHENOTYPE_SNP

li_EFFECT_ALLELE = []
li_NON_EFFECT_ALLELE = []

for i in range(len(li_PHENOTYPE_SNP)):
  condition = (df_meta.SNP == li_PHENOTYPE_SNP[i][1]) & (df_meta.PHENOTYPE == li_PHENOTYPE_SNP[i][0]) 
  
  EFFECT_ALLELE = df_meta[condition]['EFFECT_ALLELE'].values[0]
  NON_EFFECT_ALLELE = df_meta[condition]['NON_EFFECT_ALLELE'].values[0]
    
  condition2 = (df_meta.SNP == li_PHENOTYPE_SNP[i][1]) & (df_meta.PHENOTYPE == li_PHENOTYPE_SNP[i][0]) & (df_meta['NON_EFFECT_ALLELE'] == EFFECT_ALLELE) & (df_meta['EFFECT_ALLELE'] == NON_EFFECT_ALLELE)
    
  df_meta.loc[condition2, 'EFFECT_ALLELE'] = EFFECT_ALLELE
  df_meta.loc[condition2, 'NON_EFFECT_ALLELE'] = NON_EFFECT_ALLELE
  df_meta.loc[condition2, 'BETA'] *= -1
  
  li_EFFECT_ALLELE.append(EFFECT_ALLELE)
  li_NON_EFFECT_ALLELE.append(NON_EFFECT_ALLELE)

# Calculation - Weighted average of the effect sizes 
# li_BETA_META : List of Weighted average of the effect sizes corresponding to li_PHENOTYPE_SNP
# li_STD_BETA_META : List of Standard Error of Weighted average of the effect sizes corresponding to li_PHENOTYPE_SNP

li_BETA_META = []
li_STD_BETA_META = []

for i in range(len(li_PHENOTYPE_SNP)):
  
  condition = (df_meta.SNP == li_PHENOTYPE_SNP[i][1]) & (df_meta.PHENOTYPE == li_PHENOTYPE_SNP[i][0]) 
  sum_w_i = 0
  sum_w_i_beta_i = 0
  
  if len(df_meta[condition]) > 1:
    for idx, row in df_meta[condition].iterrows():
      w_i = row['BETA_SE']**(-2)
      sum_w_i += w_i
      sum_w_i_beta_i += (w_i * row['BETA'])
    li_BETA_META.append(sum_w_i_beta_i/sum_w_i)
    li_STD_BETA_META.append(sum_w_i**-0.5)
  
  elif len(df_meta[condition]) == 1:
    li_BETA_META.append(df_meta[condition]['BETA'].values[0])
    li_STD_BETA_META.append(df_meta[condition]['BETA_SE'].values[0])
 
# Calculation - Cochran's Q statistic
# li_Q : List of Cochran's Q statistic corresponding to li_PHENOTYPE_SNP

li_Q = []

for i in range(len(li_PHENOTYPE_SNP)):
  
  condition = (df_meta.SNP == li_PHENOTYPE_SNP[i][1]) & (df_meta.PHENOTYPE == li_PHENOTYPE_SNP[i][0]) 
  Q = 0
  
  if len(df_meta[condition]) > 1:
    for idx, row in df_meta[condition].iterrows():
      
      w_i = row['BETA_SE']**(-2)
      delta_beta = (row['BETA'] - li_BETA_META[i])    
      Q += (w_i * (delta_beta**2))
    li_Q.append(Q)
  
  elif len(df_meta[condition]) == 1:
    li_Q.append(0)

# Calculation - Higgin's heterogeneity metric
# li_I_square : List of Higgin's heterogeneity metric corresponding to li_PHENOTYPE_SNP

li_I_square = []

for i in range(len(li_PHENOTYPE_SNP)):
  
  condition = (df_meta.SNP == li_PHENOTYPE_SNP[i][1]) & (df_meta.PHENOTYPE == li_PHENOTYPE_SNP[i][0]) 
  count = 0
  
  if (len(df_meta[condition]) > 1) & (li_Q[i] != 0):
    for idx, row in df_meta[condition].iterrows():
      count += 1
    
    I_square = 100 * (li_Q[i] - (count - 1)) / (li_Q[i])
    li_I_square.append(max(I_square,0))
    
  else: 
    li_I_square.append('No Meta')

# Calculation - Weighted average of the effect sizes Modification - Random Effect Model
# li_BETA_META : List of Weighted average of the effect sizes corrected by a Random Effect Model
# li_STD_BETA_META : List of Standard Error of Weighted average of the effect sizes corrected by a Random Effect Model

for i in range(len(li_PHENOTYPE_SNP)):
  
  if type(li_I_square[i]) == float: 
    condition = (df_meta.SNP == li_PHENOTYPE_SNP[i][1]) & (df_meta.PHENOTYPE == li_PHENOTYPE_SNP[i][0]) 
    
    if (li_I_square[i] >= 50) & (len(df_meta[condition]) > 1):
      sum_w_i = 0
      sum_w_i_square = 0
      count = 0
      
      for idx, row in df_meta[condition].iterrows():
        w_i = row['BETA_SE']**(-2)
        sum_w_i += w_i
        sum_w_i_square += w_i**2
        count += 1
      
      tau_square = (li_Q[i]-count+1) / (sum_w_i - sum_w_i_square / sum_w_i)                   
      tau_square = max(tau_square, 0)
      
      sum_w_i_R_beta_i = 0
      sum_w_i_R = 0
      for idx2, row2 in df_meta[condition].iterrows():
        w_i = row2['BETA_SE']**(-2)
        w_i_R = 1/ (1/w_i + tau_square)
        sum_w_i_R_beta_i += w_i_R * row2['BETA']
        sum_w_i_R += w_i_R
        
      li_BETA_META[i] = sum_w_i_R_beta_i / sum_w_i_R
      li_STD_BETA_META[i] = sum_w_i_R**-0.5

# Calculation - P-value 
# li_p_value : List of P-value corresponding to li_PHENOTYPE_SNP

li_p_value = []

for i in range(len(li_PHENOTYPE_SNP)):
  
  Z = li_BETA_META[i] / li_STD_BETA_META[i]
  p_value = 2 * norm.cdf(-abs(Z))
  
  li_p_value.append(p_value)

# Output file - Meta Analysis
# df_meta_output : Data Frame of Meta-analysis Result Ouput File

df_meta_output = pd.DataFrame(columns = ['PHENOTYPE', 'SNP', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE','BETA', 'BETA_SE', 'P_VAL', 'BH_P_VAL', 'I_SQUARE', 'Q_HET'])

for i in range(len(li_PHENOTYPE_SNP)):
  df_meta_output.loc[i, 'PHENOTYPE'] = li_PHENOTYPE_SNP[i][0]
  df_meta_output.loc[i, 'SNP'] = li_PHENOTYPE_SNP[i][1]

df_meta_output['BETA'] = li_BETA_META
df_meta_output['BETA_SE'] = li_STD_BETA_META
df_meta_output['P_VAL'] = li_p_value
df_meta_output['EFFECT_ALLELE'] = li_EFFECT_ALLELE
df_meta_output['NON_EFFECT_ALLELE'] = li_NON_EFFECT_ALLELE
df_meta_output['I_SQUARE'] = li_I_square
df_meta_output['Q_HET'] = li_Q  

# P-value Correction - BH adjustment
# df_meta_output : Data Frame of Meta-analysis Result Ouput File

for idx, row in df_meta_output.iterrows():
  condition = df_meta_output.PHENOTYPE == row['PHENOTYPE']

  ascending_rank_p_val = df_meta_output.groupby('PHENOTYPE')['P_VAL'].rank(method='min', ascending=True).values[idx]
  df_meta_output.loc[idx, 'BH_P_VAL'] = df_meta_output.loc[idx, 'P_VAL'] * len(df_meta_output[condition]) / ascending_rank_p_val
  
path_meta_output = os.path.dirname(os.path.abspath(__file__)) + "/output/meta_output.xlsx" 
df_meta_output.to_excel(path_meta_output)

  
