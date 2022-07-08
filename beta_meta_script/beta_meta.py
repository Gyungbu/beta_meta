import os
import pandas as pd
import numpy as np
from scipy.stats import norm
 
# Function - Check the direction of effect between meta-studies

def sign_effect_direction(effect_allele_study1, non_effect_allele_study1, effect_allele_study2, non_effect_allele_study2):
  """
  Check the direction of effect between meta-studies
  
      Args:
          effect_allele_study1 (str): Effect allele of the study 1 for the specific phenotype and SNP
          non_effect_allele_study1 (str): Non effect allele of the study 1 for the specific phenotype and SNP
          effect_allele_study2 (str): Effect allele of the study 2 for the specific phenotype and SNP
          non_effect_allele_study2 (str): Non effect allele of the study 2 for the specific phenotype and SNP
      
      Returns:
          result (int): Same direction (result = 1) / Opposite direction (result = -1) / Effect alleles and non-effect alleles do not match between meta-studies (result = 0)
  """
  result = 0
  
  if set([effect_allele_study1, non_effect_allele_study1, effect_allele_study2, non_effect_allele_study2]).issubset(set(['A', 'G', 'T', 'C'])):
    set_allele_1 = set([effect_allele_study1, non_effect_allele_study1])
    set_allele_2 = set([effect_allele_study2, non_effect_allele_study2])
  
    if len(set_allele_1.difference(set_allele_2)) == 1:
      result = 0
                   
    elif len(set_allele_1.difference(set_allele_2)) == 0:
      if effect_allele_study1 == effect_allele_study2:
        result = 1
      
      else:
        result = -1
          
    elif len(set_allele_1.difference(set_allele_2)) == 2:
      if (set_allele_2 == set(['A', 'T'])) | (set_allele_2 == set(['G', 'C'])):
        result = 0
         
      elif (set([effect_allele_study1, effect_allele_study2]) == set(['A', 'T'])) | (set([effect_allele_study1, effect_allele_study2]) == set(['G', 'C'])):
        result = 1
       
      else:
        result = -1
                     
  return result
  
"""
# Test the Function - "sign_effect_direction"

print(sign_effect_direction('T', 'C', 'T', 'A'), 'result:0')

print(sign_effect_direction('A', 'G', 'A', 'G'), 'result:1')
print(sign_effect_direction('A', 'G', 'T', 'C'), 'result:1')
print(sign_effect_direction('A', 'G', 'C', 'T'), 'result:-1')

print(sign_effect_direction('A', 'T', 'G', 'C'), 'result:0')
print(sign_effect_direction('A', 'T', 'A', 'T'), 'result:1')
print(sign_effect_direction('A', 'T', 'T', 'A'), 'result:-1')

print(sign_effect_direction('A', 'C', 'T', 'G'), 'result:1')
print(sign_effect_direction('A', 'C', 'G', 'T'), 'result:-1')

print(sign_effect_direction('G', 'C', 'G', 'C'), 'result:1')
print(sign_effect_direction('G', 'C', 'C', 'G'), 'result:-1')

print(sign_effect_direction('-', 'A', 'G', 'A'), 'result:0')
"""

# Input Folder 
# file_list : List of File name in the Input Folder

path_meta_data_dir = os.path.dirname(os.path.abspath(__file__)) + "/input/"
file_list = os.listdir(path_meta_data_dir)

# Column Extraction / Concat the Dataframes 
# li_column : List of column name in the Input File
# df_meta_input : Data frame of Input Files to be Meta-Analyzed

li_column = ['PHENOTYPE', 'SNP', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'BETA', 'BETA_SE', 'OR', 'OR_95%CI_LOWER', 'OR_95%CI_UPPER', 'P_VAL']
df_meta_input = pd.DataFrame(columns = li_column)

for file in file_list:
  path_meta_data = path_meta_data_dir + file
  df_meta_input_data = pd.read_excel(path_meta_data)
  
  if set(li_column).issubset(set(list(df_meta_input_data.columns))):
    df_meta_input_data = df_meta_input_data.loc[:, li_column]
    df_meta_input_data = df_meta_input_data.replace({'EFFECT_ALLELE':{1:'A', 2:'C', 3:'G', 4:'T'},'NON_EFFECT_ALLELE':{1:'A', 2:'C', 3:'G', 4:'T'}})
  
    df_meta_input = pd.concat([df_meta_input, df_meta_input_data])
    
  else:
    print('Please check the columns of the input file! (Ex, PHENOTYPE, SNP, EFFECT_ALLELE, NON_EFFECT_ALLELE, BETA, BETA_SE, OR, OR_95%CI_LOWER, OR_95%CI_UPPER, P_VAL)')  

# Remove Spaces / Deduplication - df_meta_input

df_meta_input['PHENOTYPE'] = df_meta_input['PHENOTYPE'].str.strip()
df_meta_input['SNP'] = df_meta_input['SNP'].str.strip()
df_meta_input['EFFECT_ALLELE'] = df_meta_input['EFFECT_ALLELE'].str.strip()
df_meta_input['NON_EFFECT_ALLELE'] = df_meta_input['NON_EFFECT_ALLELE'].str.strip()

df_meta_input = df_meta_input.drop_duplicates(li_column)

# Calculation - BETA & BETA_SE from Odds ratio & 95% CI
# li_BETA : List of Beta corresponding to df_meta_input
# li_BETA_SE : List of Beta Standard Error corresponding to df_meta_input

li_BETA = []
li_BETA_SE = []

for idx, row in df_meta_input.iterrows():
  if (np.isnan(row['OR']) == False) & (np.isnan(row['OR_95%CI_LOWER']) == False) & (np.isnan(row['OR_95%CI_UPPER']) == False):
    Beta = np.log(row['OR'])
    Beta_SE = (np.log(row['OR_95%CI_UPPER']) - np.log(row['OR_95%CI_LOWER']))/3.92 
    
    li_BETA.append(Beta)
    li_BETA_SE.append(Beta_SE)
    
  else:
    li_BETA.append(row['BETA'])
    li_BETA_SE.append(row['BETA_SE'])
    
df_meta_input['BETA'] = li_BETA
df_meta_input['BETA_SE'] = li_BETA_SE

# Remove the NaN Data & Unnecessary Columns - df_meta_input

df_meta_input = df_meta_input.drop(['OR', 'OR_95%CI_LOWER', 'OR_95%CI_UPPER'], axis=1)
df_meta_input = df_meta_input.dropna()

# Create the list [PHENOTYPE, SNP]
# li_PHENOTYPE_SNP : List of Phenotype & SNP 

li_PHENOTYPE_SNP = []

for idx, row in df_meta_input.iterrows(): 
  if [row['PHENOTYPE'], row['SNP']] not in li_PHENOTYPE_SNP:
    li_PHENOTYPE_SNP.append([row['PHENOTYPE'], row['SNP']])

# Effect Direction Correction
# li_EFFECT_ALLELE : List of Effect Allele corresponding to li_PHENOTYPE_SNP
# li_NON_EFFECT_ALLELE : List of Non-Effect Allele corresponding to li_PHENOTYPE_SNP

df_meta_input['P_VAL'] = df_meta_input['P_VAL'].astype('float64')

li_EFFECT_ALLELE = []
li_NON_EFFECT_ALLELE = []

for j in range(len(li_PHENOTYPE_SNP)):
  condition = (df_meta_input.SNP == li_PHENOTYPE_SNP[j][1]) & (df_meta_input.PHENOTYPE == li_PHENOTYPE_SNP[j][0]) 
  
  idx_min = df_meta_input[condition]['P_VAL'].idxmin()
  EFFECT_ALLELE = df_meta_input.loc[idx_min]['EFFECT_ALLELE']
  NON_EFFECT_ALLELE = df_meta_input.loc[idx_min]['NON_EFFECT_ALLELE']
  
  for idx, row in df_meta_input[condition].iterrows():
    if sign_effect_direction(EFFECT_ALLELE, NON_EFFECT_ALLELE, row['EFFECT_ALLELE'], row['NON_EFFECT_ALLELE']) == -1:
      df_meta_input.loc[idx, 'EFFECT_ALLELE'] = EFFECT_ALLELE
      df_meta_input.loc[idx, 'NON_EFFECT_ALLELE'] = NON_EFFECT_ALLELE
      df_meta_input.loc[idx, 'BETA'] *= -1
      
    elif sign_effect_direction(EFFECT_ALLELE, NON_EFFECT_ALLELE, row['EFFECT_ALLELE'], row['NON_EFFECT_ALLELE']) == 1:
      df_meta_input.loc[idx, 'EFFECT_ALLELE'] = EFFECT_ALLELE
      df_meta_input.loc[idx, 'NON_EFFECT_ALLELE'] = NON_EFFECT_ALLELE
      df_meta_input.loc[idx, 'BETA'] *= 1      
      
    elif sign_effect_direction(EFFECT_ALLELE, NON_EFFECT_ALLELE, row['EFFECT_ALLELE'], row['NON_EFFECT_ALLELE']) == 0:
      df_meta_input = df_meta_input.drop(idx)
  
  li_EFFECT_ALLELE.append(EFFECT_ALLELE)
  li_NON_EFFECT_ALLELE.append(NON_EFFECT_ALLELE)

# Calculation - Weighted average of the effect sizes 
# li_BETA_META : List of Weighted average of the effect sizes corresponding to li_PHENOTYPE_SNP
# li_STD_BETA_META : List of Standard Error of Weighted average of the effect sizes corresponding to li_PHENOTYPE_SNP

li_BETA_META = []
li_STD_BETA_META = []

for j in range(len(li_PHENOTYPE_SNP)):
  condition = (df_meta_input.SNP == li_PHENOTYPE_SNP[j][1]) & (df_meta_input.PHENOTYPE == li_PHENOTYPE_SNP[j][0]) 
  sum_w_i = 0
  sum_w_i_beta_i = 0
  
  if len(df_meta_input[condition]) > 1:
    for idx, row in df_meta_input[condition].iterrows():
      w_i = row['BETA_SE']**(-2)
      sum_w_i += w_i
      sum_w_i_beta_i += (w_i * row['BETA'])
    
    li_BETA_META.append(sum_w_i_beta_i/sum_w_i)
    li_STD_BETA_META.append(sum_w_i**-0.5)
  
  elif len(df_meta_input[condition]) == 1:
    li_BETA_META.append(df_meta_input[condition]['BETA'].values[0])
    li_STD_BETA_META.append(df_meta_input[condition]['BETA_SE'].values[0])

# Heterogeniety Test 
# Calculation - Cochran's Q statistic
# li_Q : List of Cochran's Q statistic corresponding to li_PHENOTYPE_SNP

li_Q = []

for j in range(len(li_PHENOTYPE_SNP)):
  condition = (df_meta_input.SNP == li_PHENOTYPE_SNP[j][1]) & (df_meta_input.PHENOTYPE == li_PHENOTYPE_SNP[j][0]) 
  Q = 0
  
  if len(df_meta_input[condition]) > 1:
    for idx, row in df_meta_input[condition].iterrows():
      w_i = row['BETA_SE']**(-2)
      delta_beta = (row['BETA'] - li_BETA_META[j])    
      Q += (w_i * (delta_beta**2))
    li_Q.append(Q)
  
  elif len(df_meta_input[condition]) == 1:
    li_Q.append('No Meta')

# Calculation - Higgin's heterogeneity metric
# li_I_square : List of Higgin's heterogeneity metric corresponding to li_PHENOTYPE_SNP

li_I_square = []

for j in range(len(li_PHENOTYPE_SNP)):
  condition = (df_meta_input.SNP == li_PHENOTYPE_SNP[j][1]) & (df_meta_input.PHENOTYPE == li_PHENOTYPE_SNP[j][0]) 
  count = 0
  
  if (len(df_meta_input[condition]) > 1) & (li_Q[j] != 0):
    for idx, row in df_meta_input[condition].iterrows():
      count += 1
    
    I_square = 100 * (li_Q[j] - (count - 1)) / (li_Q[j])
    li_I_square.append(max(I_square,0))
    
  else: 
    li_I_square.append('No Meta')

# Calculation - Weighted average of the effect sizes Modification - Random Effect Model
# li_BETA_META : List of Weighted average of the effect sizes corrected by a Random Effect Model
# li_STD_BETA_META : List of Standard Error of Weighted average of the effect sizes corrected by a Random Effect Model

for j in range(len(li_PHENOTYPE_SNP)):
  if li_I_square[j] != 'No Meta': 
    condition = (df_meta_input.SNP == li_PHENOTYPE_SNP[j][1]) & (df_meta_input.PHENOTYPE == li_PHENOTYPE_SNP[j][0]) 
    
    if (li_I_square[j] >= 50) & (len(df_meta_input[condition]) > 1):
      sum_w_i = 0
      sum_w_i_square = 0
      count = 0
      
      for idx, row in df_meta_input[condition].iterrows():
        w_i = row['BETA_SE']**(-2)
        sum_w_i += w_i
        sum_w_i_square += w_i**2
        count += 1
      
      tau_square = (li_Q[j]-count+1) / (sum_w_i - sum_w_i_square / sum_w_i)                   
      tau_square = max(tau_square, 0)
      
      sum_w_i_R_beta_i = 0
      sum_w_i_R = 0
      
      for idx2, row2 in df_meta_input[condition].iterrows():
        w_i = row2['BETA_SE']**(-2)
        w_i_R = 1/(1/w_i + tau_square)
        sum_w_i_R_beta_i += w_i_R * row2['BETA']
        sum_w_i_R += w_i_R
        
      li_BETA_META[j] = sum_w_i_R_beta_i / sum_w_i_R
      li_STD_BETA_META[j] = sum_w_i_R**-0.5

# Calculation - Integrated P-value 
# li_p_value : List of P-value corresponding to li_PHENOTYPE_SNP

li_p_value = []

for j in range(len(li_PHENOTYPE_SNP)):
  if li_I_square[j] != 'No Meta': 
    Z = li_BETA_META[j] / li_STD_BETA_META[j]
    p_value = 2 * norm.cdf(-abs(Z))
    
    li_p_value.append(p_value)
  
  else:
    condition = (df_meta_input.SNP == li_PHENOTYPE_SNP[j][1]) & (df_meta_input.PHENOTYPE == li_PHENOTYPE_SNP[j][0]) 
    li_p_value.append(df_meta_input[condition]['P_VAL'].values[0])

# Output file - Meta Analysis
# df_meta_output : Data Frame of Meta-analysis Result Ouput File

df_meta_output = pd.DataFrame(columns = ['PHENOTYPE', 'SNP', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'BETA', 'BETA_SE', 'P_VAL', 'BH_P_VAL', 'I_SQUARE', 'Q_HET'])

for j in range(len(li_PHENOTYPE_SNP)):
  df_meta_output.loc[j, 'PHENOTYPE'] = li_PHENOTYPE_SNP[j][0]
  df_meta_output.loc[j, 'SNP'] = li_PHENOTYPE_SNP[j][1]

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

  
