import os, datetime
import pandas as pd
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import zepid
import sys
from zepid.graphics import EffectMeasurePlot
from scipy.stats import norm

#-------------------------------------------------------
# Common Function
#-------------------------------------------------------
def WriteLog(functionname, msg, type='INFO', fplog=None):
    #strmsg = "[%s][%s][%s] %s\n" % (datetime.datetime.now(), type, functionname, msg)
    #if( DEBUGMODE ): 
    #    print(strmsg)
    #else:
    #    if( fplog != None ):
    #        fplog.write(strmsg)
    
    head = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    writestr = f"[{head}][{functionname}] {msg}\n"
    #if( DEBUGMODE ):
    if( True ):
        #print(message)
        writestr = f"[{functionname}] {msg}\n"
        print(writestr)
        
    if( fplog != None ):
        fplog.write(writestr)
        fplog.flush()

def CalculateSignEffectDirection(effect_allele_study1, non_effect_allele_study1, effect_allele_study2, non_effect_allele_study2):
    """
    Calculate the direction sign of effect between meta-studies
    
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
# Test the Function - "CalculateSignEffectDirection"

print(CalculateSignEffectDirection('T', 'C', 'T', 'A'), 'result:0')   #If only one EA and NEA of study1 and study2 are the same and the rest are different, the two studies are not combined.

print(CalculateSignEffectDirection('A', 'G', 'A', 'G'), 'result:1')   #If the EA and NEA of study1 and study2 are the same, the beta sign of each study does not change. 
print(CalculateSignEffectDirection('A', 'G', 'T', 'C'), 'result:1')   #If EA and NEA of study1 and study2 are identical in the complementary sequence relationship, the beta sign of each study does not change.
print(CalculateSignEffectDirection('A', 'G', 'C', 'T'), 'result:-1')  #If the EA and NEA in study 1 and study 2 are opposite in the complementary sequence relationship, the beta sign changes.

print(CalculateSignEffectDirection('A', 'T', 'G', 'C'), 'result:0')   #If the EA and NEA of study1 and study2 are totally different, the two studies are not combined.
print(CalculateSignEffectDirection('A', 'T', 'A', 'T'), 'result:1')   #If the EA and NEA of study1 and study2 are the same, the beta sign of each study does not change. 
print(CalculateSignEffectDirection('A', 'T', 'T', 'A'), 'result:-1')  #If the EA and NEA in study 1 and study 2 are opposite, the beta sign changes.

print(CalculateSignEffectDirection('A', 'C', 'T', 'G'), 'result:1')   #If EA and NEA of study1 and study2 are identical in the complementary sequence relationship, the beta sign of each study does not change.
print(CalculateSignEffectDirection('A', 'C', 'G', 'T'), 'result:-1')  #If the EA and NEA in study 1 and study 2 are opposite in the complementary sequence relationship, the beta sign changes.

print(CalculateSignEffectDirection('G', 'C', 'G', 'C'), 'result:1')   #If the EA and NEA of study1 and study2 are the same, the beta sign of each study does not change.
print(CalculateSignEffectDirection('G', 'C', 'C', 'G'), 'result:-1')  #If the EA and NEA in study 1 and study 2 are opposite, the beta sign changes.

print(CalculateSignEffectDirection('-', 'A', 'G', 'A'), 'result:0')   #If the value of EA or NEA is not one of "A, G, T, C", the two studies are not combined.
"""

###################################
# MainClass
###################################

class BetaMeta:
    def __init__(self, fplog=None):
        """
        Initializes a BetaMeta object.
        """                
        
        self.__fplog=fplog
        
        curdir = os.path.abspath('')
        self.path_meta_data_dir = f"{curdir}/input/"   
        self.path_meta_output = f"{curdir}/output/meta_output.xlsx" 
        self.path_meta_forestplot_output = f"{curdir}/output/meta_forestplot.png" 
        
        self.df_meta_input = None
        self.df_meta_output = None
        
        self.file_list = None
        self.li_column = None
        self.li_PHENOTYPE_SNP = None
        self.li_EFFECT_ALLELE = None
        self.li_NON_EFFECT_ALLELE = None
        self.li_BETA_META = None
        self.li_STD_BETA_META = None
        self.li_Q = None
        self.li_I_SQUARE = None
        self.li_P_VALUE = None
        
    def ConcatData(self): 
        """
        Column Extraction & Concat the Dataframes 

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # file_list : List of File name in the Input Folder   
            self.file_list = os.listdir(self.path_meta_data_dir)           
            
            # li_column : List of column name in the Input File
            self.li_column = ['PHENOTYPE', 'SNP', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'BETA', 'BETA_SE', 'OR', 'OR_95%CI_LOWER', 'OR_95%CI_UPPER', 'P_VAL']
            # df_meta_input : Data frame of Input Files to be Meta-Analyzed
            self.df_meta_input = pd.DataFrame(columns = self.li_column)
            
            for file in self.file_list:
                path_meta_data = self.path_meta_data_dir + file
                df_meta_input_data = pd.read_excel(path_meta_data)

                if set(self.li_column).issubset(set(list(df_meta_input_data.columns))):
                    df_meta_input_data = df_meta_input_data.loc[:, self.li_column]
                    df_meta_input_data = df_meta_input_data.replace({'EFFECT_ALLELE':{'a':'A', 'c':'C', 'g':'G', 't':'T'},'NON_EFFECT_ALLELE':{'a':'A', 'c':'C', 'g':'G', 't':'T'}})

                    self.df_meta_input = pd.concat([self.df_meta_input, df_meta_input_data])

                else:
                    print('Please check the columns of the input file! (Ex, PHENOTYPE, SNP, EFFECT_ALLELE, NON_EFFECT_ALLELE, BETA, BETA_SE, OR, OR_95%CI_LOWER, OR_95%CI_UPPER, P_VAL)')  
                    sys.exit()  
        
        
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg        
        
    def DeduplicateData(self): 
        """
        Remove Spaces & Deduplication - df_meta_input

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.df_meta_input['PHENOTYPE'] = self.df_meta_input['PHENOTYPE'].str.strip()
            self.df_meta_input['SNP'] = self.df_meta_input['SNP'].str.strip()
            self.df_meta_input['EFFECT_ALLELE'] = self.df_meta_input['EFFECT_ALLELE'].str.strip()
            self.df_meta_input['NON_EFFECT_ALLELE'] = self.df_meta_input['NON_EFFECT_ALLELE'].str.strip()

            self.df_meta_input = self.df_meta_input.drop_duplicates(self.li_column)

            if not (set(list(self.df_meta_input['EFFECT_ALLELE'])).issubset(set(['A', 'G', 'T', 'C'])) & set(list(self.df_meta_input['NON_EFFECT_ALLELE'])).issubset(set(['A', 'G', 'T', 'C']))):
                print('Check the values of EFFECT_ALLELE and NON_EFFECT_ALLELE. The values should be in the form of (A, G, T, C)!')
        
        
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg           

    def CalculateBeta(self): 
        """
        Calculate BETA & BETA_SE from Odds ratio & 95% CI

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # li_BETA : List of Beta corresponding to df_meta_input
            # li_BETA_SE : List of Beta Standard Error corresponding to df_meta_input

            li_BETA = []
            li_BETA_SE = []

            for idx, row in self.df_meta_input.iterrows():
                try:
                    if (np.isnan(row['OR']) == False) & (np.isnan(row['OR_95%CI_LOWER']) == False) & (np.isnan(row['OR_95%CI_UPPER']) == False):
                        Beta = np.log(row['OR'])
                        Beta_SE = (np.log(row['OR_95%CI_UPPER']) - np.log(row['OR_95%CI_LOWER']))/3.92 

                        li_BETA.append(Beta)
                        li_BETA_SE.append(Beta_SE)

                    else:
                        li_BETA.append(row['BETA'])
                        li_BETA_SE.append(row['BETA_SE'])

                except:
                    print("Either (BETA, BETA_SE) or (OR, OR_95%CI_LOWER, OR_95%CI_UPPER) values must also be written as numbers!") 
                    sys.exit()

            self.df_meta_input['BETA'] = li_BETA
            self.df_meta_input['BETA_SE'] = li_BETA_SE
        
            # Remove the NaN Data & Unnecessary Columns - df_meta_input
            self.df_meta_input = self.df_meta_input.drop(['OR', 'OR_95%CI_LOWER', 'OR_95%CI_UPPER'], axis=1)
            self.df_meta_input = self.df_meta_input.dropna()      
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg   

    def CorrectEffectDirection(self): 
        """
        Correct the Effect Direction 

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # Create the list [PHENOTYPE, SNP]
            # li_PHENOTYPE_SNP : List of Phenotype & SNP 

            self.li_PHENOTYPE_SNP = []

            for idx, row in self.df_meta_input.iterrows(): 
                if [row['PHENOTYPE'], row['SNP']] not in self.li_PHENOTYPE_SNP:
                    self.li_PHENOTYPE_SNP.append([row['PHENOTYPE'], row['SNP']])

            # Effect Direction Correction
            # li_EFFECT_ALLELE : List of Effect Allele corresponding to li_PHENOTYPE_SNP
            # li_NON_EFFECT_ALLELE : List of Non-Effect Allele corresponding to li_PHENOTYPE_SNP

            self.df_meta_input['P_VAL'] = self.df_meta_input['P_VAL'].astype('float64')

            self.li_EFFECT_ALLELE = []
            self.li_NON_EFFECT_ALLELE = []

            for j in range(len(self.li_PHENOTYPE_SNP)):
                condition = (self.df_meta_input.SNP == self.li_PHENOTYPE_SNP[j][1]) & (self.df_meta_input.PHENOTYPE == self.li_PHENOTYPE_SNP[j][0]) 

                idx_min = self.df_meta_input[condition]['P_VAL'].idxmin()                 # For EA and NEA, the smallest p-value is the standard.
                EFFECT_ALLELE = self.df_meta_input.loc[idx_min]['EFFECT_ALLELE']    
                NON_EFFECT_ALLELE = self.df_meta_input.loc[idx_min]['NON_EFFECT_ALLELE']

                for idx, row in self.df_meta_input[condition].iterrows():
                    if CalculateSignEffectDirection(EFFECT_ALLELE, NON_EFFECT_ALLELE, row['EFFECT_ALLELE'], row['NON_EFFECT_ALLELE']) == -1:
                        self.df_meta_input.loc[idx, 'EFFECT_ALLELE'] = EFFECT_ALLELE
                        self.df_meta_input.loc[idx, 'NON_EFFECT_ALLELE'] = NON_EFFECT_ALLELE
                        self.df_meta_input.loc[idx, 'BETA'] *= -1

                    elif CalculateSignEffectDirection(EFFECT_ALLELE, NON_EFFECT_ALLELE, row['EFFECT_ALLELE'], row['NON_EFFECT_ALLELE']) == 1:
                        self.df_meta_input.loc[idx, 'EFFECT_ALLELE'] = EFFECT_ALLELE
                        self.df_meta_input.loc[idx, 'NON_EFFECT_ALLELE'] = NON_EFFECT_ALLELE
                        self.df_meta_input.loc[idx, 'BETA'] *= 1      

                    elif CalculateSignEffectDirection(EFFECT_ALLELE, NON_EFFECT_ALLELE, row['EFFECT_ALLELE'], row['NON_EFFECT_ALLELE']) == 0:
                        self.df_meta_input = self.df_meta_input.drop(idx)

                self.li_EFFECT_ALLELE.append(EFFECT_ALLELE)
                self.li_NON_EFFECT_ALLELE.append(NON_EFFECT_ALLELE)    
    
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg  

    def CalculateWeightedAverageEffectSize(self): 
        """
        Calculate the Weighted average of the effect sizes 

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # Calculation - Weighted average of the effect sizes 
            # li_BETA_META : List of Weighted average of the effect sizes corresponding to li_PHENOTYPE_SNP
            # li_STD_BETA_META : List of Standard Error of Weighted average of the effect sizes corresponding to li_PHENOTYPE_SNP

            self.li_BETA_META = []
            self.li_STD_BETA_META = []

            for j in range(len(self.li_PHENOTYPE_SNP)):
                condition = (self.df_meta_input.SNP == self.li_PHENOTYPE_SNP[j][1]) & (self.df_meta_input.PHENOTYPE == self.li_PHENOTYPE_SNP[j][0]) 
                sum_w_i = 0
                sum_w_i_beta_i = 0

                if len(self.df_meta_input[condition]) > 1:
                    for idx, row in self.df_meta_input[condition].iterrows():
                        w_i = row['BETA_SE']**(-2)
                        sum_w_i += w_i
                        sum_w_i_beta_i += (w_i * row['BETA'])

                    self.li_BETA_META.append(sum_w_i_beta_i/sum_w_i)
                    self.li_STD_BETA_META.append(sum_w_i**-0.5)

                elif len(self.df_meta_input[condition]) == 1:
                    self.li_BETA_META.append(self.df_meta_input[condition]['BETA'].values[0])
                    self.li_STD_BETA_META.append(self.df_meta_input[condition]['BETA_SE'].values[0])
    
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg  

    def CalculateQstatistic(self): 
        """
        Heterogeniety Test - Calculate the Cochran's Q statistic

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # Heterogeniety Test 
            # Calculation - Cochran's Q statistic
            # li_Q : List of Cochran's Q statistic corresponding to li_PHENOTYPE_SNP

            self.li_Q = []

            for j in range(len(self.li_PHENOTYPE_SNP)):
                condition = (self.df_meta_input.SNP == self.li_PHENOTYPE_SNP[j][1]) & (self.df_meta_input.PHENOTYPE == self.li_PHENOTYPE_SNP[j][0]) 
                Q = 0

                if len(self.df_meta_input[condition]) > 1:
                    for idx, row in self.df_meta_input[condition].iterrows():
                        w_i = row['BETA_SE']**(-2)
                        delta_beta = (row['BETA'] - self.li_BETA_META[j])    
                        Q += (w_i * (delta_beta**2))
                    self.li_Q.append(Q)

                elif len(self.df_meta_input[condition]) <= 1:
                    self.li_Q.append('Unprocessed')
    
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg  

    def CalculateHeterogeneityMetric(self): 
        """
        Heterogeniety Test - Calculate the Higgin's heterogeneity metric

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # Calculation - Higgin's heterogeneity metric
            # li_I_SQUARE : List of Higgin's heterogeneity metric corresponding to li_PHENOTYPE_SNP

            self.li_I_SQUARE = []

            for j in range(len(self.li_PHENOTYPE_SNP)):
                condition = (self.df_meta_input.SNP == self.li_PHENOTYPE_SNP[j][1]) & (self.df_meta_input.PHENOTYPE == self.li_PHENOTYPE_SNP[j][0]) 
                count = 0

                if (len(self.df_meta_input[condition]) > 1) & (self.li_Q[j] != 0):
                    for idx, row in self.df_meta_input[condition].iterrows():
                        count += 1

                    I_square = 100 * (self.li_Q[j] - (count - 1)) / (self.li_Q[j])
                    self.li_I_SQUARE.append(max(I_square,0))

                else:
                    self.li_I_SQUARE.append('Unprocessed')
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg  

    def CalculateRandomEffectModel(self): 
        """
        Random Effect Model - Calculate the Weighted average of the effect sizes 

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # Calculation - Weighted average of the effect sizes Modification - Random Effect Model
            # li_BETA_META : List of Weighted average of the effect sizes corrected by a Random Effect Model
            # li_STD_BETA_META : List of Standard Error of Weighted average of the effect sizes corrected by a Random Effect Model

            for j in range(len(self.li_PHENOTYPE_SNP)):
                if self.li_I_SQUARE[j] != 'Unprocessed': 
                    condition = (self.df_meta_input.SNP == self.li_PHENOTYPE_SNP[j][1]) & (self.df_meta_input.PHENOTYPE == self.li_PHENOTYPE_SNP[j][0]) 

                if (self.li_I_SQUARE[j] >= 50):
                    sum_w_i = 0
                    sum_w_i_square = 0
                    count = 0

                    for idx, row in self.df_meta_input[condition].iterrows():
                        w_i = row['BETA_SE']**(-2)
                        sum_w_i += w_i
                        sum_w_i_square += w_i**2
                        count += 1

                    tau_square = (self.li_Q[j]-count+1) / (sum_w_i - sum_w_i_square / sum_w_i)                   
                    tau_square = max(tau_square, 0)

                    sum_w_i_R_beta_i = 0
                    sum_w_i_R = 0

                    for idx2, row2 in self.df_meta_input[condition].iterrows():
                        w_i = row2['BETA_SE']**(-2)
                        w_i_R = 1/(1/w_i + tau_square)
                        sum_w_i_R_beta_i += w_i_R * row2['BETA']
                        sum_w_i_R += w_i_R

                    self.li_BETA_META[j] = sum_w_i_R_beta_i / sum_w_i_R
                    self.li_STD_BETA_META[j] = sum_w_i_R**-0.5

            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg  

    def CalculateIntegratedPvalue(self): 
        """
        Calculate the Integrated P-value 

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # Calculation - Integrated P-value 
            # li_P_VALUE : List of P-value corresponding to li_PHENOTYPE_SNP

            self.li_P_VALUE = []

            for j in range(len(self.li_PHENOTYPE_SNP)):
                if self.li_I_SQUARE[j] != 'Unprocessed': 
                    Z = self.li_BETA_META[j] / self.li_STD_BETA_META[j]
                    p_value = 2 * norm.cdf(-abs(Z))

                    self.li_P_VALUE.append(p_value)

                else:
                    condition = (self.df_meta_input.SNP == self.li_PHENOTYPE_SNP[j][1]) & (self.df_meta_input.PHENOTYPE == self.li_PHENOTYPE_SNP[j][0]) 
                    self.li_P_VALUE.append(self.df_meta_input[condition]['P_VAL'].values[0])

            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg  

    def CreateOutputDataFrame(self): 
        """
        Create an Output DataFrame - df_meta_output

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # Output file - Meta Analysis
            # df_meta_output : Data Frame of Meta-analysis Result Ouput File

            self.df_meta_output = pd.DataFrame(columns = ['PHENOTYPE', 'SNP', 'EFFECT_ALLELE', 'NON_EFFECT_ALLELE', 'BETA', 'BETA_SE', 'P_VAL', 'BH_P_VAL', 'I_SQUARE', 'Q_HET'])

            for j in range(len(self.li_PHENOTYPE_SNP)):
                self.df_meta_output.loc[j, 'PHENOTYPE'] = self.li_PHENOTYPE_SNP[j][0]
                self.df_meta_output.loc[j, 'SNP'] = self.li_PHENOTYPE_SNP[j][1]

            self.df_meta_output['BETA'] = self.li_BETA_META
            self.df_meta_output['BETA_SE'] = self.li_STD_BETA_META
            self.df_meta_output['P_VAL'] = self.li_P_VALUE
            self.df_meta_output['EFFECT_ALLELE'] = self.li_EFFECT_ALLELE
            self.df_meta_output['NON_EFFECT_ALLELE'] = self.li_NON_EFFECT_ALLELE
            self.df_meta_output['I_SQUARE'] = self.li_I_SQUARE
            self.df_meta_output['Q_HET'] = self.li_Q  

            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg 

    def CorrectPvalue(self): 
        """
        BH adjustment - Correct the P-value

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # P-value Correction - BH adjustment
            # df_meta_output : Data Frame of Meta-analysis Result Ouput File

            for idx, row in self.df_meta_output.iterrows():
                condition = self.df_meta_output.PHENOTYPE == row['PHENOTYPE']

                ascending_rank_p_val = self.df_meta_output.groupby('PHENOTYPE')['P_VAL'].rank(method='min', ascending=True).values[idx]
                self.df_meta_output.loc[idx, 'BH_P_VAL'] = self.df_meta_output.loc[idx, 'P_VAL'] * len(self.df_meta_output[condition]) / ascending_rank_p_val

            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg 

    def SaveOutputFiles(self): 
        """
        Save the Output Excel File & Forest Plot

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            # Save the Output Excel File             
            self.df_meta_output.to_excel(self.path_meta_output)

            # Save the Output Forest Plot
            labs = []
            measure = []
            lower = []
            upper = []
            for j in range(len(self.li_PHENOTYPE_SNP)):
                if self.li_I_SQUARE[j] != 'Unprocessed':
                    text = self.li_PHENOTYPE_SNP[j][0] + ' - ' + self.li_PHENOTYPE_SNP[j][1]
                    beta_lower = self.li_BETA_META[j] - 1.96 * self.li_STD_BETA_META[j]
                    beta_upper = self.li_BETA_META[j] + 1.96 * self.li_STD_BETA_META[j]

                    labs.append(text)
                    measure.append(self.li_BETA_META[j])
                    lower.append(beta_lower)
                    upper.append(beta_upper)

            p = EffectMeasurePlot(label=labs, effect_measure=measure, lcl=lower, ucl=upper)
            p.labels(effectmeasure='Beta')
            p.colors(pointshape="D")
            ax=p.plot(figsize=(10,5), t_adjuster=0.02)
            plt.suptitle("Beta-Meta Forest Plot",x=0.35,y=0.98)
            ax.set_xticks([0])
            ax.set_xticklabels(['$0$'])
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False) 
            ax.spines['bottom'].set_visible(True)
            ax.spines['right'].set_position(('data', 0)) 
            plt.savefig(self.path_meta_forestplot_output,bbox_inches='tight')

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            sys.exit()
    
        return rv, rvmsg     

####################################
# main
####################################
if __name__ == '__main__':
    
    betameta = BetaMeta()
    betameta.ConcatData()
    betameta.DeduplicateData()   
    betameta.CalculateBeta()    
    betameta.CorrectEffectDirection() 
    betameta.CalculateWeightedAverageEffectSize()    
    betameta.CalculateQstatistic()
    betameta.CalculateHeterogeneityMetric()
    betameta.CalculateRandomEffectModel()  
    betameta.CalculateIntegratedPvalue()  
    betameta.CreateOutputDataFrame()  
    betameta.CorrectPvalue()  
    betameta.SaveOutputFiles()  
    
    print("Meta-Analysis Complete")

    
    