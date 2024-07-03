import synapseclient
import pandas as pd
import matplotlib.pyplot as plt
syn = synapseclient.Synapse()
syn.login(authToken="dummy")
entity = syn.get("syn55234671")
entity2 = syn.get("syn55234672")
entity3 = syn.get("syn55234673")
entity4 = syn.get("syn55234796")
entity5 = syn.get("syn55234797")

#
df_CNA = pd.read_csv(entity.path, delimiter="\t")#, nrows=2)


# df_CNA = pd.read_csv('CNA.csv')
cols = df_CNA.columns
colsMask = ((cols == 'Hugo_Symbol') |
            cols.str.contains('-DFCI-') |
            cols.str.contains('-MSK-') |
            cols.str.contains('-VICC-'))
colSubset = cols[colsMask]
df_CNA = df_CNA[colSubset]



#df_CNA.to_csv('CNA.csv')
# df_patient.to_csv('patient_data.csv')
# df_clinical.to_csv('tumor_sample.csv')


df_CNA_T = df_CNA.T
df_CNA_T.columns = df_CNA_T.iloc[0]
df_CNA_T = df_CNA_T.drop('Hugo_Symbol')
del df_CNA

df_patient = pd.read_csv(entity2.path, delimiter="\t", skiprows=4)

df_clinical = pd.read_csv(entity3.path, delimiter="\t", skiprows=4)

#Let's get SampleID, Sex, Age, Cancer Type
df_patient_clinical = pd.merge(df_patient[['PATIENT_ID', 'SEX']], df_clinical[['PATIENT_ID', 'SAMPLE_ID', 'AGE_AT_SEQ_REPORT', 'CANCER_TYPE']])

df_clinical_merge = df_patient_clinical[['SAMPLE_ID', 'SEX', 'AGE_AT_SEQ_REPORT', 'CANCER_TYPE']].set_index('SAMPLE_ID')

df_CNA_T = df_CNA_T.add_suffix('_CNA')

df_CNA_clinical = pd.merge(df_CNA_T, df_clinical_merge, left_index=True, right_index=True)
df_CNA_clinical = df_CNA_clinical.dropna(axis=1) #Drop empty columns

del df_CNA_T
del df_patient
del df_clinical
del df_patient_clinical
del df_clinical_merge
# df_CNA_clinical.to_csv('mergedClinicalCNA.csv')


df_sv = pd.read_csv(entity5.path, delimiter="\t")#, nrows=2)
#df_mutations = pd.read_csv(entity4.path, delimiter="\t", nrows=1000)
instList = ['DFCI', 'MSK', 'VICC']
mask = df_sv['Center'].isin(instList)
df_sv = df_sv[mask]



#From the structural variant data, we want
        #For each gene, the
        #total count of a somatic mutation (that is, SNV and indels) was encoded
        #as a positive integer feature
    #Subset: SV_Status == "Somatic"
    #Rows: ID
    #Column: Class
    #Cell: Category (Insertion/Deletion/Duplication/Inversion)

somatic_mask = df_sv['SV_Status'] == 'SOMATIC' #Just somatic variants
df_sv = df_sv.loc[somatic_mask, ]
df_svCounts = df_sv.groupby(['Sample_Id','Site1_Hugo_Symbol'])['Class'].count().unstack() 
df_svCounts = df_svCounts.fillna(0) #Fill missing counts with 0

df_svCounts = df_svCounts.add_suffix('_SV')

df_CNA_SV = pd.merge(df_svCounts, df_CNA_clinical,left_index=True, right_index=True) #on=['INSTITUTION','SAMPLE'])

del df_CNA_clinical
del df_sv
del df_svCounts

#Let's also remove age <18, >89
import pandas as pd
# df = pd.read_csv('CNA_SV.csv').set_index('Unnamed: 0')
# df_CNA_SV = df
ageMask = ((df_CNA_SV['AGE_AT_SEQ_REPORT'] == '<18') | (df_CNA_SV['AGE_AT_SEQ_REPORT'] == '>89') | (df_CNA_SV['AGE_AT_SEQ_REPORT'] == 'Unknown'))
df_CNA_SV = df_CNA_SV[~ageMask]
sexMask = df_CNA_SV['SEX'] != 'Unknown'
df_CNA_SV = df_CNA_SV[sexMask]
df_CNA_SV.to_csv('CNA_SV.csv')


#From the mutation data, we want
    #Rows: ID
    #Column: Mutation gene
    #Cell: Polyphen score
#somatic_mask2 = df_mutations['SV_Status'] == 'SOMATIC'  #Do we need this?
# df_mutations