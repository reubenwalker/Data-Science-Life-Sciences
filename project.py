import synapseclient
import pandas as pd
import matplotlib.pyplot as plt
syn = synapseclient.Synapse()
syn.login(authToken="dummy")
# entity = syn.get("syn55234671")
entity2 = syn.get("syn55234672")
entity3 = syn.get("syn55234673")
entity4 = syn.get("syn55234796")
entity5 = syn.get("syn55234797")


df_sv = pd.read_csv(entity5.path, delimiter="\t", nrows=2)
df_mutations = pd.read_csv(entity4.path, delimiter="\t", nrows=2)

# df_CNA = pd.read_csv(entity.path, delimiter="\t")#, nrows=2)
df_patient = pd.read_csv(entity2.path, delimiter="\t", skiprows=4)

df_clinical = pd.read_csv(entity3.path, delimiter="\t", skiprows=4)

df_CNA = pd.read_csv('CNA.csv')
cols = df_CNA.columns
colsMask = ((cols == 'Hugo_Symbol') |
            # cols.str.contains('-DFCI-') |
            # cols.str.contains('-MSK-') |
            cols.str.contains('-VICC-'))
colSubset = cols[colsMask]
df_CNA = df_CNA[colSubset]

#df_CNA.to_csv('CNA.csv')
df_patient.to_csv('patient_data.csv')
df_clinical.to_csv('tumor_sample.csv')


df_CNA_T = df_CNA.T
df_CNA_T.columns = df_CNA_T.iloc[0]
df_CNA_T = df_CNA_T.drop('Hugo_Symbol')

df_clinical_merge = df_clinical[['SAMPLE_ID', 'CANCER_TYPE']].set_index('SAMPLE_ID')

df_CNA_clinical = pd.merge(df_CNA_T, df_clinical_merge, left_index=True, right_index=True)#')

df_CNA_clinical.to_csv('mergedClinicalCNA.csv')

patientDF = pd.merge(df_patient, df_clinical[['PATIENT_ID', 'SAMPLE_ID', 'AGE_AT_SEQ_REPORT', 'CANCER_TYPE']], on='PATIENT_ID')
#For Kaplan Meier, we need 
    #days from tumor sample
    #dead or censored
    #group
    
#We have:
    #Age at tumor sample
    #'Interval in days from DOB to date of last contact'
    #'Interval in days from DOB to DOD'
    #'Year of last contact'
    #'Vital Status'
    #'Year of death'
    
patientDF.groupby('AGE_AT_SEQ_REPORT')['PATIENT_ID'].count()
# Out[12]:
# AGE_AT_SEQ_REPORT
# 17            5
# 18         3146
# 19          428
# 20          387
# 21          403
           # ...
# 87          654
# 88          515
# <18        5737
# >89        1517
# Unknown    1234
# Name: PATIENT_ID, Length: 75, dtype: int64    
    
###So we have <18, >89, and Unknown as ages.
excludeMask = ((patientDF['AGE_AT_SEQ_REPORT'] != '<18') & (patientDF['AGE_AT_SEQ_REPORT'] != '>89') & (patientDF['AGE_AT_SEQ_REPORT'] != 'Unknown'))
kaplanDF = patientDF[excludeMask] 
kaplanDF['AGE_DAYS'] = kaplanDF['AGE_AT_SEQ_REPORT'].astype(int)*365
kaplanDF.groupby('INT_CONTACT')['PATIENT_ID'].count()
# Out[44]:
# INT_CONTACT
# 10000               1
# 10001              14
# 10002               1
# 10003               1
# 10004               2
                 # ...
# <6570            4041
# >32485           2006
# Not Collected    2142
# Not Released      480
# Unknown          2149
# Name: PATIENT_ID, Length: 23327, dtype: int64

###So we need to remove <6570, >32485, Not Collected, Not Released, an Unknown
excludeMask2 = pd.to_numeric(kaplanDF['INT_CONTACT'], errors='coerce').notnull()
kaplanDF = kaplanDF[excludeMask2]

kaplanDF['INT_SINCE_SAMPLE'] = kaplanDF['INT_CONTACT'].astype(int) - kaplanDF['AGE_DAYS']

kaplanDF.groupby('DEAD')['PATIENT_ID'].count()
# Out[47]:
# DEAD
# FALSE      54454
# False      47913
# TRUE       39769
# True       31665
# Unknown     1098
# Name: PATIENT_ID, dtype: int64

def binaryConvert(cell):
    if cell.lower() == 'true':
        return 1
    else:
        return 0
kaplanDF['DEAD_'] = kaplanDF['DEAD'].apply(lambda x: binaryConvert(x))

def groupCancer(cell):
    if cell == 'Cancer of Unknown Primary':
        return 'CUP'
    #if cell == 'Melanoma':
    #    return 'Melanoma'
    else:
        return 'Other'

kaplanDF['CANCER_TYPE_'] = kaplanDF['CANCER_TYPE'].apply(lambda x: groupCancer(x))

import kaplanmeier as km
# Compute Survival
results = km.fit(kaplanDF['INT_SINCE_SAMPLE'], kaplanDF['DEAD_'], kaplanDF['CANCER_TYPE_'])

# Plot
# plt.figure()
km.plot(results)
plt.savefig('kaplanmeier.png')
# plt.show()

#We have 2310 samples with negative times from last contact to sample date.
#Let's toss them for now.
kaplanDF = kaplanDF[kaplanDF['INT_SINCE_SAMPLE'] > 0]
#Let's also just look up to 1825 days (5 years)?
# Compute Survival
results = km.fit(kaplanDF['INT_SINCE_SAMPLE'], kaplanDF['DEAD_'], kaplanDF['CANCER_TYPE_'])

# Plot
# plt.figure()
km.plot(results)
plt.tight_layout()
plt.savefig('kaplanmeier2.png')


kaplanDF.groupby('CANCER_TYPE')['PATIENT_ID'].count().sort_values(ascending=False)
def groupCancer2(cell):
    if cell == 'Cancer of Unknown Primary':
        return 'CUP'
    if cell == 'Melanoma':
       return 'Melanoma'
    if cell == 'Non-Small Cell Lung Cancer':
       return 'Non-Small Cell Lung Cancer'
    if cell == 'Breast Cancer':
       return 'Breast Cancer'    
    else:
        return 'Other'

kaplanDF['CANCER_TYPE_'] = kaplanDF['CANCER_TYPE'].apply(lambda x: groupCancer2(x))
mask = kaplanDF['CANCER_TYPE_'] != 'Other'
results = km.fit(kaplanDF.loc[mask,'INT_SINCE_SAMPLE'], kaplanDF.loc[mask,'DEAD_'], kaplanDF.loc[mask,'CANCER_TYPE_'])

# Plot
# plt.figure()
km.plot(results)
# plt.tight_layout()
plt.savefig('kaplanmeier3.png')

# df_CNA = pd.read_csv(entity.path, delimiter="\t")#, nrows=2)
patient_cols = pd.read_csv(entity2.path, delimiter="\t", nrows=4)

clinical_cols = pd.read_csv(entity3.path, delimiter="\t", nrows=4)
