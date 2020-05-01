#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd
import numpy as np


# In[2]:


pd.__version__


# In[3]:


from platform import python_version
print(python_version())


# In[24]:


df_DS8_low = pd.read_pickle('PC9-DS8_param-scan_twoState_lowDivsMoreDIP.pkl')
df_DS8_lowMed = pd.read_pickle('PC9-DS8_param-scan_twoState_lowMedDivsMoreDIP.pkl')
df_DS8_medHigh = pd.read_pickle('PC9-DS8_param-scan_twoState_medHighDivsMoreDIP.pkl')
df_DS8_high = pd.read_pickle('PC9-DS8_param-scan_twoState_highDivsMoreDIP.pkl')

# df_DS8 = pd.concat([df_DS8_low, df_DS8_lowMed, df_DS8_medHigh, df_DS8_high])

# df_DS8 = df_DS8.round({'death rate': 5, 'division rate': 5, 'DIP rate': 5})
df_DS8_low['cell.line'] = 'PC9.DS8'
df_DS8_lowMed['cell.line'] = 'PC9.DS8'
df_DS8_medHigh['cell.line'] = 'PC9.DS8'
df_DS8_high['cell.line'] = 'PC9.DS8'

print(df_DS8_low)
print(type(df_DS8_low))



# In[29]:


DS8_dipdf_low = pd.DataFrame(df_DS8_low['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
DS8_divdf_low = pd.DataFrame(df_DS8_low['division rate'].values.tolist(), columns=['div1','div2'])
DS8_dthdf_low = pd.DataFrame(df_DS8_low['death rate'].values.tolist(), columns=['dth1','dth2'])
DS8_rest_low = df_DS8_low[['p-value', 'cell.line']]
DS8_rest_low.reset_index(drop=True)

DS8_dipdf_lowMed = pd.DataFrame(df_DS8_lowMed['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
DS8_divdf_lowMed = pd.DataFrame(df_DS8_lowMed['division rate'].values.tolist(), columns=['div1','div2'])
DS8_dthdf_lowMed = pd.DataFrame(df_DS8_lowMed['death rate'].values.tolist(), columns=['dth1','dth2'])
DS8_rest_lowMed = df_DS8_lowMed[['p-value', 'cell.line']]
DS8_rest_lowMed.reset_index(drop=True)

DS8_dipdf_medHigh = pd.DataFrame(df_DS8_medHigh['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
DS8_divdf_medHigh = pd.DataFrame(df_DS8_medHigh['division rate'].values.tolist(), columns=['div1','div2'])
DS8_dthdf_medHigh = pd.DataFrame(df_DS8_medHigh['death rate'].values.tolist(), columns=['dth1','dth2'])
DS8_rest_medHigh = df_DS8_medHigh[['p-value', 'cell.line']]
DS8_rest_medHigh.reset_index(drop=True)

DS8_dipdf_high = pd.DataFrame(df_DS8_high['DIP rate'].values.tolist(), columns=['DIP1','DIP2'])
DS8_divdf_high = pd.DataFrame(df_DS8_high['division rate'].values.tolist(), columns=['div1','div2'])
DS8_dthdf_high = pd.DataFrame(df_DS8_high['death rate'].values.tolist(), columns=['dth1','dth2'])
DS8_rest_high = df_DS8_high[['p-value', 'cell.line']]
DS8_rest_high.reset_index(drop=True)

print(DS8_dipdf_low)
print(DS8_divdf_low)
print(DS8_dthdf_low)
print(DS8_rest_low)

DS8_dipdf_low.to_csv("DS8_twoState_lowDiv_dips.csv")
DS8_divdf_low.to_csv("DS8_twoState_lowDiv_divs.csv")
DS8_dthdf_low.to_csv("DS8_twoState_lowDiv_dths.csv")
DS8_rest_low.to_csv("DS8_twoState_lowDiv_rest.csv")


# In[48]:


DS8_low = pd.concat([DS8_rest_low.reset_index(drop=True), 
                     DS8_dipdf_low.reset_index(drop=True), 
                     DS8_divdf_low.reset_index(drop=True),
                     DS8_dthdf_low], axis = 1)

DS8_lowMed = pd.concat([DS8_rest_lowMed.reset_index(drop=True), 
                     DS8_dipdf_lowMed.reset_index(drop=True), 
                     DS8_divdf_lowMed.reset_index(drop=True),
                     DS8_dthdf_lowMed], axis = 1)

DS8_medHigh = pd.concat([DS8_rest_medHigh.reset_index(drop=True), 
                     DS8_dipdf_medHigh.reset_index(drop=True), 
                     DS8_divdf_medHigh.reset_index(drop=True),
                     DS8_dthdf_medHigh], axis = 1)

DS8_high = pd.concat([DS8_rest_high.reset_index(drop=True), 
                     DS8_dipdf_high.reset_index(drop=True), 
                     DS8_divdf_high.reset_index(drop=True),
                     DS8_dthdf_high], axis = 1)
print(DS8_low)
print(DS8_lowMed)
print(DS8_medHigh)
print(DS8_high)

# DS8_low.to_csv("DS8_twoState_lowDivs.csv")
# DS8_lowMed.to_csv("DS8_twoState_lowMedDivs.csv")
# DS8_medHigh.to_csv("DS8_twoState_medHighDivs.csv")
# DS8_high.to_csv("DS8_twoState_highDivs.csv")


# In[50]:


DS8_all = pd.concat([DS8_low, DS8_lowMed,
                      DS8_medHigh, DS8_high], ignore_index = True)

print(DS8_all)
DS8_all.to_csv('DS8_all.csv')

DS8 = DS8_all


# In[51]:


DS8_wholeRange = pd.read_csv("DS8_wholeRange.csv")
np.shape(DS8_wholeRange)


# In[52]:


DS8['cell.line'] = np.where(DS8['p-value']>0.1, 'PC9-DS8', 'not.assigned')
DS8['param.pair'] = range(DS8.shape[0])
print(DS8)


# In[53]:


DS8_sig1 = DS8[['p-value', 'cell.line', 'DIP1', 'div1', 'dth1', 'param.pair']]
DS8_sig2 = DS8[['p-value', 'cell.line', 'DIP2', 'div2', 'dth2', 'param.pair']]

DS8_sig1.rename(columns={'DIP1': 'DIP Rate',
                         'div1': 'Division Rate',
                         'dth1': 'Death Rate'},
                 inplace=True)
DS8_sig2.rename(columns={'DIP2': 'DIP Rate',
                         'div2': 'Division Rate',
                         'dth2': 'Death Rate'},
                inplace=True)

print(DS8_sig1)
print(DS8_sig2)


# In[54]:


DS8_sig1['Cell Line'] = np.where(DS8_sig1['cell.line'] == "PC9-DS8", 'PC9-DS8.1', 'not.assigned')
DS8_sig2['Cell Line'] = np.where(DS8_sig2['cell.line'] == "PC9-DS8", 'PC9-DS8.2', 'not.assigned')

print(DS8_sig1)
print(DS8_sig2)



# In[55]:


DS8_sig_all = pd.concat([DS8_sig1, DS8_sig2])
print(DS8_sig_all)


# In[56]:


DS8_sig_all.to_csv('DS8_twoState_tile_wholeRange.csv')


# In[ ]:





# In[ ]:




