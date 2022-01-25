#!/usr/bin/env python
# coding: utf-8

# In[31]:


import pandas as pd
import seaborn as sns


# In[3]:


data = pd.read_csv("ampir_predictions.csv")


# In[27]:


data = data.sort_values(by=['AMP Probability'], ascending=False).dropna()


# In[29]:


data.to_csv("ampir_predictions_nonan.csv")


# In[38]:


data['AMP Probability'][46313]


# In[49]:


displot_all = sns.displot(data, x='AMP Probability', binwidth=0.1)
fig = displot_all.fig
fig.savefig("2022_01_25_displot_allProteins_Ampir.pdf") 


# In[46]:


#now for anything with probability over 0.5
data_50 = data[data['AMP Probability']>= 0.5]


# In[50]:


displot_all = sns.displot(data_50, x='AMP Probability', binwidth=0.1)
fig = displot_all.fig
fig.savefig("2022_01_25_displot_over0.5_Ampir.pdf")


# In[ ]:




