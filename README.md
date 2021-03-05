**Some data and codes for the study of Wang et al., 2021, Individual Uniqueness in the Neonatal Functional Connectome**

Data: 
1.	Basic information of the neonates, including participant ID, sex, age, and head motion parameter.
Basic_info.mat
Note: In the third column, ‘-1’ and ‘1’ denote female and male, respectively. The fourth and fifth columns denote the gestational age at birth and postmenstrual age at scan, respectively. The sixth and seventh columns denote mFD of Section 1 and Section 2, respectively. mFD, mean framewise displacement.
2.	Individual identification based on preprocessed fMRI data with GSR
	   indiv_GSR.mat
3.	Individual identification based on preprocessed fMRI data without GSR
	   indiv_NGSR.mat
4.	Individual identification based on scrubbed fMRI data
	   indiv_scrubbed.mat
5.	Individual identification based on only strong connections
	   indiv_strong.mat 

Codes:
1.	Individual identification (Finn et al. 2015)
	   ID_predict.m
2.	Calculate group consistency and differential power values of edges (Finn et al. 2015)
	   Cal_GC_DP.m
3.	Generate a mask denoting strong connections in common for both sections
	   strong_mask.m
4.	Examine age, sex, and their interaction effects
	   glm.m
5.	Calculate the inter-subject variability at the connection and ROI levels
	   variability.m

References

Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang J, Chun MM, Papademetris X, Constable RT. 2015. Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity. Nat Neurosci. 18:1664-1671.

