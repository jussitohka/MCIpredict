# MCIpredict
Code for the paper 
"Predicting the future development of mild cognitive impairment in the cognitively healthy elderly" 
by  
Bryan A. Strange, Linda Zhang, Christopher J. Long, Alba Sierra-Marcos, Eva Alfayate, Jussi Tohka, Miguel Medina  

Utilizes Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N.
http://www.stanford.edu/~hastie/glmnet_matlab/ package which can be obtained from above mentioned URL.

Some functions of the GLMNET package have been modified to allow for multiple runs of CV and AUC computation by pooling. 
The modified code is provided in this repository. The modified code is not designed to be used for any other purpose than this. 

prepare_vallecas_data_visitmatched_ds.m is the main data preparation script
vallecas_glmnet_visitmatchedconverters_ds.m is the main data analysis and result collection script

24 Oct 2020: A bug correction in knn_imputation file. This caused the method to use mean imputation instead of the knn imputation.   
