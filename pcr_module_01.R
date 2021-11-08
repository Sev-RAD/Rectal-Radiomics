#### Radiomics Source v1.0
# Written by Sungwon Kim, MD, PhD, Jaeseung Shin, MD, PhD, Kyunghwa Han, PhD. 
# Department of Radiology, Yonsei University College of Medicine
# 
# The radiomics-related R codes that are important for performing this study are described here.
# Radiomic features were extracted separately for the T2WIs and ADC maps using PyRadiomics 
# open-source Python package (version 2.1.2; https://pyradiomics.readthedocs.io/). 
# All radiomic features implemented in PyRadiomics were extracted in the original image and 
# filtered images, including wavelet and Laplacian of Gaussian. Before the feature extraction,
# Z-score normalization of the MRI signal intensities only for T2WIs, gray-level discretization 
# with fixed bin width values of 3 (T2WIs) and 20 (ADC maps), and voxel size resampling by 1*1*1mm
# were performed using PyRadiomics. Consequently, 1132 features were obtained for each T2WI and ADC map. 
# Feature selection and classification model building. Hierarchical feature clustering was performed
# using Spearman correlation coefficient to reduce redundancy among radiomic features. Subsequently,
# the least absolute shrinkage and selection operator (LASSO) method was used to select the most useful
# predictive features from the training set. A radiomics score (RAD score) was calculated for each 
# patient as a linear combination of the selected features weighted by their respective coefficients.
# To find an optimal regulation weight (λ) in LASSO logistic regression, ten-fold cross validation
# with minimum criteria was employed, where the final value of λ yielded minimum binomial deviance. 
# 
# For related questions, please contact the following e-mail: dinbe@yuhs.ac (Sungwon Kim)



#### MODULE 1: initialization

rm(list=ls())
setwd(WORKING_DIRECTORY) # set the working directory
SEED = 1234 
library(caret)

# The radiomic features were loaded in the variable 'doc1' according to the MRI sequences using command line console.
# Radiomic features with low ICC values were pre-removed and prepared.



#### MODULE 2: hierarchical clustering

CLUSTER_THR = 0.95
CUT_DIST = (1.0 - CLUSTER_THR)

dist = as.dist(1 - abs(cor(d, method = "spearman"))) # Variable 'd' contains patient-specific radiomic features.
fit = hclust(dist)
plot(fit, main = "Dendrogram") # draw a dendrogram.

groups = cutree(fit, h = CUT_DIST) # cut the tree with appropriate value.
rect.hclust(fit, h = CUT_DIST, border = "red") # plotting with red box.
num_group = max(groups)
# By examining the element features for each clustered group, representative features for each group were selected.



#### MODULE 3: cross-validated lasso using cv.glmnet library
# Variable 'sdoc' contains patient-specific radiomic features.

library(DMwR)
library(glmnet)
library(c060)

x = model.matrix(pCR ~., sdoc)[,-1]
y = sdoc$pCR

set.seed(SEED)
cvres = cv.glmnet(x, y, alpha = 1, family = "binomial") 
print(cvres$lambda.min)
res = cvres$glmnet.fit
plot(cvres)

cof = coef(res, s = cvres$lambda.min)
active.index = which(cof != 0)
ftcofs = cof[active.index]
ftnames = rownames(cof)[active.index]
Plot.coef.glmnet(cvfit = cvres, betas = ftnames[-1]) # -intercept.



#### MODULE 4: calculate radscore

doc$radscore = ftcofs[1] # (intercept)
for (i in 2:length(ftnames)) {
    doc$radscore = doc$radscore + (doc[, ftnames[i]] * ftcofs[i])
}



#### MODULE 5: validation set
# Variable 'bestcutoff' contains the best cutoff value derived from the training process.
# All parameters related to the model were determined during the training process.

vdoc$decision = ifelse(vdoc$radscore >= bestcutoff, "YES", "NO") 
vdoc$decision = as.factor(vdoc$decision)
ans = data.frame (
    Pred = vdoc$decision, 
    Ref = vdoc$pCR # reference standard values
)
print(confusionMatrix(ans$Pred, ans$Ref, "YES")) # calculate confusion matrix.




