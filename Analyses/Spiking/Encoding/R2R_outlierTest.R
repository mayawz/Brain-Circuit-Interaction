# outlier test

setwd("/Volumes/GoogleDrive/My Drive/01Data/00 w_Anatomy/2020final/R2R/")
library(foreign)
library(car) # for Boxplot() function
source("showmelm.R")
# shift+command+C comment / uncomment

# betas_in = [OFFERin.md(mdN).bev1(:,binN),...
#             OFFERin.md(mdN).bev2(:,binN),...
#             OFFERin.md(17).blev(:,453),...
#             OFFERin.md(17).brev(:,453)];
# betas_out = [OFFERout.md(mdN).bev1(:,binN),...
#              OFFERout.md(mdN).bev2(:,binN),...
#              OFFERout.md(17).blev(:,453),...
#              OFFERout.md(17).brev(:,453)];
# betas_pcc = [OFFERpcc.md(mdN).bev1(:,binN),...
#              OFFERpcc.md(mdN).bev2(:,binN),...
#              OFFERpcc.md(17).blev(:,453),...
             # OFFERpcc.md(17).brev(:,453)];

# LOAD DATA####
# Data = read.csv("betas_in.csv", header=FALSE, sep=",")
# Data = read.csv("betas_out.csv", header=FALSE, sep=",")
Data = read.csv("betas_pcc.csv", header=FALSE, sep=",")
head(Data)

Boxplot(Data$V1) 
Boxplot(Data$V2) 
Boxplot(Data$V3) 
Boxplot(Data$V4) 


# mymodel1 = lm(V1 ~ V2, data = Data)
mymodel1 = lm(V2 ~ V1, data = Data)
mymodel1
# create object with model diagnostics
mymodel1.diag = ls.diag(mymodel1)
names(mymodel1.diag) # available diagnostics
Data$cooks1 = mymodel1.diag$cooks
Boxplot(Data$cooks1, ylab="Cook's D value", cex = 0.8)
text(1:nrow(Data), Data$cooks1, cex = 0.3, pos = 3)
# title(main="cOFCm  - ev1 v ev2")
# title(main="cOFCl  - ev1 v ev2")
title(main="PCC  - ev1 v ev2")

mymodel2 = lm(V4 ~ V3, data = Data)
mymodel2
# create object with model diagnostics
mymodel2.diag = ls.diag(mymodel2)
Data$cooks2 = mymodel2.diag$cooks
Boxplot(Data$cooks2, ylab="Cook's D value", cex = 0.8)
text(1:nrow(Data), Data$cooks2, cex = 0.3, pos = 3)
# title(main="cOFCm  - evL v evR")
# title(main="cOFCl  - evL v evR")
title(main="PCC  - evL v evR")




