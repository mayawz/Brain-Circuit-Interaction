# analysis for OFC-PCC paper @ NCOMM
############################## some examples ############################## 
prop.test(x, n, p = NULL,
          alternative = c("two.sided", "less", "greater"),
          conf.level = 0.95, correct = TRUE)

# prob of 5 or less heads of fair coin out of 10 flips
# rbinom(5, 10, .5)
rbinom(n, size, prob)
# 1-sample test for equality of proportions with continuity correction
prop.test(x = c(9, 3), n = c(12, 12), correct = TRUE, alternative ="greater")

# comparing 47% and 30%: one big N --> p-value = 0.05522/ NO CORRECTION p-value = 9.382e-05
prop.test(x = c(75, 117.5), n = c(250, 250), correct = FALSE, alternative ="two.sided")

# comparing 47% and 30%: one small N --> p-value = 0.0001367/ NO CORRECTION p-value = 0.01605
prop.test(x = c(28.5, 44.65), n = c(95, 95), correct = FALSE, alternative ="two.sided")

# comparing 47% and 30%: 2 different N --> p-value = 3.089e-08
prop.test(x = c(109.5, 415.95), n = c(365, 885), correct = FALSE, alternative ="two.sided")
############################## end of examples ############################## 

# percent correct 
prop.test(x = c(26951, 31699*0.5667), n = c(31699, 31699), correct = TRUE, alternative ="two.sided")

# E CORRECT VS D CORRECT
prop.test(x = c(14051, 12900), n = c(15914, 15785), correct = TRUE, alternative ="greater")

# choseOffer vs optimal 14047
prop.test(x = c(14047, 31699*0.4), n = c(31699, 31699), correct = TRUE, alternative ="two.sided")
# choOE vs choOD
prop.test(x = c(6761, 7286), n = c(15914, 15785), correct = TRUE, alternative ="two.sided")

binom.test(15, 125, 0.05, alternative="two.sided")
# EffectSize: Relative Risk
0.12/0.05
prop.test(x = c(8, 7), n = c(15, 15), correct = TRUE, alternative ="greater")

binom.test(21, 125, 0.05, alternative="two.sided")
# EffectSize: Relative Risk
0.168/0.05
prop.test(x = c(14, 7), n = c(21, 21), correct = TRUE, alternative ="greater")

binom.test(12, 125, 0.05, alternative="two.sided")
# EffectSize: Relative Risk
0.096/0.05
prop.test(x = c(9, 3), n = c(12, 12), correct = TRUE, alternative ="greater")

binom.test(16, 125, 0.05, alternative="two.sided")
# EffectSize: Relative Risk
0.128/0.05
prop.test(x = c(9, 3), n = c(12, 12), correct = TRUE, alternative ="greater")

prop.test(x = c(15, 21), n = c(125, 125), correct = TRUE, alternative ="two.sided")

prop.test(x = c(21, 16), n = c(125, 125), correct = TRUE, alternative ="two.sided")

# p.adjust(p, method = p.adjust.methods, n = length(p))
# 
# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none")

library("pwr")
pwr.r.test(n = NULL , r = 0.33, sig.level = 0.05, power = 0.60)

pwr.r.test(n = 44 , r = NULL, sig.level = 0.05, power = 0.60)

pwr.r.test(n = 54 , r = NULL, sig.level = 0.05, power = 0.60)

pwr.r.test(n = NULL, r = 0.2, sig.level = 0.05, power = 0.60)


pwr.r.test(n = NULL, r = -0.36, sig.level = 0.05, power = 0.65, 
           alternative = c("two.sided", "less","greater"))

pwr.r.test(n = 125, r = NULL, sig.level = 0.05, power = 0.7)




# prop of cell
prop.test(x = c(3, 9), n = c(12, 12), correct = TRUE, alternative ="two.sided")

prop.test(x = c(0.168*125, 0.128*125), n = c(125, 125), correct = TRUE, alternative ="greater")

# NN_decoder chance level 625*0.2=125 
prop.test(x = c(625*0.2452, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.2620, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.2732, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.3960, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.5136, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.3196, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.2204, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.2156, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.1736, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.2060, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.1840, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.2324, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.5136, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.1712, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

# SVM_replication
prop.test(x = c(625*0.2400, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.2852, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.2844, 125), n = c(625, 625), correct = TRUE, alternative ="greater")

prop.test(x = c(625*0.4384, 125), n = c(625, 625), correct = TRUE, alternative ="greater")


# SVM_Modality chance level 250*0.5=125 
prop.test(x = c(250*0.2400, 125), n = c(250, 250), correct = TRUE, alternative ="greater")

# outlier test

setwd("/Users/mwang/Dropbox/01Data/01 Simulation/00Analysis_Simulation/R")
library(foreign)
library(car) # for Boxplot() function
source("showmelm.R")

# LOAD DATA####
Data = read.csv("betas.csv", header=FALSE, sep=",")
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
title(main="Described offer - outcome")

mymodel2 = lm(V4 ~ V3, data = Data)
mymodel2
# create object with model diagnostics
mymodel2.diag = ls.diag(mymodel2)
Data$cooks2 = mymodel2.diag$cooks
Boxplot(Data$cooks2, ylab="Cook's D value", cex = 0.8)
text(1:nrow(Data), Data$cooks2, cex = 0.3, pos = 3)
title(main="Experienced offer - outcome")

mymodel3 = lm(V3 ~ V1, data = Data)
mymodel3
# create object with model diagnostics
mymodel3.diag = ls.diag(mymodel3)
Data$cooks3 = mymodel3.diag$cooks
Boxplot(Data$cooks3, ylab="Cook's D value", cex = 0.8)
text(1:nrow(Data), Data$cooks3, cex = 0.3, pos = 3)
title(main="Described offer - experienced offer")

mymodel4 = lm(V4 ~ V2, data = Data)
mymodel4
# create object with model diagnostics
mymodel4.diag = ls.diag(mymodel4)
Data$cooks4 = mymodel4.diag$cooks
Boxplot(Data$cooks4, ylab="Cook's D value", cex = 0.8)
text(1:nrow(Data), Data$cooks4, cex = 0.3, pos = 3)
title(main="Described outcome - experienced outcome",cex.lab=0.75)




