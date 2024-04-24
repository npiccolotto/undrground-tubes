library(pwr)
library(MASS)
library(ggplot2)
library(car)
library(psych)
library(plyr)
library(moments)
library(rstatix)
library(emmeans)
library(ggpubr)

setwd('/Users/npiccolotto/Projects/cvast/bssvis/ensemble-set-rendering/eval/')

# load file
qual <- read.csv('qualitative_feedback.csv',colClasses = c("numeric", "factor", "factor", "numeric", "factor"),)

## Confidence

###  Task 1
qual_subset <- qual[qual$task==1 & qual$metric=='confidence',]
describeBy(qual_subset[,c('rating','alpha')], group=qual_subset$alpha)
k <- kruskal.test(qual_subset$rating ~ qual_subset$alpha)
print(k)

d <-  dunn_test(rating ~ alpha, data = qual_subset, p.adjust.method = "bonferroni")
print(d)

'0 significantly better than 1 and 5, 5 better than 1'

###  Task 2
qual_subset <- qual[qual$task==2 & qual$metric=='confidence',]
describeBy(qual_subset[,c('rating','alpha')], group=qual_subset$alpha)
k <- kruskal.test(qual_subset$rating ~ qual_subset$alpha)
print(k)

d <-  dunn_test(rating ~ alpha, data = qual_subset, p.adjust.method = "bonferroni")
print(d)

'1 significantly better than 0 and 5, not so much difference between 0 and 5'

###  Task 3
qual_subset <- qual[qual$task==3 & qual$metric=='confidence',]
describeBy(qual_subset[,c('rating','alpha')], group=qual_subset$alpha)
k <- kruskal.test(qual_subset$rating ~ qual_subset$alpha)
print(k)

d <-  dunn_test(rating ~ alpha, data = qual_subset, p.adjust.method = "bonferroni")
print(d)

'0 significantly worse than 1 and 5, not so much difference between 1 and 5'

###  Task 4
qual_subset <- qual[qual$task==4 & qual$metric=='confidence',]
describeBy(qual_subset[,c('rating','alpha')], group=qual_subset$alpha)
k <- kruskal.test(qual_subset$rating ~ qual_subset$alpha)
print(k)

d <-  dunn_test(rating ~ alpha, data = qual_subset, p.adjust.method = "bonferroni")
print(d)

'5 better than 1, but no difference otherwise' 

## Preference 
###  Task 1
qual_subset <- qual[qual$task==1 & qual$metric=='preference',]
describeBy(qual_subset[,c('rating','alpha')], group=qual_subset$alpha)
k <- kruskal.test(qual_subset$rating ~ qual_subset$alpha)
print(k)

d <-  dunn_test(rating ~ alpha, data = qual_subset, p.adjust.method = "bonferroni")
print(d)

'0 better than 5 better than 1, all very significant'

###  Task 2
qual_subset <- qual[qual$task==2 & qual$metric=='preference',]
describeBy(qual_subset[,c('rating','alpha')], group=qual_subset$alpha)
k <- kruskal.test(qual_subset$rating ~ qual_subset$alpha)
print(k)

d <-  dunn_test(rating ~ alpha, data = qual_subset, p.adjust.method = "bonferroni")
print(d)

'0 significantly worse than 1 and 5, 1 slightly better than 5'

###  Task 3
qual_subset <- qual[qual$task==3 & qual$metric=='preference',]
describeBy(qual_subset[,c('rating','alpha')], group=qual_subset$alpha)
k <- kruskal.test(qual_subset$rating ~ qual_subset$alpha)
print(k)

d <-  dunn_test(rating ~ alpha, data = qual_subset, p.adjust.method = "bonferroni")
print(d)

'0 significantly worse than 1 and 5, 1 slightly better than 5'

###  Task 4
qual_subset <- qual[qual$task==4 & qual$metric=='preference',]
describeBy(qual_subset[,c('rating','alpha')], group=qual_subset$alpha)
k <- kruskal.test(qual_subset$rating ~ qual_subset$alpha)
print(k)

d <-  dunn_test(rating ~ alpha, data = qual_subset, p.adjust.method = "bonferroni")
print(d)

'5 better than 1 ***, 0 better than 1 **, no difference between 0 and 5'

# Check correlation between confidence and accuracy
acc <- read.csv('survey_long.csv',colClasses = c("numeric", "numeric", "factor", "factor", "factor", "factor"),)

qual_subset <- qual[qual$metric == 'confidence',]
qual_subset <- dplyr::arrange(qual_subset, participant, task, alpha)

acc_subset <- acc[acc$metric=='a' & acc$pipeline=='h',]
acc_subset <- dplyr::arrange(acc_subset, participant, task, alpha)

cor(qual_subset$rating, acc_subset$value, method='spearman')


# Check correlation between confidence and preference - are these the same answers?
qual_subset1 <- qual[qual$metric == 'preference',]
qual_subset1 <- dplyr::arrange(qual_subset1, participant, task, alpha)

qual_subset2 <- qual[qual$metric == 'confidence',]
qual_subset2 <- dplyr::arrange(qual_subset2, participant, task, alpha)

cor(qual_subset1$rating, qual_subset2$rating, method='spearman')

'kind of'
