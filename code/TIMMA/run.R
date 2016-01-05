install.packages("QCA")
install.packages('timma_1.2.1.tar.gz', repo=NULL, type='source')
library(timma)

rm(list = ls())

root_folder = '/home/shahin/Dropbox/Dream/'
setwd(paste(root_folder, 'experiment/code/TIMMA', sep=''))

Dream_interaction_binary = read.table('../../Dream_interaction_binary.csv', sep=',', header=T)
Interactions = Dream_interaction_binary[, 2:dim(Dream_interaction_binary)[2]];
row.names(Interactions) = Dream_interaction_binary[, 1]

Dream_sensitivity = read.table('../../Dream_sensitivity.csv', sep=',', header=T)
Sensitivity = Dream_sensitivity[, 2:dim(Dream_sensitivity)[2]];
row.names(Sensitivity) = Dream_sensitivity[, 1]


median_sensitivity<-Sensitivity[, 1]
results<-timma(as.matrix(Interactions), median_sensitivity)

Syn = read.csv('predictedDrugScoring.csv')