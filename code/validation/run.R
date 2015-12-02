source('ch2scoring_fc.R')
conf_pred = '../../output/predictions/leaderBoard/confidence_matrix.csv'

syn_pred = '../../output/predictions/leaderBoard/synergy_matrix.csv'
getPrecision_ch2(syn_pred, threshold=30)