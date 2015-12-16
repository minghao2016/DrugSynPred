source('ch2scoring_fc.R')
syn_pred = '../../output/predictions/leaderBoard/synergy_matrix.csv'
getPrecision_ch2(syn_pred, threshold=20)
getPrecision_ch2(syn_pred, threshold=30)
getPrecision_ch2(syn_pred, threshold=40)
getGlobalScore_ch2(syn_pred)
# getOneDimScore_ch2(syn_pred)
