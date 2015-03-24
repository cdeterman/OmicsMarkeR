
# binary class feature ranking
svmrfeFeatureRanking(x = vars,
                     y = groups, 
                     c = 0.1,
                     num.rem = 10)

# multiclass
svmrfeFeatureRankingForMulticlass(x = vars,
                                  y = groups, 
                                  c = 0.1,
                                  num.rem = 10)
