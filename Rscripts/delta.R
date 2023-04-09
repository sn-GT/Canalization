

# Delta = right tail deviation - left tail deviation 

# right tail deviation : Difference in mean prevalence between high and low risk envionemt at PRS > 2sd above the mean
# left tail deviation : Difference in mean prevalence between high and low risk envionemt at PRS < -2sd above the mean


# PP is the the prevalence percentile dataframe with columns:
# PGS = percentile, Prev, GROUP = 2 env groups, meanPRS = mean PRS at each percentile  

#   PGS      Prev   meanPRS    n GROUP 
# 1   1 0.5687204 -2.664593 1055 White 
# 2   2 0.4743833 -2.178900 1054 White 
# 3   3 0.4743833 -1.969699 1054 White 
# 4   4 1.1385199 -1.818655 1054 White 

# 2 env groups: 
l = names(table(PP$GROUP))

# prev-perc for high-risk env
PP_high <- PP[PP$GROUP == paste0(l[1]),]
# get threshold at meanPRS +/- 2sd 
rightcut <- mean(PP_high$meanPRS) + 2*sd(PP_high$meanPRS)
leftcut <- mean(PP_high$meanPRS) - 2*sd(PP_high$meanPRS)
# Mean prev for high-risk curve at left and right tail
high_lefttail <- mean(PP_high$Prev[PP_high$meanPRS < leftcut])
high_righttail <- mean(PP_high$Prev[PP_high$meanPRS > rightcut])

# prev-perc for low-risk env
PP_low <- PP[PP$GROUP == paste0(l[2]),]
# get threshold at meanPRS +/- 2sd 
rightcut <- mean(PP_low$meanPRS) + 2*sd(PP_low$meanPRS)
leftcut <- mean(PP_low$meanPRS) - 2*sd(PP_low$meanPRS)
# Mean prev for low-risk curve at left and right tail
low_lefttail <- mean(PP_low$Prev[PP_low$meanPRS < leftcut])
low_righttail <- mean(PP_low$Prev[PP_low$meanPRS > rightcut])

# left and right tail deviation
righttail_dev <- high_righttail - low_righttail
lefttail_dev <- high_lefttail - low_lefttail
delta <- righttail_dev - lefttail_dev

