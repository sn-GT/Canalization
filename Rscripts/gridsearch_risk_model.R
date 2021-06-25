setwd("~/Dropbox (GaTech)/IBD/Canalization/CAD/Data/")
library(data.table)
library(optparse)

option_list = list(
  make_option(c("-f", "--p"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

########### Real Data Stats #############################
field = "f.1070.0.0"
# Prevalence vs percentile data for all 153 exposures and 1 trait
ppdf <- read.table("PPdf_1e03.txt", header= T, sep = "\t")
ppdf <- ppdf[ppdf$field == field,]
ppdf$GROUP <- as.character(ppdf$GROUP)

# title/description of exposure
desname <- ppdf$desc[ppdf$field == field][1]
g <- names(table(ppdf$GROUP))

# 1 for high exposure and 2 for low exposure
g_index = 2

# For real data, get the mean, median, min, max, derivative-mean and derivative-sd: 
# Prev-perc stats for high/low group in real data
real_mean = mean(ppdf$Prev[ppdf$GROUP == g[g_index]])
real_min = min(ppdf$Prev[ppdf$GROUP == g[g_index]])
real_max = max(ppdf$Prev[ppdf$GROUP == g[g_index]])
real_median = median(ppdf$Prev[ppdf$GROUP == g[g_index]])

# Derivative of Prev-perc curve for real data
sub <- ppdf[ppdf$GROUP == g[g_index],]
derivative=diff(predict(lm(Prev~poly(PGS,3, raw = T), data = sub)))
real_mean_deriv <- mean(derivative)
real_sd_deriv <- sd(derivative)

# Get the PRS mean and sd in low and high exposure for that trait
prs_stats <- fread("Canalization_Index/results/for_high_low_plots/CAD_separate.txt", header= T, sep = "\t")
prs_stats <- prs_stats[prs_stats$FIELD  == field,]

# PRS stats for low/high exposure 
real_prs_mean = prs_stats$mean_PRS[g_index]
real_prs_sd = prs_stats$std_PRS[g_index]

##################################################

# Simulate Data
# Disease Risk Model = SP + MBL* exp(PRS/CF) + Env
# where SP non-modifiable adjustment factor to ensure median risk is as observed​
# MBL: modifiable baseline prevalence interacting with PRS in multiplicative manner​
# PRS: with mean 0 and observed sd​
# CF: calibration factor adjusting effect of PRS on disease​
# Env: Individual stochastic risk dist (0, 2)

# GRID SEARCH
# Define grid on which to perform search to find the best optimal parameters for simulated data.
# Best optimal parameters are obtained by minimizing the error between mean, median, range, derivative dist between real and simulated data.
Env = seq(-4,6, by = 0.5)
BL = seq(0,1, by = 0.1)
SF = seq(1,3, by = 0.2)

# get total number of combinations to perform search on
print(length(BL)*length(Env)*length(SF))

# Simulate data for 10000 individuals
totsamp=10000
# creating IDs for 10k individuals
ID <- seq(10001,20000, by = 1 )

# PRS mean = 0 and sd = observed sd in high/low exposure
prs_norm <- qnorm(runif(totsamp),0, 1)
df <- data.frame("ID"=ID, "PRS_norm"= prs_norm)
df <- df[order(df$PRS_norm, decreasing = T),]
df$PRS <- (df$PRS_norm/(1/real_prs_sd))

# Env mean = 0, sd = 2
df$Ind <- qnorm(runif(df$PRS),0, 2)


param_df <- data.frame()
prev_eqn_df <- data.frame()
# 1. Iterating over all combinations SP, MBL, CF as provided in the grid above 
# 2. Computing Risk from model Risk = SP + MBL*exp(PRS/CF) + Env
# 3. Compute prevalence vs percentile for simulated data: x axis: Percentile PRS, y axis = prevalence as the average risk within each percentile bin
# 4. Compute Error = sum of square deviations between mean, median, min, max, deivative mean, sd between real data stats and simulated data
# 5. Pick the combination of SP, MBL, CF having minimum error
# 6. Repeat this 10 times and get best optimal parameter combinations (shown in fig S5)
for(b in BL){
  for(e in Env){
    for(s in SF){
      print(b)
      # DISEASE RISK MODEL
      df$Risk <- e + df$Ind + (b*exp(df$PRS/s))
      
      # Prevalence vs percentile with disease risk model
      # x axis is percentile of PRS, y axis is Prevalence as the average Risk (obtained from above disease risk eqn) within that PRS percentile
      n=100
      prevdf <- data.frame()
      for(i in 0:(n-1)){
        if(i == 0){
          dx = df[which(df$PRS >= quantile(df$PRS, i/n) & df$PRS <= quantile(df$PRS, (i+1)/n)),] 
          prev = mean(dx$Risk)
        }
        else {
          dx = df[which(df$PRS > quantile(df$PRS, i/n) & df$PRS <= quantile(df$PRS, (i+1)/n)),]
          prev = mean(dx$Risk)
        }
        
        if(is.na(prev)){prev = 0}
        pdf = data.frame("PGS" = i+1, "Prev" = prev) 
        prevdf <- rbind(prevdf, pdf)
      }

      # Get the mean, median, min, max, derivative of simulated data
      meanP <- mean(prevdf$Prev)
      minP <- min(prevdf$Prev)
      maxP <- max(prevdf$Prev)
      medianP <- median(prevdf$Prev)
      
      derivativeP <- diff(predict(lm(Prev~poly(PGS,3, raw = T), data = prevdf)))
      mean_derivP <- mean(derivativeP)
      sd_derivP <- sd(derivativeP)
      
      # Compute Error as sum of squared differences between real and simulated data 
      err <- (real_mean - meanP)^2 + (real_min - minP)^2 + (real_max - maxP)^2 + (real_median - medianP)^2 + 
        (real_mean_deriv - mean_derivP)^2 + (real_sd_deriv - sd_derivP)^2 
      
      prevdx <- prevdf
      prevdx$EQN <- paste0("SP",e,", MBL",b,", SF", s, ", Err", round(err,2))
      prevdx$SP <- e
      prevdx$MBL <- b
      prevdx$SF <- s
      prevdx$Err <- err
      prevdx$realmean <- real_mean
      prevdx$computedmean <- meanP
      prevdx$realmedian <- real_median
      prevdx$computedmedian <- medianP
      prevdx$realmin <- real_min
      prevdx$computedmin <- minP
      prevdx$realmax <- real_max
      prevdx$computedmax <- maxP
      prevdx$realderivmean <- real_mean_deriv
      prevdx$computed_derivmean <- mean_derivP
      prevdx$realderivsd <- real_sd_deriv
      prevdx$computed_derivsd <- sd_derivP
      
      prev_eqn_df <- rbind(prev_eqn_df, prevdx)
      
    }
  }
}
#prev_eqn_df$iter <- opt$p
prev_eqn_df$field <- field
prev_eqn_df$desc <- desname
# Pick the parameter comibnation with minimum error
finap_prev <- prev_eqn_df[prev_eqn_df$Err == min(prev_eqn_df$Err),]

#write.table(finap_prev, paste0("Canalization_Index/Parameter_search/multi_iters_param/CAD/",field,"_low.txt"), 
 #           row.names = F, col.names = T, quote = F, sep = "\t")


###################################################################


