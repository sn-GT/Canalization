# Liability threshold model to compute expected prevalences at each percentile PGS, under null expectation of no PGSxE 
# For details on the model, please see the Models section on our Rshiny: https://canalization-gibsonlab.shinyapps.io/rshiny/
# ui = t-invcdf(1-Pi) = a meanPRSi + b Env(0/1) + c
# ui stochasiticity 10 iterations, N(ui,sd), sd = residual error of linear model fitted above

# For each disease and exposure 
# PP (prev-percentilePGS) is the dataframe with observed prevalence of disease vs percentile PGS 
# for high vs low of one environmental exposure (for example high vs low for SES). 
# This data will have 200 points: 
# For HIGH env exposure: PGS 1:100 with correspoding observed prevalence
# For LOW env exposure: PGS 1:100 with correspoding observed prevalence
PP <- read.table("example_PrevPerc.txt", header = T, sep = "\t")

# Get names ofhigh vs low group for environmental exposure
l <- names(table(PP$GROUP))

# Overall prevalence of disease
overallprev <- 0.056
# overall threshold t from overall population prevalence
t_overall = qnorm(1-overallprev) #inverse cdf to get t 

# Setting up dataframe to fit regression model: 

# Env (1/0) based on high vs low group of environmental exposure
PP$ENV[PP$GROUP == l[1]] <- 1
PP$ENV[PP$GROUP == l[2]] <- 0

# Compute threshold t at each percentile PGS based on observed prevalence at each PGS.
PP$t_PGS <- qnorm(1-PP$Prev)
PP$overallP <- overallprev
PP$t <- t_overall

# Compute underlying mean ui of the liability model thresholds  
PP$ui <- PP$t - PP$t_PGS

# Fit linear regression model, ui = t-invcdf(1-Pi) = a meanPRSi + b Env(0/1) + c
model <- lm(data = PP, ui ~ meanPRS + ENV)
modfit <- summary(model)

# get regression coefficients
# mean PRS
a1= model$coefficients[2]
# env 
a2 = model$coefficients[3]
# intercept
a3 = model$coefficients[1]

# estimated ui from regression coeff
PP$Estimated_ui <- a1*PP$meanPRS + a2*PP$ENV + a3
# sd = residual error from linear regression model
modsigma <- (modfit$sigma)
  
if(nrow(PP[PP$Prev == 0,]) == 0){
	df <- data.frame()
	# 10 iterations to add stochasticity to underlying estimated ui
	# Computed expected prevalence from estimated ui:  1-cdf(N(ui,1), toverall)
	for(it in 1:10){
		exprev <- c()
		for(estui in PP$Estimated_ui){
			estui_se <- rnorm(1, estui, modsigma)
			# expected prevalence from estimated ui
			expecprev_sub <- (1 - pnorm(t, mean = estui_se, sd = 1))
			exprev <- c(exprev, expecprev_sub)
		}
		PP$ExpectedPrev <- exprev
		PP$iteration <- it
		df <- rbind(df, PP) 
	}	
}
