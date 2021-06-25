path=""
setwd(path)

library(data.table)
library(ggplot2)
library(cowplot)
library(scales)
library(plyr)
require(gridExtra)
require(grid)
library(ggpmisc)

# prevalence vs percentile for all exposures by high and low categoires
PP_df <- fread("PPdf.txt", header = T, sep = "\t")
# prevalence vs percentile for all exposures, randomly taken high and low categories, keeping sample size the same as real data
rand_iter <- fread("Random_iterations.txt", header = T, sep = "\t")
# list of all 153 exposures field ids
cont <- names(table(PP_df$field))

cubic_curves <- data.frame()
derivative_cubicfit <- data.frame()
coefdf <- data.frame()
canalization_linear <- data.frame()
canalization_derivative <- data.frame()
for(k in cont){
      fld=k
      pdf_males <- PP_df[PP_df$field == fld,]

      # Prevalence vs percentile divided by high and low exposure group
      p1 <- ggplot(pdf, aes(x=PGS, y=Prev, col = GROUP, group = GROUP)) + 
        geom_point( size = 2.3) + theme_bw() + xlab("Percentile of PRS-BMI") + ylab("Prevalance of Obesity (%)") + 
        ggtitle("GRS p < 1e-03", subtitle = paste0("Red=", pdf$HIGH[1], ", Blue=",pdf$LOW[1])) + 
        geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1, se = F) + theme(legend.position="none", plot.title = element_text(face="bold"))

      ############################################################################################
      
      # Assesment of Canalization H - L curves for real data and random iterations  + Linear slope of H-L
      # - Plotting H-L curves and see that they are deviated from H-L curves of 100 random iterations
      # - Linear slope of H-L for real data is highest compared to 100 random iterations.
      
      # High - Low difference of prevalence over percentile PRS
      lev <- table(pdf_males$GROUP)
      lev <- lev[lev !=0]
      set1 <- pdf_males[pdf_males$GROUP == names(lev)[1],]
      set2 <- pdf_males[pdf_males$GROUP == names(lev)[2],]
      hminusl <- merge(set1, set2, by ="PGS")
      # H-L difference between prevalence over percentile of PRS 
      hminusl$diff <- hminusl$Prev.x - hminusl$Prev.y
      hminusl$diff_scaled <- hminusl$scaled.x - hminusl$scaled.y
      names(hminusl)[5] <- "FIELD"
      names(hminusl)[6] <- "desc"
      names(hminusl)[7] <- "HIGH"
      names(hminusl)[8] <- "LOW"
      hminusl <- hminusl[,c("PGS", "diff", "diff_scaled","FIELD", "desc", "HIGH", "LOW")]
      hminusl$iter <- "REAL_DATA"
      
      # H - L from random iterations
      iter_df <- rand_iter[rand_iter$FIELD ==k,]
      
      # combine real data and random iterations
      fiter <- rbind(hminusl, iter_df)
      
      # setting alpha for faded lines
      sub <- fiter
      sub$ALPHA <- "FALSE"
      sub$ALPHA[sub$iter =="REAL_DATA"] <- "TRUE"
      alphas1 <- ifelse(sub$ALPHA, 1, 0.85)
      
      # Plot H-L for real data and random iterations 
      p2x <- ggplot(sub, aes(x=PGS, y=diff, col = iter)) + theme_bw() + xlab("Percentile PRS-T2D") + 
        ylab(" ") + 
        ggtitle("HIGH - LOW") + 
        geom_line(stat="smooth",method = "lm", formula = y ~ poly(x, 3, raw = T),size = 1.5, aes(alpha = alphas1)) + 
        scale_alpha(guide=FALSE) + scale_color_discrete(name = "Iteration", guide=F)
      p2 <- p2x + scale_color_viridis_d(direction = -1, option = "plasma", guide=F) 
      
      # Plot H-L SCALED for real data and random iterations 
      px <- ggplot(sub, aes(x=PGS, y=diff_scaled, col = iter)) + theme_bw() + xlab(" ") + 
        ylab(" ") + 
        ggtitle("HIGH - LOW: Mean Centered") + 
        geom_line(stat="smooth",method = "lm", formula = y ~ poly(x, 3, raw = T),size = 1.5, aes(alpha = alphas1)) + 
        scale_alpha(guide=FALSE) + scale_color_discrete(name = "Iteration", guide=F) 
      p3 <- px + scale_color_viridis_d(direction = -1, option = "plasma", guide=F) 


      # Linear slope of H-L for real data must be greater than slope of H-L for 100 random iterations, to qualify for canalization.
      linear <- ddply(sub, ~iter, summarise, linslope=lm((predict(lm(diff_scaled ~ PGS))) ~ PGS)[[1]][2] )
      linear$abs <- abs(linear$linslope)
  
      if(linear[which.max(linear$abs), "iter"] == "REAL_DATA"){
        canln <- data.frame("FIELD" = fld, "LnSlope" = linear$abs[linear$iter == "REAL_DATA"],"LnSlopeDIR" = linear$linslope[linear$iter == "REAL_DATA"], "Canalization_linear" = "Yes")
      }
      else {
        canln <- data.frame("FIELD" = fld, "LnSlope" = max(linear$abs), "LnSlopeDIR" = max(linear$linslope), "Canalization_linear" = "No")
      }
      canalization_linear <- rbind(canalization_linear, canln)

      # storing data: H-L real data and random iteration for all exposures
      cubic_curves <- rbind(cubic_curves, sub)
      # storing coefficients of cubic fit of H-L curves for real and random iterations
      coef_cub <- ddply(sub,~ iter,summarise, x=lm(diff~poly(PGS,3, raw =T))$coefficients[2], x2=lm(diff~poly(PGS,3, raw =T))$coefficients[3], x3=lm(diff~poly(PGS,3, raw =T))$coefficients[4])
      coef_cub$FIELD <- fld
      coefdf <- rbind(coefdf, coef_cub)

      ###############################################################################################

      # Assessment of canalization: 
      # - Derivative of cubic fit to see if H-L curve is monotonically increase/decreasing
      # - Canalization if derivative distribution of H-L for real data is  > 0 or < 0 (shifted relative to derivative dist of 100 random iterations),
      #   the median of derivative for real data must be the highest.

      # Derivative of cubic fit
      sumx <- ddply(sub,~ iter,summarise, deriv=diff(predict(lm(diff~poly(PGS,3)))))
      sumx$FIELD <- sub$FIELD[1]
      sumx$desc <- sub$desc[1]
      sumx$GROUP<- "RANDOM_ITERATION"
      sumx$GROUP[sumx$iter == "REAL_DATA"] <- "REAL_DATA"
      
      # setting alpha for random iteration vs real data boxplots
      sumx$ALPHA <- "FALSE"
      sumx$ALPHA[sumx$iter =="REAL_DATA"] <- "TRUE"
      alphas <- ifelse(sumx$ALPHA, 0.9, 0.4)
      
      derivative_cubicfit <- rbind(derivative_cubicfit, sumx)
      
      # boxplots showing the distribution of derivative of H-L cubic curves for real data vs random iterations.
      # For exposures showing canalization, the H-L is monotonically increasing/decreasing, 
      # so the median of the real data is shifted relative to random iterations. 
      phist <- ggplot(data=sumx, aes(x=reorder(iter, deriv, FUN=median), y=deriv, group=iter, fill=iter, color = iter)) +
      geom_boxplot() + gghighlight(iter=="REAL_DATA") + theme_classic() + ggtitle(" ") +
      ylab("Derivative of Cubic Fit") + xlab("Iterations")
      p4 <- phist + scale_fill_viridis_d(direction = 1, option = "plasma") +
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
        theme(legend.position="top", legend.title = element_blank(), legend.text = element_text(face="bold", size = 12))  

      # Median of the derivative of H-L cubic curve for real data must be greater than 100 random iterations to qualify as canalization  
      mx2 <- ddply(sumx,~ iter,summarise, medians=median(deriv))
      mx2$abs_median <- abs(mx2$medians)
      
      # Get max median and check if it is for real data or one of random iteration
      if(mx2[which.max(mx2$abs_median),"iter"] =="REAL_DATA"){
        candx <- data.frame("Median" = mx2$medians[mx2$iter == "REAL_DATA"], "Realdatamedian" = mx2$abs_median[mx2$iter == "REAL_DATA"], "Altmedian" = min(mx2$abs_median), "Canalization_deriv" = "Yes", "Canalization_Index" = ((mx2$abs_median[mx2$iter == "REAL_DATA"])/min(mx2$abs_median)))
      }
      else if(mx2[which.min(mx2$abs_median),"iter"] =="REAL_DATA"){
        candx <- data.frame("Median" = mx2$medians[mx2$iter == "REAL_DATA"],"Realdatamedian" = mx2$abs_median[mx2$iter == "REAL_DATA"], "Altmedian" = max(mx2$abs_median), "Canalization_deriv" = "Yes", "Canalization_Index" = ((max(mx2$abs_median))/(mx2$abs_median[mx2$iter == "REAL_DATA"])))
      }
      else { candx <- data.frame("Median" = mx2$medians[mx2$iter == "REAL_DATA"],"Realdatamedian" = max(mx2$abs_median) , "Altmedian" = min(mx2$abs_median), "Canalization_deriv" = "No", "Canalization_Index"=0)}
      
      candx$FIELD <- sumx$FIELD[1]
      candx$desc <- sumx$desc[1]
      canalization_derivative <- rbind(canalization_derivative, candx)    
}


# To qualify "Yes" for Canalization, 2 conditions must be met
# 1. Highest Linear Slope: Slope of H-L for real data is higher than 100 random iterations
# 2. Highest Median Derivative: Median of deivative of H_L for real data is higher than 100 random iterations

canalization  <- merge(canalization_derivative, canalization_linear, by="FIELD")
# final canalization table with field name of exposure, description, coefficeint of cubic fit x, x2, x3, median of derivative for iterations, real data median,
# linear slope, absolute linear slope, canalization_derivative, canalization_linear
canalization <- canalization[,c("FIELD", "desc","x", "x2", "x3", "Median", "Realdatamedian","LnSlopeDIR", "LnSlope", "Canalization_deriv", "Canalization_linear")]
canalization$CANALIZATION <- "No"
# canalization is "yes" if both canalization_deriv and canalization_linear is True
canalization$CANALIZATION[canalization$Higest_median_derivative == "Yes" & canalization$Higest_linear_slope == "Yes"] <- "Yes"

# save final table
write.table(canalization, "Canalization_BMI_153exposures.txt", row.names = F, col.names = T, quote = F, sep = "\t")

