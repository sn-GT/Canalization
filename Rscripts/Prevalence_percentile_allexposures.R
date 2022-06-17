# Compute prevalence vs percentile for all exposures
path="~/Dropbox (GaTech)/IBD/Canalization/"

library(data.table)
library(ggplot2)
library(cowplot)
library(scales)

# PRS 
prs <- data.frame(fread("PRS_1e03_1MBr2.sscore", header = T))
names(prs)[1] <- "FID"

# Trait: For example, BMI.
# define Obesity case-control
bmi <- fread("bmi_all_individuals.txt")
names(bmi)[1] <- "FID"
bmi <- na.omit(bmi)
bmi$STATUS[bmi$f.21001.0.0 >= 30] <- 1
bmi$STATUS[bmi$f.21001.0.0 < 30] <- 0

# Environmental exposures field name and its title description file
field_dict <- read.csv(paste0(path, "UKB_exposures.csv"))
cont <- field_dict$feid

PP_df <- data.frame()
for(k in cont){
    print(k)
    field=k

    # reading environmental exposure file with high and low categories defined for each exposure
    # columns: Individual IID, phenotype value from UKB, GROUP with high vs low group category
    pheno_sub <- fread(paste0(path,"/UKB_Categories/Categories/",k,".txt"), header = T, sep = "\t")

    # merged df having PRS, obesity case-control and high and low groups for exposures
    pre_df <- merge(prs_pheno_df, pheno_sub, by="FID")
    names(pre_df)[names(pre_df) == "SCORE1_SUM"] <- "PRS"
    pre_df$NUM_HIGH <- table(pre_df$GROUP)[[1]]
    pre_df$NUM_LOW <- table(pre_df$GROUP)[[2]]
    
    # computing prevalence of obesity over percentile bins of PRS
    PP <- data.frame()
    # computing for HIGH and LOW groups separately
    l <- names(table(pre_df$GROUP))
    for(t in l){
        n=100
        df <- data.frame()
        for(i in 0:(n-1)){
            print(i)
            if(i == 0){
                dx <- br_grs[br_grs$GROUP == t,]
                dx1 <- dx[which(dx$invGRS >= quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
            }
            else {
                dx <- br_grs[br_grs$GROUP == t,]
                dx1 <- dx[which(dx$invGRS > quantile(dx$invGRS, i/n) & dx$invGRS <= quantile(dx$invGRS, (i+1)/n)),]
            }
            #prev = data.frame(table(dx1$STATUS)/nrow(dx1))[2,2] * 100
            prev=table(dx1$STATUS)[2]/nrow(dx1)*100
            if(is.na(prev)){prev = 0}
            df2 = data.frame("PGS" = i+1, "Prev" = prev) 
            df <- rbind(df, df2)
        }
        df$GROUP <- paste0(t)
        df$scaled <- scale(df$Prev, scale = F)
        PP <- rbind(PP, df)   
    }
    PP$field <- field
    PP$desc <- field_dict$Desc[field_dict$id == field]
    PP$HIGH <- table(pre_df$GROUP)[[1]]
    PP$LOW <- table(pre_df$GROUP)[[2]]
    PP_df <- rbind(PP_df, PP)
  
    # Plot prevalence vs percentile divided by high and low exposure group 
    p1 <- ggplot(PP, aes(x=PGS, y=Prev, col = GROUP, group = GROUP)) + 
      geom_point( size = 2.3) + theme_bw() + xlab("Percentile of PRS-BMI") + ylab("Prevalance of Obesity (%)") + 
       geom_smooth(method = "lm", formula = y ~ poly(x, 3, raw = T), size = 1, se = F) + 
      # adding field name of the exposure, it title description and number of samples in high and low exposure groups in the plot
      ggtitle(paste0(field, ": ", field_dict$Desc[field_dict$id == field]), 
              subtitle = paste0("HIGH=", table(pre_df$GROUP)[[1]], ", LOW=",table(pre_df$GROUP)[[2]])) 
}

# Save prevalence vs percentile for all exposures
write.table(PP_df, "PrevPerc_BMI_PRS1e03_allexposures.txt", row.names = F, col.names = T, quote =F, sep = "\t")
