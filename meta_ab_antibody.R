library(ggplot2)
library (meta)
library(dplyr)
sessionInfo()

R.version.string
getwd()
setwd("/Users/workingdirectory")

# Data import####
df1 <- read.csv("meta_ab_antibody_data.csv", header = TRUE)

# Calculate the mean_cdrsb ####
# of the baseline age, mmse, apoe4 carrier percentage, and amyloid beta burden on PET
df1$age_mean <- with(df1, (age_e*Ne_c + age_c*Nc_c) / (Ne_c+Nc_c))
df1$mmse_mean <- with(df1, (mmse_e*Ne_c + mmse_c*Nc_c) / (Ne_c+Nc_c))
df1$apoe4_mean <- with(df1, (apoe4_e*Ne_c + apoe4_c*Nc_c) / (Ne_c+Nc_c))
df1$ab_centi_mean <- with(df1, (ab_centi_e*Ne_c + ab_centi_c*Nc_c) / (Ne_c+Nc_c))

#mean_adascog
df1$age_mean_a <- with(df1, (age_e*Ne_a + age_c*Nc_a) / (Ne_a+Nc_a))
#df1 <- df1 %>% mutate (age_mean2 = if_else(Number == 10|Number==11, 73.660, age_mean))
df1$mmse_mean_a <- with(df1, (mmse_e*Ne_a + mmse_c*Nc_a) / (Ne_a+Nc_a))
df1$apoe4_mean_a <- with(df1, (apoe4_e*Ne_a + apoe4_c*Nc_a) / (Ne_a+Nc_a))
df1$ab_centi_mean_a <- with(df1, (ab_centi_e*Ne_a + ab_centi_c*Nc_a) / (Ne_a+Nc_a))
df1$adascog_mean_a<-with(df1, (Me_pre_a*Ne_a+Mc_pre_a*Nc_a)/ (Ne_a+Nc_a))

#mean_ariae
df1$age_mean_ariae <- with(df1, (age_e*Ne_ariae + age_c*Nc_ariae) / (Ne_ariae+Nc_ariae))
df1$mmse_mean_ariae <- with(df1, (mmse_e*Ne_ariae + mmse_c*Nc_ariae) / (Ne_ariae+Nc_ariae))
df1$apoe4_mean_ariae <- with(df1, (apoe4_e*Ne_ariae + apoe4_c*Nc_ariae) / (Ne_ariae+Nc_ariae))
df1$ab_centi_mean_ariae <- with(df1, (ab_centi_e*Ne_ariae + ab_centi_c*Nc_ariae) / (Ne_ariae+Nc_ariae))

#mean_SAE
df1$age_mean_SAE <- with(df1, (age_e*Ne_SAE + age_c*Nc_SAE) / (Ne_SAE+Nc_SAE))
df1$mmse_mean_SAE <- with(df1, (mmse_e*Ne_SAE + mmse_c*Nc_SAE) / (Ne_SAE+Nc_SAE))
df1$apoe4_mean_SAE <- with(df1, (apoe4_e*Ne_SAE + apoe4_c*Nc_SAE) / (Ne_SAE+Nc_SAE))
df1$ab_centi_mean_SAE <- with(df1, (ab_centi_e*Ne_SAE + ab_centi_c*Nc_SAE) / (Ne_SAE+Nc_SAE))


#mean_death
df1$age_mean_death <- with(df1, (age_e*Ne_death + age_c*Nc_death) / (Ne_death+Nc_death))
df1$mmse_mean_death <- with(df1, (mmse_e*Ne_death + mmse_c*Nc_death) / (Ne_death+Nc_death))
df1$apoe4_mean_death <- with(df1, (apoe4_e*Ne_death + apoe4_c*Nc_death) / (Ne_death+Nc_death))
df1$ab_centi_mean_death <- with(df1, (ab_centi_e*Ne_death + ab_centi_c*Nc_death) / (Ne_death+Nc_death))

#mean_Head
df1$age_mean_head <- with(df1, (age_e*Ne_Head + age_c*Nc_Head) / (Ne_Head+Nc_Head))
df1$mmse_mean_head<- with(df1, (mmse_e*Ne_Head + mmse_c*Nc_Head) / (Ne_Head+Nc_Head))
df1$apoe4_mean_head <- with(df1, (apoe4_e*Ne_Head + apoe4_c*Nc_Head) / (Ne_Head+Nc_Head))
df1$ab_centi_mean_head <- with(df1, (ab_centi_e*Ne_Head + ab_centi_c*Nc_Head) / (Ne_Head+Nc_Head))

#mean_fall
df1$age_mean_fall <- with(df1, (age_e*Ne_fall + age_c*Nc_fall) / (Ne_fall+Nc_fall))
df1$mmse_mean_fall<- with(df1, (mmse_e*Ne_fall + mmse_c*Nc_fall) / (Ne_fall+Nc_fall))
df1$apoe4_mean_fall <- with(df1, (apoe4_e*Ne_fall + apoe4_c*Nc_fall) / (Ne_fall+Nc_fall))
df1$ab_centi_mean_fall <- with(df1, (ab_centi_e*Ne_fall + ab_centi_c*Nc_fall) / (Ne_fall+Nc_fall))

#mean_dizziness
df1$age_mean_dizzi <- with(df1, (age_e*Ne_dizzi + age_c*Nc_dizzi) / (Ne_dizzi+Nc_dizzi))
df1$mmse_mean_dizzi<- with(df1, (mmse_e*Ne_dizzi + mmse_c*Nc_dizzi) / (Ne_dizzi+Nc_dizzi))
df1$apoe4_mean_dizzi <- with(df1, (apoe4_e*Ne_dizzi + apoe4_c*Nc_dizzi) / (Ne_dizzi+Nc_dizzi))
df1$ab_centi_mean_dizzi <- with(df1, (ab_centi_e*Ne_dizzi + ab_centi_c*Nc_dizzi) / (Ne_dizzi+Nc_dizzi))


#CDR-SB####
#Calculation of SD and Mean change of the experimental group :CDR-SB
df1 <- df1 %>%
  mutate(
    SE_e_cdr_calc=ifelse (is.na(SE_e_c),(Me_u_c-Me_l_c)/(2*1.96),SE_e_c), #1 SE, #2 Calculation of SE from 95%CI
    SE_e_cdr_calc=ifelse(is.na(SE_e_cdr_calc), sqrt(((diff_u_c-diff_l_c)/(2*1.96))^2-(SE_c_c)^2),SE_e_cdr_calc), #3 Calculation of SE(experimental group) from 95%CI of difference 
    SE_c_cdr_calc=ifelse (is.na(SE_c_c),(Mc_u_c-Mc_l_c)/(2*1.96),SE_c_c),#1 SE, #2 Calculation of SE from 95%CI
    SD_e_cdr_calc=ifelse (is.na(SDe_c),SE_e_cdr_calc*sqrt(Ne_c),SDe_c), #Calculation of SD from SE
    SD_c_cdr_calc=ifelse (is.na (SDc_c),SE_c_cdr_calc*sqrt(Nc_c), SDc_c), #Calculation of SD from SE
    Me_calc_cdr=ifelse (is.na(Me_c), Mc_c+diff_m_c, Me_c),#Calculation of Me(mean change of the experimental group) by adding difference to the Mc_c (mean change of the control group)
    SE_diff_cdr=ifelse(is.na(diff_l_c), sqrt((SDe_c^2/Ne_c)+(SDc_c^2/Nc_c)), (diff_u_c-diff_l_c)/(2*1.96))
  )

#Combine low-dose and high-dose: CDR-SB#### 
df1 <- df1 %>%
  group_by(paper) %>%
  mutate(
    N1 = sum(Ne_c[low1_high2_other0 == 1]), # N of low dose group #low1_high2_ohter0: low dose group==1, high dose group==2, other(one dose only in the paper)==0
    N2 = sum(Ne_c[low1_high2_other0 == 2]), # N of high dose group
    M1 = sum(Me_calc_cdr[low1_high2_other0 == 1]), # Mean change of low dose group
    M2 = sum(Me_calc_cdr[low1_high2_other0 == 2]), # Mean change of high dose group
    SD1 = sum(SD_e_cdr_calc[low1_high2_other0 == 1]), # SD of low dose group
    SD2 = sum(SD_e_cdr_calc[low1_high2_other0 == 2]), #SD of high dose group
  ) %>%
  mutate(
    Ne_Combined1_cdr = sum(Ne_c[low1_high2_other0 %in% c(1, 2)]), #Sample size of the combined groups
    Me_Combined1_cdr = (N1 * M1 + N2 * M2) / (N1 + N2),# Mean change of combined groups 
    Me_Combined2_cdr=ifelse(is.na(Me_Combined1_cdr), Me_calc_cdr,Me_Combined1_cdr),
    SDe_Combined1_cdr = sqrt(((N1 - 1) * SD1^2 + (N2 - 1) * SD2^2 + (N1 * N2 / (N1 + N2)) * (M1^2 + M2^2 - 2 * M1 * M2)) / (N1 + N2 - 1)), #  Combined SD
  ) %>%
  ungroup() # ungroup

df1 <- df1 %>%
  mutate(
    DIFF_combined_cdr = ifelse(low1_high2_other0==0, diff_m_c, Me_Combined1_cdr-Mc_c),
    DIFFSE_combined_cdr=ifelse(low1_high2_other0==0, SE_diff_cdr, sqrt((SDe_Combined1_cdr^2/Ne_Combined1_cdr)+(SD_c_cdr_calc^2/Nc_c)))
  )


#Exclude Low dose（currently, the same entries are included for both low dose and high dose, so select the high dose as representative）
#Note that control(placebo) groups are the same between the low and the high-dose groups
df2<-subset(df1, low1_high2_other0 %in% c(0,2))

#sensitivity analysis
#df2 <- filter(df2, antibody == "1_Lecanemab" |antibody=="2_Gantenerumab" | antibody == "3_Donanemab" |antibody=="5_Solanezumab"|antibody=="6_Bapineuzumab") #modify as needed 
df2$antibody <- factor(df2$antibody, levels = c( "1_Lecanemab","2_Gantenerumab", "3_Donanemab", "4_Aducanumab", "5_Solanezumab", "6_Bapineuzumab"  ))

# Metagen_CDRSB####
metagen_cdrsb <- metagen(
  TE = df2$DIFF_combined_cdr,         # Mean Difference
  seTE = df2$DIFFSE_combined_cdr,    # SE of Mean Difference
  data = df2,        
  studlab = df2$Trial_3, 
  sm = "MD",
  subgroup=df2$antibody,
  method.tau="PM",
  prediction=TRUE
)
#forest cdrsb####
tiff("forest_cdrsb_com.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metagen_cdrsb, xlim=c(-1.0,1.0), label.right="Placebo better", label.left="Antibody better", col.subgroup="black") #Forrestplot
dev.off()

#funnel####
funnel(metagen_cdrsb, ylim=c(0.3,0),xlim=c(-0.8,0.4),xlab="Mean Difference")
metabias(metagen_cdrsb, method="linreg")


##Function to run meta-regression analysis####
run_meta_reg<-function(metagen_result,modifier){
  meta_reg_result <- metareg(metagen_result, as.formula(paste("~", modifier)))
  return(summary(meta_reg_result))
}

#Set the levels as needed
df2$AD_stage<-factor(df2$AD_stage, levels=c("mildmoderate","early" ))
df2$antibody <- factor(df2$antibody, levels = c("6_Bapineuzumab" , "5_Solanezumab", "4_Aducanumab", "3_Donanemab","2_Gantenerumab","1_Lecanemab"  ))
df2$anti_type<-factor(df2$anti_type, levels=c("humanized", "human"))
df2$mechanism<-factor(df2$mechanism, levels=c("monomers","oligoaggre","aggregates" ))
df2$Biological<-factor(df2$Biological, levels=c("non","effective" ))

#meta_regression_CDRSB
run_meta_reg(metagen_cdrsb, "age_mean")
run_meta_reg(metagen_cdrsb, "mmse_mean")
run_meta_reg(metagen_cdrsb, "apoe4_mean")
run_meta_reg(metagen_cdrsb, "ab_centi_mean")

#run_meta_reg(metagen_cdrsb, "AD_stage") #AD stage (early vs. mildmoderate)
run_meta_reg(metagen_cdrsb, "antibody") #Drug
run_meta_reg(metagen_cdrsb, "anti_type") #Antibody type (human vs. humanized)
run_meta_reg(metagen_cdrsb, "mechanism") #Binding mechanism 
#run_meta_reg(metagen_cdrsb, "Biological") #Biological effect (Non vs. effective)

## bubble plots: Continuous variables#####
plot_data <- with(df2, data.frame(
  ab_centi_mean,
  yi = df2$DIFF_combined_cdr, # Change this formula if needed to reflect your specific effect size calculation
  vi = df2$DIFFSE_combined_cdr
))

# Now use 'plot_data' in ggplot
ggplot(plot_data, aes(x = ab_centi_mean, y = yi, size = 1/vi)) +
  geom_point(alpha = 0.6) +
  scale_size_area(max_size = 10) +
  geom_smooth(method = "lm", se = FALSE, color = "black", show.legend=FALSE) +
  #scale_x_continuous(breaks = seq(1, 2, by = 1)) + 
  xlab("Mean ab_centi") +
  ylab("Mean difference") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 16),  # Y axis title font size
    axis.text.y = element_text(size = 14),   # Y axis title font size
    axis.text.x = element_text(size = 14),   # Y axis title font size
  )


## Create_bubble_plot Function: Categorical variables ####
create_bubble_plot <- function(df, modifier_col, levels, yi_col, vi_col) {
  # Convert Categorical variables to numeric variables
  plot_data <- df %>%
    dplyr::mutate(modifier_numeric = dplyr::case_when(
      .data[[modifier_col]] == levels[1] ~ 1,
      .data[[modifier_col]] == levels[2] ~ 2,
      .data[[modifier_col]] == levels[3] ~ 3,
      .data[[modifier_col]] == levels[4] ~ 4,
      .data[[modifier_col]] == levels[5] ~ 5,
      .data[[modifier_col]] == levels[6] ~ 6,
      TRUE ~ NA_real_  # NA , if the level is not specified
    ),
    yi = .data[[yi_col]],
    vi = .data[[vi_col]])
  
  # ggplot
  ggplot(plot_data, aes(x = modifier_numeric, y = yi, size = 1/vi)) +
    geom_point(alpha = 0.6) +
    scale_size_area(max_size = 10) +
    geom_smooth(method = "lm", se = FALSE, color = "black", show.legend = FALSE) +  # Add regression line
    scale_x_continuous(breaks = 1:length(levels), labels = levels, minor_breaks=NULL,limits = c(0.5,2.5)) +  # X axis label
    ylim(-0.3,0.1)+
    xlab(modifier_col) + 
    ylab("Standardized mean difference") +
    theme_minimal() +
    theme(
      axis.title.y = element_text(size = 16),  # Y axis title font size
      axis.text.y = element_text(size = 14),   # Y axis title font size
      axis.text.x = element_text(size = 14),   # Y axis title font size
    )
}

# Calling the function
create_bubble_plot(df2, "AD_stage", c ("early", "mildmoderate"), "DIFF_combined_cdr", "DIFFSE_combined_cdr")
create_bubble_plot(df2, "anti_type", c("human", "humanized"), "DIFF_combined_cdr", "DIFFSE_combined_cdr")
create_bubble_plot(df2, "Biological", c("non", "effective"), "DIFF_combined_cdr", "DIFFSE_combined_cdr")
create_bubble_plot(df2, "antibody", c("1_Lecanemab", "2_Gantenerumab", "3_Donanemab", "4_Aducanumab", "5_Solanezumab", "6_Bapineuzumab"),"DIFF_combined_cdr", "DIFFSE_combined_cdr")
create_bubble_plot(df2, "mechanism", c("monomers", "oligoaggre", "aggregates"), "DIFF_combined_cdr", "DIFFSE_combined_cdr")

# Calculate 'yi' and 'vi' for the bubble plot
plot_data <- with(df2, data.frame(
  ab_centi_mean,
  yi = df2$DIFF_combined_cdr, # Change this formula if needed to reflect your specific effect size calculation
  vi = df2$DIFFSE_combined_cdr
))

# Now use 'plot_data' in ggplot####
ggplot(plot_data, aes(x = ab_centi_mean, y = yi, size = vi)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black", show.legend = FALSE) +
  xlab("Mean Age") +
  ylab("Mean difference") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 16),  # Y axis title font size
    axis.text.y = element_text(size = 14),   # Y axis title font size
    axis.text.x = element_text(size = 14),   # Y axis title font size
  )
###Combined result_ADASCOG####
#ADAS-COG
df1 <- df1 %>%
  group_by(paper) %>%
  mutate(
    Me_a_cal = ifelse(is.na(Me_a), Mc_a + diff_m_a, Me_a),
    SDe_a_cal = ifelse(!is.na(SDe_a), SDe_a, (SE_e_a * sqrt(Ne_a))), #Calculate SD from SE
    SDc_a_cal = ifelse(!is.na(SDc_a), SDc_a, (SE_c_a * sqrt(Nc_a))), #Calculate SD from SE
    SDe_a_cal= ifelse (!is.na(SDe_a_cal), SDe_a_cal, ((Me_u_a-Me_l_a)/(2^1.96))*sqrt(Ne_a)), 
    SDc_a_cal= ifelse (!is.na(SDc_a_cal), SDc_a_cal, ((Mc_u_a-Mc_l_a)/(2^1.96))*sqrt(Nc_a)),
    SDe_a_cal= ifelse (!is.na(SDe_a_cal), SDe_a_cal, sqrt(((diff_u_a-diff_l_a)/(2*1.96))^2-(SE_c_a)^2)),
    diff_m_a_cal=ifelse(!is.na(diff_m_a), diff_m_a, Me_a_cal-Mc_a),
    SE_diff_a=ifelse(is.na(diff_l_a), sqrt((SDe_a_cal^2/Ne_a)+(SDc_a_cal^2/Nc_a)), (diff_u_a-diff_l_a)/(2*1.96)),
    SD_pooled_raw=sqrt((Ne_a*SDe_a_cal^2+(Nc_a-1)*SDc_a_cal^2)/(Ne_a+Nc_a-2)),
    SMD_a_raw=diff_m_a_cal/SD_pooled_raw,
    SE_SMD_a_raw=SE_diff_a/SD_pooled_raw,
    N1a= sum(Ne_a[low1_high2_other0 == 1]),
    N2a = sum(Ne_a[low1_high2_other0 == 2]), # high dose  N
    M1a = sum(Me_a_cal[low1_high2_other0 == 1]), # low dose  Mean change
    M2a = sum(Me_a_cal[low1_high2_other0 == 2]), # high dose  Mean change
    SD1a = sum(SDe_a_cal[low1_high2_other0 == 1]), # low dose  SD
    SD2a = sum(SDe_a_cal[low1_high2_other0 == 2])
  )


df1 <- df1 %>%
  group_by(paper) %>%
  mutate(
    Ne_Combined_a2 = ifelse(low1_high2_other0==0, Ne_a, sum(Ne_a[low1_high2_other0 %in% c(1, 2)])),
    Me_Combined_a2 = ifelse(low1_high2_other0==0, Me_a_cal, (N1a * M1a + N2a * M2a) / (N1a + N2a)), 
    SDe_Combined_a2 = ifelse(low1_high2_other0==0, SDe_a_cal, sqrt(((N1a - 1) * SD1a^2 + (N2a - 1) * SD2a^2 + (N1a * N2a / (N1a + N2a)) * (M1a^2 + M2a^2 - 2 * M1a * M2a)) / (N1a + N2a - 1))), 
    DIFF_combined_a = ifelse(low1_high2_other0==0, diff_m_a_cal, Me_Combined_a2-Mc_a),
    DIFFSE_combined_a=ifelse(low1_high2_other0==0, SE_diff_a, sqrt((SDe_Combined_a2^2/Ne_Combined_a2)+(SDc_a_cal^2/Nc_a))),
    SD_pooled=sqrt((Ne_Combined_a2*SDe_Combined_a2^2+(Nc_a-1)*SDc_a_cal^2)/(Ne_Combined_a2+Nc_a-2)),
    SMD_a=DIFF_combined_a/SD_pooled,
    SE_SMD_a=DIFFSE_combined_a/SD_pooled
  )

#Exclude Low dose（currently, the same entries are included for both low dose and high dose, so select the high dose as representative）
#Note that control(placebo) groups are the same between the low and the high-dose groups
df3<-subset(df1, low1_high2_other0 %in% c(0,2))
df3$antibody <- factor(df3$antibody, levels = c( "1_Lecanemab" ,"2_Gantenerumab","3_Donanemab","4_Aducanumab","5_Solanezumab", "6_Bapineuzumab" ))

#Sensitivity analysis
#df3 <- filter(df3, antibody == "1_Lecanemab" | antibody=="2_Gantenerumab"| antibody == "3_Donanemab"|antibody=="5_Solanezumab"|antibody=="6_Bapineuzumab") #modify as needed

# Metagen_ADASCog
metagen_adascog <- metagen(
  TE = df3$SMD_a,         # Mean Difference
  seTE = df3$SE_SMD_a,    # SE of Mean Difference
  data = df3,        
  studlab = df3$Trial_3, 
  sm = "SMD",
  subgroup=df3$antibody,
  method.tau="PM",
  prediction=TRUE
)

#forrest_adascog####
tiff("forest_adascog_com.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metagen_adascog, xlim=c(-0.4,0.4), label.right="Placebo better", label.left="Antibody better", col.subgroup="black") #Forrestplot
dev.off()

#funnel####
funnel(metagen_adascog, ylim=c(0.08,0),xlim=c(-0.3,0.1),xlab="Standardized Mean Difference")
metabias(metagen_adascog, method="linreg") #Eger's test Comment subgroup in metagen when performing this

#meta_regression_ADASCOG
run_meta_reg(metagen_adascog, "age_mean_a")
run_meta_reg(metagen_adascog, "mmse_mean_a")
run_meta_reg(metagen_adascog, "apoe4_mean_a")
run_meta_reg(metagen_adascog, "ab_centi_mean_a")

#Set the levels as needed
df3$AD_stage<-factor(df3$AD_stage, levels=c("mildmoderate","early" ))
df3$antibody <- factor(df3$antibody, levels = c("3_Donanemab","4_Aducanumab","2_Gantenerumab","5_Solanezumab","6_Bapineuzumab", "1_Lecanemab"))
df3$anti_type<-factor(df3$anti_type, levels=c("humanized","human"))
df3$mechanism<-factor(df3$mechanism, levels=c("monomers","oligoaggre","aggregates" ))
df3$Biological<-factor(df3$Biological, levels=c("non","effective" ))

run_meta_reg(metagen_adascog, "AD_stage") #AD stage (early vs. mild vs. mildmoderate)
run_meta_reg(metagen_adascog, "antibody") #Drug
run_meta_reg(metagen_adascog, "anti_type") #Antibody type (human vs. humanized)
run_meta_reg(metagen_adascog, "mechanism") #Binding mechanism 
run_meta_reg(metagen_adascog, "Biological") #Biological effect (Non vs. effective)

# Calling the function Bubble plot. ADASCOG SMD
create_bubble_plot(df3, "AD_stage", c ("early","mildmoderate"), "SMD_a", "SE_SMD_a")
create_bubble_plot(df3, "anti_type", c("human", "humanized"), "SMD_a", "SE_SMD_a")
create_bubble_plot(df3, "Biological", c("non", "effective"),"SMD_a", "SE_SMD_a")
create_bubble_plot(df3, "antibody", c("1_Lecanemab", "2_Gantenerumab", "3_Donanemab", "4_Aducanumab", "5_Solanezumab", "6_Bapineuzumab"),"SMD_a", "SE_SMD_a")
create_bubble_plot(df3, "mechanism", c("monomers", "oligoaggre", "aggregates"), "SMD_a", "SE_SMD_a")


# Calculate 'yi' and 'vi' for the bubble plot
plot_data <- with(df3, data.frame(
  ab_centi_mean_a,
  yi = df3$SMD_a, # Change this formula if needed to reflect your specific effect size calculation
  vi = df3$SE_SMD_a
))

# Now use 'plot_data' in ggplot####
ggplot(plot_data, aes(x = ab_centi_mean_a, y = yi, size = 1/vi)) +
  geom_point(alpha = 0.6) +
  scale_size_area(max_size = 10) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", show.legend = FALSE) +
  xlab("Mean MMSE") +
  ylab("Standardized Mean difference") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 16),  # Y axis title font size
    axis.text.y = element_text(size = 14),   # Y axis title font size
    axis.text.x = element_text(size = 14),   # Y axis title font size
  )

#Adverse_effect combined ####
df1 <- df1 %>%
  group_by(paper) %>%
  mutate(
    Ne_Combined_h2 = ifelse(low1_high2_other0==0, Ne_ariah, sum(Ne_ariah[low1_high2_other0 %in% c(1, 2)])),
    ARIAH_Combined2= ifelse(low1_high2_other0==0, ARIA_H_e, sum(ARIA_H_e[low1_high2_other0 %in% c(1, 2)])), 
    Ne_Combined_e2 = ifelse(low1_high2_other0==0, Ne_ariae, sum(Ne_ariae[low1_high2_other0 %in% c(1, 2)])),
    ARIAE_Combined2= ifelse(low1_high2_other0==0, ARIA_E_e, sum(ARIA_E_e[low1_high2_other0 %in% c(1, 2)])), 
    Ne_Combined_death2 = ifelse(low1_high2_other0==0, Ne_death, sum(Ne_death[low1_high2_other0 %in% c(1, 2)])),
    Death_Combined2= ifelse(low1_high2_other0==0, death_e, sum(death_e[low1_high2_other0 %in% c(1, 2)])),
    Ne_Combined_Head2 = ifelse(low1_high2_other0==0, Ne_Head, sum(Ne_Head[low1_high2_other0 %in% c(1, 2)])),
    Head_Combined2= ifelse(low1_high2_other0==0, Headache_e, sum(Headache_e[low1_high2_other0 %in% c(1, 2)])), 
    Ne_Combined_Fall2 = ifelse(low1_high2_other0==0, Ne_fall, sum(Ne_fall[low1_high2_other0 %in% c(1, 2)])),
    Fall_Combined2= ifelse(low1_high2_other0==0, Fall_e, sum(Fall_e[low1_high2_other0 %in% c(1, 2)])), 
    Ne_Combined_dizzi2 = ifelse(low1_high2_other0==0, Ne_dizzi, sum(Ne_dizzi[low1_high2_other0 %in% c(1, 2)])),
    dizz_Combined2= ifelse(low1_high2_other0==0, Dizziness_e, sum(Dizziness_e[low1_high2_other0 %in% c(1, 2)])), 
    Ne_Combined_ch2 = ifelse(low1_high2_other0==0, Ne_ch, sum(Ne_ch[low1_high2_other0 %in% c(1, 2)])),
    ch_Combined2= ifelse(low1_high2_other0==0, ch_e, sum(ch_e[low1_high2_other0 %in% c(1, 2)])), 
    Ne_Combined_sae2 = ifelse(low1_high2_other0==0, Ne_SAE, sum(Ne_SAE[low1_high2_other0 %in% c(1, 2)])),
    sae_Combined2= ifelse(low1_high2_other0==0, SAE_E, sum(SAE_E[low1_high2_other0 %in% c(1, 2)])),
  ) 

#Remove Low-dose (Currently, both high-dose and low-dose contain the same thing, so choose high-dose as a representative. )
df4<-subset(df1, low1_high2_other0 %in% c(0,2))

#Sensitivity analysis
#df4 <- filter(df4, antibody == "1_Lecanemab" | antibody == "2_Gantenerumab" | antibody == "3_Donanemab")

df4$antibody<- factor(df4$antibody, levels = c( "1_Lecanemab" ,"2_Gantenerumab","3_Donanemab","4_Aducanumab","5_Solanezumab", "6_Bapineuzumab" ))
#metabin####
metabin_death=metabin(Death_Combined2, Ne_Combined_death2, death_c, Nc_death,data=df4, studlab=df4$Trial_2,method.tau="PM", subset = complete.cases(df4$Ne_Combined_death2))
metabin_sae=metabin(sae_Combined2, Ne_Combined_sae2, SAE_c, Nc_SAE,data=df4, studlab=df4$Trial_2,method.tau="PM",subset = complete.cases(df4$Ne_Combined_sae2))
metabin_ariae=metabin(ARIAE_Combined2, Ne_Combined_e2, ARIA_E_c, Nc_ariae, data=df4, studlab=df4$Trial_2,method.tau="PM", subset = complete.cases(df4$Ne_Combined_e2))
metabin_ariah=metabin(ARIAH_Combined2, Ne_Combined_h2, ARIA_H_c, Nc_ariah,data=df4, studlab=df4$Trial_2,method.tau="PM", subset = complete.cases(df4$Ne_Combined_h2))
metabin_head=metabin(Head_Combined2, Ne_Combined_Head2, Headache_c, Nc_Head,data=df4, studlab=df4$Trial_2,method.tau="PM", subset = complete.cases(df4$Ne_Combined_Head2))
metabin_fall=metabin(Fall_Combined2, Ne_Combined_Fall2, Fall_c, Nc_fall,data=df4, studlab=df4$Trial_2,method.tau="PM", subset = complete.cases(df4$Ne_Combined_Fall2))
metabin_dizzi=metabin(dizz_Combined2, Ne_Combined_dizzi2, Dizziness_c, Nc_dizzi,data=df4, studlab=df4$Trial_2,method.tau="PM", subset = complete.cases(df4$Ne_Combined_dizzi2&df4$dizz_Combined2))
metabin_ch=metabin(ch_Combined2, Ne_Combined_ch2,ch_c, Nc_ch, data=df4, studlab=df4$Trial_2, method.tau="PM", subset=complete.cases(df4$Ne_Combined_ch2))

#forestplot####
tiff("forest_death.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metabin_death,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()

tiff("forest_sae.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metabin_sae,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()

tiff("forest_ariae.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metabin_ariae,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()

tiff("forest_ariah.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metabin_ariah,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()
tiff("forest_head.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metabin_head,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()
tiff("forest_fall.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metabin_fall,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()
tiff("forest_dizzi.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metabin_dizzi,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()

tiff("forest_ch.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metabin_ch,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()

#meta_reg <- metareg(metabin21, ~ age_mean)
run_meta_reg2<-function(metabin_result,modifier){
  meta_reg_summary <- metareg(metabin_result, as.formula(paste("~", modifier)))
  meta_reg_df <- data.frame(   # Create a data frame to store the results
    Parameter = rownames(meta_reg_summary$b),
    RiskRatio = exp(meta_reg_summary$b[, 1]), # Calculating Risk Ratio
    `P value` = meta_reg_summary$pval,
    `95% CI Lower Bound` = exp(meta_reg_summary$ci.lb),
    `95% CI Upper Bound` = exp(meta_reg_summary$ci.ub)
  )
  file_name <- paste0("meta_reg_rr_summary_", modifier, ".csv")
  write.csv(meta_reg_df, file_name, row.names = FALSE) # Save the results to a CSV file
  return(meta_reg_summary) # Return the meta regression summary
}

run_meta_reg2(metabin_sae, "age_mean_SAE") #Change the metabin and the age_mean
run_meta_reg2(metabin_sae, "mmse_mean_SAE")
run_meta_reg2(metabin_sae, "apoe4_mean_SAE")
run_meta_reg2(metabin_sae, "ab_centi_mean_SAE")

#Set the levels as needed
df4$AD_stage<-factor(df4$AD_stage, levels=c("early","mildmoderate" ))
df4$antibody <- factor(df4$antibody, levels = c("2_Gantenerumab","3_Donanemab","4_Aducanumab","5_Solanezumab","1_Lecanemab","6_Bapineuzumab"))
df4$anti_type<-factor(df4$anti_type, levels=c("human","humanized"))
df4$mechanism<-factor(df4$mechanism, levels=c("monomers","oligoaggre","aggregates" ))
df4$Biological<-factor(df4$Biological, levels=c("effective","non" ))

run_meta_reg2(metabin_dizzi, "AD_stage")
run_meta_reg2(metabin_dizzi, "antibody")
run_meta_reg2(metabin_dizzi, "anti_type")
run_meta_reg2(metabin_ariae, "mechanism")
run_meta_reg2(metabin_dizzi, "Biological")

# Bubble plot Function2####
create_bubble_plot2 <- function(metabin_result, modifier) {
  meta_reg_summary <- metareg(metabin_result, as.formula(paste("~", modifier))) #Perform meta-regression
  
  logRR <- meta_reg_summary$yi.f #estimates=logRR
  vi <- meta_reg_summary$vi.f
  modifier_numeric <- as.numeric(factor(meta_reg_summary$data[[modifier]]))
  
  meta_reg_df <- data.frame(             #Generate the data frame
    logRR = logRR,                                       # Effect size
    vi = vi,                                             # variance
    modifier = factor(meta_reg_summary$data[[modifier]]),  # Modifier
    RiskRatio = exp(logRR)                             # Change log RR to Risk Ratio
  )
  plot <- ggplot(meta_reg_df, aes(x = modifier_numeric, y = RiskRatio, size = 1/vi)) + #Bubble plot
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "black",show.legend = FALSE) +  # Liner regression line. Delete if you don't need.
    scale_x_continuous(breaks = 1:length(levels(factor(meta_reg_summary$data[[modifier]]))),
                       labels = levels(factor(meta_reg_summary$data[[modifier]])),
                       limits = c(0.5,3.5),#Change the limit depending on the number of categories
                       minor_breaks = NULL) +  
    scale_y_continuous(limits = c(0.8, NA)) +
    scale_size_area(max_size = 15) +  # Size of the bubbles
    xlab(modifier) +
    ylab("Risk Ratio") +
    theme_minimal() +
    theme(
      axis.title.y = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1)
    )
  print(plot)
}

create_bubble_plot2(metabin_ariae, "AD_stage")
create_bubble_plot2(metabin_ariae, "anti_type")
create_bubble_plot2(metabin_ariae, "Biological")

create_bubble_plot2(metabin_ariae, "antibody")
create_bubble_plot2(metabin_ariae, "mechanism")

#meta-reg for continuous modifier ####
create_bubble_plot3 <- function(metabin_result, modifier) {
  meta_reg_summary <- metareg(metabin_result, as.formula(paste("~", modifier))) 
  logRR <- meta_reg_summary$yi.f  # logRR
  vi <- meta_reg_summary$vi.f  # variance
  modifier_values <- meta_reg_summary$data[[modifier]]
  meta_reg_df <- data.frame(
    logRR = logRR,  #effect size（logRR）
    vi = vi,  # variation
    modifier = modifier_values,  
    RiskRatio = exp(logRR)  # convert from log RR to RR
    )
  
  # Create bubble plot
  plot <- ggplot(meta_reg_df, aes(x = modifier, y = RiskRatio, size = 1 / vi)) + 
    geom_point(alpha = 0.6) + 
    geom_smooth(method = "lm", se = FALSE, color = "black", show.legend = FALSE) +  #linear regression
    scale_y_continuous(limits = c(0.8, NA)) +  # Y axis acale, change if you need
    scale_size_area(max_size = 15) +  
    xlab(modifier) +  # X axis label
    ylab("Risk Ratio") +  # Y axis label
    theme_minimal() + 
    theme(
      axis.title.y = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.text.x = element_text(size = 14)
    )
  
  print(plot)  #Print the plot
}

create_bubble_plot3(metabin_ariae, "age_mean_ariae")
create_bubble_plot3(metabin_ariae, "mmse_mean_ariae")
create_bubble_plot3(metabin_ariae, "apoe4_mean_ariae")
create_bubble_plot3(metabin_ariae, "ab_centi_mean_ariae")

#funnel categorical outcomes
par(cex.axis = 1.5,  # Size of the labels in x and y axis
    cex.lab = 1.5)   
funnel(metabin_death)
metabias(metabin_death, method="linreg")

#Sensitivity analysis low-dose only ####
#df5: exclude high-dose　(Only low dose)
df5<-subset(df1, low1_high2_other0 %in% c(0,2)) #change to c(0,2) when you perform analysis for "only high dose group".
df5$antibody<- factor(df5$antibody, levels = c( "1_Lecanemab" ,"2_Gantenerumab","3_Donanemab","4_Aducanumab","5_Solanezumab", "6_Bapineuzumab" ))
# Metagen_CDRSB
metagen_cdrsb <- metagen(
  TE = df5$diff_m_c,         # Mean Difference
  seTE = df5$SE_diff_cdr,    # SE of Mean Difference
  data = df5,        
  studlab = df5$Trial_3, 
  sm = "MD",
  subgroup=df5$antibody,
  method.tau="PM",
  prediction=TRUE
)
#forest cdrsb####
tiff("forest_cdrsb_high.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metagen_cdrsb, xlim=c(-1.0,1.0), label.right="Placebo better", label.left="Antibody better", col.subgroup="black") #Forrestplot
dev.off()

# Metagen_ADASCog
metagen_adascog <- metagen(
  TE = df5$SMD_a_raw,         # Mean Difference
  seTE = df5$SE_SMD_a_raw,    # SE of Mean Difference
  data = df5,        
  studlab = df5$Trial_3, 
  sm = "SMD",
  subgroup=df5$antibody,
  method.tau="PM",
  prediction=TRUE
)
#forrest_adascog####
tiff("forest_adascog_high.tiff",width = 10, height = 10, units = 'in', res = 300)
forest(metagen_adascog, xlim=c(-0.4,0.4), label.right="Placebo better", label.left="Antibody better", col.subgroup="black") #Forrestplot
dev.off()

#adverse_event####
metabin_death=metabin(death_e, Ne_death, death_c, Nc_death,data=df5, studlab=df5$Trial_2,method.tau="PM", subset = complete.cases(df5$Ne_death))
metabin_sae=metabin(SAE_E, Ne_SAE, SAE_c, Nc_SAE,data=df5, studlab=df5$Trial_2,method.tau="PM",subset = complete.cases(df5$Ne_SAE))
metabin_ariae=metabin(ARIA_E_e, Ne_ariae, ARIA_E_c, Nc_ariae, data=df5, studlab=df5$Trial_2,method.tau="PM", subset = complete.cases(df5$Ne_ariae))
metabin_ariah=metabin(ARIA_H_e, Ne_ariah, ARIA_H_c, Nc_ariah,data=df5, studlab=df5$Trial_2,method.tau="PM", subset = complete.cases(df5$Ne_ariah))
metabin_head=metabin(Headache_e, Ne_Head, Headache_c, Nc_Head,data=df5, studlab=df5$Trial_2,method.tau="PM", subset = complete.cases(df5$Ne_Head))
metabin_fall=metabin(Fall_e, Ne_fall, Fall_c, Nc_fall,data=df5, studlab=df5$Trial_2,method.tau="PM", subset = complete.cases(df5$Ne_fall))
metabin_dizzi=metabin(Dizziness_e, Ne_dizzi, Dizziness_c, Nc_dizzi,data=df5, studlab=df5$Trial_2,method.tau="PM", subset = complete.cases(df5$Ne_dizzi&df5$Dizziness_e))
metabin_ch=metabin(ch_e, Ne_ch, ch_c, Nc_ch, data=df5, studlab=df5$Trial_2, method.tau="PM",subset=complete.cases(df5$Ne_ch))

tiff("forest_death.tiff",width = 10, height = 5, units = 'in', res = 300)
forest(metabin_death,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()

tiff("forest_sae.tiff",width = 10, height = 5, units = 'in', res = 300)
forest(metabin_sae,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()

tiff("forest_ariae.tiff",width = 10, height = 5, units = 'in', res = 300)
forest(metabin_ariae,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()

tiff("forest_ariah.tiff",width = 10, height = 5, units = 'in', res = 300)
forest(metabin_ariah,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()
tiff("forest_head.tiff",width = 10, height = 5, units = 'in', res = 300)
forest(metabin_head,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()
tiff("forest_fall.tiff",width = 10, height = 5, units = 'in', res = 300)
forest(metabin_fall,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()
tiff("forest_dizzi.tiff",width = 10, height = 5, units = 'in', res = 300)
forest(metabin_dizzi,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()
tiff("forest_ch.tiff",width = 10, height = 5, units = 'in', res = 300)
forest(metabin_ch,  prediction=TRUE, label.right="Favors Placebo", label.left="Favors Antibody",col.subgroup="black")
dev.off()

#NNH Categorical outcome (secondary outcome)####
#Calculate Risk Difference
metabin_ariae2 =metabin(ARIAE_Combined2, Ne_Combined_e2, ARIA_E_c, Nc_ariae, data = df4, studlab = df4$Trial_2, method.tau = "PM", sm = "RD", subset = complete.cases(df4$Ne_Combined_e2))
metabin_ariah2=metabin(ARIAH_Combined2, Ne_Combined_h2, ARIA_H_c, Nc_ariah,data=df4, studlab=df4$Trial_2,method.tau="PM", sm = "RD",subset = complete.cases(df4$Ne_Combined_h2))
metabin_head2=metabin(Head_Combined2, Ne_Combined_Head2, Headache_c, Nc_Head,data=df4, studlab=df4$Trial_2,method.tau="PM", sm = "RD",subset = complete.cases(df4$Ne_Combined_Head2))

summary(metabin_head2)

#Calculate NNH
  NNH <- 1 / abs(metabin_head2$TE.random) #Calculation of Number need to harm (NNH)
  NNH_ariae <- 1 / abs(metabin_ariae2$TE.random) #Calculation of Number need to harm (NNH)
  NNH_ariah <- 1 / abs(metabin_ariah2$TE.random) #Calculation of Number need to harm (NNH)
  print(paste("NNH: ", round(NNH, 2)))
  print(paste("NNH: ", round(NNH_ariae, 2)))
  print(paste("NNH: ", round(NNH_ariah, 2)))

#NNT Continuous outcome (primary outcome)####
#CDR-SB
df2<-filter(df2, antibody!="6_Bapineuzumab"& antibody!="5_Solanezumab"& antibody!="4_Aducanumab")
  df3<-filter(df3, antibody!="6_Bapineuzumab"& antibody!="5_Solanezumab"& antibody!="4_Aducanumab")
  df4<-filter(df4, antibody!="6_Bapineuzumab"& antibody!="5_Solanezumab"& antibody!="4_Aducanumab")
  group1 <- df2$Me_Combined2_cdr #Mean change of the experimental group (combined)　#for ADAS-Cog, it is df3$Me_Combined_a2
  group2 <- df2$Mc_c #Mean change of the placebo group #for ADAS-Cog, it is df3$Mc_a
  
#ADAS-Cog
  
  # Normalize the ADAS-Cog scores based on the version used
  # Example: Assuming df3 has a column 'ADAS_version' that indicates which version is used
  
  df3$Me_Combined_a2_normalized <- ifelse(df3$adascog_ver == 11, df3$Me_Combined_a2/ 70,
                                          ifelse(df3$adascog_ver == 13, df3$Me_Combined_a2 / 85,
                                                 df3$Me_Combined_a2 / 90))
  
  df3$Mc_a_normalized <- ifelse(df3$adascog_ver == 11, df3$Mc_a / 70,
                                ifelse(df3$adascog_ver == 13, df3$Mc_a / 85,
                                       df3$Mc_a / 90))
  
  # Use the normalized scores in your NNT analysis
  group1 <- df3$Me_Combined_a2_normalized
  group2 <- df3$Mc_a_normalized
  
  
  combined_data <- c(group1, group2)
  
  # Rank the combined data
  ranks <- rank(combined_data)
  
  # Separate the ranks for group1 and group2
  ranks_group1 <- ranks[1:length(group1)]
  ranks_group2 <- ranks[(length(group1) + 1):length(combined_data)]
  
  # Calculate the sum of ranks for group1 (which corresponds to W)
  W_group1 <- sum(ranks_group1)
  # Calculate U value
  N1 <- length(group1)
  N2 <- length(group2)
  U1 <- W_group1 - (N1 * (N1 + 1)) / 2  # U-value for group 1
  U2 <- N1 * N2 - U1             # U-value for group 2
  # Take the smaller of U1 and U2
  U_value <- min(U1, U2)

AUC <- U_value / (N1 * N2)
RD<-(2*AUC-1)
Q1 <- AUC / (2 - AUC)
Q2 <- 2 * AUC^2 / (1 + AUC)
SE_AUC <- sqrt((AUC * (1 - AUC) + (N1 - 1) * (Q1 - AUC^2) + (N2 - 1) * (Q2 - AUC^2)) / (N1 * N2))
SE_RD <- 2 * SE_AUC
# Create a dataframe with RD and SE_RD for meta-analysis
df <- data.frame(RD = RD, SE_RD = SE_RD)
meta_rd<-metagen(TE=df$RD, 
                 seTE=df$SE_RD, 
                 sm="RD",
                 method.tau="PM")
summary(meta_rd)
# Extract the pooled Risk Difference (RD) from the meta-analysis
pooled_RD <- meta_rd$TE.random  # This is the pooled RD from the random-effects model

# Calculate NNT using the pooled RD
NNT <- 1 / abs(pooled_RD)

# Print the NNT
cat("Number Needed to Treat (NNT):", NNT, "\n")

