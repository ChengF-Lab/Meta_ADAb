library (meta)
library(dplyr)
library(metafor)
sessionInfo()

R.version.string
getwd()
setwd("/Users/workingdirectory")

#Data import####
df1 <- read.csv("meta_ab_antibody_data.csv", header = TRUE)

#CDR-SB####
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

#Combine low-dose and high-dose: CDR-SB
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
df2$antibody <- factor(df2$antibody, levels = c( "1_Lecanemab","2_Gantenerumab", "3_Donanemab", "4_Aducanumab", "5_Solanezumab", "6_Bapineuzumab" ))

# Case study 1####
## Meta_Metagen_CDRSB####
metagen_cdrsb <- metagen(
  TE = df2$DIFF_combined_cdr,         # Mean Difference
  seTE = df2$DIFFSE_combined_cdr,    # SE of Mean Difference
  data = df2,        
  studlab = df2$Trial_3, 
  sm = "MD",
  subgroup=df2$antibody, #The results don't change with or without defining subgroups
  method.tau="PM",
  prediction=TRUE
)
summary(metagen_cdrsb)

## Metafor_CDRSB#### 
res <- rma(
  yi = df2$DIFF_combined_cdr,                   
  sei = df2$DIFFSE_combined_cdr,               
  method = "PM",           
  data = df2,               
  slab = df2$Trial_3          
)
summary(res)

#Case study 2####
##Meta####
run_meta_reg<-function(metagen_result,modifier){
  meta_reg_result <- metareg(metagen_result, as.formula(paste("~", modifier)))
  return(summary(meta_reg_result))
}

df2$antibody <- factor(df2$antibody, levels = c("1_Lecanemab","2_Gantenerumab","3_Donanemab","4_Aducanumab","5_Solanezumab","6_Bapineuzumab"))
run_meta_reg(metagen_cdrsb, "antibody") #Drug

## Metafor_CDRSB#### 
res <- rma(
  yi = df2$DIFF_combined_cdr,                   
  sei = df2$DIFFSE_combined_cdr,               
  mods = ~ antibody,     #Modifier antibody    
  method = "PM",  
  data = df2,               
  slab = df2$Trial_3          
)
summary(res)

#Case study 3 ####
#Adverse_effect combined 
df1 <- df1 %>%
  group_by(paper) %>%
  mutate(
    Ne_Combined_e2 = ifelse(low1_high2_other0==0, Ne_ariae, sum(Ne_ariae[low1_high2_other0 %in% c(1, 2)])),
    ARIAE_Combined2= ifelse(low1_high2_other0==0, ARIA_E_e, sum(ARIA_E_e[low1_high2_other0 %in% c(1, 2)])), 
  ) 

#Remove Low-dose (Currently, both high-dose and low-dose contain the same thing, so choose high-dose as a representative. )
df4<-subset(df1, low1_high2_other0 %in% c(0,2))
df4$ARIA_E_c <- as.numeric(df4$ARIA_E_c)
df4_subset <- df4[complete.cases(df4$Ne_Combined_e2), ]
##Meta####
metabin_ariae=metabin(ARIAE_Combined2, Ne_Combined_e2, ARIA_E_c, Nc_ariae, data=df4_subset, studlab=df4_subset$Trial_2,method.tau="PM")
summary(metabin_ariae)

##Metafor####
res <- rma(
  ai = ARIAE_Combined2,         
  n1i = Ne_Combined_e2,         
  ci = ARIA_E_c,               
  n2i = Nc_ariae,              
  measure = "RR",               
  method = "PM",              
  data = df4_subset,                   
  slab = Trial_2               
)
summary(res)

#Case study 4####
##Meta####
run_meta_reg2<-function(metabin_result,modifier){
  meta_reg_summary <- metareg(metabin_result, as.formula(paste("~", modifier)))
  meta_reg_df <- data.frame(   # Create a data frame to store the results
    Parameter = rownames(meta_reg_summary$b),
    RiskRatio = exp(meta_reg_summary$b[, 1]), # Calculating Risk Ratio
    `P value` = meta_reg_summary$pval,
    `95% CI Lower Bound` = exp(meta_reg_summary$ci.lb),
    `95% CI Upper Bound` = exp(meta_reg_summary$ci.ub)
  )
  return(meta_reg_summary) # Return the meta regression summary
}
df4$antibody <- factor(df4$antibody, levels = c("2_Gantenerumab","3_Donanemab","4_Aducanumab","5_Solanezumab","1_Lecanemab","6_Bapineuzumab"))
run_meta_reg2(metabin_ariae, "antibody")

##Metafor####
res <- rma(
  ai = ARIAE_Combined2,         
  n1i = Ne_Combined_e2,        
  ci = ARIA_E_c,                
  n2i = Nc_ariae,              
  measure = "RR",               
  method = "PM",                
  mods = ~ antibody, #modifier: antibody
  data = df4_subset,            
  slab = Trial_2             
)
summary(res)

