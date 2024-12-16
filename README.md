# **Systematic review and meta-analysis for passive immunotherapies targeting amyloid beta in Alzheimer's disease**

## This repository contains the data and R script for the meta-analysis of passive immunotheapies targeting amyloid beta in Alzheimer's disease.

- meta_ab_antibody_search_results.xlsx   - Contains the search results from pubmed, embase and clinical trials.gov, and reasons for inclusion/exclusion for the meta-anslysis.
- meta_ab_antibody_data.csv  - Includes trial-level data collected from each trial
- meta_ab_antibody_definition.xlsx -Includes definitions for each variable in the "meta_ab_antibody_data.csv". 
- meta_ab_antibody.R   - R script used in the study. This script can be run on the meta_ab_antibody_data.csv dataset.
- s_meta_ab_antibody_data_additional.csv -Includes trial-level data collected from trials, including those halted due to futility that ended with fewer than 200 patients in each arm.
- s_meta_ab_antibody_TextF.R. - R script used in the comparason between meta and metafor packages (TextF in S1 Supporting information)