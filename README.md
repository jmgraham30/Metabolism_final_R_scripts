# Title of Dataset: Pervasive effects of Wolbachia on host body mass and metabolic rate
---

Dataset includes data files for estimates of Drosophila host body mass, metabolic rate, and locomotor activity.

## Description of the data and file structure

metabolism_summaryData.csv: File includes measures of host body mass and metabolic rate corresponding to Figures 2, 3, and all Supplemental Figures and Tables. Each row corresponds to a measurement of standard metabolic rate (SMR) from an individual fly.
Date: Date of MAVEn run
Time: Start time of SMR measurement
Chamber: Chamber (1-16) on MAVEn that each fly was randomly assigned to
runNum: Sequential count of SMR recordings for the day
Replicate: Replicate of the SMR recording for an individual fly
SMRstart.sec: Start time of SMR measurement converted to seconds
SMRend.sec: End time of SMR measurement converted to seconds
aveSMR: Estimate of SMR (nmol min-1)
aveCO2ppm: Average CO2 production during SMR recording (ppm)
aveAjdCO2ppm: Average CO2 production during SMR recording (ppm) after correcting for drift
aveFRC_mlmin: Average flow rate during SMR recording (ml min-1)
aveWVppt: Average water vapor during SMR recording (ppt)
aveTempC: Average temperature during SMR recording (C)
aveRH_pct: Average relative humidity during SMR recording (%)
aveLight_lux: Average light intensity during SMR recording (lux)
aveActivity: Average fly locomotor activity during SMR recording (ADS)
flyID: Unique individual identifier for each fly
Genotype: Genotype of fly
Host: Species of Drosophila fly
Infected: Infection status of individual fly (U = uninfected, I = infected)
Sex: Sex of fly
Age_days: Age of adult fly in days
Weight_mg: Wet weight of fly (mg)
Notes: Notes from run

activity_summaryData.csv: File includes measures of host locomotor activity over a 3-hour period corresponding to Supplemental Figure S5 (Hague et al., 2021). Data obtained from https://doi.org/10.5061/dryad.6t1g1jwxv. See linked Dryad site for additional details.

phylogram_wb.tre: Bayesian phylogram .tre file of A- and B-group Wolbachia corresponding to Figure 1. 

phylogram_wb_wNo.tre: Bayesian phylogram .tre file of A- and B-group Wolbachia excluding the B-group strain wNo.

phylogram_2_host.tre: Bayesian phylogram .tre file of Drosophila host species corresponding to Figure 2.