library(caper)
library(phytools)
library(ape)
library(geiger)
library(egg)
library(emmeans)
library(lme4)
library(Hmisc)
library(ggtree)
library(performance)
library(tidyverse)
theme_set(theme_bw())
options(rgl.useNULL = TRUE)
library(geomorph)
library(gt)

#---------Load data---------------------------------------------
#load wolbachia phylogeny
wb.tree <- read.nexus("data/phylogram_wb.tre")
wb_wNo.tree <- read.nexus("data/phylogram_wb_wNo.tre")

# Load file with summary data
hague_df <- read.csv(file = "data/metabolism_summaryData.csv")

#load discrete weight data
wb.weight <- tribble(
                    ~RowName,          ~state_change,
                    "Mel_CS_berk",      "yes",    
                    "mel_FFD25",        "no",   
                    "mel_pc75",         "yes",    
                    "wYak",             "yes",    
                    "wSan_Cooper2019",  "yes",   
                    "wTei",             "no",    
                    "wRi",              "yes",   
                    "wSuz_Soizios2013", "no",   
                    "wAura_Turelli2018","no",    
                    "wHa",              "yes",    
                    "LD15",             "yes",    
                    "wMau",             "yes",   
                  ) |>
                    column_to_rownames("RowName")

#load discrete metabolism data
wb.metabolism <- tribble(
                  ~RowName,          ~noWeight, ~Weight,
                  "Mel_CS_berk",      "yes",    "no",
                  "mel_FFD25",        "yes",    "yes",
                  "mel_pc75",         "yes",    "no",
                  "wYak",             "no",     "no",
                  "wSan_Cooper2019",  "yes",    "yes",
                  "wTei",             "no",     "no",
                  "wRi",              "yes",    "no",
                  "wSuz_Soizios2013", "yes",    "no",
                  "wAura_Turelli2018","no",     "no",
                  "wHa",              "yes",    "yes",
                  "LD15",             "no",     "no",
                  "wMau",             "no",     "no"
                ) |>
                  column_to_rownames("RowName")

#list wb variants to include in phylogeny (remove wNo because it's co-infected)
wb.included <- c("Mel_CS_berk", "mel_FFD25", "mel_pc75", "wYak", "wSan_Cooper2019", "wTei",
                 "wRi", "wSuz_Soizios2013", "wAura_Turelli2018", "wHa", "LD15", "wMau")

#Plot tree with wNo with ggtree
p1 <- ggtree(wb_wNo.tree) +
  geom_tiplab(align=TRUE) +
  scale_x_continuous(expand=expansion(0.5)) + # make more room for the labels
  geom_treescale() 
p1

p1 <- flip(p1, 4, 19)
p1 <- flip(p1, 3, 2)
p1 <- flip(p1, 10, 11)
p1

ggsave("output/phylogram_wb_wNo.pdf", plot = p1, width = 15, height = 6, dpi=300, useDingbats=FALSE)

#enter number of permutations for D simulations
D.perms <- 1000

#-----------Test for phylogenetic signal of weight data----------------------------------
#-----------Phylogenetic signal of discrete traits using Fritz and Purvis' D-------------
#https://duke.atlassian.net/wiki/spaces/AnthroTree/pages/15926294/5.5.1+Fritz+and+Purvis+D+in+R
#D is expected to be close to 1 if traits are random with regard to phylogeny 
#      (i.e. Probability of E(D) resulting from no (random) phylogenetic structure, is D significantly different from 1)
#close to 0 if binary traits are determined by an underlying continuous trait following a Brownian motion random walk model of evolution 
#      (i.e. Probability of E(D) resulting from Brownian phylogenetic structure, is D significantly different from 0) 
#and <0 if traits are highly structured in the phylogeny.

wb.weight$strains <- row.names(wb.weight)
wb.comparative <- comparative.data(phy = wb.tree, data = wb.weight, names.col = strains, 
                                   vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

#D is 1 if the distribution of the binary trait is random with respect to phylogeny, 
#and greater than 1 if the distribution of the trait is more overdispersed than the random expectation.
set.seed(72)
D.result <- phylo.d(data=wb.comparative, binvar = state_change, permut = 1000)
D.result
plot(D.result)

p2 <- ggplot(as.data.frame(D.result$Permutations), aes(x=random)) +
  geom_density(color = "blue") +
  geom_density(aes(x=brownian), color= "red") +
  geom_vline(xintercept = D.result$Parameters$MeanRandom, color="blue") +
  geom_vline(xintercept = D.result$Parameters$MeanBrownian, color="red") +
  geom_vline(xintercept = D.result$Parameters$Observed, linetype="dotted") +
  scale_x_continuous(sec.axis = sec_axis(~(. - mean(D.result$Permutations$brownian))/(mean(D.result$Permutations$random) - mean(D.result$Permutations$brownian)),
                                         name = "D", breaks = c(0, as.numeric(D.result$DEstimate), 1))) +
  xlab("Sum of charecter change") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2

#-----------Test for phylogenetic signal of metabolic data-------------------------------
#-----------Phylogenetic signal of discrete traits using Fritz and Purvis' D-------------

wb.metabolism$strains <- row.names(wb.metabolism)
wb.comparative <- comparative.data(phy = wb.tree, data = wb.metabolism, names.col = strains, 
                                   vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

# Calculate D statistic for SMR analysis with no weight
D.result <- phylo.d(data=wb.comparative, binvar = noWeight, permut = 1000)
D.result
plot(D.result)

p3 <- ggplot(as.data.frame(D.result$Permutations), aes(x=random)) +
  geom_density(color = "blue") +
  geom_density(aes(x=brownian), color= "red") +
  geom_vline(xintercept = D.result$Parameters$MeanRandom, color="blue") +
  geom_vline(xintercept = D.result$Parameters$MeanBrownian, color="red") +
  geom_vline(xintercept = D.result$Parameters$Observed, linetype="dotted") +
  scale_x_continuous(sec.axis = sec_axis(~(. - mean(D.result$Permutations$brownian))/(mean(D.result$Permutations$random) - mean(D.result$Permutations$brownian)),
                                         name = "D", breaks = c(0, as.numeric(D.result$DEstimate), 1))) +
  xlab("Sum of charecter change") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3

# Calculate D statistic for SMR analysis with adjusted weight
D.result <- phylo.d(data=wb.comparative, binvar = Weight, permut = 1000)
D.result
plot(D.result)

p4 <- ggplot(as.data.frame(D.result$Permutations), aes(x=random)) +
  geom_density(color = "blue") +
  geom_density(aes(x=brownian), color= "red") +
  geom_vline(xintercept = D.result$Parameters$MeanRandom, color="blue") +
  geom_vline(xintercept = D.result$Parameters$MeanBrownian, color="red") +
  geom_vline(xintercept = D.result$Parameters$Observed, linetype="dotted") +
  scale_x_continuous(sec.axis = sec_axis(~(. - mean(D.result$Permutations$brownian))/(mean(D.result$Permutations$random) - mean(D.result$Permutations$brownian)),
                                         name = "D", breaks = c(0, as.numeric(D.result$DEstimate), 1))) +
  xlab("Sum of charecter change") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4


#----------------------Visualize data-------------------------
#-----------Build tree-------------------------------------
#Plot tree with ggtree
p5 <- ggtree(wb.tree) +
  geom_tiplab(align=TRUE) +
  scale_x_continuous(expand=expansion(0.5)) + # make more room for the labels
  geom_treescale()
p5

p5 <- flip(p5, 9, 23)
p5 <- flip(p5, 6, 8)
p5 <- flip(p5, 1, 2)
p5

#---------Clean up weight data and calculate means---------------------
# filter out rows with missing flyID or Weight_mg, and rows with Notes
hague_df <- hague_df |>
  filter(!is.na(flyID) & !is.na(Weight_mg)) |>
  filter(is.na(Notes)) |>
  mutate(log_aveSMR = log10(aveSMR)) |>
  mutate(log_Weight = log10(Weight_mg))

hague_df <- hague_df |>
  mutate(Infected = relevel(factor(Infected),ref="U"),
         SMRstart.sec = SMRstart.sec/10000)

# reorder Genotype levels
hague_df <- hague_df |>
  mutate(Genotype = factor(Genotype,levels=c("CSBerk", "FFD25", "PC75", "yakB13L5", "san", "teiB13L11", 
                                             "R84", "suz", "auraL2", "Car5", "sech", "mau31", "sim198")))

# create a reduced dataset with mean Weight_mg per flyID
hague_df_weight_reduced <- hague_df |>
  dplyr::group_by(flyID) |>
  dplyr::summarize(Weight_mg=mean(Weight_mg))

glimpse(hague_df_weight_reduced)

# extract Infected and Sex per ID
hague_df_id_info <- hague_df |>
  select(flyID,Infected,Sex,Genotype,Host) |>
  distinct()

glimpse(hague_df_id_info)

# merge mean weight info with hague_df_id_info
hague_df_weight_reduced <- hague_df_weight_reduced |>
  left_join(hague_df_id_info,by="flyID") 

glimpse(hague_df_weight_reduced)

hague_df_weight_reduced |>
  group_by(Genotype) |>
  summarise(n=n())

nrow(hague_df_weight_reduced) #total number of flies assayed

weight_means_melt <- hague_df_weight_reduced %>%
  group_by(Genotype, Sex, Infected) %>%
  summarise(as_tibble(as.list(smean.cl.normal(Weight_mg, conf.int=0.95))))

#--------------summarize all weight data and generate figure--------------------------------------
## Plot mean data
dodge <- position_dodge(width = .5)

weight_means_melt <- weight_means_melt %>%
  mutate(group = interaction(Sex, Infected, sep = "_"),
         group = factor(group, levels = c("F_U","F_I","M_U","M_I")))

p7 <- ggplot(weight_means_melt, aes(Genotype, y=Mean, group=group)) +
  geom_errorbar(aes(ymin=Lower,
                    ymax=Upper),
                width=0.3, alpha=1,
                position = dodge) +
  geom_point(aes(fill = Infected, shape = Sex), colour = "black", size = 2, position = dodge) +
  scale_fill_manual(values = c("U" = "gray80", "I" = "brown1")) +
  scale_shape_manual(values = c("F" = 21, "M" = 22)) +

  scale_y_continuous(breaks=c(0.5,1,1.5)) +
  ylab("mean weight (mg)")
p7

#-------Calculate LS mean data for SMR WITHOUT weight----------------------------------
fit_model_by_species <- function(Genotype_val,data_df = hague_df){
  
  df <- data_df |>
    filter(Genotype == Genotype_val)
  
  lmer_fit <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + 
                     aveTempC + aveLight_Lux + SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                     (1 | flyID), 
                     data=df)
  
  emm_sex_infected <- emmeans(lmer_fit,specs = pairwise ~ Sex:Infected) 
  emm_sex_infected_df <- as.data.frame(emm_sex_infected$emmeans)
  contrasts_sex_infected_df <- as.data.frame(emm_sex_infected$contrasts)
  
  return(tibble(Genotype = Genotype_val,
                fitted_model = list(lmer_fit),
                model_summary = list(broom.mixed::tidy(lmer_fit)),
                emm_sex_infected = list(emm_sex_infected_df),
                contrasts_sex_infected = list(contrasts_sex_infected_df)))
  
}

# apply the function to each unique Genotype and combine results
metabolic_activity_models <- map(unique(hague_df$Genotype),
                                 fit_model_by_species) |>
  bind_rows()

glimpse(metabolic_activity_models)

### LS mean tables for models without weight
lsmean_table <- metabolic_activity_models |>
  select(Genotype, emm_sex_infected) |>
  unnest(emm_sex_infected) |>
  
  # back-transform LS means and confidence intervals
  mutate(emmean = 10^(emmean),
         SE = 10^(SE),
         lower.CL = 10^(lower.CL),
         upper.CL = 10^(upper.CL)) |>
  
  # Create a combined factor for columns
  unite(
    col = "sex_infection",
    Sex, Infected,
    sep = "_"
  ) |>
  
  # Keep only LS means
  select(Genotype, sex_infection, emmean, lower.CL, upper.CL, SE) |>
  
  pivot_wider(
    names_from  = sex_infection,
    values_from = c(emmean, lower.CL, upper.CL, SE)
  )

lsmean_table$No_weight_model <- "No_weight"

lsmean_table <- lsmean_table |>
  mutate(Genotype = factor(Genotype, levels=c("sim198", "mau31", "sech", "Car5","auraL2", "suz", "R84", "teiB13L11", "san", "yakB13L5", "PC75", "FFD25", "CSBerk"))) |>
  arrange(desc(Genotype))

glimpse(lsmean_table)

#-----summarize all SMR WITHOUT weight and generate figure------------------------------------------------------
lsmean_table

lsmean_long <- lsmean_table |>
  pivot_longer(
    cols           = -c(Genotype, No_weight_model),
    names_to       = c(".value", "group"),
    names_pattern  = "(emmean|lower\\.CL|upper\\.CL|SE)_([FM]_[UI])"
  ) |>
  separate(group, into = c("Sex", "Infected"), sep = "_")

lsmean_long$Genotype <- factor(lsmean_long$Genotype,levels=c("CSBerk", "FFD25", "PC75", "yakB13L5", "san", "teiB13L11", 
                                                  "R84", "suz", "auraL2", "Car5", "sech", "mau31", "sim198"))
lsmean_long$Sex <- factor(lsmean_long$Sex, levels = c("F", "M"))

lsmean_long <- lsmean_long %>%
  mutate(group = interaction(Sex, Infected, sep = "_"),
         group = factor(group, levels = c("F_U", "F_I", "M_U", "M_I")))

p11 <- ggplot(lsmean_long, aes(Genotype, y=emmean, group=group)) +
  geom_errorbar(aes(ymin=lower.CL,
                    ymax=upper.CL),
                width=0.3, alpha=1,
                position = dodge) +
  geom_point(aes(fill = Infected, shape = Sex), colour = "black", size = 2, position = dodge) +
  scale_fill_manual(values = c("U" = "gray80", "I" = "brown1")) +
  scale_shape_manual(values = c("F" = 21, "M" = 22)) +
  guides(color = guide_legend(reverse=TRUE)) +
  ylab("LS mean SMR (nmol min-1)")
p11

#-------Calculate LS mean data for SMR WITH weight----------------------------------
fit_model_by_species_weight <- function(Genotype_val,data_df = hague_df){
  
  df <- data_df |>
    filter(Genotype == Genotype_val)
  
  lmer_fit <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + log_Weight +
                     aveTempC + aveLight_Lux + SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                     (1 | flyID), 
                   data=df)
  
  emm_sex_infected <- emmeans(lmer_fit,specs = pairwise ~ Sex:Infected) 
  emm_sex_infected_df <- as.data.frame(emm_sex_infected$emmeans)
  contrasts_sex_infected_df <- as.data.frame(emm_sex_infected$contrasts)
  
  return(tibble(Genotype = Genotype_val,
                fitted_model = list(lmer_fit),
                model_summary = list(broom.mixed::tidy(lmer_fit)),
                emm_sex_infected = list(emm_sex_infected_df),
                contrasts_sex_infected = list(contrasts_sex_infected_df)))
  
}

# apply the function to each unique Genotype and combine results
metabolic_activity_models_weight <- map(unique(hague_df$Genotype),
                                 fit_model_by_species_weight) |>
  bind_rows()

glimpse(metabolic_activity_models_weight)

### LS mean tables for models WITH weight
lsmean_table_weight <- metabolic_activity_models_weight |>
  select(Genotype, emm_sex_infected) |>
  unnest(emm_sex_infected) |>
  
  # back-transform LS means and confidence intervals
  mutate(emmean = 10^(emmean),
         SE = 10^(SE),
         lower.CL = 10^(lower.CL),
         upper.CL = 10^(upper.CL)) |>
  
  # Create a combined factor for columns
  unite(
    col = "sex_infection",
    Sex, Infected,
    sep = "_"
  ) |>
  
  # Keep only LS means
  select(Genotype, sex_infection, emmean, lower.CL, upper.CL, SE) |>
  
  pivot_wider(
    names_from  = sex_infection,
    values_from = c(emmean, lower.CL, upper.CL, SE)
  )

lsmean_table_weight$weight_model <- "weight"

lsmean_table_weight <- lsmean_table_weight |>
  mutate(Genotype = factor(Genotype, levels=c("sim198", "mau31",  "sech", "Car5","auraL2", "suz", "R84", "teiB13L11", "san", "yakB13L5", "PC75", "FFD25", "CSBerk"))) |>
  arrange(desc(Genotype))

glimpse(lsmean_table_weight)

#-----summarize all SMR WITH weight and generate figure------------------------------------------------------
lsmean_table_weight

lsmean_long_weight <- lsmean_table_weight |>
  pivot_longer(
    cols           = -c(Genotype, weight_model),
    names_to       = c(".value", "group"),
    names_pattern  = "(emmean|lower\\.CL|upper\\.CL|SE)_([FM]_[UI])"
  ) |>
  separate(group, into = c("Sex", "Infected"), sep = "_")

lsmean_long_weight$Genotype <- factor(lsmean_long_weight$Genotype,levels=c("CSBerk", "FFD25", "PC75", "yakB13L5", "san", "teiB13L11", 
                                                             "R84", "suz", "auraL2", "Car5", "sech", "mau31", "sim198"))

lsmean_long_weight$Sex <- factor(lsmean_long_weight$Sex, levels = c("F", "M"))

lsmean_long_weight <- lsmean_long_weight %>%
  mutate(group = interaction(Sex, Infected, sep = "_"),
         group = factor(group, levels = c("F_U", "F_I", "M_U", "M_I")))

p13 <- ggplot(lsmean_long_weight, aes(Genotype, y=emmean, group=group)) +
  geom_errorbar(aes(ymin=lower.CL,
                    ymax=upper.CL),
                width=0.3, alpha=1,
                position = dodge) +
  geom_point(aes(fill = Infected, shape = Sex), colour = "black", size = 2, position = dodge) +
  scale_fill_manual(values = c("U" = "gray80", "I" = "brown1")) +
  scale_shape_manual(values = c("F" = 21, "M" = 22)) +
  guides(color = guide_legend(reverse=TRUE)) +
  ylab("Mass-adjusted LS mean SMR (nmol min-1)")
p13

combined_plots <- ggarrange(p7, p11, p13, nrow=3)
ggsave("output/phylo_comparisons_metabolism.wb.pdf", plot = combined_plots, width = 10, height = 8, dpi=300, useDingbats=FALSE)

###########################################################
##############Allometry analysis##########################
###########################################################
#Examine just uninfected males from the sim198 genotype
# View(subset(hague_df, Genotype == "sim198" & Sex == "M" & Infected =="U"))
#Fly 3152 has an unusually large body mass value
subset(hague_df, flyID == "3152")
#Remove fly 3152
# hague_df <- subset(hague_df, flyID != "3152")

# Add log10_Weight and Group columns needed for allometric analysis
hague_df <- hague_df |>
  mutate(
    Group = factor(
      paste(Sex, Infected, sep = "_"),
      levels = c("F_U", "F_I", "M_U", "M_I")
    )
  )

# One row per fly: weight is constant across repeated SMR measurements
hague_df_fly <- hague_df |>
  group_by(flyID, Genotype, Sex, Infected, Group) |>
  summarise(Weight_mg = mean(Weight_mg, na.rm = TRUE), .groups = "drop")

mass_range_summary <- hague_df_fly |>
  group_by(Genotype, Sex, Infected, Group) |>
  summarise(
    n           = n(),
    mean_weight = mean(Weight_mg, na.rm = TRUE),
    min_weight  = min(Weight_mg,  na.rm = TRUE),
    max_weight  = max(Weight_mg,  na.rm = TRUE),
    # log10(max/min): same metric as Chown et al. (2007) Table 1 "size range"
    log10_range = log10(max_weight / min_weight),
    .groups     = "drop"
  )

mass_range_summary |>
  arrange(Genotype, Group) |>
  gt(groupname_col = "Genotype") |>
  tab_header(
    title    = md("**Body mass range per genotype and Sex x Infection group**"),
    subtitle = "Contextualizes reliability of allometric slope estimates"
  ) |>
  cols_label(
    Sex         = "Sex",
    Infected    = "Infected",
    Group       = "Group",
    n           = md("*n*"),
    mean_weight = "Mean (mg)",
    min_weight  = "Min (mg)",
    max_weight  = "Max (mg)",
    log10_range = md("log10 range")
  ) |>
  fmt_number(
    columns  = c(mean_weight, min_weight, max_weight, log10_range),
    decimals = 3
  ) |>
  tab_source_note(md(
    "log10 range = log10(max weight / min weight). Values > 0.5 indicate
    adequate size variation for reliable slope estimation."
  )) |>
  opt_table_outline() |>
  opt_all_caps()

# fit_allometric_model() fits the covariate-adjusted LMM for a single
# genotype and returns a list-column tibble with the fitted model object,
# emtrends()-derived slopes for the four Sex x Infected groups, and a
# tidy model summary for diagnostics.
fit_allometric_model <- function(genotype_val, data_df = hague_df) {
  
  df <- data_df |> filter(Genotype == genotype_val)
  
  # Covariate-adjusted model — mirrors the nuisance-covariate structure used
  # in the existing manuscript SMR analyses.
  # The three-way interaction gives a separate allometric slope for each of
  # the four Sex x Infected groups.
  # (1 | flyID) accounts for repeated SMR measures per fly (~2 per fly).
  model <- lmer(
    log_aveSMR ~
      log_Weight * Sex * Infected +
      aveTempC     +   # chamber temperature (degrees C)
      aveLight_Lux +   # light intensity (lux)
      SMRstart.sec +   # time of day at start of SMR recording (seconds)
      aveFRC_mlmin +   # flow rate (ml/min)
      aveWVppt     +   # water vapor pressure (ppt)
      aveActivity  +   # activity index (ADS)
      (1 | flyID),
    data = df,
    REML = TRUE
  )
  
  # Extract slopes via emtrends() — reports the slope of log10_Weight at
  # each combination of Sex and Infected with correct SE and 95% CI,
  # accounting for all other terms in the model.
  slopes_df <- emtrends(
    model,
    specs = ~ Sex * Infected,
    var   = "log_Weight"
  ) |>
    as.data.frame() |>
    rename(
      slope    = log_Weight.trend,
      SE       = SE,
      df_resid = df,
      ci_lower = lower.CL,
      ci_upper = upper.CL
    ) |>
    mutate(Genotype = genotype_val)
  
  tibble(
    Genotype      = genotype_val,
    fitted_model  = list(model),
    slopes        = list(slopes_df),
    model_summary = list(broom.mixed::tidy(model))
  )
}

# Apply to all 13 genotypes
allometric_models <- map(
  unique(hague_df$Genotype),
  fit_allometric_model
) |>
  bind_rows()

# The nested slopes tibbles already contain a Genotype column added inside
# fit_allometric_model(), so we drop the outer Genotype before unnesting to
# avoid a duplicate-column error.
slopes_all <- allometric_models |>
  select(slopes) |>          # drop outer Genotype; inner one is retained
  unnest(slopes) |>
  mutate(
    Group = factor(
      paste(Sex, Infected, sep = "_"),
      levels = c("F_U", "F_I", "M_U", "M_I")
    )
  ) |>
  left_join(
    mass_range_summary |>
      select(Genotype, Group, n, min_weight, max_weight, log10_range),
    by = c("Genotype", "Group")
  ) |>
  mutate(
    # Does the 95% CI exclude zero?
    sig_nonzero  = !(ci_lower <= 0 & ci_upper >= 0),
    # Does the 95% CI include the isometry benchmark (b = 1.0)?
    # TRUE = consistent with isometric scaling; FALSE = hypometric scaling.
    includes_100 = ci_lower <= 1.0 & ci_upper >= 1.0
  ) |>
  select(
    Genotype, Sex, Infected, Group,
    n, min_weight, max_weight, log10_range,
    slope, SE, df_resid, ci_lower, ci_upper,
    sig_nonzero, includes_100
  ) |>
  arrange(Genotype, Group)

slopes_all <- slopes_all |>
  mutate(Genotype = factor(Genotype,levels=c("CSBerk", "FFD25", "PC75", "yakB13L5", "san", "teiB13L11",
                                             "R84", "suz", "auraL2", "Car5", "sech", "mau31", "sim198")))

slopes_all$Group <- factor(slopes_all$Group, levels = c("F_U", "F_I", "M_U", "M_I") )

p15 <- ggplot(slopes_all, aes(Genotype, y=slope, group=Group)) +
  geom_errorbar(aes(ymin=ci_lower,
                    ymax=ci_upper),
                width=0.3, alpha=1,
                position = dodge) +
  geom_point(aes(fill = Infected, shape = Sex), colour = "black", size = 2, position = dodge) +
  geom_hline(yintercept = 0,   linetype = "dashed", color = "grey50", linewidth = 0.5) + # Zero reference line
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey50", linewidth = 0.5) + # Isometry benchmark line (b = 1.0)
  scale_fill_manual(values = c("U" = "gray80", "I" = "brown1")) +
  scale_shape_manual(values = c("F" = 21, "M" = 22)) +
  guides(color = guide_legend(reverse=TRUE)) +
  ylab("Mass scaling exponent")
p15

ggsave("output/phylo_comparisons_scaling.wb.pdf", plot = p15, width = 10, height = 2.6, dpi=300, useDingbats=FALSE)

slopes_all$Group <- factor(slopes_all$Group, levels = c("F_U", "F_I", "M_U", "M_I") )

group_colors <- c(
  "F_U" = "#E69F00",  # blue       — Female Uninfected
  "F_I" = "#D55E00",  # sky blue   — Female Infected
  "M_U" = "#0072B2",  # vermilion  — Male Uninfected
  "M_I" = "#56B4E9"   # orange     — Male Infected
)

group_labels <- c(
  "F_U" = "Female Uninfected",
  "F_I" = "Female Infected",
  "M_U" = "Male Uninfected",
  "M_I" = "Male Infected"
)

p16 <- hague_df |>
        mutate(
          Genotype_plot = factor(Genotype,levels=c("CSBerk", "FFD25", "PC75", "yakB13L5", "san", "teiB13L11",
                                                   "R84", "suz", "auraL2", "Car5", "sech", "sim198", "mau31"))
        ) |>
        ggplot(aes(
          x     = log_Weight,
          y     = log_aveSMR,
          color = Group,
          fill  = Group
        )) +
        geom_point(alpha = 0.12, size = 0.75, show.legend = FALSE) +
        geom_smooth(
          method    = "lm",
          formula   = y ~ x,
          se        = TRUE,
          alpha     = 0.12,
          linewidth = 0.65
        ) +
        facet_wrap(~ Genotype_plot, scales = "fixed", ncol = 7) +
        scale_color_manual(values = group_colors, labels = group_labels, name = NULL) +
        scale_fill_manual( values = group_colors, labels = group_labels, name = NULL) +
        labs(
          x = expression(log[10] * "(Weight, mg)"),
          y = expression(log[10] * "(SMR, nmol " * min^{-1} * ")")
        ) +
        theme(
          strip.text      = element_text(size = 8, face = "italic"),
          axis.text       = element_text(size = 7),
          legend.position = "right",
          legend.text     = element_text(size = 8.5)
        )

p16

ggsave("output/allData_allometry.pdf", plot=p16, width=18, height=5, useDingbats=FALSE)

slopes_all |>
  arrange(Genotype, Group) |>
  mutate(CI_95 = sprintf("[%.3f, %.3f]", ci_lower, ci_upper)) |>
  select(Genotype, Group, n, log10_range, slope, SE, CI_95,
         sig_nonzero, includes_100) |>
  gt(groupname_col = "Genotype") |>
  tab_header(
    title    = md("**Allometric scaling exponents — covariate-adjusted model**"),
    subtitle = md(paste(
      "log10(aveSMR) ~ log10(Weight) x Sex x Infected",
      "+ covariates + (1 | flyID)"
    ))
  ) |>
  cols_label(
    Group        = "Group",
    n            = md("*n*"),
    log10_range  = md("Mass range\n(log10)"),
    slope        = "Slope",
    SE           = "SE",
    CI_95        = "95% CI",
    sig_nonzero  = md("Slope \u2260 0"),
    includes_100 = md("CI includes\nb = 1.0")
  ) |>
  fmt_number(columns = c(slope, SE, log10_range), decimals = 3) |>
  fmt_number(columns = n, decimals = 0) |>
  cols_align(
    align   = "center",
    columns = c(n, log10_range, slope, SE, CI_95, sig_nonzero, includes_100)
  ) |>
  tab_spanner(label = "Allometric slope",
              columns = c(slope, SE, CI_95)) |>
  tab_spanner(label = "Isometry check",
              columns = c(sig_nonzero, includes_100)) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(columns = slope, rows = sig_nonzero == TRUE)
  ) |>
  tab_style(
    style     = cell_fill(color = "#DAEEF8"),
    locations = cells_body(rows = Group %in% c("F_U", "F_I"))
  ) |>
  tab_style(
    style     = cell_fill(color = "#FCE8D0"),
    locations = cells_body(rows = Group %in% c("M_U", "M_I"))
  ) |>
  tab_source_note(md(paste(
    "Groups: F_U = Female Uninfected, F_I = Female Infected,",
    "M_U = Male Uninfected, M_I = Male Infected.",
    "**Bold slope**: 95% CI excludes zero.",
    "CI includes b = 1.0: TRUE indicates the slope is consistent with",
    "isometric scaling; FALSE indicates hypometric scaling.",
    "Mass range = log10(max/min weight) within group."
  ))) |>
  opt_table_outline() |>
  opt_all_caps()

# Helper: round all double columns to 3 decimal places for cleaner export
round_numerics <- function(df, digits = 3) {
  df |> mutate(across(where(is.double), \(x) round(x, digits)))
}

# CSV files
write_csv(
  slopes_all |> round_numerics(),
  "output/allometric_slopes_covadjusted.csv"
)
###########################################################
# Plot activity data from Hague et al 2021 Biol Lett
###########################################################
# Load file with summary data
hague_df_activity <- read.csv(file = "data/activity_summaryData.csv")

# Remove row where activity recorder seems to not have worked
hague_df_activity <- hague_df_activity[hague_df_activity$hours123_rawActivity_count > 0,]
hague_df_activity <- subset(hague_df_activity, is.na(flyID) == FALSE)

# Log transformation
hague_df_activity$log_hours123_meanActivity <- log10(hague_df_activity$hours123_meanActivity)

#-------Calculate LS mean data for activity----------------------------------
fit_model_by_species_activity <- function(Genotype_val,data_df = hague_df_activity){
  
  df <- data_df |>
    filter(Genotype == Genotype_val)
  
  lm_fit <- lm(log_hours123_meanActivity ~ Infected*Sex + 
                   Chamber + hours123_startTime + hours123_meanWVppt + hours123_meanRH_pct + hours123_meanTempC + hours123_meanLight_Lux, 
                   data=df)
  
  emm_sex_infected <- emmeans(lm_fit,specs = pairwise ~ Sex:Infected) 
  emm_sex_infected_df <- as.data.frame(emm_sex_infected$emmeans)
  contrasts_sex_infected_df <- as.data.frame(emm_sex_infected$contrasts)
  
  return(tibble(Genotype = Genotype_val,
                fitted_model = list(lm_fit),
                model_summary = list(broom.mixed::tidy(lm_fit)),
                emm_sex_infected = list(emm_sex_infected_df),
                contrasts_sex_infected = list(contrasts_sex_infected_df)))
  
}

# apply the function to each unique Genotype and combine results
activity_models <- map(unique(hague_df_activity$Genotype),
                                  fit_model_by_species_activity) |>
  bind_rows()

glimpse(activity_models)

### LS mean tables for models WITH weight
lsmean_table_activity <- activity_models |>
  select(Genotype, emm_sex_infected) |>
  unnest(emm_sex_infected) |>
  
  # back-transform LS means and confidence intervals
  mutate(emmean = 10^(emmean),
         SE = 10^(SE),
         lower.CL = 10^(lower.CL),
         upper.CL = 10^(upper.CL)) |>
  
  # Create a combined factor for columns
  unite(
    col = "sex_infection",
    Sex, Infected,
    sep = "_"
  ) |>
  
  # Keep only LS means
  select(Genotype, sex_infection, emmean, lower.CL, upper.CL, SE) |>
  
  pivot_wider(
    names_from  = sex_infection,
    values_from = c(emmean, lower.CL, upper.CL, SE)
  )

lsmean_table_activity$weight_model <- "activity"

lsmean_table_activity <- lsmean_table_activity |>
  mutate(Genotype = factor(Genotype, levels=c("sim198", "mau31",  "sech", "Car5","auraL2", "suz", "R84", "teiB13L11", "san", "yakB13L5", "PC75", "FFD25", "CSBerk"))) |>
  arrange(desc(Genotype))

glimpse(lsmean_table_activity)

#-----summarize all SMR WITH weight and generate figure------------------------------------------------------
lsmean_table_activity

lsmean_long_activity <- lsmean_table_activity |>
  pivot_longer(
    cols           = -c(Genotype, weight_model),
    names_to       = c(".value", "group"),
    names_pattern  = "(emmean|lower\\.CL|upper\\.CL|SE)_([FM]_[UI])"
  ) |>
  separate(group, into = c("Sex", "Infected"), sep = "_")

lsmean_long_activity$Genotype <- factor(lsmean_long_activity$Genotype,levels=c("CSBerk", "FFD25", "PC75", "yakB13L5", "san", "teiB13L11", 
                                                                           "R84", "suz", "auraL2", "Car5", "sech", "mau31", "sim198"))

lsmean_long_activity$Sex <- factor(lsmean_long_activity$Sex, levels = c("F", "M"))

lsmean_long_activity <- lsmean_long_activity %>%
  mutate(group = interaction(Sex, Infected, sep = "_"),
         group = factor(group, levels = c("F_U", "F_I", "M_U", "M_I")))

p17 <- ggplot(lsmean_long_activity, aes(Genotype, y=emmean, group=group)) +
  geom_errorbar(aes(ymin=lower.CL,
                    ymax=upper.CL),
                width=0.3, alpha=1,
                position = dodge) +
  geom_point(aes(fill = Infected, shape = Sex), colour = "black", size = 2, position = dodge) +
  scale_fill_manual(values = c("U" = "gray80", "I" = "brown1")) +
  scale_shape_manual(values = c("F" = 21, "M" = 22)) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_y_log10() +
  annotation_logticks(sides="l") +
  ylab("LS mean activity (ADS)")
p17

combined_plots <- ggarrange(p7, p11, p13, p15, p17, nrow=5)
ggsave("output/phylo_comparisons_scaling+activity.wb.pdf", plot = combined_plots, width = 10, height = 13, dpi=300, useDingbats=FALSE)


