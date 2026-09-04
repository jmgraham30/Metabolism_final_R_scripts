# load essential packages
library(tidyverse) # for data manipulation and visualization
library(lmerTest) # for linear mixed-effects models
library(car) # for Anova() Type III F-tests
library(pbkrtest) # for Kenward-Roger denominator df in F-tests
library(modelsummary) # for data and model summaries
library(marginaleffects) # for model interrogation and interpretation 
library(emmeans) # for estimated marginal means
library(patchwork) # for combining ggplots
library(ggthemes) # for additional color themes
library(ape) # for phylogenetic analysis 
library(ggtree) # for phylogenetic tree visualization
library(gt) # for creating summary tables

# set a default theme for ggplot2
theme_set(theme_bw())

# read in the dataset and phylo tree
hague_df <- read_csv("data/metabolism_summaryData.csv")

# filter out rows with missing flyID or Weight_mg, and rows with Notes
hague_df <- hague_df |>
  filter(!is.na(flyID) & !is.na(Weight_mg)) |>
  filter(is.na(Notes)) |>
  mutate(aveSMRpl = log10(aveSMR + abs(min(aveSMR) + 1))) |>
  mutate(log_aveSMR = log10(aveSMR)) |>
  mutate(log_Weight = log10(Weight_mg))

# scaling: log10(aveSMR + abs(min(aveSMR) + 1))
shift <- abs(min(hague_df$aveSMR, na.rm = TRUE) + 1)

# re-level Infected factor so U is the reference level
hague_df <- hague_df |>
  mutate(Infected = relevel(factor(Infected),ref="U"),
         SMRstart.sec = SMRstart.sec/10000)

#Examine just uninfected males from the sim198 genotype
# View(subset(hague_df, Genotype == "sim198" & Sex == "M" & Infected =="U"))
#Fly 3152 has an unusually large body mass value
subset(hague_df, flyID == "3152")
#Remove fly 3152
# hague_df <- subset(hague_df, flyID != "3152")

# reorder Genotype levels
hague_df <- hague_df |>
  mutate(Genotype = factor(Genotype,levels=c("CSBerk", "FFD25", "PC75", "yakB13L5", "san", "teiB13L11", "R84", "suz", "auraL2", "Car5", "sech", "sim198", "mau31")))

# plot histogram of aveActivity, faceted by Genotype
p0 <- hague_df |>
          ggplot(aes(x=aveActivity)) + 
          geom_histogram(fill="steelblue",color="white") +
          facet_wrap(~Genotype, nrow=2) + 
          labs(x="Average Activity (ADS)", 
               y="Count") 

ggsave("output/allData_aveActivity.pdf", plot=p0, width=15, height=5, useDingbats=FALSE)

# peak at the data
glimpse(hague_df)

dodge <- position_dodge(width = .5)
min <- min(hague_df$aveSMR)
max <- max(hague_df$aveSMR)

p1 <- ggplot(data=hague_df, aes(x=Sex, y=aveSMR, fill=Infected, group=interaction(Sex, Infected))) + 
  facet_wrap(. ~ Genotype, nrow=2) +
  geom_violin(color=NA, position=dodge, aes(fill=Infected), width=1) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, position=dodge, fill="white") +
  scale_fill_manual(values=c("gray80", "brown1")) +
  scale_y_continuous(trans='log10', limits = c(min, max)) +
  annotation_logticks(sides="l") +
  labs(y = "SMR")

p1

ggsave("output/allData_SMR.pdf", plot=p1, width=9, height=5, useDingbats=FALSE)

###############################################################################
############ Count unique flyIDs per genotype per analysis ####################
###############################################################################

smr_common_vars <- c("aveSMRpl", "Sex", "Infected",
                     "aveTempC", "aveLight_Lux", "SMRstart.sec",
                     "aveFRC_mlmin", "aveWVppt", "aveActivity")

n_smr <- hague_df |>
  drop_na(all_of(c(smr_common_vars, "Weight_mg"))) |>
  group_by(Genotype) |>
  summarise(n = n_distinct(flyID), .groups = "drop")

################################################################################
################## Fit Linear Mixed Models WITHOUT weight ######################
################################################################################

# function to fit linear mixed models
fit_model_by_species <- function(Genotype_val,data_df = hague_df){
  
  df <- data_df |>
    filter(Genotype == Genotype_val)
  
  lmer_fit <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + aveTempC + aveLight_Lux +
                     SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                     (1 | flyID), 
                     data=df)
  
  emm_sex_infected <- emmeans(lmer_fit,specs = pairwise ~ Sex:Infected) 
  emm_sex_infected_df <- as.data.frame(emm_sex_infected$emmeans)
  contrasts_sex_infected_df <- as.data.frame(emm_sex_infected$contrasts)
  
  # Type III ANOVA with F-tests and Kenward-Roger denominator degrees of freedom.
  # https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
  lmer_fit_sum <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + aveTempC + aveLight_Lux +
                         SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                         (1 | flyID), 
                       data = df,
                       contrasts = list(Sex = contr.sum, Infected = contr.sum))
  
  anova_F <- car::Anova(lmer_fit_sum, type = "III", test.statistic = "F")
  anova_df <- anova_F |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    dplyr::rename(num_df = Df, den_df = Df.res, F_value = `F`, p.value = `Pr(>F)`)
  
  return(tibble(Genotype = Genotype_val,
                fitted_model = list(lmer_fit),
                model_summary = list(broom.mixed::tidy(lmer_fit)),
                anova_F_results = list(anova_df),
                emm_sex_infected = list(emm_sex_infected_df),
                contrasts_sex_infected = list(contrasts_sex_infected_df)))
  
}

# test the function on one species
test_model <- fit_model_by_species("R84")
glimpse(test_model)
test_model |> unnest(model_summary)
test_model |>
  unnest(anova_F_results)
test_model |>
  unnest(emm_sex_infected)
test_model |>
  unnest(contrasts_sex_infected)

# apply the function to each unique Genotype and combine results
metabolic_activity_models <- map(unique(hague_df$Genotype),
                                        fit_model_by_species) |>
  bind_rows()


glimpse(metabolic_activity_models)

# unnest model summaries for all species
metabolic_activity_models_summary <- metabolic_activity_models |>
  unnest(model_summary)
glimpse(metabolic_activity_models_summary)

# plot estimates and confidence intervals for SexM, InfectedI, SexM:InfectedI per species
effects_to_plot <- metabolic_activity_models_summary |>
  filter(term %in% c("SexM", "InfectedI", "SexM:InfectedI")) |>
  mutate(significant = ifelse(p.value < 0.05, "Yes", "No"))

glimpse(effects_to_plot)

# create a ggplot of the effects
p2 <- effects_to_plot |>
  ggplot(aes(x=Genotype, y=estimate, color=term, shape=significant)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_errorbar(aes(ymin=estimate - 1.96 * std.error,
                    ymax=estimate + 1.96 * std.error),
                position=position_dodge(width=0.5), width=0, alpha=1) +
  labs(title="Effects of Sex and Infection Status on Metabolic Rate by Species",
       x="Genotype",
       y="Estimate (with 95% CI)",
       color="Term",
       shape="p < 0.05") +
  scale_color_manual(values=c("#D55E00", "#E69f00", "#56B4E9")) +
  scale_shape_manual(
    values = c("Yes" = 16, "No" = 1)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
p2

ggsave("output/wb.SMR_noWeight.pdf", plot = p2, width = 6, height = 4, dpi=300, useDingbats=FALSE)

# unnest model summaries for all species
metabolic_activity_models_summary <- metabolic_activity_models |>
  unnest(model_summary) |>
  select(Genotype, term, estimate, statistic, p.value) |>
  filter(term %in% c("aveFRC_mlmin","aveLight_Lux", "aveTempC", "aveWVppt", "aveActivity", "InfectedI",
                     "SMRstart.sec","SexM","SexM:InfectedI")) |>
  # add significance column
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

glimpse(metabolic_activity_models_summary)

# define predictor order for column grouping
terms <- c("SexM", "InfectedI", "SexM:InfectedI",
                  "aveTempC", "aveLight_Lux", "SMRstart.sec",
                  "aveFRC_mlmin", "aveWVppt", "aveActivity")
suffixes <- c("_estimate", "_statistic", "_p.value", "_significance")
col_order <- c("Genotype", "n",
                      as.vector(t(outer(terms, suffixes, paste0))))

csv_table_wide <- metabolic_activity_models_summary |>
  mutate(
    estimate  = round(estimate, 3),
    statistic = round(statistic, 3)
  ) |>
  pivot_wider(
    names_from = term,
    values_from = c(estimate, statistic, p.value, significance),
    names_glue = "{term}_{.value}"
  ) |>
  left_join(n_smr, by = "Genotype") |>
  select(all_of(col_order)) |>
  arrange(Genotype)

# save as CSV
write_csv(csv_table_wide, "output/wb.SMR_noWeight.csv")

#-------------------------------------------------------------------------------
# Type III ANOVA tables
#-------------------------------------------------------------------------------

# unnest the per-Genotype ANOVA results
anova_table <- metabolic_activity_models |>
  select(Genotype, anova_F_results) |>
  unnest(anova_F_results) |>
  # drop the intercept row
  filter(term != "(Intercept)") |>
  # add significance stars
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

glimpse(anova_table)

# define predictor order for column grouping
anova_terms <- c("Sex", "Infected", "Sex:Infected",
                 "aveTempC", "aveLight_Lux", "SMRstart.sec",
                 "aveFRC_mlmin", "aveWVppt", "aveActivity")
anova_suffixes <- c("_F_value", "_num_df", "_den_df", "_p.value", "_significance")
anova_col_order <- c("Genotype", "n",
                     as.vector(t(outer(anova_terms, anova_suffixes, paste0))))

anova_table_wide <- anova_table |>
  mutate(
    F_value = round(F_value, 3),
    den_df  = round(den_df, 1),
    p.value = round(p.value, 4)
  ) |>
  pivot_wider(
    names_from  = term,
    values_from = c(F_value, num_df, den_df, p.value, significance),
    names_glue  = "{term}_{.value}"
  ) |>
  left_join(n_smr, by = "Genotype") |>
  select(any_of(anova_col_order)) |>
  arrange(Genotype)

glimpse(anova_table_wide)

# save as CSV
write_csv(anova_table_wide, "output/wb.SMR_noWeight_ANOVA.csv")

PC75_test <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + aveTempC + aveLight_Lux +
                   SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                   (1 | flyID), 
                   data=subset(hague_df, Genotype == "PC75"))

summary(PC75_test)

PC75_test_ANOVA <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + aveTempC + aveLight_Lux +
                       SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                       (1 | flyID), 
                       data = subset(hague_df, Genotype == "PC75"),
                       contrasts = list(Sex = contr.sum, Infected = contr.sum))

car::Anova(PC75_test_ANOVA, type = "III", test.statistic = "F")

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
  # filter(Genotype != "sim198") |>
  mutate(Genotype = factor(Genotype, levels=c("mau31", "sim198", "sech", "Car5","auraL2", "suz", "R84", "teiB13L11", "san", "yakB13L5", "PC75", "FFD25", "CSBerk"))) |>
  arrange(desc(Genotype))

glimpse(lsmean_table)

################################################################################
################## Fit Linear Mixed Models WITH weight #########################
################################################################################

# function to fit linear mixed models
fit_model_by_species_weight <- function(Genotype_val,data_df = hague_df){
  
  df <- data_df |>
    filter(Genotype == Genotype_val)
  
  lmer_fit <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + aveTempC + aveLight_Lux +
                     SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                     log_Weight +
                     (1 | flyID), 
                     data=df)
  
  emm_sex_infected <- emmeans(lmer_fit,specs = pairwise ~ Sex:Infected) 
  emm_sex_infected_df <- as.data.frame(emm_sex_infected$emmeans)
  contrasts_sex_infected_df <- as.data.frame(emm_sex_infected$contrasts)
  
  # Type III ANOVA with F-tests and Kenward-Roger denominator degrees of freedom.
  # https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
  lmer_fit_sum <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + aveTempC + aveLight_Lux +
                         SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                         log_Weight +
                         (1 | flyID), 
                       data = df,
                       contrasts = list(Sex = contr.sum, Infected = contr.sum))
  
  anova_F <- car::Anova(lmer_fit_sum, type = "III", test.statistic = "F")
  anova_df <- anova_F |>
    as.data.frame() |>
    tibble::rownames_to_column("term") |>
    dplyr::rename(num_df = Df, den_df = Df.res, F_value = `F`, p.value = `Pr(>F)`)
  
  return(tibble(Genotype = Genotype_val,
                fitted_model = list(lmer_fit),
                model_summary = list(broom.mixed::tidy(lmer_fit)),
                anova_F_results = list(anova_df),
                emm_sex_infected = list(emm_sex_infected_df),
                contrasts_sex_infected = list(contrasts_sex_infected_df)))
  
}

# test the function on one species
test_model <- fit_model_by_species_weight("R84")
glimpse(test_model)
test_model |> unnest(model_summary)
test_model |>
  unnest(anova_F_results)
test_model |>
  unnest(emm_sex_infected)
test_model |>
  unnest(contrasts_sex_infected)

# apply the function to each unique Genotype and combine results
metabolic_activity_models_weight <- map(unique(hague_df$Genotype),
                                        fit_model_by_species_weight) |>
  bind_rows()


glimpse(metabolic_activity_models_weight)

# unnest model summaries for all species
metabolic_activity_models_summary_weight <- metabolic_activity_models_weight |>
  unnest(model_summary)
glimpse(metabolic_activity_models_summary_weight)

# plot estimates and confidence intervals for SexM, InfectedI, SexM:InfectedI, and log_Weight per species
effects_to_plot_weight <- metabolic_activity_models_summary_weight |>
  filter(term %in% c("SexM", "InfectedI", "SexM:InfectedI", "log_Weight")) |>
  mutate(significant = ifelse(p.value < 0.05, "Yes", "No"))

glimpse(effects_to_plot_weight)

# create a ggplot of the effects
p3 <- effects_to_plot_weight |>
  ggplot(aes(x=Genotype, y=estimate, color=term, shape=significant)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_errorbar(aes(ymin=estimate - 1.96 * std.error,
                    ymax=estimate + 1.96 * std.error),
                position=position_dodge(width=0.5), width=0, alpha=1) +
  labs(title="Effects of Sex, Infection Status, and Weight on Metabolic Rate by Species",
       x="Genotype",
       y="Estimate (with 95% CI)",
       color="Term",
       shape="p < 0.05") +
  scale_color_manual(values=c("#D55E00", "#009E73", "#E69f00", "#56B4E9")) +
  scale_shape_manual(
    values = c("Yes" = 16, "No" = 1)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
p3

ggsave("output/wb.SMR_Weight.pdf", plot = p3, width = 6, height = 4, dpi=300, useDingbats=FALSE)

# unnest model summaries for all species
metabolic_activity_models_summary_weight <- metabolic_activity_models_weight |>
  unnest(model_summary) |>
  select(Genotype, term, estimate, statistic, p.value) |>
  filter(term %in% c("aveFRC_mlmin","aveLight_Lux", "aveTempC", "aveWVppt", "aveActivity", "InfectedI",
                     "SMRstart.sec","log_Weight","SexM","SexM:InfectedI")) |>
  # add significance column
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

glimpse(metabolic_activity_models_summary_weight)

# define predictor order for column grouping
terms_weight <- c("SexM", "InfectedI", "SexM:InfectedI", "log_Weight",
                  "aveTempC", "aveLight_Lux", "SMRstart.sec",
                  "aveFRC_mlmin", "aveWVppt", "aveActivity")
col_order_weight <- c("Genotype", "n",
                      as.vector(t(outer(terms_weight, suffixes, paste0))))

csv_table_wide_weight <- metabolic_activity_models_summary_weight |>
  mutate(
    estimate  = round(estimate, 3),
    statistic = round(statistic, 3)
  ) |>
  pivot_wider(
    names_from = term,
    values_from = c(estimate, statistic, p.value, significance),
    names_glue = "{term}_{.value}"
  ) |>
  left_join(n_smr, by = "Genotype") |>
  select(all_of(col_order_weight)) |>
  arrange(Genotype)

# save as CSV
write_csv(csv_table_wide_weight, "output/wb.SMR_Weight.csv")

#-------------------------------------------------------------------------------
# Type III ANOVA tables
#-------------------------------------------------------------------------------

# unnest the per-genotype ANOVA results
anova_table_weight <- metabolic_activity_models_weight |>
  select(Genotype, anova_F_results) |>
  unnest(anova_F_results) |>
  # drop the intercept row
  filter(term != "(Intercept)") |>
  # add significance stars
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

glimpse(anova_table_weight)

# define predictor order for column grouping
anova_terms_weight <- c("Sex", "Infected", "Sex:Infected", "log_Weight",
                        "aveTempC", "aveLight_Lux", "SMRstart.sec",
                        "aveFRC_mlmin", "aveWVppt", "aveActivity")
anova_suffixes <- c("_F_value", "_num_df", "_den_df", "_p.value", "_significance")
anova_col_order_weight <- c("Genotype", "n",
                            as.vector(t(outer(anova_terms_weight, anova_suffixes, paste0))))

anova_table_wide_weight <- anova_table_weight |>
  mutate(
    F_value = round(F_value, 3),
    den_df  = round(den_df, 1),
    p.value = round(p.value, 4)
  ) |>
  pivot_wider(
    names_from  = term,
    values_from = c(F_value, num_df, den_df, p.value, significance),
    names_glue  = "{term}_{.value}"
  ) |>
  left_join(n_smr, by = "Genotype") |>
  select(any_of(anova_col_order_weight)) |>
  arrange(Genotype)

glimpse(anova_table_wide_weight)

# save as CSV
write_csv(anova_table_wide_weight, "output/wb.SMR_Weight_ANOVA.csv")

PC75_test <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + aveTempC + aveLight_Lux +
                    SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                    log_Weight +
                    (1 | flyID), 
                  data=subset(hague_df, Genotype == "PC75"))

summary(PC75_test)

PC75_test_ANOVA <- lmer(log_aveSMR  ~ 1 + Sex + Infected + Sex:Infected + aveTempC + aveLight_Lux +
                          SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity + 
                          log_Weight +
                          (1 | flyID), 
                        data = subset(hague_df, Genotype == "PC75"),
                        contrasts = list(Sex = contr.sum, Infected = contr.sum))

car::Anova(PC75_test_ANOVA, type = "III", test.statistic = "F")

### LS mean tables for Weight Only Models
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

lsmean_table_weight$weight_model <- "Weight_Only"

lsmean_table_weight <- lsmean_table_weight |>
  mutate(Genotype = factor(Genotype, levels=c("mau31", "sim198", "sech", "Car5","auraL2", "suz", "R84", "teiB13L11", "san", "yakB13L5", "PC75", "FFD25", "CSBerk"))) |>
  arrange(desc(Genotype))

glimpse(lsmean_table_weight)



