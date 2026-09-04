# load essential packages
library(tidyverse) # for data manipulation and visualization
library(lme4) # for linear mixed-effects models
library(modelsummary) # for data and model summaries
library(emmeans) # for estimated marginal means
library(marginaleffects) # for model interrogation and interpretation 
library(ggthemes) # for additional color themes
library(broom) # for tidying models 
library(car) # for Anova() Type III F-tests
library(gt) # for creating summary tables
library(Hmisc) # for 95% CIs

# set a default theme for ggplot2
theme_set(theme_bw())

# read in the dataset
hague_df <- read_csv("data/metabolism_summaryData.csv")

# filter out rows with missing flyID or Weight_mg, and rows with Notes
hague_df <- hague_df |>
  filter(!is.na(flyID) & !is.na(Weight_mg)) |>
  filter(is.na(Notes)) |>
  mutate(aveSMRpl = log10(aveSMR + abs(min(aveSMR) + 1))) |>
  mutate(log_aveSMR = log10(aveSMR)) |>
  mutate(log_Weight = log10(Weight_mg)) 

# re-level Infected factor so U is the reference level
hague_df <- hague_df |>
  mutate(Infected = relevel(factor(Infected),ref="U"))

# reorder Genotype levels
hague_df <- hague_df |>
  mutate(Genotype = factor(Genotype,levels=c("CSBerk", "FFD25", "PC75", "yakB13L5", "san", "teiB13L11", "R84", "suz", "auraL2", "Car5", "sech", "sim198", "mau31")))

# peak at the data
glimpse(hague_df)

# create a reduced dataset with mean Weight_mg per flyID
hague_df_weight_reduced <- hague_df |>
  dplyr::group_by(flyID) |>
  dplyr::summarize(Weight_mg=mean(Weight_mg))

glimpse(hague_df_weight_reduced)

# extract Infected and Sex per ID
hague_df_id_info <- hague_df |>
  select(flyID,Infected,Sex,Genotype) |>
  distinct()

glimpse(hague_df_id_info)

# merge mean weight info with hague_df_id_info
hague_df_weight_reduced <- hague_df_weight_reduced |>
  left_join(hague_df_id_info,by="flyID") 

glimpse(hague_df_weight_reduced)

weight_means_melt <- hague_df_weight_reduced %>%
  group_by(Genotype, Sex, Infected) %>%
  summarise(
    n = n(),
    as_tibble(as.list(smean.cl.normal(Weight_mg, conf.int=0.95))),
    .groups = "drop")

weight_means_melt

#Examine just uninfected males from the sim198 genotype
# View(subset(hague_df_weight_reduced, Genotype == "sim198" & Sex == "M" & Infected =="U"))
#Fly 3152 has an unusually large body mass value
subset(hague_df_weight_reduced, flyID == "3152")
#Remove fly 3152
# hague_df_weight_reduced <- subset(hague_df_weight_reduced, flyID != "3152")

################################################################################
################## Basic Plots ###############################################
################################################################################
# histogram of Weight_mg, faceted by Genotype
hague_df_weight_reduced |>
  ggplot(aes(x=Weight_mg)) +
  geom_histogram(fill="steelblue",color="white") +
  facet_wrap(~Genotype)

# plot boxplot of Weight_mg by Infected status and Genotype
hague_df_weight_reduced |>
  ggplot(aes(x=Genotype,y=Weight_mg,fill=Infected)) +
  geom_boxplot() +
  facet_wrap(~Sex) +
  scale_fill_colorblind() +
  # rotate x-axis labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################
################## BY SPECIES MODELS ###########################################
################################################################################

species_counts <- hague_df_weight_reduced |>
  group_by(Genotype) |>
  summarise(n=n())

# function to fit linear model and extract residuals
lm_model_by_species <- hague_df_weight_reduced |>
  group_by(Genotype) |>
  summarise(
    model = list(lm(Weight_mg ~ Sex + Infected + Sex:Infected, data = pick(everything()))),

    model_sum = list(lm(Weight_mg ~ Sex + Infected + Sex:Infected,
                        data = pick(everything()),
                        contrasts = list(Sex = contr.sum, Infected = contr.sum)))
  ) |>
  mutate(
    mod_resids = map(model, residuals),
    tidy_output = map(model, tidy),        
    glance_output = map(model, glance),     
    
    anova_output = map(
      model_sum,
      ~ car::Anova(.x, type = "III") |>
        as.data.frame() |>
        tibble::rownames_to_column("term") |>
        dplyr::rename(sum_sq = `Sum Sq`, num_df = Df,
                      F_value = `F value`, p.value = `Pr(>F)`)
    ),
    emmeans_output = map(
      model,
      ~ emmeans(.x, ~ Sex * Infected)
    ),
    emmeans_summary = map(
      emmeans_output,
      ~ summary(.x, infer = TRUE) 
    )
  )

glimpse(lm_model_by_species)

lm_model_by_species_tidy <- lm_model_by_species |>
  unnest(tidy_output)

glimpse(lm_model_by_species_tidy)

lm_model_by_species_glance <- lm_model_by_species |>
  unnest(glance_output)

glimpse(lm_model_by_species_glance)

p1 <- lm_model_by_species_tidy |>
        filter(term %in% c("SexM","InfectedI","SexM:InfectedI")) |>
        select(Genotype,term,estimate,std.error,p.value) |>
        left_join(species_counts,by="Genotype") |>
        mutate(significant = ifelse(p.value < 0.05, "yes", "no")) |>
        mutate(lwr_ci = estimate - 1.96 * std.error,
               upr_ci = estimate + 1.96 * std.error) |>
        ggplot(aes(x=Genotype,y=estimate,shape=significant,color=term)) +
        geom_point(position=position_dodge(width=0.5), size=3) +
        geom_errorbar(aes(ymin=lwr_ci,ymax=upr_ci),width=0,alpha=1,
                      position=position_dodge(width=0.5)) +
        scale_color_manual(values=c("#D55E00", "#E69f00", "#56B4E9")) +
        scale_shape_manual(
          values = c("yes" = 16, "no" = 1)) +
        theme(axis.text.x = element_text(angle=45,hjust=1)) +
        labs(y="Estimate with 95% CI",color="Trem")
p1

ggsave("output/wb.weight.pdf", plot = p1, width = 6, height = 4, dpi=300, useDingbats=FALSE)

#-------------------------------------------------------------------------------
# CSV 1: Regression coefficients
#-------------------------------------------------------------------------------

coef_table <- lm_model_by_species_tidy |>
  select(Genotype, term, estimate, std.error, statistic, p.value) |>
  filter(term %in% c("SexM", "InfectedI", "SexM:InfectedI")) |>
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

coef_terms <- c("SexM", "InfectedI", "SexM:InfectedI")
coef_suffixes <- c("_estimate", "_std.error", "_statistic", "_p.value", "_significance")
coef_col_order <- c("Genotype", "n",
                    as.vector(t(outer(coef_terms, coef_suffixes, paste0))))

coef_table_wide <- coef_table |>
  mutate(
    estimate  = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 3),
    p.value   = round(p.value, 4)
  ) |>
  pivot_wider(
    id_cols     = Genotype,
    names_from  = term,
    values_from = c(estimate, std.error, statistic, p.value, significance),
    names_glue  = "{term}_{.value}"
  ) |>
  left_join(species_counts, by = "Genotype") |>
  select(any_of(coef_col_order)) |>
  arrange(Genotype)

glimpse(coef_table_wide)

# save as CSV
write_csv(coef_table_wide, "output/wb.weight_coefficients.csv")

#-------------------------------------------------------------------------------
# CSV 2: Type III F-test ANOVA
#-------------------------------------------------------------------------------

anova_table <- lm_model_by_species |>
  select(Genotype, anova_output) |>
  unnest(anova_output) |>
  filter(!(term %in% c("(Intercept)", "Residuals"))) |>
  mutate(significance = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

glimpse(anova_table)

anova_resid_df <- lm_model_by_species |>
  select(Genotype, anova_output) |>
  unnest(anova_output) |>
  filter(term == "Residuals") |>
  select(Genotype, den_df = num_df)

anova_terms <- c("Sex", "Infected", "Sex:Infected")
anova_suffixes <- c("_F_value", "_num_df", "_p.value", "_significance")
anova_col_order <- c("Genotype", "n", "den_df",
                     as.vector(t(outer(anova_terms, anova_suffixes, paste0))))

anova_table_wide <- anova_table |>
  mutate(
    F_value = round(F_value, 3),
    p.value = round(p.value, 4)
  ) |>
  pivot_wider(
    id_cols     = Genotype,
    names_from  = term,
    values_from = c(F_value, num_df, p.value, significance),
    names_glue  = "{term}_{.value}"
  ) |>
  left_join(species_counts, by = "Genotype") |>
  left_join(anova_resid_df, by = "Genotype") |>
  select(any_of(anova_col_order)) |>
  arrange(Genotype)

glimpse(anova_table_wide)

# save as CSV
write_csv(anova_table_wide, "output/wb.weight_ANOVA.csv")

sech_test <- lm(Weight_mg ~ Sex + Infected + Sex:Infected,
                data = subset(hague_df_weight_reduced, Genotype == "sech"),
                contrasts = list(Sex = contr.sum, Infected = contr.sum))

car::Anova(sech_test, type = "III")

dodge <- position_dodge(width = .5)
min <- min(hague_df_weight_reduced$Weight_mg)
max <- max(hague_df_weight_reduced$Weight_mg)

p2 <- ggplot(data=hague_df_weight_reduced, aes(x=Sex, y=Weight_mg, fill=Infected, group=interaction(Sex, Infected))) + 
  facet_wrap(. ~ Genotype, nrow=2) +
  geom_violin(color=NA, position=dodge, aes(fill=Infected), width=1) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, position=dodge, fill="white") +
  scale_fill_manual(values=c("gray80", "brown1")) +
  scale_y_continuous(trans='log10', limits = c(min, max)) +
  annotation_logticks(sides="l") +
  labs(y = "Body mass (mg)")

p2

ggsave("output/allData_weight.pdf", plot=p2, width=9, height=5, useDingbats=FALSE)


