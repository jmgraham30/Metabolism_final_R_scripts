library(caper)
library(phytools)
library(ape)
library(geiger)
library(ggplot2); theme_set(theme_bw())
library(pmc)
library(dplyr)
library(tibble)
library(reshape)
library(TreeSim)
library(egg)
library(emmeans)
library(lme4)
library(Hmisc)
library(ggtree)
library(readr)
library(performance)
library(tidyr)
options(rgl.useNULL = TRUE)
library(geomorph)
library(tidyverse)

#---------Load data---------------------------------------------
#load host phylogeny
host.tree <- read.nexus("data/phylogram_2_host.tre")

host_map <- c(
  "CSBerk"    = "DMel_CSBerk",
  "FFD25"     = "DMel_FFD25",
  "PC75"      = "DMel_PC75",
  "yakB13L5"  = "DYak",
  "san"       = "san-Quija630-39",
  "teiB13L11" = "teis-cascade_4_2",
  "R84"       = "DSim",
  "suz"       = "DSuz",
  "auraL2"    = "auraria_bronski",
  "Car5"      = "DSim",
  "sech"      = "DSech",
  "sim198"    = "DSim",
  "mau31"     = "mauritiana"
)

#---------------------------------------------------------------------------------
#---------------------Analysis of weight data-------------------------------------
#---------------------------------------------------------------------------------
# import weight data
hague_df <- read_csv("data/metabolism_summaryData.csv")

# filter out rows with missing flyID or Weight_mg, and rows with Notes
hague_df <- hague_df |>
  filter(!is.na(flyID) & !is.na(Weight_mg)) |>
  filter(is.na(Notes)) |>
  mutate(log_aveSMR = log10(aveSMR)) |>
  mutate(log_Weight = log10(Weight_mg))

# re-level Infected factor so U is the reference level
hague_df <- hague_df |>
  mutate(Infected = relevel(factor(Infected),ref="U"),
         SMRstart.sec = SMRstart.sec/10000)

# reorder Genotype levels
hague_df <- hague_df |>
                mutate(Genotype = factor(Genotype,levels=c("CSBerk", "FFD25", "PC75", "yakB13L5", "san", "teiB13L11", 
                                                           "R84", "suz", "auraL2", "Car5", "sech", "sim198", "mau31")))

# calculate mean and std. of temperature from each SMR measurement
mean(hague_df$aveTempC)
sd(hague_df$aveTempC)

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

hague_df_weight_reduced_U <- subset(hague_df_weight_reduced, Infected == "U")

nrow(hague_df_weight_reduced_U) #total number of uninfected flies assayed

weightFIT <- lm(Weight_mg ~ Genotype + Sex, data=hague_df_weight_reduced_U) #Using all Dsim genotypes in analysis

summary(weightFIT)

#Check assumptions of model for linear model
check_model(weightFIT, check = c("normality", "linearity", "homogeneity", "outliers", "vif", "qq"))

#Type III SS
car::Anova(weightFIT, type="III") 

#--------------E-PLGS linear model to test if weight differs by species------------
# Adams & Collyer (2024) Methods in Ecology and Evolution
# Glynne & Adams (2026) Evolution
# https://cran.r-project.org/web/packages/geomorph/geomorph.pdf

# Estimate phylogenetic variance-covariance matrix (not needed for E-PGLS)
C <- vcv.phylo(host.tree)

e_pgls_model <- extended.pgls(
  f1 = Weight_mg ~ Genotype + Sex,  
  data = hague_df_weight_reduced_U, #Include all Dsim in analysis
  species = "Host",
  phy = host.tree)

anova(e_pgls_model)

#--------Pagel's lambda to test for phylogenetic signal-------------------
glimpse(hague_df_weight_reduced_U)

weight_means_melt <- hague_df_weight_reduced_U %>%
                        group_by(Genotype, Sex) %>%
                        summarise(as_tibble(as.list(smean.cl.normal(Weight_mg, conf.int=0.95))))

weight_means <- weight_means_melt %>%
                  ungroup() %>%
                  pivot_wider(
                    id_cols     = Genotype,
                    names_from  = Sex,
                    values_from = c(Mean, Lower, Upper)
                  )

weight_means$Host <- host_map[as.character(weight_means$Genotype)]

glimpse(weight_means)

weight_means <- subset(weight_means, !(Genotype %in% c("Car5", "sim198"))) #Use R84 for Dsim data
# weight_means <- subset(weight_means, !(Genotype %in% c("R84", "sim198"))) #Use Car5 for Dsim data
# weight_means <- subset(weight_means, !(Genotype %in% c("R84", "Car5"))) #Use sim198 for Dsim data

weight_means <- column_to_rownames(weight_means, var="Host")

#--------------------------------------------------------------
#Plot and prune host phylogeny
#http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html
plotTree(host.tree, type="phylogram")

#Check for concordance between tree and trait data
name.check(host.tree, weight_means)

#break polytomies randomly
set.seed(42)
host.tree <- multi2di(host.tree, random=TRUE) #for subsequent analyses
plotTree(host.tree, type="phylogram")

#----------Run analysis for uninfected FEMALES------------------------------------------
host.weight_phylo <- as.data.frame(weight_means$Mean_F)
row.names(host.weight_phylo) <- row.names(weight_means)
host.weight_geiger <- treedata(host.tree, host.weight_phylo)

#Using phytools
phylosig(host.weight_geiger$phy, host.weight_geiger$data, method = "lambda", test = TRUE)

#----------Run analysis for uninfected MALES------------------------------------------
host.weight_phylo <- as.data.frame(weight_means$Mean_M)
row.names(host.weight_phylo) <- row.names(weight_means)
host.weight_geiger <- treedata(host.tree, host.weight_phylo)

#Using phytools
phylosig(host.weight_geiger$phy, host.weight_geiger$data, method = "lambda", test = TRUE)

#--------------------------------------------------------------------------------
#---------------------Analysis of SMR data WITHOUT body mass---------------------
#--------------------------------------------------------------------------------

#----------Pagel's lambda test for phylogenetic signal with SMR-----------------

# function to fit genotype-specific LMMs to extract emmeans
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

# apply the function to each unique genotype and combine results
metabolic_activity_models <- map(unique(hague_df$Genotype),
                                        fit_model_by_species) |>
  bind_rows()

glimpse(metabolic_activity_models)

### LS mean tables for no weight models
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

lsmean_table$weight_model <- "No_Weight"

lsmean_table <- lsmean_table |>
  mutate(Genotype = factor(Genotype, levels=c("mau31", "sim198", "sech", "Car5","auraL2", "suz", "R84", "teiB13L11", "san", "yakB13L5", "PC75", "FFD25", "CSBerk"))) |>
  arrange(desc(Genotype))

glimpse(lsmean_table)

lsmean_table$Host <- host_map[as.character(lsmean_table$Genotype)]
lsmean_table <- as.data.frame(lsmean_table)

glimpse(lsmean_table)

lsmean_table <- subset(lsmean_table, !(Genotype %in% c("Car5", "sim198"))) #Use R84 for Dsim data
# lsmean_table <- subset(lsmean_table, !(Genotype %in% c("R84", "sim198"))) #Use Car5 for Dsim data
# lsmean_table <- subset(lsmean_table, !(Genotype %in% c("R84", "Car5"))) #Use sim198 for Dsim data

rownames(lsmean_table) <- lsmean_table$Host

#--------------------------------------------------------------
#Plot and prune host phylogeny
#http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html
plotTree(host.tree, type="phylogram")

#Check for concordance between tree and trait data
name.check(host.tree, lsmean_table)

#break polytomies randomly
# host.tree <- multi2di(host.tree, random=TRUE) #for subsequent analyses
# plotTree(host.tree, type="phylogram")

#----------Run analysis for uninfected FEMALES------------------------------------------
host.metabolism_phylo <- as.data.frame(lsmean_table$emmean_F_U)
row.names(host.metabolism_phylo) <- row.names(lsmean_table)
host.metabolism_geiger <- treedata(host.tree, host.metabolism_phylo)

#Using phytools
phylosig(host.metabolism_geiger$phy, host.metabolism_geiger$data, method = "lambda", test = TRUE)

#----------Run analysis for uninfected MALES------------------------------------------
host.metabolism_phylo <- as.data.frame(lsmean_table$emmean_M_U)
row.names(host.metabolism_phylo) <- row.names(lsmean_table)
host.metabolism_geiger <- treedata(host.tree, host.metabolism_phylo)

#Using phytools
phylosig(host.metabolism_geiger$phy, host.metabolism_geiger$data, method = "lambda", test = TRUE)

#---------------------------------------------------------------------------------
#---------------------Analysis of SMR data WITH body mass-------------------------
#---------------------------------------------------------------------------------
hague_df_SMR_U <- subset(hague_df, Infected == "U") #select only uninfected flies

# peak at the data
glimpse(hague_df_SMR_U)

hague_df %>% 
  group_by(Host, Genotype) %>% 
  summarise(n_flies = n_distinct(flyID))

length(unique(hague_df_SMR_U$flyID)) #total number of flies assayed

SMR_fit <- lmer(log_aveSMR ~ 1 + Genotype + Sex + log_Weight + Sex*log_Weight + 
                + aveTempC + aveLight_Lux + SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity +
                (1 | flyID), 
                data=hague_df_SMR_U)

summary(SMR_fit)

#Check assumptions of model for linear model
check_model(SMR_fit, check = c("normality", "linearity", "homogeneity", "outliers", "vif", "qq"))

#-----Linear mixed model: test for effect of individual variables----------------------------------
#http://glmm.wikidot.com/
#----Type III ANOVA with F-tests and Kenward-Roger denominator degrees of freedom----------------------------
# Type III SS
car::Anova(SMR_fit, type = "III", test.statistic = "F")

#--------------E-PLGS linear model to test if SMR differs by species------------
# Adams & Collyer (2024) Methods in Ecology and Evolution
# Glynne & Adams (2026) Evolution
# https://cran.r-project.org/web/packages/geomorph/geomorph.pdf

# subset data to only the first replicate because E-PGLS can't accommodate random effects/repeated measures
hague_df_SMR_U_reduced <- subset(hague_df_SMR_U, hague_df_SMR_U$Replicate == 1)

# Estimate phylogenetic variance-covariance matrix (not needed for E-PGLS)
C <- vcv.phylo(host.tree)

e_pgls_model <- extended.pgls(
  f1 = log_aveSMR ~ Genotype + Sex + log_Weight + Sex*log_Weight +
          aveTempC + aveLight_Lux + SMRstart.sec + aveFRC_mlmin + aveWVppt + aveActivity,  
          data = hague_df_SMR_U_reduced, #Include all Dsim in analysis
          species = "Host",
          phy = host.tree)

anova(e_pgls_model)

#---------------------------------------------------------------------------------
#--- Extract mass scaling exponent slopes and 95% CIs ----------------------------
#---------------------------------------------------------------------------------
iso_slope <- 1   # isometric expectation to test against

# ---- lmer: pooled slope from the existing SMR_fit ----
summary(SMR_fit)

lmm_pooled_b  <- fixef(SMR_fit)["log_Weight"]

lmm_pooled_se <- sqrt(vcov(SMR_fit)["log_Weight", "log_Weight"])
lmm_pooled_ci <- lmm_pooled_b + c(-1, 1) * 1.96 * lmm_pooled_se
cat("LMM pooled slope:", round(lmm_pooled_b, 4),
    "  95% CI: [", round(lmm_pooled_ci[1], 4), ",", round(lmm_pooled_ci[2], 4), "]\n")

lmm_sex_trends    <- emtrends(SMR_fit, ~ Sex, var = "log_Weight")
lmm_sex_trends_df <- as.data.frame(summary(lmm_sex_trends, infer = TRUE, level = 0.95))
lmm_sex_trends_df

geno_means_lmm <- hague_df_SMR_U |>
  group_by(Genotype, Sex) |>
  summarise(
    x_mean = mean(log_Weight, na.rm = TRUE),
    x_sd   = sd(log_Weight,   na.rm = TRUE),
    y_mean = mean(log_aveSMR, na.rm = TRUE),
    y_sd   = sd(log_aveSMR,   na.rm = TRUE),
    .groups = "drop")

allom_lmm_plot <- ggplot(geno_means_lmm,
                         aes(x = x_mean, y = y_mean, color = Sex)) +
  geom_errorbar(aes(ymin = y_mean - y_sd, ymax = y_mean + y_sd), width = 0) +
  geom_errorbarh(aes(xmin = x_mean - x_sd, xmax = x_mean + x_sd), height = 0) +
  geom_point(shape = 18, size = 3) +
  geom_smooth(method = "lm", se = FALSE,, linewidth = 1) +
  geom_abline(slope = 1,
              intercept = mean(geno_means_lmm$y_mean) - 1 * mean(geno_means_lmm$x_mean),
              linetype = "dashed", linewidth = 0.6, color = "black") +
  scale_color_manual(values = c("F" = "#E69F00", "M" = "#0072B2")) +
  labs(x = expression(Log[10]~"(Body mass, mg)"),
       y = expression(Log[10]~"(SMR)"),
       color = "Sex") +
  theme_bw()
allom_lmm_plot

ggsave("output/SMR_allometry_LMM.pdf", plot = allom_lmm_plot,
       width = 7, height = 5, dpi = 300, useDingbats = FALSE)

#----------Pagel's lambda test for phylogenetic signal with mass-adjusted SMR-----------------

# function to fit genotype-specific LMMs to extract emmeans
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

# apply the function to each unique genotype and combine results
metabolic_activity_models_weight <- map(unique(hague_df$Genotype),
                                        fit_model_by_species_weight) |>
  bind_rows()

glimpse(metabolic_activity_models_weight)

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

lsmean_table_weight$Host <- host_map[as.character(lsmean_table_weight$Genotype)]
lsmean_table_weight <- as.data.frame(lsmean_table_weight)

glimpse(lsmean_table_weight)

lsmean_table_weight <- subset(lsmean_table_weight, !(Genotype %in% c("Car5", "sim198"))) #Use R84 for Dsim data
# lsmean_table_weight <- subset(lsmean_table_weight, !(Genotype %in% c("R84", "sim198"))) #Use Car5 for Dsim data
# lsmean_table_weight <- subset(lsmean_table_weight, !(Genotype %in% c("R84", "Car5"))) #Use sim198 for Dsim data

rownames(lsmean_table_weight) <- lsmean_table_weight$Host

#--------------------------------------------------------------
#Plot and prune host phylogeny
#http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html
plotTree(host.tree, type="phylogram")

#Check for concordance between tree and trait data
name.check(host.tree, lsmean_table_weight)

#break polytomies randomly
# host.tree <- multi2di(host.tree, random=TRUE) #for subsequent analyses
# plotTree(host.tree, type="phylogram")

#----------Run analysis for uninfected FEMALES------------------------------------------
host.metabolism_phylo <- as.data.frame(lsmean_table_weight$emmean_F_U)
row.names(host.metabolism_phylo) <- row.names(lsmean_table_weight)
host.metabolism_geiger <- treedata(host.tree, host.metabolism_phylo)

#Using phytools
phylosig(host.metabolism_geiger$phy, host.metabolism_geiger$data, method = "lambda", test = TRUE)

#----------Assess power using Boettiger et al. (2012)------------------------------------------------------
# https://github.com/cboettig/pmc
# http://www2.uaem.mx/r-mirror/web/packages/pmc/vignettes/pmc_tutorial.pdf
# Generate ML estimates for tree
# bm_v_lambda_obs_F_U <- pmc(host.metabolism_geiger$phy, host.metabolism_geiger$data, "BM", "lambda", nboot = 1000)
# save(bm_v_lambda_obs_F_U, file="bm_v_lambda_obs_SMR_FU.Rda")
load("output/bm_v_lambda_obs_SMR_FU.Rda")

lambdas_obs <- bm_v_lambda_obs_F_U$par_dists %>% filter(comparison=="BB", parameter=="lambda")

est_obs <- coef(bm_v_lambda_obs_F_U[["B"]])[["lambda"]]
est_obs

mean(lambdas_obs$value)
cast(lambdas_obs, comparison ~ parameter, function(x) quantile(x, c(0.025, 0.975)), value = c("lower", "upper"))

p1 <- ggplot(lambdas_obs) +
  geom_histogram(aes(value), bins=50) +
  geom_vline(xintercept=est_obs, linetype="dashed") +
  coord_cartesian(xlim=c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x="Estimated lambda", y="Count")
p1

# Simulate larger trees with observed lambda value to assess power
# n = 25 taxa
# simTree_25 <- sim.bdtree(n=25)
# simData_25 <- sim.char(phytools::rescale(simTree_25, "lambda", est_obs), 1)[,1,]
# bm_v_lambda_sim25_F_U <- pmc(simTree_25, simData_25, "BM", "lambda", nboot = 1000)
# save(bm_v_lambda_sim25_F_U, file="bm_v_lambda_sim25_F_U.Rda")
load("output/bm_v_lambda_sim25_F_U.Rda")

lambdas_sim25 <- bm_v_lambda_sim25_F_U$par_dists %>% filter(comparison=="BB", parameter=="lambda")

est_sim25 <- coef(bm_v_lambda_sim25_F_U[["B"]])[["lambda"]]
est_sim25

mean(lambdas_sim25$value)
cast(lambdas_sim25, comparison ~ parameter, function(x) quantile(x, c(0.025, 0.975)), value = c("lower", "upper"))

p2 <- ggplot(lambdas_sim25) +
  geom_histogram(aes(value), bins=50) +
  geom_vline(xintercept=est_sim25, linetype="dashed") +
  coord_cartesian(xlim=c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  labs(x="Estimated lambda", y="Count")
p2

# simTree_25_geiger <- treedata(simTree_25, simData_25)
# phylosig(simTree_25_geiger$phy, simTree_25_geiger$data, method = "lambda", test = TRUE)

# n = 50 taxa
# simTree_50 <- sim.bdtree(n=50)
# simData_50 <- sim.char(phytools::rescale(simTree_50, "lambda", est_obs), 1)[,1,]
# bm_v_lambda_sim50_F_U <- pmc(simTree_50, simData_50, "BM", "lambda", nboot = 1000)
# save(bm_v_lambda_sim50_F_U, file="bm_v_lambda_sim50F_U.Rda")
load("output/bm_v_lambda_sim50F_U.Rda")

lambdas_sim50 <- bm_v_lambda_sim50_F_U$par_dists %>% filter(comparison=="BB", parameter=="lambda")

est_sim50 <- coef(bm_v_lambda_sim50_F_U[["B"]])[["lambda"]]
est_sim50

mean(lambdas_sim50$value)
cast(lambdas_sim50, comparison ~ parameter, function(x) quantile(x, c(0.025, 0.975)), value = c("lower", "upper"))

p3 <- ggplot(lambdas_sim50) +
  geom_histogram(aes(value), bins=50) +
  geom_vline(xintercept=est_sim50, linetype="dashed") +
  coord_cartesian(xlim=c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  labs(x="Estimated lambda", y="Count")
p3

# simTree_50_geiger <- treedata(simTree_50, simData_50)
# phylosig(simTree_50_geiger$phy, simTree_50_geiger$data, method = "lambda", test = TRUE)

finalPlot <- ggarrange(p1, p2, p3, ncol=3)

ggsave("output/pmc_simulations_F_U.pdf", finalPlot, width=9, height=3, dpi=300, useDingbats=FALSE)

#----------Run analysis for uninfected MALES------------------------------------------
host.metabolism_phylo <- as.data.frame(lsmean_table_weight$emmean_M_U)
row.names(host.metabolism_phylo) <- row.names(lsmean_table_weight)
host.metabolism_geiger <- treedata(host.tree, host.metabolism_phylo)

#Using phytools
phylosig(host.metabolism_geiger$phy, host.metabolism_geiger$data, method = "lambda", test = TRUE)

#---------------------------------------------------------------------------------
#---------------------Plot summary data-------------------------------------------
#---------------------------------------------------------------------------------

#-----------Build tree-------------------------------------
host.tree <- read.nexus("data/phylogram_2_host.tre")

p5 <- ggtree(host.tree) +
        geom_tiplab(align=TRUE) +
        scale_x_continuous(expand=expansion(0.5)) + # make more room for the labels
        geom_treescale()
      # geom_text(aes(label=node))
p5

#make sure tips are in the correct order
p5 <- flip(p5,4,6)
p5 <- flip(p5,7,8)
p5

#-------Plot weight data------------------
dodge <- position_dodge(width = .5)
pd <- position_dodge(.3)

weight_means_melt

host.order <- c("CSBerk", "FFD25", "PC75", "R84","Car5", "sim198", "mau31", "sech",  
                   "yakB13L5", "san", "teiB13L11", "suz", "auraL2")
weight_means_melt$Genotype <- factor(weight_means_melt$Genotype, levels=host.order)
weight_means_melt$Sex <- factor(weight_means_melt$Sex, levels = c("M", "F"))

p6 <- ggplot(subset(weight_means_melt, !(Genotype %in% c("Car5", "sim198"))), 
             aes(Genotype, y=Mean, group=interaction(Genotype, Sex))) +
  geom_errorbar(aes(ymin=Lower,
                    ymax=Upper),
                width=0.3, alpha=1,
                position = dodge) +
  geom_point(aes(shape=Sex),
    colour="black", fill="gray80", size=2, position=dodge) +
  scale_shape_manual(values = c("F" = 21, "M" = 22)) +
  scale_x_discrete(limits=rev) + #reverse order of crosses
  scale_y_continuous(breaks=c(0.5,1,1.5)) +
  coord_flip() +
  guides(color = guide_legend(reverse=TRUE)) +
  ylab("mean weight (mg)")
p6

#-------Plot SMR data------------------
lsmean_table

lsmean_long <- lsmean_table |>
  pivot_longer(
    cols           = -c(Genotype, weight_model, Host),
    names_to       = c(".value", "group"),
    names_pattern  = "(emmean|lower\\.CL|upper\\.CL|SE)_([FM]_[UI])"
  ) |>
  separate(group, into = c("Sex", "Infected"), sep = "_")

lsmean_long$Genotype <- factor(lsmean_long$Genotype, levels=host.order)
lsmean_long_U <- subset(lsmean_long, Infected=="U")
lsmean_long_U$Sex <- factor(lsmean_long_U$Sex, levels = c("M", "F"))

p7 <- ggplot(subset(lsmean_long_U, !(Genotype %in% c("Car5", "sim198"))), 
             aes(Genotype, y=emmean, group=interaction(Genotype, Sex))) +
  geom_errorbar(aes(ymin=lower.CL,
                    ymax=upper.CL),
                width=0.3, alpha=1,
                position = dodge) +
  geom_point(aes(shape=Sex),
             colour="black", fill="gray80", size=2, position=dodge) +
  scale_shape_manual(values = c("F" = 21, "M" = 22)) +
  scale_x_discrete(limits=rev) + #reverse order of crosses
  coord_flip() +
  guides(color = guide_legend(reverse=TRUE)) +
  ylab("LS mean SMR (nmol min-1)")
p7


#-------Plot mass-adjusted SMR data----
lsmean_table_weight

lsmean_long_weight <- lsmean_table_weight |>
  pivot_longer(
    cols           = -c(Genotype, weight_model, Host),
    names_to       = c(".value", "group"),
    names_pattern  = "(emmean|lower\\.CL|upper\\.CL|SE)_([FM]_[UI])"
  ) |>
  separate(group, into = c("Sex", "Infected"), sep = "_")

lsmean_long_weight$Genotype <- factor(lsmean_long_weight$Genotype, levels=host.order)
lsmean_long_weight_U <- subset(lsmean_long_weight, Infected=="U")
lsmean_long_weight_U$Sex <- factor(lsmean_long_weight_U$Sex, levels = c("M", "F"))

p8 <- ggplot(subset(lsmean_long_weight_U, !(Genotype %in% c("Car5", "sim198"))), 
             aes(Genotype, y=emmean, group=interaction(Genotype, Sex))) +
  geom_errorbar(aes(ymin=lower.CL,
                    ymax=upper.CL),
                width=0.3, alpha=1,
                position = dodge) +
  geom_point(aes(shape=Sex),
             colour="black", fill="gray80", size=2, position=dodge) +
  scale_shape_manual(values = c("F" = 21, "M" = 22)) +
  scale_x_discrete(limits=rev) + #reverse order of crosses
  coord_flip() +
  guides(color = guide_legend(reverse=TRUE)) +
  ylab("mass-adjusted LS mean SMR (nmol min-1)")
p8

combined_plots <- ggarrange(p5, p6, p7, p8, nrow=1, widths = c(1,0.2,0.2,0.2))

ggsave("output/phylo_comparisons_metabolism.host.pdf", plot = combined_plots, width = 24, height = 7, dpi=300, useDingbats=FALSE)


