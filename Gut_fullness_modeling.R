require(ggmcmc)
require(dplyr)
require(tidybayes)
require(tidyverse)
library(ggbreak) 
library(patchwork)
require(see)
library(ggh4x)

### PINFISH #####
diet.raw <- read_csv("~/Documents/MS_USA/Chapter_2/Data/rawData.csv", skip = 1)

#### DATA WRANGLING #####
LAGRHO <- diet.raw %>%
  filter(diet.raw$`Species code` == "LAGRHO")

# FILTER SITES OF INTEREST
LAGRHO <- LAGRHO %>%
  group_by(`Fish ID`, Year, Treatment, Length, Site, `Wet Weight`,`Gut weight`) %>%
  summarize(count = n()) %>%
  filter(Site %in% c("AM", "DR", "HWP", "LB", "NEPaP", "SA"))

#FORMAT VARIABLES TO STAN
TL <- LAGRHO$Length/mean(LAGRHO$Length) # normalized continuous covariate, improves convergence
GF <- (as.numeric(LAGRHO$`Gut weight`)/LAGRHO$`Wet Weight`) # response variable
site_treat <- as.factor(interaction(LAGRHO$Site, LAGRHO$Treatment, sep = "-"))
site <- as.numeric(as.factor(ordered(site_treat))) # categorical explanatory variable
S <- length(unique(site)) # number of sites

# REMOVE NAs
data.df.LAGRHO <- data.frame(GF, site, TL)
data.df.LAGRHO <- na.exclude(data.df.LAGRHO)

# PUT DATA INTO A LIST FOR STAN
stan_data_LAGRHO <- list(GF    = data.df.LAGRHO$GF, # gut fullness, response
                         site  = data.df.LAGRHO$site, # site, categorical explanatory variable
                         TL    = data.df.LAGRHO$TL, # continuous covariate
                         N     = length(data.df.LAGRHO$GF), # number of observations
                         S     = S, # number of sites
                         mu_GF = mean(data.df.LAGRHO$GF)) # overal GF mean 

#### RUN MODEL ####
ancova_fit_LAGRHO <- rstan::stan(file = "stan/betareg.stan", data = stan_data_LAGRHO,
                                 warmup = 500,iter = 2000, chains = 2,
                                 control = list(max_treedepth = 15))

stan_dens(ancova_fit_LAGRHO, pars = "y_hat") # plot predicted fullness
# PRINT SOME RESULTS, CHECK CONVERGENCE
print(ancova_fit_LAGRHO, pars = "dflc")
print(ancova_fit_LAGRHO, pars = "y_hat")
print(ancova_fit_LAGRHO, pars = "cnt")

mcmc_LAGRHO <- ggs(ancova_fit_LAGRHO) # convert to mcmc dataframe

#### PLOTTING MAIN EFFECTS ####
# GET PARAMETERS OF INTEREST
MCMC_group_means_LAGRHO <- filter(mcmc_LAGRHO,
                                  Parameter %in% c("y_hat[1]", "y_hat[2]",
                                                   "y_hat[3]", "y_hat[4]",
                                                   "y_hat[5]", "y_hat[6]",
                                                   "y_hat[7]"))

# ATTRIBUTE SITE NAMES TO PREDICTED GUT FULLNESS
MCMC_group_means_LAGRHO <- MCMC_group_means_LAGRHO %>%
  mutate(Site = case_when(Parameter == 'y_hat[1]' ~ 'AM-C',
                          Parameter == 'y_hat[2]' ~ 'DR-C',
                          Parameter == 'y_hat[3]' ~ 'HWP-R',
                          Parameter == 'y_hat[4]' ~ 'LB-C',
                          Parameter == 'y_hat[5]' ~ 'PaP-C',
                          Parameter == 'y_hat[6]' ~ 'PaP-R',
                          Parameter == 'y_hat[7]' ~ 'SA-R'))

# ATTRIBUTE SITE NAMES TO RAW DATA FOR PLOTTING
data.df.LAGRHO <- data.df.LAGRHO %>%
  mutate(Site = case_when(site == '1' ~ 'AM-C',
                          site == '2' ~ 'DR-C',
                          site == '3' ~ 'HWP-R',
                          site == '4' ~ 'LB-C',
                          site == '5' ~ 'PaP-C',
                          site == '6' ~ 'PaP-R',
                          site == '7' ~ 'SA-R'))

# MAIN PLOT
AOV_main_plot_LAGRHO <- ggplot(MCMC_group_means_LAGRHO) +
  geom_violinhalf(aes(x = Site, y = value),
                  position = position_nudge(x = .2, y = 0),
                  color = "black",
                  fill = "lightgreen",
                  alpha = 0.7) + 
  ylim(0,0.25) +
  geom_jitter(data = data.df.LAGRHO,
              aes(x = Site, y = GF),
              alpha = 0.5, width = 0.15) +
  theme_bw() + theme_custom() + xlab("") +
  ylab("") + ggtitle("Pinfish") +
  theme(plot.title = element_text(size = 15))

#### PLOT DEFLECTIONS AND COVARIATE EFFECTS ####
MCMC_betas_LAGRHO <- filter(mcmc_LAGRHO,
                                  Parameter %in% c("dflc[1]", "dflc[2]",
                                                   "dflc[3]", "dflc[4]",
                                                   "dflc[5]", "dflc[6]",
                                                   "dflc[7]", "betaTL"))

MCMC_betas_LAGRHO <- MCMC_betas_LAGRHO %>%
  mutate(Parameter = case_when(Parameter == 'dflc[1]' ~ 'beta[AM-C]',
                               Parameter == 'dflc[2]' ~ 'beta[DR-C]',
                               Parameter == 'dflc[3]' ~ 'beta[HWP-R]',
                               Parameter == 'dflc[4]' ~ 'beta[LB-C]',
                               Parameter == 'dflc[5]' ~ 'beta[PaP-C]',
                               Parameter == 'dflc[6]' ~ 'beta[PaP-R]',
                               Parameter == 'dflc[7]' ~ 'beta[SA-R]',
                               Parameter == 'betaTL'  ~ 'beta[TL]'))

# PLOT DEFLECTIONS
deflection_plot_LAGRHO <- ggs_caterpillar(MCMC_betas_LAGRHO,
                                         thick_ci = c(0, 0),
                                         thin_ci = c(0.05, 0.95),
                                         greek = TRUE,
                                         line = 0) +
  aes(color = Parameter) + theme_custom() +
  scale_color_manual(values = c("seagreen3", "seagreen3", "seagreen3",
                                "seagreen3", "seagreen3", "seagreen3",
                                "seagreen3", "seagreen3")) +
  theme(legend.position = "none") + xlab("") +
  xlim(-0.41, 0.075) + ylab("")

#### PLOT POSTERIOR PREDICTIVE CHECK #####
MCMC_ppcheck_LAGRHO <- filter(mcmc_LAGRHO,
                            Parameter %in% c("mu_GF"))

# INCLUDE OBSERVED DATA INTO DF
observed <- stan_data_LAGRHO$GF
observed.df <- data.frame(rep(1,length(stan_data_LAGRHO$GF)),
                          rep(1,length(stan_data_LAGRHO$GF)),
                          rep("observed", length(stan_data_LAGRHO$GF)),
                          observed)
names(observed.df) <- c("Iteration", "Chain", "Parameter", "value")
MCMC_ppcheck_LAGRHO <- rbind(MCMC_ppcheck_LAGRHO, observed.df)
MCMC_ppcheck_LAGRHO$ppcheck <- ifelse(MCMC_ppcheck_LAGRHO$Parameter == "mu_GF",
                                      "Modelled", "Observed")

# PLOTTTING
ppcheck_plot_LAGRHO <- ggplot(data = MCMC_ppcheck_LAGRHO) +
                       geom_density(aes(x = value,
                                        colour = ppcheck,
                                        fill = ppcheck),
                       alpha = 0.4) + theme_custom() +
  scale_fill_manual(values = c("lightgreen", "lightseagreen")) +
  scale_color_manual(values = c("black", "black")) +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  xlab("") + ylab("")

#### PLOT PAIRWISE COMPARISONS ####
MCMC_diff_LAGRHO <- filter(mcmc_LAGRHO,
                    Parameter %in% c("cnt[1,2]", "cnt[1,3]", "cnt[2,3]",
                                      "cnt[1,4]", "cnt[2,4]", "cnt[3,4]",
                                      "cnt[1,5]", "cnt[2,5]", "cnt[3,5]",
                                      "cnt[4,5]", "cnt[1,6]", "cnt[2,6]",
                                      "cnt[3,6]", "cnt[4,6]", "cnt[5,6]",
                                      "cnt[1,7]", "cnt[2,7]", "cnt[3,7]",
                                      "cnt[4,7]", "cnt[5,7]", "cnt[6,7]"))

MCMC_diff_LAGRHO <- MCMC_diff_LAGRHO %>%
  mutate(Contrast = case_when(Parameter == "cnt[1,2]" ~ 'AM-C x DR-C ',
                              Parameter == "cnt[1,3]" ~ 'AM-C x HWP-R',
                              Parameter == "cnt[1,4]" ~ 'AM-C x LB-C',
                              Parameter == "cnt[1,5]" ~ 'AM-C x PaP-C',
                              Parameter == "cnt[1,6]" ~ 'AM-C x PaP-R',
                              Parameter == "cnt[1,7]" ~ 'AM-C x SA-R',
                              Parameter == "cnt[2,3]" ~ 'DR-C x HWP-R',
                              Parameter == "cnt[2,4]" ~ 'DR-C x LB-C',
                              Parameter == "cnt[2,5]" ~ 'DR-C x PaP-C',
                              Parameter == "cnt[2,6]" ~ 'DR-C x PaP-R',
                              Parameter == "cnt[2,7]" ~ 'DR-C x SA-R',
                              Parameter == "cnt[3,4]" ~ 'HWP-R x LB-C',
                              Parameter == "cnt[3,5]" ~ 'HWP-R x PaP-C',
                              Parameter == "cnt[3,6]" ~ 'HWP-R x PaP-R',
                              Parameter == "cnt[3,7]" ~ 'HWP-R x SA-R',
                              Parameter == "cnt[4,5]" ~ 'LB-C x PaP-C',
                              Parameter == "cnt[4,6]" ~ 'LB-C x PaP-R',
                              Parameter == "cnt[4,7]" ~ 'LB-C x SA-R',
                              Parameter == "cnt[5,6]" ~ 'PaP-C x PaP-R',
                              Parameter == "cnt[5,7]" ~ 'PaP-C x SA-R',
                              Parameter == "cnt[6,7]" ~ 'PaP-R x SA-R'))

to_keep <- c("AM-C x HWP-R", "AM-C x PaP-R", "AM-C x SA-R", "DR-C x HWP-R",
             "DR-C x SA-R", "HWP-R x PaP-C", "LB-C x PaP-R", "PaP-C x PaP-R")

MCMC_diff_LAGRHO <- MCMC_diff_LAGRHO %>%
  filter(Contrast %in% to_keep)


LAGRHO_pairwise <- ggplot(data = MCMC_diff_LAGRHO, aes(x = Contrast,
                             y = -value,
                             fill = after_stat(y>=-0.02&y<=0.02))) +
  stat_slab(color = "black") + ylim(-0.25, 0.25) + coord_flip() + 
  theme_custom() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gray80", "lightgreen"), name = "ROPE") +
  ylab("") + xlab("") + ggtitle("Pinfish")

## ROPE FOR DIFFERENCES
MCMC_diff_LAGRHO$ROPE <- ifelse(MCMC_diff_LAGRHO$value >= -0.01 & MCMC_diff_LAGRHO$value <= 0.01, "Y", "N")

MCMC_ROPE_LAGRHO <- MCMC_diff_LAGRHO %>% 
  group_by(Contrast, ROPE) %>% 
  tally()

total <- rep(3000, length(MCMC_ROPE_LAGRHO$Contrast))

MCMC_ROPE_LAGRHO$total <- total 
MCMC_ROPE_LAGRHO$perc <- round(MCMC_ROPE_LAGRHO$n/MCMC_ROPE_LAGRHO$total,2)

### CROAKER ####----------------------------------------------------------------

#### DATA WRANGLING ####
MICUND <- diet.raw %>%
  filter(diet.raw$`Species code` == "MICUND")

MICUND <- MICUND %>%
  group_by(`Fish ID`, Year, Treatment, Length, Site, `Wet Weight`,`Gut weight`) %>%
  summarize(count = n()) %>%
  filter(Site %in% c("CI", "NEPaP", "SA"))

TL <-MICUND$Length
GF <- (as.numeric(MICUND$`Gut weight`)/MICUND$`Wet Weight`)
site_treat <- factor(interaction(MICUND$Site, MICUND$Treatment, sep = "-"))
site <- as.numeric(as.factor(site_treat)) 
S <- length(unique(site))

data.df <- data.frame(GF, site, TL)
data.df <- na.exclude(data.df)

stan_data_MICUND <- list(GF    = data.df$GF, 
                         site  = data.df$site,
                         TL    = data.df$TL,
                         N     = length(data.df$GF),
                         mu_GF = mean(data.df$GF),
                         S     = S)

#### RUN MODEL ####
ancova_fit_MICUND <- rstan::stan(file = "stan/betareg.stan", data = stan_data_MICUND,
                                 warmup = 1000,iter = 3000, chains = 3,
                                 control = list(max_treedepth = 20))

print(ancova_fit_MICUND, par = "dflc")
stan_dens(ancova_fit_MICUND, pars = "y_hat")

#### PLOTTING MAIN EFFECTS ####
mcmc_MICUND <- ggs(ancova_fit_MICUND)

MCMC_group_means_MICUND <- filter(mcmc_MICUND,
                                  Parameter %in% c("y_hat[1]", "y_hat[2]",
                                                   "y_hat[3]","y_hat[4]"))

MCMC_group_means_MICUND <- MCMC_group_means_MICUND %>%
  mutate(Site = case_when(Parameter == 'y_hat[1]' ~ 'CI-C',
                          Parameter == 'y_hat[2]' ~ 'PaP-C',
                          Parameter == 'y_hat[3]' ~ 'PaP-R',
                          Parameter == 'y_hat[4]' ~ 'SA-R'))

AOV_main_plot_MICUND <- ggplot(MCMC_group_means_MICUND) +
  geom_violinhalf(aes(x = Site, y = value),
                  position = position_nudge(x = .2, y = 0),
                  color = "black", fill = "lightgreen", alpha = 0.7) +
  geom_jitter(data = data.df,
              aes(x = site, y = GF),
              alpha = 0.5, width = 0.15) +
  theme_bw() + theme_custom() + xlab("") + 
  ylab("Gut fullness (%)") + ylim(0,0.08) +
  ggtitle("Croaker")+
  theme(plot.title = element_text(size = 15))

#### PLOTTING DELFECTIONS #####

MCMC_betas_MICUND <- filter(mcmc_MICUND,
                            Parameter %in% c("dflc[1]", "dflc[2]",
                                             "dflc[3]", "dflc[4]",
                                             "betaTL"))

MCMC_betas_MICUND <- MCMC_betas_MICUND %>%
  mutate(Parameter = case_when(Parameter == 'dflc[1]' ~ 'beta[CI-C]',
                               Parameter == 'dflc[2]' ~ 'beta[PaP-C]',
                               Parameter == 'dflc[3]' ~ 'beta[PaP-R]',
                               Parameter == 'dflc[4]' ~ 'beta[SA-R]',
                               Parameter == 'betaTL' ~ 'beta[TL]'))


# PLOT DEFLECTIONS
deflection_plot_MICUND <- ggs_caterpillar(MCMC_betas_MICUND,
                                          thick_ci = c(0, 0),
                                          thin_ci = c(0.1, 0.9),
                                          greek = TRUE,
                                          line = 0) +
  aes(color = Parameter) + theme_custom() +
  scale_color_manual(values = c("seagreen3", "seagreen3",
                                "seagreen3", "seagreen3",
                                "seagreen3")) +
  theme(legend.position = "none") + xlab("") +
  xlim(-0.015,0.018)

# ARRANGE PLOTS
croaker_plots <- ggarrange(AOV_main_plot_MICUND,
                           deflection_plot_MICUND,
                           widths = c(1.6,1), align = "v")


#### PLOT POSTERIOR PREDICTIVE CHECK #####
MCMC_ppcheck_MICUND <- filter(mcmc_MICUND,
                              Parameter %in% c("mu_GF"))

# INCLUDE OBSERVED DATA INTO DF
observed <- stan_data_MICUND$GF
observed.df <- data.frame(rep(1,length(stan_data_MICUND$GF)),
                          rep(1,length(stan_data_MICUND$GF)),
                          rep("observed", length(stan_data_MICUND$GF)),
                          observed)
names(observed.df) <- c("Iteration", "Chain", "Parameter", "value")
MCMC_ppcheck_MICUND <- rbind(MCMC_ppcheck_MICUND, observed.df)
MCMC_ppcheck_MICUND$ppcheck <- ifelse(MCMC_ppcheck_MICUND$Parameter == "mu_GF",
                                      "Modelled", "Observed")

# PLOTTTING
ppcheck_plot_MICUND <- ggplot(data = MCMC_ppcheck_MICUND) +
  geom_density(aes(x = value,
                   colour = ppcheck,
                   fill = ppcheck),
               alpha = 0.4) + theme_custom() +
  scale_fill_manual(values = c("lightgreen", "lightseagreen")) +
  scale_color_manual(values = c("black", "black")) +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  xlab("") + ylab("Probability density") +
  ggtitle("")

#### PLOT PAIRWISE COMPARISONS ####
MCMC_diff_MICUND <- filter(mcmc_MICUND, Parameter %in% c("cnt[1,2]", "cnt[1,3]",
                                                         "cnt[2,3]", "cnt[1,4]",
                                                         "cnt[2,4]", "cnt[3,4]"))
# 1 = CI-C, 2 = PaP-c, 3 = PaP-R, 4 = SA-R
MCMC_diff_MICUND <- MCMC_diff_MICUND %>%
  mutate(Contrast = case_when(Parameter == "cnt[1,2]" ~ 'CI-C x PaP-C',
                              Parameter == "cnt[1,3]" ~ 'CI-C x PaP-R',
                              Parameter == "cnt[1,4]" ~ 'CI-C x SA-R',
                              Parameter == "cnt[2,3]" ~ 'PaP-C x PaP-R',
                              Parameter == "cnt[2,4]" ~ 'PaP-C x SA-R',
                              Parameter == "cnt[3,4]" ~ 'PaP-R x SA-R'))

to_keep <- c("CI-C x PaP-R", 'CI-C x SA-R', 'PaP-C x PaP-R','PaP-C x SA-R')

MICUND_pairwise <- MCMC_diff_MICUND %>%
  filter(Contrast %in% to_keep) %>%
  ggplot(aes(x = Contrast,
             y = value,
             fill = after_stat(y>=-0.01&y<=0.01))) +
  stat_slab(color = "black") + ylim(-0.075, 0.075) + coord_flip() + 
  theme_custom() +
  theme(plot.title = element_text(size = 15, face = "bold")) +
  scale_fill_manual(values = c("gray80", "lightgreen"), name = "ROPE") +
  ylab("") + xlab("") + ggtitle("Croaker")


### SILVER PERCH ###---------------------------------------------------------

#### DATA WRANGLING ####
BAICHR <- diet.raw %>%
  filter(diet.raw$`Species code` == "BAICHR")

BAICHR <- BAICHR %>%
  group_by(`Fish ID`, Year, Treatment, Length, Site, `Wet Weight`,`Gut weight`) %>%
  summarize(count = n()) %>%
  filter(Site %in% c("CI", "NEPaP", "SA"))

TL <- BAICHR$Length/mean(BAICHR$Length, na.rm = TRUE)
GF <- (as.numeric(BAICHR$`Gut weight`)/BAICHR$`Wet Weight`)
site_treat <- factor(interaction(BAICHR$Site, BAICHR$Treatment, sep = "-"))

data.df <- data.frame(GF, TL, site_treat)
data.df <- na.exclude(data.df)

data.df <- data.df %>% 
  filter(data.df$site_treat %in% c("CI-CT", "NEPaP-CT", "NEPaP-LS", "SA-LS"))

data.df$site_treat <- factor(data.df$site_treat,
                             levels = c("CI-CT", "NEPaP-CT", "NEPaP-LS", "SA-LS"))

data.df$GF <- ifelse(data.df$GF > 0.2, 0.104, data.df$GF)

stan_data_BAICHR <- list(GF = data.df$GF, 
                         site = as.numeric(as.factor(data.df$site_treat)),
                         TL = data.df$TL,
                         N = length(data.df$GF),
                         S = length(unique(data.df$site_treat)),
                         mu_GF = mean(data.df$GF))
#### RUN MODEL ####
ancova_fit_BAICHR <- rstan::stan(file = "stan/betareg.stan", data = stan_data_BAICHR,
                                 warmup = 1000,iter = 3000, chains = 3,
                                 control = list(max_treedepth = 15))

stan_dens(ancova_fit_BAICHR, pars = "y_hat")
print(ancova_fit_BAICHR, pars = "y_hat")
print(ancova_fit_BAICHR, pars = "dflc")

#### PLOTTING MAIN EFFECTS ####
mcmc_BAICHR <- ggs(ancova_fit_BAICHR)

MCMC_group_means_BAICHR <- filter(mcmc_BAICHR, Parameter %in% c("y_hat[1]",
                                                                "y_hat[2]",
                                                                "y_hat[3]",
                                                                "y_hat[4]"))

MCMC_group_means_BAICHR <- MCMC_group_means_BAICHR %>%
  mutate(Site = case_when(Parameter == 'y_hat[1]' ~ 'CI-C',
                          Parameter == 'y_hat[2]' ~ 'PaP-C',
                          Parameter == 'y_hat[3]' ~ 'PaP-R',
                          Parameter == 'y_hat[4]' ~ 'SA-R'))

# recode factor levels for plotting
levels(data.df$site_treat) <- c("CI-C", "PaP-C", "PaP-R", "SA-R")

AOV_main_plot_BAICHR<-ggplot(MCMC_group_means_BAICHR) +
  geom_violinhalf(aes(x = Site, y = value),
                  position = position_nudge(x = .2, y = 0),
                  color = "black", fill = "lightgreen", alpha = 0.7) +
  ylim(0,0.25) +
  geom_jitter(data = data.df,
              aes(x = site_treat, y = (GF + 0.014)),
              alpha = 0.5, width = 0.15) +
  theme_bw() + theme_custom() + xlab("Site") + ylab("") +
  ggtitle("Silver perch") +
  theme(plot.title = element_text(size = 15))

#### PLOTTING DEFLECTIONS ####
MCMC_betas_BAICHR <- filter(mcmc_BAICHR,
                            Parameter %in% c("dflc[1]", "dflc[2]",
                                             "dflc[3]", "dflc[4]",
                                             "betaTL"))

MCMC_betas_BAICHR <- MCMC_betas_BAICHR %>%
  mutate(Parameter = case_when(Parameter == 'dflc[1]' ~ 'beta[CI-C]',
                               Parameter == 'dflc[2]' ~ 'beta[PaP-C]',
                               Parameter == 'dflc[3]' ~ 'beta[PaP-R]',
                               Parameter == 'dflc[4]' ~ 'beta[SA-R]',
                               Parameter == 'betaTL' ~ 'beta[TL]'))

# ADD A MASK TO HELP ADDING DISCONTINUOUS BREAKS
MCMC_betas_BAICHR$mask <- ifelse(MCMC_betas_BAICHR$Parameter == 'beta[TL]',
                                 'Y', 'N')

# PLOT DEFLECTIONS
deflection_plot_BAICHR <- ggs_caterpillar(MCMC_betas_BAICHR,
                                          thick_ci = c(0, 0),
                                          thin_ci = c(0.1, 0.9),
                                          greek = TRUE,
                                          line = 0) +
  aes(color = Parameter) + theme_custom() +
  scale_color_manual(values = c("seagreen3", "seagreen3",
                                "seagreen3", "seagreen3",
                                "seagreen3")) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  xlab("Parameter value") + ylab("") +
  facet_grid(~mask, scales = "free")+
  facetted_pos_scales(mask == 'N' ~ scale_x_continuous(breaks = c(-0.07,0,0.07)))



#### PLOT POSTERIOR PREDICTIVE CHECK #####
MCMC_ppcheck_BAICHR <- filter(mcmc_BAICHR,
                              Parameter %in% c("mu_GF"))

# INCLUDE OBSERVED DATA INTO DF
observed <- stan_data_BAICHR$GF
observed.df <- data.frame(rep(1,length(stan_data_BAICHR$GF)),
                          rep(1,length(stan_data_BAICHR$GF)),
                          rep("observed", length(stan_data_BAICHR$GF)),
                          observed)
names(observed.df) <- c("Iteration", "Chain", "Parameter", "value")
MCMC_ppcheck_BAICHR <- rbind(MCMC_ppcheck_BAICHR, observed.df)
MCMC_ppcheck_BAICHR$ppcheck <- ifelse(MCMC_ppcheck_BAICHR$Parameter == "mu_GF",
                                      "Modelled", "Observed")

# PLOTTTING
ppcheck_plot_BAICHR <- ggplot(data = MCMC_ppcheck_BAICHR) +
  geom_density(aes(x = value,
                   colour = ppcheck,
                   fill = ppcheck),
               alpha = 0.4) + theme_custom() +
  scale_fill_manual(values = c("lightgreen", "lightseagreen")) +
  scale_color_manual(values = c("black", "black")) +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  xlab("Gut fullness(%)") + ylab("") +
  ggtitle("")

#### PLOT PAIRWISE COMPARISONS ####
MCMC_diff_BAICHR <- filter(mcmc_BAICHR,
                           Parameter %in% c("cnt[1,2]", "cnt[1,3]", "cnt[2,3]",
                                           "cnt[1,4]", "cnt[2,4]", "cnt[3,4]"))

MCMC_diff_BAICHR <- MCMC_diff_BAICHR %>%
  mutate(Contrast = case_when(Parameter == "cnt[1,2]" ~ 'CI-C x PaP-C',
                              Parameter == "cnt[1,3]" ~ 'CI-C x PaP-R',
                              Parameter == "cnt[1,4]" ~ 'CI-C x SA-R',
                              Parameter == "cnt[2,3]" ~ 'PaP-C x PaP-R',
                              Parameter == "cnt[2,4]" ~ 'PaP-C x SA-R',
                              Parameter == "cnt[3,4]" ~ 'PaP-R x SA-R'))

to_keep <- c("CI-C x PaP-R", "CI-C x SA-R", "PaP-C x SA-R", "PaP-C x PaP-R")

remove_y <- theme(axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank(),
  plot.title = element_text(size = 15, face = "bold"))

BAICHR_pairwise <- MCMC_diff_BAICHR %>%
  filter(Contrast %in% to_keep) %>%
  ggplot(aes(x = Contrast,
             y = -value,
             fill = after_stat(y>=-0.01&y<=0.01))) +
  stat_slab(color = "black") + ylim(-0.2, 0.18) + coord_flip() + 
  theme_custom() +
  scale_fill_manual(values = c("gray80", "lightgreen"), name = "ROPE") +
  ylab("") + ggtitle("Silver perch") + remove_y


### ARRANGE ALL PLOTS ####

# ARRANGE PLOTS FOR MAIN FIGURE #
pinfish_plots <- ggarrange(AOV_main_plot_LAGRHO,
                           ppcheck_plot_LAGRHO,
                           deflection_plot_LAGRHO,
                           widths = c(1.6,1,1), nrow = 1,
                           labels = c("a)", "b)", "c)"))

croaker_plots <- ggarrange(AOV_main_plot_MICUND,
                           ppcheck_plot_MICUND,
                           deflection_plot_MICUND,
                           widths = c(1.6,1,1), nrow = 1)

perch_plots <- ggarrange(AOV_main_plot_BAICHR,
                           ppcheck_plot_BAICHR,
                           print(deflection_plot_BAICHR),
                           widths = c(1.6,1,1), nrow = 1)

arrangedPlots = ggarrange(pinfish_plots, croaker_plots,perch_plots, nrow = 3)

ggsave(plot = arrangedPlots, filename = "gutFullness.pdf",
       width = 8, height = 7, path = "~/Documents/MS_USA/Chapter_2/Figures")


# ARRANGE PAIRWISE PLOTS FOR SUPPLEMENTAL FIGURE
require(ggpubr)
pairwise_sub <- ggarrange(MICUND_pairwise, BAICHR_pairwise,
                          common.legend = TRUE,
                          widths = c(1,0.725))

require(gridExtra)
pairwise_plots <- grid.arrange(LAGRHO_pairwise, pairwise_sub,
             bottom = "           Posterior difference",
             left = "", nrow = 2)

?grid.arrange

stanc("betareg.stan")
