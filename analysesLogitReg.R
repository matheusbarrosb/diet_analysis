# PACKAGES
library(readr)
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(forcats)
library(flextable)
library(tibble)

setwd("~/Documents/MS_USA/Chapter_2/Code")

# LOADING DATA
rawData <- read_csv("~/Documents/MS_USA/Chapter_2/Data/rawData.csv", skip = 1)

# LOAD FUNCTIONS
source("functions/theme_custom.R")
source("functions/plotOntoShifts.R")

### 1. PINFISH ###--------------------------------------------------------------
#### DATA WRANGLING --------------------------------------------------------------
LAGRHO = rawData %>%
  filter(`Species code` == "LAGRHO") %>%
  group_by(`Fish ID`,Year, Treatment, Group, Length, Site) %>%
  filter(Group %in% c('SAV', 'Amphipod', 'Crustacean',
                      'Polychaete', 'Isopod', 'Fish', 'Tanaidacea')) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = Group,
              values_from = count,
              values_fill = 0)

## GET PREY MATRIX
preyMat = as.matrix(LAGRHO[,6:ncol(LAGRHO)])
for (j in 1:ncol(preyMat)) {
  for (i in 1:nrow(preyMat))
    preyMat[i,j] = ifelse(preyMat[i,j] == 0, 0, 1)  
}

# order is Amphipod, Polychaete, Crustacean, Fish, Tanaidacea, SAV, Isopod

## GET CATEGORICAL COVARIATE
# Levels: AM-CT DR-CT LB-CT NEPaP-CT HWP-LS NEPaP-LS SA-LS
# Parameters to extract: p[1], p[3], p[6], p[7], p[10], p[11], p[12] 
LAGRHO$siteTreat = factor(interaction(LAGRHO$Site, LAGRHO$Treatment, sep = "-"))


# GET MEAN SIZE PER SITE
meanTLs_LAGRHO = LAGRHO %>%
  group_by(siteTreat) %>%
  summarise(meanTLs = mean(Length, na.rm = TRUE))


## MAKE STAN DATA
LAGRHO_amphipod = list(
  N       = length(LAGRHO$`Fish ID`), # number of observations
  S       = length(unique(LAGRHO$siteTreat)),
  site    = as.numeric(as.factor(LAGRHO$siteTreat)),
  TL      = as.numeric(LAGRHO$Length),
  y       = as.numeric(preyMat[,1]),
  meanTLs = meanTLs_LAGRHO$meanTLs
)

LAGRHO_polychaete = list(
  N    = length(LAGRHO$`Fish ID`), # number of observations
  S    = length(unique(LAGRHO$siteTreat)),
  site = as.numeric(as.factor(LAGRHO$siteTreat)),
  TL   = as.numeric(LAGRHO$Length),
  y    = as.numeric(preyMat[,2]),
  meanTLs = meanTLs_LAGRHO$meanTLs
)

LAGRHO_crustacean = list(
  N    = length(LAGRHO$`Fish ID`), # number of observations
  S    = length(unique(LAGRHO$siteTreat)),
  site = as.numeric(as.factor(LAGRHO$siteTreat)),
  TL   = as.numeric(LAGRHO$Length),
  y    = as.numeric(preyMat[,3]),
  meanTLs = meanTLs_LAGRHO$meanTLs
)

LAGRHO_fish = list(
  N    = length(LAGRHO$`Fish ID`), # number of observations
  S    = length(unique(LAGRHO$siteTreat)),
  site = as.numeric(as.factor(LAGRHO$siteTreat)),
  TL   = as.numeric(LAGRHO$Length),
  y    = as.numeric(preyMat[,4]),
  meanTLs = meanTLs_LAGRHO$meanTLs
)

LAGRHO_tanaid = list(
  N    = length(LAGRHO$`Fish ID`), # number of observations
  S    = length(unique(LAGRHO$siteTreat)),
  site = as.numeric(as.factor(LAGRHO$siteTreat)),
  TL   = as.numeric(LAGRHO$Length),
  y    = as.numeric(preyMat[,5]),
  meanTLs = meanTLs_LAGRHO$meanTLs
)

LAGRHO_SAV = list(
  N    = length(LAGRHO$`Fish ID`), # number of observations
  S    = length(unique(LAGRHO$siteTreat)),
  site = as.numeric(as.factor(LAGRHO$siteTreat)),
  TL   = as.numeric(LAGRHO$Length),
  y    = as.numeric(preyMat[,6]),
  meanTLs = meanTLs_LAGRHO$meanTLs
)

LAGRHO_isopod = list(
  N    = length(LAGRHO$`Fish ID`), # number of observations
  S    = length(unique(LAGRHO$siteTreat)),
  site = as.numeric(as.factor(LAGRHO$siteTreat)),
  TL   = as.numeric(LAGRHO$Length),
  y    = as.numeric(preyMat[,7]),
  meanTLs = meanTLs_LAGRHO$meanTLs
)

#### RUNNING MODELS ##------------------------------------------------------------
LAGRHO_amphipod_res = rstan::stan(file   = "stan/logitReg.stan",
                                  data   = LAGRHO_amphipod,
                                  warmup = 1000,
                                  iter   = 10000,
                                  chains = 3)

LAGRHO_polychaete_res = rstan::stan(file   = "stan/logitReg.stan",
                                    data   = LAGRHO_polychaete,
                                    warmup = 1000,
                                    iter   = 10000,
                                    chains = 3)

LAGRHO_crustacean_res = rstan::stan(file   = "stan/logitReg.stan",
                                    data   = LAGRHO_crustacean,
                                    warmup = 1000,
                                    iter   = 10000,
                                    chains = 3)

LAGRHO_fish_res = rstan::stan(file   = "stan/logitReg.stan",
                              data   = LAGRHO_fish,
                              warmup = 1000,
                              iter   = 10000,
                              chains = 3)

LAGRHO_tanaid_res = rstan::stan(file   = "stan/logitReg.stan",
                                data   = LAGRHO_tanaid,
                                warmup = 1000,
                                iter   = 10000,
                                chains = 3)

LAGRHO_SAV_res = rstan::stan(file   = "stan/logitReg.stan",
                             data   = LAGRHO_SAV,
                             warmup = 1000,
                             iter   = 10000,
                             chains = 3)

LAGRHO_isopod_res = rstan::stan(file   = "stan/logitReg.stan",
                                data   = LAGRHO_isopod,
                                warmup = 1000,
                                iter   = 10000,
                                chains = 3)

#### EXTRACTING PARAMETERS -----------------------------------------------------

require(ggmcmc)

LAGRHO_amphipod_params = ggs(LAGRHO_amphipod_res)
LAGRHO_amphipod_params = LAGRHO_amphipod_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[3]", "p[6]",
                          "p[7]", "p[10]", "p[11]", "p[12]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Amphipod")


LAGRHO_polychaete_params = ggs(LAGRHO_polychaete_res)
LAGRHO_polychaete_params = LAGRHO_polychaete_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[3]", "p[6]",
                          "p[7]", "p[10]", "p[11]", "p[12]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Polychaete")


LAGRHO_crustacean_params = ggs(LAGRHO_crustacean_res)
LAGRHO_crustacean_params = LAGRHO_crustacean_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[3]", "p[6]",
                          "p[7]", "p[10]", "p[11]", "p[12]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Crustacean")


LAGRHO_fish_params = ggs(LAGRHO_fish_res)
LAGRHO_fish_params = LAGRHO_fish_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[3]", "p[6]",
                          "p[7]", "p[10]", "p[11]", "p[12]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Fish")


LAGRHO_tanaid_params = ggs(LAGRHO_tanaid_res)
LAGRHO_tanaid_params = LAGRHO_tanaid_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[3]", "p[6]",
                          "p[7]", "p[10]", "p[11]", "p[12]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Tanaidacea")


LAGRHO_SAV_params = ggs(LAGRHO_SAV_res)
LAGRHO_SAV_params = LAGRHO_SAV_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[3]", "p[6]",
                          "p[7]", "p[10]", "p[11]", "p[12]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "SAV")


LAGRHO_isopod_params = ggs(LAGRHO_isopod_res)
LAGRHO_isopod_params = LAGRHO_isopod_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[3]", "p[6]",
                          "p[7]", "p[10]", "p[11]", "p[12]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Isopod")


LAGRHO_params = rbind(LAGRHO_amphipod_params,
                      LAGRHO_crustacean_params,
                      LAGRHO_fish_params,
                      LAGRHO_isopod_params,
                      LAGRHO_polychaete_params,
                      LAGRHO_SAV_params,
                      LAGRHO_tanaid_params)


#### PLOTTING PROBABILITIES OF ENCOUNTER ---------------------------------------
LAGRHO_PoE_PLOT = LAGRHO_params %>%
  filter(Parameter %in% c("p[1]", "p[3]", "p[6]",
                          "p[7]", "p[10]", "p[11]", "p[12]")) %>%
  mutate(Site = case_when(Parameter == "p[1]" ~ "AM-C",
                          Parameter == "p[3]" ~ "DR-C", 
                          Parameter =="p[6]" ~ "LB-C",
                          Parameter =="p[7]" ~ "PaP-C", 
                          Parameter =="p[10]" ~ "HWP-R",
                          Parameter =="p[11]" ~ "PaP-R",
                          Parameter =="p[12]" ~ "SA-R")) %>%
  mutate(Status = case_when(Site == "AM-C" ~ "Control",
                            Site == "DR-C" ~ "Control",
                            Site == "LB-C" ~ "Control",
                            Site == "PaP-C" ~ "Control",
                            Site == "HWP-R" ~ "Restored",
                            Site == "PaP-R" ~ "Restored",
                            Site == "SA-R" ~ "Restored")) %>%
  ggplot() +
  geom_errorbar(aes(x = Site,
                    ymin = mean - sd,
                    ymax = mean + sd), 
                width = 0,
                color = "black") +
  geom_point(aes(x = Site,
                 y = mean),
             color = "green4") +
  facet_wrap(~Prey,
             ncol = 2,
             scales = "free_y") +
  xlab("") +
  ylab("PoE (%)") +
  theme_custom() +
  theme(axis.text.x = 
          element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Pinfish")+
  scale_y_continuous(n.breaks = 4)

LAGRHO_PoE_PLOT

#### PLOTTING SIZE RELATIONSHIPS -----------------------------------------------
ontoShift_LAGRHO_PLOT = plotOntoShifts(
  outputs = list(LAGRHO_amphipod_res,
                 LAGRHO_polychaete_res,
                 LAGRHO_SAV_res,
                 LAGRHO_isopod_res),
  prey = list("Amphipods", "Polychaetes", "SAV", "Isopods"),
  minTL = min(LAGRHO$Length), maxTL = max(LAGRHO$Length),
  plotNcols = 4
)

ontoShift_LAGRHO_PLOT

#### PLOTTING PRIORS VS POSTERIORS ----------------------------------------------
resList   = list(LAGRHO_amphipod_res, LAGRHO_polychaete_res, LAGRHO_crustacean_res,
                 LAGRHO_fish_res, LAGRHO_tanaid_res, LAGRHO_SAV_res, LAGRHO_isopod_res)
priorList = c("alpha_hat", "beta_hat")
postList  = c("alpha", "beta")
parNames  = postList
preyNames = c("Amphipod", "Polychaete", "Crustacean", "Fish", "Tanaids", "SAV", "Isopod")

LAGRHO_postPrior_plot = plotPostPrior_logit(resList    = resList,
                                             priorList = priorList,
                                             postList  = postList,
                                             parNames  = parNames,
                                             preyNames = preyNames,
                                             plotTitle = "Pinfish")

ggsave(plot = LAGRHO_postPrior_plot, filename = "postPriorPlot_logit_LAGRHO.pdf",
       width = 4.1, height = 7.6,
       path = "~/Documents/MS_USA/Chapter_2/Figures/Supplemental")

### 2. CROAKER -----------------------------------------------------------------
#### DATA WRANGLING --------------------------------------------------------------
MICUND = rawData %>%
  filter(`Species code` == "MICUND") %>%
  group_by(`Fish ID`,Year, Treatment, Group, Length, Site) %>%
  filter(Group %in% c('Amphipod', 'Crustacean',
                      'Polychaete', 'Isopod', 'Fish', 'Tanaidacea')) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = Group,
              values_from = count,
              values_fill = 0)

MICUND = na.exclude(MICUND)

## GET PREY MATRIX
preyMat = as.matrix(MICUND[,6:ncol(MICUND)])
for (j in 1:ncol(preyMat)) {
  for (i in 1:nrow(preyMat))
    preyMat[i,j] = ifelse(preyMat[i,j] == 0, 0, 1)  
}

head(preyMat)

# order is Amphipod, Crustacean, Polychaete, Fish, Isopod, Tanaidacea

## GET CATEGORICAL COVARIATE
# Levels: CI-C, PaP-C, PaP-R, SA-R
# To extract: p[1], p[4], p[7], p[8]
MICUND$siteTreat = factor(interaction(MICUND$Site, MICUND$Treatment, sep = "-"))


# GET MEAN SIZE PER SITE
meanTLs_MICUND = MICUND %>%
  group_by(siteTreat) %>%
  summarise(meanTLs = mean(Length, na.rm = TRUE))

## MAKE STAN DATA
MICUND_amphipod = list(
  N       = length(MICUND$`Fish ID`), # number of observations
  S       = length(unique(MICUND$siteTreat)),
  site    = as.numeric(as.factor(MICUND$siteTreat)),
  TL      = as.numeric(MICUND$Length),
  y       = as.numeric(preyMat[,2]),
  meanTLs = meanTLs_MICUND$meanTLs
)

MICUND_crustacean = list(
  N       = length(MICUND$`Fish ID`), # number of observations
  S       = length(unique(MICUND$siteTreat)),
  site    = as.numeric(as.factor(MICUND$siteTreat)),
  TL      = as.numeric(MICUND$Length),
  y       = as.numeric(preyMat[,3]),
  meanTLs = meanTLs_MICUND$meanTLs
)

MICUND_polychaete = list(
  N    = length(MICUND$`Fish ID`), # number of observations
  S    = length(unique(MICUND$siteTreat)),
  site = as.numeric(as.factor(MICUND$siteTreat)),
  TL   = as.numeric(MICUND$Length),
  y    = as.numeric(preyMat[,1]),
  meanTLs = meanTLs_MICUND$meanTLs
)

MICUND_fish = list(
  N    = length(MICUND$`Fish ID`), # number of observations
  S    = length(unique(MICUND$siteTreat)),
  site = as.numeric(as.factor(MICUND$siteTreat)),
  TL   = as.numeric(MICUND$Length),
  y    = as.numeric(preyMat[,5]),
  meanTLs = meanTLs_MICUND$meanTLs
)

MICUND_isopod = list(
  N    = length(MICUND$`Fish ID`), # number of observations
  S    = length(unique(MICUND$siteTreat)),
  site = as.numeric(as.factor(MICUND$siteTreat)),
  TL   = as.numeric(MICUND$Length),
  y    = as.numeric(preyMat[,6]),
  meanTLs = meanTLs_MICUND$meanTLs
)

MICUND_tanaid = list(
  N    = length(MICUND$`Fish ID`), # number of observations
  S    = length(unique(MICUND$siteTreat)),
  site = as.numeric(as.factor(MICUND$siteTreat)),
  TL   = as.numeric(MICUND$Length),
  y    = as.numeric(preyMat[,4]),
  meanTLs = meanTLs_MICUND$meanTLs
)

#### RUNNING MODELS ##------------------------------------------------------------
MICUND_amphipod_res = rstan::stan(file   = "stan/logitReg.stan",
                                  data   = MICUND_amphipod,
                                  warmup = 1000,
                                  iter   = 10000,
                                  chains = 2)

MICUND_crustacean_res = rstan::stan(file   = "stan/logitReg.stan",
                                    data   = MICUND_crustacean,
                                    warmup = 1000,
                                    iter   = 10000,
                                    chains = 2)

MICUND_polychaete_res = rstan::stan(file   = "stan/logitReg.stan",
                                    data   = MICUND_polychaete,
                                    warmup = 1000,
                                    iter   = 10000,
                                    chains = 2)

MICUND_fish_res = rstan::stan(file   = "stan/logitReg.stan",
                              data   = MICUND_fish,
                              warmup = 1000,
                              iter   = 10000,
                              chains = 2)

MICUND_isopod_res = rstan::stan(file   = "stan/logitReg.stan",
                                data   = MICUND_isopod,
                                warmup = 1000,
                                iter   = 10000,
                                chains = 2)

MICUND_tanaid_res = rstan::stan(file   = "stan/logitReg.stan",
                                data   = MICUND_tanaid,
                                warmup = 1000,
                                iter   = 10000,
                                chains = 2)

#### EXTRACTING PARAMETERS -----------------------------------------------------

require(ggmcmc)

MICUND_amphipod_params = ggs(MICUND_amphipod_res)
MICUND_amphipod_params = MICUND_amphipod_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[4]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Amphipod")


MICUND_polychaete_params = ggs(MICUND_polychaete_res)
MICUND_polychaete_params = MICUND_polychaete_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta",  "p[1]", "p[4]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Polychaete")


MICUND_crustacean_params = ggs(MICUND_crustacean_res)
MICUND_crustacean_params = MICUND_crustacean_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[4]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Crustacean")


MICUND_fish_params = ggs(MICUND_fish_res)
MICUND_fish_params = MICUND_fish_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[4]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Fish")


MICUND_tanaid_params = ggs(MICUND_tanaid_res)
MICUND_tanaid_params = MICUND_tanaid_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[4]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Tanaidacea")


MICUND_isopod_params = ggs(MICUND_isopod_res)
MICUND_isopod_params = MICUND_isopod_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[1]", "p[4]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Isopod")


MICUND_params = rbind(MICUND_amphipod_params,
                      MICUND_crustacean_params,
                      MICUND_fish_params,
                      MICUND_isopod_params,
                      MICUND_polychaete_params,
                      MICUND_tanaid_params)


#### PLOTTING PROBABILITIES OF ENCOUNTER ---------------------------------------
# Levels: CI-C, PaP-C, PaP-R, SA-R

MICUND_PoE_PLOT = MICUND_params %>%
  filter(Parameter %in% c("p[1]", "p[4]", "p[7]",
                          "p[8]")) %>%
  mutate(Site = case_when(Parameter == "p[1]" ~ "CI-C",
                          Parameter == "p[4]" ~ "PaP-C", 
                          Parameter == "p[7]" ~ "PaP-R", 
                          Parameter == "p[8]" ~ "SA-R")) %>%
  ggplot() +
  geom_errorbar(aes(x = Site,
                    ymin = mean - sd,
                    ymax = mean + sd), width = 0, color = "black") +
  geom_point(aes(x = Site, y = mean), color = "green4") +
  facet_wrap(~Prey, ncol = 2, scales = "free_y") +
  xlab("") +
  ylab("") +
  theme_custom()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Croaker")+
  scale_y_continuous(n.breaks = 4)

MICUND_PoE_PLOT

#### PLOTTING SIZE RELATIONSHIPS -----------------------------------------------
micundPrey = list(
  "Amphipods", "Polychaetes",
  "Tanaids", "Other crustaceans"
  )

micundOutputs = list(
  MICUND_amphipod_res,
  MICUND_polychaete_res,
  MICUND_tanaid_res,
  MICUND_crustacean_res)

ontoShift_MICUND_PLOT = plotOntoShifts(
  outputs = micundOutputs,
  prey    = micundPrey,
  minTL   = min(MICUND$Length), maxTL = max(MICUND$Length),
  plotNcols = 4
)

ontoShift_MICUND_PLOT

#### PLOTTING PRIORS VS POSTERIORS ----------------------------------------------
resList   = list(MICUND_amphipod_res, MICUND_polychaete_res, MICUND_crustacean_res,
                 MICUND_fish_res, MICUND_tanaid_res, MICUND_isopod_res)
priorList = c("alpha_hat", "beta_hat")
postList  = c("alpha", "beta")
parNames  = postList
preyNames = c("Amphipod", "Polychaete", "Crustacean", "Fish", "Tanaids", "Isopod")

MICUND_postPrior_plot = plotPostPrior_logit(resList   = resList,
                                            priorList = priorList,
                                            postList  = postList,
                                            parNames  = parNames,
                                            preyNames = preyNames,
                                            plotTitle = "Croaker")

ggsave(plot = MICUND_postPrior_plot, filename = "postPriorPlot_logit_MICUND.pdf",
       width = 4.1, height = 7.6,
       path = "~/Documents/MS_USA/Chapter_2/Figures/Supplemental")

### 3. SILVER PERCH ------------------------------------------------------------
#### DATA WRANGLING ------------------------------------------------------------
# Levels: CI-C, PaP-C, PaP-R, SA-R

BAICHR = rawData %>%
  filter(`Species code` == "BAICHR") %>%
  group_by(`Fish ID`,Year, Treatment, Group, Length, Site) %>%
  filter(Group %in% c('Amphipod', 'Crustacean',
                      'Polychaete', 'Isopod', 'Fish', 'Tanaidacea')) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = Group,
              values_from = count,
              values_fill = 0)

BAICHR = na.exclude(BAICHR)

# order is Crustacean Amphipod Isopod Fish Polychaete Tanaidacea

## GET CATEGORICAL COVARIATE
# Levels: DR-CT NEPaP-CT HWP-LS NEPaP-LS SA-LS
# to extract: 2, 5, 7, 8
BAICHR$siteTreat = factor(interaction(BAICHR$Site, BAICHR$Treatment, sep = "-"))

BAICHR = droplevels(BAICHR[!BAICHR$siteTreat == "SA-DG",])
BAICHR = droplevels(BAICHR[!BAICHR$siteTreat == "CI-LS",])

# GET MEAN SIZE PER SITE
meanTLs_BAICHR = BAICHR %>%
  group_by(siteTreat) %>%
  summarise(meanTLs = mean(Length, na.rm = TRUE))


## GET PREY MATRIX
preyMat = as.matrix(BAICHR[,6:ncol(BAICHR)])
for (j in 1:ncol(preyMat)) {
  for (i in 1:nrow(preyMat))
    preyMat[i,j] = ifelse(preyMat[i,j] == 0, 0, 1)  
}

head(preyMat)


## MAKE STAN DATA
BAICHR_crustacean = list(
  N    = length(BAICHR$`Fish ID`), # number of observations
  S    = length(unique(BAICHR$siteTreat)),
  site = as.numeric(as.factor(BAICHR$siteTreat)),
  TL   = as.numeric(BAICHR$Length),
  y    = as.numeric(preyMat[,1]),
  meanTLs = meanTLs_BAICHR$meanTLs
)

BAICHR_amphipod = list(
  N    = length(BAICHR$`Fish ID`), # number of observations
  S    = length(unique(BAICHR$siteTreat)),
  site = as.numeric(as.factor(BAICHR$siteTreat)),
  TL   = as.numeric(BAICHR$Length),
  y    = as.numeric(preyMat[,2]),
  meanTLs = meanTLs_BAICHR$meanTLs
)

BAICHR_isopod = list(
  N    = length(BAICHR$`Fish ID`), # number of observations
  S    = length(unique(BAICHR$siteTreat)),
  site = as.numeric(as.factor(BAICHR$siteTreat)),
  TL   = as.numeric(BAICHR$Length),
  y    = as.numeric(preyMat[,3]),
  meanTLs = meanTLs_BAICHR$meanTLs
)

BAICHR_fish = list(
  N    = length(BAICHR$`Fish ID`), # number of observations
  S    = length(unique(BAICHR$siteTreat)),
  site = as.numeric(as.factor(BAICHR$siteTreat)),
  TL   = as.numeric(BAICHR$Length),
  y    = as.numeric(preyMat[,4]),
  meanTLs = meanTLs_BAICHR$meanTLs
)

BAICHR_polychaete = list(
  N    = length(BAICHR$`Fish ID`), # number of observations
  S    = length(unique(BAICHR$siteTreat)),
  site = as.numeric(as.factor(BAICHR$siteTreat)),
  TL   = as.numeric(BAICHR$Length),
  y    = as.numeric(preyMat[,5]),
  meanTLs = meanTLs_BAICHR$meanTLs
)

BAICHR_tanaid = list(
  N    = length(BAICHR$`Fish ID`), # number of observations
  S    = length(unique(BAICHR$siteTreat)),
  site = as.numeric(as.factor(BAICHR$siteTreat)),
  TL   = as.numeric(BAICHR$Length),
  y    = as.numeric(preyMat[,6]),
  meanTLs = meanTLs_BAICHR$meanTLs
)

#### RUNNING MODELS ##------------------------------------------------------------
BAICHR_amphipod_res = rstan::stan(file   = "stan/logitReg.stan",
                                  data   = BAICHR_amphipod,
                                  warmup = 500,
                                  iter   = 3000,
                                  chains = 2)

BAICHR_crustacean_res = rstan::stan(file   = "stan/logitReg.stan",
                                    data   = BAICHR_crustacean,
                                    warmup = 500,
                                    iter   = 3000,
                                    chains = 2)

BAICHR_polychaete_res = rstan::stan(file   = "stan/logitReg.stan",
                                    data   = BAICHR_polychaete,
                                    warmup = 500,
                                    iter   = 3000,
                                    chains = 2)

BAICHR_fish_res = rstan::stan(file   = "stan/logitReg.stan",
                              data   = BAICHR_fish,
                              warmup = 500,
                              iter   = 3000,
                              chains = 2)

BAICHR_isopod_res = rstan::stan(file   = "stan/logitReg.stan",
                                data   = BAICHR_isopod,
                                warmup = 500,
                                iter   = 3000,
                                chains = 2)

BAICHR_tanaid_res = rstan::stan(file   = "stan/logitReg.stan",
                                data   = BAICHR_tanaid,
                                warmup = 500,
                                iter   = 3000,
                                chains = 2)

#### EXTRACTING PARAMETERS -----------------------------------------------------

require(ggmcmc)

BAICHR_amphipod_params = ggs(BAICHR_amphipod_res)
BAICHR_amphipod_params = BAICHR_amphipod_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[2]", "p[5]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Amphipod")


BAICHR_polychaete_params = ggs(BAICHR_polychaete_res)
BAICHR_polychaete_params = BAICHR_polychaete_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[2]", "p[5]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Polychaete")


BAICHR_crustacean_params = ggs(BAICHR_crustacean_res)
BAICHR_crustacean_params = BAICHR_crustacean_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[2]", "p[5]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Crustacean")


BAICHR_fish_params = ggs(BAICHR_fish_res)
BAICHR_fish_params = BAICHR_fish_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[2]", "p[5]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Fish")


BAICHR_tanaid_params = ggs(BAICHR_tanaid_res)
BAICHR_tanaid_params = BAICHR_tanaid_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[2]", "p[5]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Tanaidacea")


BAICHR_isopod_params = ggs(BAICHR_isopod_res)
BAICHR_isopod_params = BAICHR_isopod_params %>%
  group_by(Parameter) %>%
  filter(Parameter %in% c("beta", "p[2]", "p[5]", "p[7]",
                          "p[8]")) %>%
  summarise(mean = mean(value),
            sd   = sd(value)) %>% 
  mutate(Prey = "Isopod")


BAICHR_params = rbind(BAICHR_amphipod_params,
                      BAICHR_crustacean_params,
                      BAICHR_fish_params,
                      BAICHR_isopod_params,
                      BAICHR_polychaete_params,
                      BAICHR_tanaid_params)


#### PLOTTING PROBABILITIES OF ENCOUNTER ---------------------------------------
# Plot order: CI-C, PaP-C, PaP-R, SA-R

BAICHR_PoE_PLOT = BAICHR_params %>%
  filter(Parameter %in% c("p[2]", "p[5]", "p[7]",
                          "p[8]")) %>%
  mutate(Site = case_when(Parameter == "p[2]" ~ "CI-C",
                          Parameter == "p[5]" ~ "PaP-C", 
                          Parameter =="p[7]" ~ "PaP-R",
                          Parameter =="p[8]" ~ "SA-R")) %>%
  ggplot() +
  geom_errorbar(aes(x = Site,
                    ymin = mean - sd,
                    ymax = mean + sd), width = 0, color = "black") +
  geom_point(aes(x = Site, y = mean), color = "green4") +
  facet_wrap(~Prey, ncol = 2, scales = "free_y") +
  xlab("") +
  ylab("") +
  theme_custom() + ggtitle("Silver perch")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(n.breaks = 4)

BAICHR_PoE_PLOT

#### PLOTTING SIZE RELATIONSHIPS -----------------------------------------------
baichrPrey    = list("Amphipods", "Fish", "Tanaids", "Other crustaceans")
baichrOutputs = list(BAICHR_amphipod_res, BAICHR_fish_res,
                     BAICHR_tanaid_res, BAICHR_crustacean_res) 

ontoShift_BAICHR_PLOT = plotOntoShifts(
  outputs   = baichrOutputs,
  prey      = baichrPrey,
  minTL     = min(BAICHR$Length),
  maxTL     = max(BAICHR$Length),
  plotNcols = 4
)

ontoShift_BAICHR_PLOT

#### PLOTTING PRIORS VS POSTERIORS ----------------------------------------------
resList   = list(BAICHR_amphipod_res, BAICHR_polychaete_res, BAICHR_crustacean_res,
                 BAICHR_fish_res, BAICHR_tanaid_res, BAICHR_isopod_res)
priorList = c("alpha_hat", "beta_hat")
postList  = c("alpha", "beta")
parNames  = postList
preyNames = c("Amphipod", "Polychaete", "Crustacean", "Fish", "Tanaids", "Isopod")

BAICHR_postPrior_plot = plotPostPrior_logit(resList   = resList,
                                            priorList = priorList,
                                            postList  = postList,
                                            parNames  = parNames,
                                            preyNames = preyNames,
                                            plotTitle = "Silver perch")

ggsave(plot = BAICHR_postPrior_plot, filename = "postPriorPlot_logit_BAICHR.pdf",
       width = 4.1, height = 7.6,
       path = "~/Documents/MS_USA/Chapter_2/Figures/Supplemental")

# 4. MERGING PLOTS -------------------------------------------------------------
# ontogenetic shift plots
ontoShift_merged = cowplot::plot_grid(
                   ontoShift_LAGRHO_PLOT + theme(axis.title.x = element_blank(),
                                                 axis.title.y = element_blank()),
                   ontoShift_MICUND_PLOT + theme(axis.title.x = element_blank()),
                   ontoShift_BAICHR_PLOT + theme(axis.title.y = element_blank()),
                   nrow = 3, align = "hv", labels = "auto",
                   axis = "lr", label_size = 16)

ggsave(plot = ontoShift_merged, filename = "ontoShift.pdf",
       width = 7.5, height = 6.5, path = "~/Documents/MS_USA/Chapter_2/Figures")


# PoE plots
plot01 = ggarrange(MICUND_PoE_PLOT, BAICHR_PoE_PLOT, ncol = 1)
PoE_plot = ggarrange(LAGRHO_PoE_PLOT, plot01)

ggsave(plot = PoE_plot, filename = "PoE.pdf",
       width = 8, height = 8, path = "~/Documents/MS_USA/Chapter_2/Figures")


# 5. POSTERIOR PREDICTIVE CHECK ------------------------------------------------

# Parameters to extract: p[1], p[3], p[6], p[7], p[10], p[11], p[12] 
LAGRHO_params = LAGRHO_params %>%
  filter(Parameter %in% c("p[1]", "p[3]", "p[6]",
                          "p[7]", "p[10]", "p[11]", "p[12]")) %>%
  mutate(Site = case_when(Parameter == "p[1]"  ~ "AM-C",
                          Parameter == "p[3]"  ~ "DR-C", 
                          Parameter == "p[6]"  ~ "LB-C",
                          Parameter == "p[7]"  ~ "PaP-C", 
                          Parameter == "p[10]" ~ "HWP-R",
                          Parameter == "p[11]" ~ "PaP-R",
                          Parameter == "p[12]" ~ "SA-R")) %>%
  mutate(Species = "Pinfish")

# To extract: p[1], p[4], p[7], p[8]
MICUND_params = MICUND_params %>%
  filter(Parameter %in% c("p[1]", "p[4]", "p[7]",
                          "p[8]")) %>%
  mutate(Site = case_when(Parameter == "p[1]" ~ "CI-C",
                          Parameter == "p[4]" ~ "PaP-C", 
                          Parameter == "p[7]" ~ "PaP-R", 
                          Parameter == "p[8]" ~ "SA-R")) %>%
  mutate(Species = "Croaker")

# to extract: 2, 5, 7, 8
BAICHR_params = BAICHR_params %>% 
  filter(Parameter %in% c("p[2]", "p[5]", "p[7]",
                          "p[8]")) %>%
  mutate(Site = case_when(Parameter == "p[2]" ~ "CI-C",
                          Parameter == "p[5]" ~ "PaP-C", 
                          Parameter == "p[7]" ~ "PaP-R", 
                          Parameter == "p[8]" ~ "SA-R"))%>%
  mutate(Species = "Silver perch")

PoE_DF = rbind(LAGRHO_params, MICUND_params, BAICHR_params)

# get observed PoEs

preyMat_LAGRHO = LAGRHO[,6:12]
for(i in 1:ncol(preyMat_LAGRHO)) {
  for (j in 1:nrow(preyMat_LAGRHO)) {
    preyMat_LAGRHO[j,i] = ifelse(preyMat_LAGRHO[j,i] > 1, 1, preyMat_LAGRHO[j,i])
  }
}

preyMat_MICUND = MICUND[,6:11]
for(i in 1:ncol(preyMat_MICUND)) {
  for (j in 1:nrow(preyMat_MICUND)) {
    preyMat_MICUND[j,i] = ifelse(preyMat_MICUND[j,i] > 1, 1, preyMat_MICUND[j,i])
  }
}

preyMat_BAICHR = BAICHR[,6:11]
for(i in 1:ncol(preyMat_BAICHR)) {
  for (j in 1:nrow(preyMat_BAICHR)) {
    preyMat_BAICHR[j,i] = ifelse(preyMat_BAICHR[j,i] > 1, 1, preyMat_BAICHR[j,i])
  }
}


preyMat_LAGRHO$siteTreat = LAGRHO$siteTreat
preyMat_BAICHR$siteTreat = BAICHR$siteTreat
preyMat_MICUND$siteTreat = MICUND$siteTreat

observed_LAGRHO = preyMat_LAGRHO %>%
  group_by(siteTreat) %>%
  summarise(Amphipod = sum(Amphipod)/n(),
            Polychaete = sum(Polychaete)/n(),
            Crustacean = sum(Crustacean)/n(),
            Fish = sum(Fish)/n(),
            Tanaidacea = sum(Tanaidacea)/n(),
            Isopod = sum(Isopod)/n(),
            SAV = sum(SAV)/n())%>% 
  filter(!siteTreat %in% c("CI-CT", "GB-CT",
                           "HB-CT", "WFR-CT",
                           "CI-LS")) %>%
  pivot_longer(!siteTreat, names_to = "Prey",
               values_to = "Observed") %>%
  mutate(Species = "Pinfish")

observed_MICUND = preyMat_MICUND %>%
  group_by(siteTreat) %>%
  filter(siteTreat %in% c("CI-CT", "NEPaP-CT", "NEPaP-LS", "SA-LS")) %>%
  summarise(Amphipod = sum(Amphipod)/n(),
            Polychaete = sum(Polychaete)/n(),
            Crustacean = sum(Crustacean)/n(),
            Isopod = sum(Isopod)/n(),
            Fish = sum(Fish)/n(),
            Tanaidacea = sum(Tanaidacea)/n()) %>%
  pivot_longer(!siteTreat, names_to = "Prey",
               values_to = "Observed") %>%
  mutate(Species = "Croaker")

observed_BAICHR = preyMat_BAICHR %>%
  group_by(siteTreat) %>%
  filter(siteTreat %in% c("CI-CT", "NEPaP-CT", "NEPaP-LS", "SA-LS")) %>%
summarise(Amphipod = sum(Amphipod)/n(),
          Polychaete = sum(Polychaete)/n(),
          Crustacean = sum(Crustacean)/n(),
          Fish = sum(Fish)/n(),
          Isopod = sum(Isopod)/n(),
          Tanaidacea = sum(Tanaidacea)/n()) %>%
  pivot_longer(!siteTreat, names_to = "Prey",
               values_to = "Observed") %>%
  mutate(Species = "Silver perch")

# observed_LAGRHO = LAGRHO %>%
#   group_by(siteTreat) %>%
#   filter(siteTreat %in% c("AM-CT", "DR_CT", "LB-CT", "NEPaP-CT", "HWP-LS", "NEPaP-LS", 
#                           "SA-LS", "CI-CT")) %>%
#   summarize(SAV = sum(SAV, na.rm = T)/n(),
#             Crustacean = sum(Crustacean, na.rm = T)/n(),
#             Amphipod = sum(Amphipod, na.rm = T)/n(),
#             Isopod = sum(Isopod, na.rm = T)/n(),
#             Fish = sum(Fish, na.rm = T)/n(),
#             Polychaete = sum(Polychaete, na.rm = T)/n(),
#             Tanaidacea = sum(Tanaidacea, na.rm = T)/n()) %>% 
#   pivot_longer(!siteTreat, names_to = "Prey",
#                values_to = "Observed") %>%
#   mutate(Species = "Pinfish")
# 
# 
# observed_MICUND = MICUND %>%
#   group_by(siteTreat) %>%
#   filter(siteTreat %in% c("NEPaP-CT", "NEPaP-LS", 
#                           "SA-LS", "CI-CT")) %>%
#   summarize(Crustacean = sum(Crustacean, na.rm = T)/n(),
#             Amphipod = sum(Amphipod, na.rm = T)/n(),
#             Isopod = sum(Isopod, na.rm = T)/n(),
#             Fish = sum(Fish, na.rm = T)/n(),
#             Polychaete = sum(Polychaete, na.rm = T)/n(),
#             Tanaidacea = sum(Tanaidacea, na.rm = T)/n()) %>%
#   pivot_longer(!siteTreat, names_to = "Prey",
#                values_to = "Observed") %>%
#   mutate(Species = "Croaker")
# 
# 
# observed_BAICHR = BAICHR %>%
#   group_by(siteTreat) %>%
#   filter(siteTreat %in% c("NEPaP-CT", "NEPaP-LS", 
#                           "SA-LS", "CI-CT")) %>%
#   summarize(Crustacean = sum(Crustacean, na.rm = T)/n(),
#             Amphipod = sum(Amphipod, na.rm = T)/n(),
#             Isopod = sum(Isopod, na.rm = T)/n(),
#             Fish = sum(Fish, na.rm = T)/n(),
#             Polychaete = sum(Polychaete, na.rm = T)/n(),
#             Tanaidacea = sum(Tanaidacea, na.rm = T)/n()) %>%
#   pivot_longer(!siteTreat, names_to = "Prey",
#                values_to = "Observed") %>%
#   mutate(Species = "Silver perch")


observed_merged = rbind(observed_LAGRHO, observed_MICUND, observed_BAICHR) %>% 
  mutate(site = case_when(siteTreat == "AM-CT" ~ "AM-C",
                          siteTreat == "CI-CT" ~ "CI-C",
                          siteTreat == "DR-CT" ~ "DR-C",
                          siteTreat == "LB-CT" ~ "LB-C",
                          siteTreat == "NEPaP-CT" ~ "PaP-C",
                          siteTreat == "HWP-LS" ~ "HWP-R",
                          siteTreat == "NEPaP-LS" ~ "PaP-R",
                          siteTreat == "SA-LS" ~ "SA-R")) %>% 
  na.exclude()
  
head(observed_merged)


estimated_merged = tibble(PoE_DF$Site, PoE_DF$Prey, PoE_DF$mean,
                              PoE_DF$Species, PoE_DF$sd)
names(estimated_merged) = c("site", "prey", "POE", "species", "sd", "category")
head(estimated_merged)

# ORDER AND CHECK IF TWO DATAFRAMES MATCH
estimated_merged = estimated_merged[with(estimated_merged, order(site)),]
observed_merged  = observed_merged[with(observed_merged, order(site)),]

estimated_merged = estimated_merged[with(estimated_merged, order(prey)),]
observed_merged  = observed_merged[with(observed_merged, order(Prey)),]

estimated_merged = estimated_merged[with(estimated_merged, order(species)),]
observed_merged  = observed_merged[with(observed_merged, order(Species)),]


PPCHECK_DF = data.frame(estimated_merged$site, observed_merged$site, estimated_merged$prey,
                        estimated_merged$species, estimated_merged$POE,
                        observed_merged$Observed, estimated_merged$sd)
names(PPCHECK_DF) = c("site", "site2", "prey", "species", "estimated", "observed", "sd")
PPCHECK_DF$diff = abs(PPCHECK_DF$estimated - PPCHECK_DF$observed)

# plot with prey coding to highlight deviances are due to ontogenetic shifts

PPCHECK_DF %>%
  ggplot() +
  geom_point(aes(x = observed,
                 y = estimated,
                 color = prey,
                 alpha = 0.7), size = 3) +
  facet_wrap(~species) +
  geom_abline(slope=1, linetype = "dashed", color="Red") +
  xlim(0,1) + ylim(0,1) + theme_custom() +
  xlab("Observed FO (%)") +
  ylab("PoE (%)") +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  scale_color_brewer(palette = "Set1") +
  guides(size = "none",
         alpha = "none")


ggsave("ppcheck.pdf", width = 9, height = 5, path = "~/Documents/MS_USA/Chapter_2/Figures")

# 6. FO % TABLES ---------------------------------------------------------------

# create base dataframes
LAGRHO_FO_df = PPCHECK_DF %>% 
  group_by(prey, site, species) %>%
  filter(species == "Pinfish") %>%
  tidyr::pivot_wider(id_cols = !c(site2, estimated, diff, sd),
                     names_from = 'prey',
                     values_from = 'observed')

LAGRHO_meanTLs = LAGRHO %>%
  group_by(siteTreat) %>%
  filter(siteTreat %in% c("AM-CT", "DR-CT", "HWP-LS", "LB-CT",
                          "NEPaP-CT", "NEPaP-LS", "SA-LS")) %>%
  mutate(siteTreat = factor(siteTreat,
                            levels = c("AM-CT", "DR-CT", "HWP-LS", "LB-CT",
                                       "NEPaP-CT", "NEPaP-LS", "SA-LS"))) %>%
  summarise(mean   = round(mean(Length, na.rm = T),1),
            range  = paste0(min(Length, na.rm = T), "-",
                            max(Length, na.rm = T)),
            N      = n())

LAGRHO_FO_df2 = data.frame(
  LAGRHO_FO_df$site,
  LAGRHO_meanTLs$mean,
  LAGRHO_meanTLs$range,
  LAGRHO_meanTLs$N,
  round(LAGRHO_FO_df$Amphipod,2),
  round(LAGRHO_FO_df$Crustacean,2),
  round(LAGRHO_FO_df$Fish,2),
  round(LAGRHO_FO_df$Isopod,2),
  round(LAGRHO_FO_df$Polychaete,2),
  round(LAGRHO_FO_df$SAV,2),
  round(LAGRHO_FO_df$Tanaidacea,2)
)

names(LAGRHO_FO_df2) = c("Site", "Mean_TL", "TL_range", "N", 
                        colnames(LAGRHO_FO_df[,3:ncol(LAGRHO_FO_df)]))

print(LAGRHO_FO_df2)


MICUND_FO_df = PPCHECK_DF %>% 
  group_by(prey, site, species) %>%
  filter(species == "Croaker") %>%
  tidyr::pivot_wider(id_cols = !c(site2, estimated, diff, sd),
                     names_from = 'prey',
                     values_from = 'observed')

MICUND_meanTLs = MICUND %>%
  group_by(siteTreat) %>%
  filter(siteTreat %in% c("CI-CT", "NEPaP-CT", "NEPaP-LS", "SA-LS")) %>%
  mutate(siteTreat = factor(siteTreat,
                            levels = c("CI-CT", "NEPaP-CT", "NEPaP-LS", "SA-LS"))) %>%
  summarise(mean   = round(mean(Length, na.rm = T),1),
            range  = paste0(min(Length, na.rm = T), "-",
                            max(Length, na.rm = T)),
            N      = n())

MICUND_FO_df2 = data.frame(
  MICUND_FO_df$site,
  MICUND_meanTLs$mean,
  MICUND_meanTLs$range,
  MICUND_meanTLs$N,
  round(MICUND_FO_df$Amphipod,2),
  round(MICUND_FO_df$Crustacean,2),
  round(MICUND_FO_df$Fish,2),
  round(MICUND_FO_df$Isopod,2),
  round(MICUND_FO_df$Polychaete,2),
  round(MICUND_FO_df$Tanaidacea,2)
)

names(MICUND_FO_df2) = c("Site", "Mean_TL", "TL_range", "N",
                         colnames(MICUND_FO_df[,3:ncol(MICUND_FO_df)]))

print(MICUND_FO_df2)

BAICHR_FO_df = PPCHECK_DF %>% 
  group_by(prey, site, species) %>%
  filter(species == "Silver perch") %>%
  tidyr::pivot_wider(id_cols = !c(site2, estimated, diff, sd),
                     names_from = 'prey',
                     values_from = 'observed')

BAICHR_meanTLs = BAICHR %>%
  group_by(siteTreat) %>%
  filter(siteTreat %in% c("CI-CT", "NEPaP-CT", "NEPaP-LS", "SA-LS")) %>%
  mutate(siteTreat = factor(siteTreat,
                            levels = c("CI-CT", "NEPaP-CT", "NEPaP-LS", "SA-LS"))) %>%
  summarise(mean   = round(mean(Length, na.rm = T),1),
            range  = paste0(min(Length, na.rm = T), "-",
                            max(Length, na.rm = T)),
            N      = n())

BAICHR_FO_df2 = data.frame(
  BAICHR_FO_df$site,
  BAICHR_meanTLs$mean,
  BAICHR_meanTLs$range,
  BAICHR_meanTLs$N,
  round(BAICHR_FO_df$Amphipod,2),
  round(BAICHR_FO_df$Crustacean,2),
  round(BAICHR_FO_df$Fish,2),
  round(BAICHR_FO_df$Isopod,2),
  round(BAICHR_FO_df$Polychaete,2),
  round(BAICHR_FO_df$Tanaidacea,2)
)

names(BAICHR_FO_df2) = c("Site", "Mean TL", "TL range", "N",
                         colnames(BAICHR_FO_df[,3:ncol(BAICHR_FO_df)]))

# create tables
LAGRHO_FO_table = flextable(
  data     = LAGRHO_FO_df2,
  col_keys = c("Site", "Mean", "range", "N",
               'Amphipod', 'Crustacean', 'Fish', 'Isopod',
               'Polychaete', 'SAV', 'Tanaidacea')) %>%
  bold(part = "header", bold = TRUE) %>%
  align(align = "center", part = "all")

MICUND_FO_table = flextable(
  data     = MICUND_FO_df2,
  col_keys = c("Site", "Mean", "range", "N",
               'Amphipod', 'Crustacean', 'Fish', 'Isopod',
               'Polychaete', 'Tanaidacea')) %>%
  bold(part = "header", bold = TRUE) %>%
  align(align = "center", part = "all")

BAICHR_FO_table = flextable(
  data     = BAICHR_FO_df2,
  col_keys = c("Site", "Mean", "range", "N",
               'Amphipod', 'Crustacean', 'Fish', 'Isopod',
               'Polychaete', 'Tanaidacea')) %>%
  bold(part = "header", bold = TRUE) %>%
  align(align = "center", part = "all")


save_as_docx(LAGRHO_FO_table, path = "~/Documents/MS_USA/Chapter_2/Tables/LAGRHO_FO_table.docx")
save_as_docx(MICUND_FO_table, path = "~/Documents/MS_USA/Chapter_2/Tables/MICUND_FO_table.docx")
save_as_docx(BAICHR_FO_table, path = "~/Documents/MS_USA/Chapter_2/Tables/BAICHR_FO_table.docx")

# 7. STACKED AREA CHART TO SHOW ONTOGENETIC SHIFTS -----------------------------

# PINFISH ##
LAGRHO_AREA = LAGRHO %>%
  pivot_longer(cols = 6:12, names_to = "Prey", values_to = "Ocurrence") %>%
  mutate(Ints = cut(as.numeric(Length), breaks = c(30,40,50,60,70,80,90,170,180), right = TRUE)) %>%
  group_by(Prey, Ints) %>%
  summarise(TotalOcc = sum(Ocurrence, na.rm = TRUE),
            n = n())

LAGRHO_AREA$Midpoints = rep(c(35,45,55,65,75,85,130,175,200), length(unique(LAGRHO_AREA$Prey)))

LAGRHO_areaPlot = LAGRHO_AREA %>%
  group_by(as.factor(Midpoints)) %>%
  mutate(PercOcc = TotalOcc/n) %>%
  
  ggplot(aes(x = Midpoints, y = PercOcc, fill = Prey)) +
  geom_area(alpha = 0.7) + 
  theme_custom() + 
  scale_fill_brewer(palette = "Set1") +
  xlab("Total length (mm)") +
  ylab("FO %") + xlim(35,130) +
  ggtitle("Pinfish")

# CROAKER ##
MICUND_AREA = MICUND %>%
  pivot_longer(cols = 6:11, names_to = "Prey", values_to = "Ocurrence") %>%
  mutate(Ints = cut(as.numeric(Length), breaks = c(25, 50, 75, 100, 125, 150, 200), right = TRUE)) %>%
  group_by(Prey, Ints) %>%
  summarise(TotalOcc = sum(Ocurrence, na.rm = TRUE),
            n = n())

MICUND_AREA$Midpoints = rep(c(37, 62, 87, 112, 175, 200),length(unique(MICUND_AREA$Prey)))

MICUND_areaPlot = MICUND_AREA %>%
  group_by(as.factor(Midpoints)) %>%
  mutate(PercOcc = TotalOcc/n) %>%
  
  ggplot(aes(x = Midpoints, y = PercOcc, fill = Prey)) +
  geom_area(alpha = 0.7) + 
  theme_custom() + 
  scale_fill_brewer(palette = "Set1") +
  xlab("Total length (mm)") +
  ylab("FO %") + ggtitle("Croaker")

# SILVER PERCH ##

BAICHR_AREA = BAICHR %>%
  pivot_longer(cols = 6:11, names_to = "Prey", values_to = "Ocurrence") %>%
  mutate(Ints = cut(as.numeric(Length), breaks = c(20, 40, 60, 80, 100, 120, 130), right = TRUE)) %>%
  group_by(Prey, Ints) %>%
  summarise(TotalOcc = sum(Ocurrence, na.rm = TRUE),
            n = n())

BAICHR_AREA$Midpoints = rep(c(30, 50, 70, 90, 120),length(unique(BAICHR_AREA$Prey)))

BAICHR_areaPlot = BAICHR_AREA %>%
  group_by(as.factor(Midpoints)) %>%
  mutate(PercOcc = TotalOcc/n) %>%
  
  ggplot(aes(x = Midpoints, y = PercOcc, fill = Prey)) +
  geom_area(alpha = 0.7) + 
  theme_custom() + 
  scale_fill_brewer(palette = "Set1") +
  xlab("Total length (mm)") +
  ylab("FO %") + ggtitle("Silver perch")

###### MERGING

area_merged = cowplot::plot_grid(
  LAGRHO_areaPlot + theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          legend.position = "top",
                          legend.title = element_blank(),
                          legend.text = element_text(size = 8),
                          legend.key.size = unit(0.4, 'cm')),
  MICUND_areaPlot + theme(axis.title.x = element_blank(),
                          legend.position = "none"),
  BAICHR_areaPlot + theme(axis.title.y = element_blank(),
                          legend.position = "none"),
  nrow = 3, align = "v", axis = "lr",
  rel_heights = c(0.38,0.31,0.31))

ggsave(plot = area_merged, filename = "rawOntoShift.pdf",
       width = 5.5, height = 7.58, path = "~/Documents/MS_USA/Chapter_2/Figures/Supplemental")
