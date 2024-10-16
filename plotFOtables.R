library(forcats)
library(flextable)
library(readr)

rawData <- read_csv("~/Documents/MS_USA/Chapter_2/Data/rawData.csv", skip = 1)

LAGRHO = rawData %>%
  filter(`Species code` == "LAGRHO") %>%
  group_by(`Fish ID`,Year, Treatment, Group, Length, Site) %>%
  filter(Group %in% c('SAV', 'Amphipod', 'Crustacean',
                      'Polychaete', 'Isopod', 'Fish', 'Tanaidacea')) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = Group,
              values_from = count,
              values_fill = 0)

LAGRHO$siteTreat = factor(interaction(LAGRHO$Site, LAGRHO$Treatment, sep = "-"))


LAGRHOdf = data.frame(LAGRHO$siteTreat, LAGRHO$Length,
                      LAGRHO$Amphipod, LAGRHO$Polychaete, LAGRHO$Crustacean,
                      LAGRHO$Fish, LAGRHO$Tanaidacea, LAGRHO$SAV,
                      LAGRHO$Isopod)
names(LAGRHOdf) = c('site', 'TL', 'Amphipod', 'Polychaete', 'Crustacean', 'Fish',
                      'Tanaidacea', 'SAV', 'Isopod')
LAGRHOdf = LAGRHOdf %>%
  group_by(site) %>%
  filter(site %in% c("AM-CT", "DR-CT", "LB-CT", "NEPaP-CT", "HWP-LS", "NEPaP-LS", "SA-LS")) %>%
  summarize(TL_range = paste(min(TL), max(TL), sep = "-"),
            mean_TL    = paste(round(mean(TL, na.rm = T),1)),
            Amphipod   = sum(Amphipod),
            Polychaete = sum(Polychaete),
            Crustacean = sum(Crustacean),
            Fish       = sum(Fish),
            Tanaidacae = sum(Tanaidacea),
            SAV        = sum(SAV),
            Isopod     = sum(Isopod)) %>%
  mutate(N = paste(rowSums(.[4:10])), .after = "mean_TL") %>%
  mutate(site = fct_recode(site,
                           "AM-C" = "AM-CT",
                           "DR-C" = "DR-CT",
                           "LB-C" = "LB-CT",
                           "PaP-C" = "NEPaP-CT",
                           "HWP-R" = "HWP-LS",
                           "PaP-R" = "NEPaP-LS",
                           "SA-R" = "SA-LS"
                           ))

LAGRHOrawMat = as.matrix(LAGRHOdf[,5:11])

LAGRHOmat = matrix(NA, nrow = nrow(LAGRHOrawMat), ncol = ncol(LAGRHOrawMat))
for (i in 1:(nrow(LAGRHOmat))) {
  for (j in 1:(ncol(LAGRHOmat))) {
    LAGRHOmat[i,j] = (LAGRHOrawMat[i,j]/as.numeric(LAGRHOdf$N[i]))*100
  }
}

LAGRHOdf = data.frame(cbind(LAGRHOdf[,1:4], round(LAGRHOmat,1)))
names(LAGRHOdf) = c("Site", "Size range", "Mean size", "N", 'Amphipod', 'Polychaete', 'Crustacean', 'Fish',
                    'Tanaidacea', 'SAV', 'Isopod')

LAGRHO_FOtable = flextable(
  data     = LAGRHOdf,
  col_keys = c("Site", "Size range", "Mean size", "N", 'Amphipod', 'Polychaete', 'Crustacean', 'Fish',
               'Tanaidacea', 'SAV', 'Isopod')) %>%
  colformat_num(suffix = "%") %>%
  bold(part = "header", bold = TRUE) %>%
  align(align = "center", part = "all")

LAGRHO_FOtable

save_as_docx(LAGRHO_FOtable, path = "~/Documents/MS_USA/Chapter_2/Tables/LAGRHO_FO_table.docx")


#-------------------------------------------------------------------------------

MICUND = rawData %>%
  filter(`Species code` == "MICUND") %>%
  group_by(`Fish ID`,Year, Treatment, Group, Length, Site) %>%
  filter(Group %in% c('Amphipod', 'Crustacean',
                      'Polychaete', 'Isopod', 'Fish', 'Tanaidacea')) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = Group,
              values_from = count,
              values_fill = 0)

MICUND$siteTreat = factor(interaction(MICUND$Site, MICUND$Treatment, sep = "-"))


MICUNDdf = data.frame(MICUND$siteTreat, MICUND$Length,
                      MICUND$Amphipod, MICUND$Polychaete, MICUND$Crustacean,
                      MICUND$Fish, MICUND$Tanaidacea, MICUND$Isopod)
names(MICUNDdf) = c('site', 'TL', 'Amphipod', 'Polychaete', 'Crustacean', 'Fish',
                    'Tanaidacea', 'Isopod')

MICUNDdf = MICUNDdf %>%
  group_by(site) %>%
  filter(site %in% c("CI-CT", "NEPaP-CT", "NEPaP-LS", "SA-LS")) %>%
  summarize(TL_range = paste(min(TL, na.rm = T), max(TL, na.rm = T), sep = "-"),
            mean_TL    = paste(round(mean(TL, na.rm = T),1)),
            Amphipod   = sum(Amphipod),
            Polychaete = sum(Polychaete),
            Crustacean = sum(Crustacean),
            Fish       = sum(Fish),
            Tanaidacae = sum(Tanaidacea),
            Isopod     = sum(Isopod)) %>%
  mutate(N = paste(rowSums(.[4:9])), .after = "mean_TL") %>%
  mutate(site = fct_recode(site,
                           "CI-C" = "CI-CT",
                           "PaP-C" = "NEPaP-CT",
                           "PaP-R" = "NEPaP-LS",
                           "SA-R" = "SA-LS"
  ))

MICUNDrawMat = as.matrix(MICUNDdf[,5:10])

MICUNDmat = matrix(NA, nrow = nrow(MICUNDrawMat), ncol = ncol(MICUNDrawMat))
for (i in 1:(nrow(MICUNDmat))) {
  for (j in 1:(ncol(MICUNDmat))) {
    MICUNDmat[i,j] = (MICUNDrawMat[i,j]/as.numeric(MICUNDdf$N[i]))*100
  }
}

MICUNDdf = data.frame(cbind(MICUNDdf[,1:4], round(MICUNDmat,1)))
names(MICUNDdf) = c("Site", "Size range", "Mean size", "N", 'Amphipod', 'Polychaete', 'Crustacean', 'Fish',
                    'Tanaidacea', 'Isopod')

MICUND_FOtable = flextable(
  data     = MICUNDdf,
  col_keys = c("Site", "Size range", "Mean size", "N", 'Amphipod', 'Polychaete', 'Crustacean', 'Fish',
               'Tanaidacea', 'Isopod')) %>%
  colformat_num(suffix = "%") %>%
  bold(part = "header", bold = TRUE) %>%
  align(align = "center", part = "all")

MICUND_FOtable

save_as_docx(MICUND_FOtable, path = "~/Documents/MS_USA/Chapter_2/Tables/MICUND_FO_table.docx")
