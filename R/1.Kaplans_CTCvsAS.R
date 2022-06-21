# Author:                      Job van Riet
# Date:                        21-06-2022
# Function:                    Compare performance of dichotomized CTC vs. mFAST-SeqS aneuploidy scores.

# Import libraries ----

library(dplyr)
library(ggplot2)
library(survminer)
library(gtsummary)
library(extrafont)
library(patchwork)

# Helper themes and functions.
source('R/misc_Themes.R')
source('R/misc_Functions.R')


# Import data. ----

data.Patient <- list()
data.Patient$metrics <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Sample Overview') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$clinicalData <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Clinical Characteristics') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))


# Convert and clean data. ----

data.Patient$survival <- data.Patient$clinicalData %>%
  # Convert dates.
  dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    dateCensor = ifelse(!is.na(`Date: Death`), `Date: Death`, na.omit(c(`Date: Last follow-up`, `Date: End of study`, `Date: Pre-screening`))),
    dateCensor = as.Date(dateCensor, origin = '1970-01-01'),
    daysFromPreScreeningToEnd = dateCensor - `Date: Pre-screening`,
    monthsFromPreScreeningToEnd = daysFromPreScreeningToEnd / (365.25 / 12)
  ) %>%
  dplyr::inner_join(data.Patient$metrics) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    `Dichotomized CTC count (Baseline)` = ifelse(`CTC Count (Baseline – 7.5mL)` >= 5, 'CTC Count ≥5', 'CTC Count <5'),
    `WHO status (Pooled)` = ifelse(`WHO/ECOG PS at registration` %in% c(1,2), '1-2', `WHO/ECOG PS at registration`)
  ) %>% 
  dplyr::select(`Genome-wide status (Baseline)`, `Dichotomized CTC count (Baseline)`, `WHO status (Pooled)`, monthsFromPreScreeningToEnd, Survival) %>% 
  dplyr::filter(`Genome-wide status (Baseline)` != '.') %>% 
  dplyr::mutate(combinedScores = paste(`Genome-wide status (Baseline)`, `Dichotomized CTC count (Baseline)`, sep = ' & '))


# Multivariate analysis. ----

data.Patient$survival %>% 
  dplyr::filter(!is.na(Survival)) %>% 
  survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ `Genome-wide status (Baseline)` + `Dichotomized CTC count (Baseline)` + `WHO status (Pooled)`, data = ., ties = 'breslow') %>% 
  plotHR(., withQ = T)


# Survival Analysis (Cox regression) ----

survData <- data.frame(data.Patient$survival)
plotFits <- list()

## Dichotomized mFAST-SeqS. ----
fit.mFASTSeqs <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status..Baseline., data = survData)
names(fit.mFASTSeqs$strata) <-  base::gsub('.*=', '', names(fit.mFASTSeqs$strata))
plotFits$fit.mFASTSeqs <- plotSurvival(fit.mFASTSeqs, data = survData, ylim = 45, palette = c('#648FFF', '#FE6100'))

## Dichotomized CTC Counts. ----
fit.CTC <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count..Baseline., data = survData)
names(fit.CTC$strata) <-  base::gsub('.*=', '', names(fit.CTC$strata))
plotFits$fit.CTC <- plotSurvival(fit.CTC, data = survData, ylim = 45, palette = c('#f23005', '#ffbe73'))

## Dichotomized WHO. ----
fit.WHO <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.status..Pooled., data = survData)
names(fit.WHO$strata) <-  base::gsub('.*=', 'WHO: ', names(fit.WHO$strata))
plotFits$fit.WHO <- plotSurvival(fit.WHO, data = survData, ylim = 45, palette = c('#2C3D4F', '#1ABB9A'))


## Combine plots. ----
layout <- 'ABC
DEF'

plotFits$fit.mFASTSeqs$plot +
  plotFits$fit.CTC$plot +
  plotFits$fit.WHO$plot +
  plotFits$fit.mFASTSeqs$table +
  plotFits$fit.CTC$table +
  plotFits$fit.WHO$table +
  patchwork::plot_layout(design = layout, heights = c(1, .2), guides = 'auto') +
  patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Arial'))


## Compare models. ----

fit.mFASTSeqS <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status, data = survData2)
fit.CTC <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count, data = survData2)

AIC(fit.mFASTSeqS, fit.CTC)


# Combine CTC and mFAST-SeqS ----

fit.Combined <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ combinedScores, data = survData)
names(fit.Combined$strata) <-  base::gsub('.*=', '', names(fit.Combined$strata))
plotFits$fit.Combined <- plotSurvival(fit.Combined, data = survData, ylim = 45, palette = c('#E64B35FF', '#4DBBD5FF', '#00A087FF', '#3C5488FF'))
