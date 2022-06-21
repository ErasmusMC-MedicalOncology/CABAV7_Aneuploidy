# Author:                      Job van Riet
# Date:                        21-06-2022
# Function:                    Compare performance of dichotomized CTC vs. mFAST-SeqS aneuploidy scores (CABARESC).

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

data.CABARESC <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'CABARESC', skip = 1) %>% 
  dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
  dplyr::mutate(
    `Subject Number` = as.character(`Subject Number`),
    daysFromPreScreeningToEnd = `Date: Last contact (FU)` - `Date: Registration`, 
    monthsFromPreScreeningToEnd = daysFromPreScreeningToEnd / (365.25 / 12),
    `Dichotomized CTC count` = ifelse(`CTC Count` >= 5, 'CTC Count ≥5', 'CTC Count <5'),
    `WHO status (Pooled)` = ifelse(`WHO performance status [0-5]` %in% c('WHO 1', 'WHO2'), '1-2', '0')
  ) %>% 
  dplyr::select(`Genome-wide status`, `Dichotomized CTC count`, `WHO status (Pooled)`, monthsFromPreScreeningToEnd, Survival) %>% 
  dplyr::mutate(combinedScores = paste(`Genome-wide status`, `Dichotomized CTC count`, sep = ' & '))

# Multivariate analysis. ----

data.CABARESC %>% 
  survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ `Genome-wide status` + `Dichotomized CTC count` + `WHO status (Pooled)`, data = ., ties = 'breslow') %>% 
  plotHR(., withQ = T)

# Survival Analysis (Cox regression) ----

survData <- data.frame(data.CABARESC)
plotFits <- list()

## Dichotomized mFAST-SeqS. ----
fit.mFASTSeqs <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status, data = survData)
names(fit.mFASTSeqs$strata) <-  base::gsub('.*=', '', names(fit.mFASTSeqs$strata))
plotFits$fit.mFASTSeqs <- plotSurvival(fit.mFASTSeqs, data = survData, ylim = 41, palette = c('#648FFF', '#FE6100'))

## Dichotomized CTC Counts. ----
fit.CTC <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count, data = survData)
names(fit.CTC$strata) <-  base::gsub('.*=', '', names(fit.CTC$strata))
plotFits$fit.CTC <- plotSurvival(fit.CTC, data = survData, ylim = 41, palette = c('#f23005', '#ffbe73'))

## Dichotomized WHO. ----
fit.WHO <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.status..Pooled., data = survData)
names(fit.WHO$strata) <-  base::gsub('.*=', 'WHO: ', names(fit.WHO$strata))
plotFits$fit.WHO <- plotSurvival(fit.WHO, data = survData, ylim = 41, palette = c('#2C3D4F', '#1ABB9A'))

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
