# Author:                      Job van Riet
# Date:                        09-05-2022
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
    dplyr::filter(`Genome-wide status (Baseline)` != '.')


# Multivariate analysis. ----

data.Patient$survival %>% 
    dplyr::filter(!is.na(Survival)) %>% 
    survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ `Genome-wide status (Baseline)` + `Dichotomized CTC count (Baseline)` + `WHO status (Pooled)`, data = ., ties = 'breslow') %>% 
    plotHR(., withQ = T)


# Survival Analysis (Cox regression ----

survData <- data.frame(data.Patient$survival)
plotFits <- list()

## Dichotomized mFAST-SeqS. ----
fit.mFASTSeqs <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status..Baseline., data = survData)
names(fit.mFASTSeqs$strata) <-  base::gsub('.*=', '', names(fit.mFASTSeqs$strata))
plotFits$fit.mFASTSeqs <- plotSurvival(fit.mFASTSeqs, data = survData, ylim = 45, palette = c('#648FFF', '#FE6100'))

## Dichotomized mFAST-SeqS. ----
fit.AllInClusion.CTC <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count..Baseline., data = survData)
names(fit.AllInClusion.CTC$strata) <-  base::gsub('.*=', '', names(fit.AllInClusion.CTC$strata))
plotFits$fit.AllInClusion.CTC <- plotSurvival(fit.AllInClusion.CTC, data = survData, ylim = 45, palette = c('#f23005', '#ffbe73'))

## Dichotomized WHO. ----
fit.AllInClusion.WHO <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.status..Pooled., data = survData)
names(fit.AllInClusion.WHO$strata) <-  base::gsub('.*=', 'WHO: ', names(fit.AllInClusion.WHO$strata))
plotFits$fit.AllInClusion.WHO <- plotSurvival(fit.AllInClusion.WHO, data = survData, ylim = 45, palette = c('#2C3D4F', '#1ABB9A'))


## Combine plots. ----
layout <- "ABC
DEF"

plotFits$fit.mFASTSeqs$plot +
    plotFits$fit.AllInClusion.CTC$plot +
    plotFits$fit.AllInClusion.WHO$plot +
    plotFits$fit.mFASTSeqs$table +
    plotFits$fit.AllInClusion.CTC$table +
    plotFits$fit.AllInClusion.WHO$table +
    patchwork::plot_layout(design = layout, heights = c(1, .2), guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Arial'))
