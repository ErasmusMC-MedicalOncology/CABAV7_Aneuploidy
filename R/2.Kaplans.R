# Author:                      Job van Riet
# Date:                        24-03-2023
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
data.Patient$CABAV7.Overview <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Overview (CABAV7)') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$CABARESC.Overview <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Overview (CABARESC)') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))

data.Patient$CABAV7.clinical <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Clinical (CABAV7)') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$CABARESC.clinical <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Clinical (CABARESC)') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))


# Convert and clean data. ----

## CABA-V7 ----
data.Patient$CABAV7.OS <- data.Patient$CABAV7.clinical %>%
  # Convert dates.
  dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    dateCensor = ifelse(!is.na(`Date: Death`), `Date: Death`, na.omit(c(`Date: Last follow-up`, `Date: End of study`, `Date: Pre-screening`))),
    dateCensor = as.Date(dateCensor, origin = '1970-01-01'),
    daysFromPreScreeningToEnd = dateCensor - `Date: Pre-screening`,
    monthsFromPreScreeningToEnd = daysFromPreScreeningToEnd / (365.25 / 12)
  ) %>%
  dplyr::inner_join(data.Patient$CABAV7.Overview) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    `Dichotomized CTC count (Baseline)` = ifelse(`CTC Count (Baseline – 7.5mL)` >= 5, 'CTC Count ≥5', 'CTC Count <5'),
    `WHO status (Pooled)` = ifelse(`WHO/ECOG PS at registration` %in% c(1,2), '1-2', `WHO/ECOG PS at registration`)
  ) %>% 
  dplyr::select(`Genome-wide status (Baseline)`, `Dichotomized CTC count (Baseline)`, `WHO status (Pooled)`, monthsFromPreScreeningToEnd, Survival) %>% 
  dplyr::filter(`Genome-wide status (Baseline)` != '.') %>% 
  dplyr::mutate(
    `Genome-wide status (Baseline)` = gsub('Genome-wide Z-', 'Aneuploidy ', `Genome-wide status (Baseline)`),
    `Combined Scores` = paste(`Genome-wide status (Baseline)`, `Dichotomized CTC count (Baseline)`, sep = ' & ')
  )


## CABARESC ----

data.Patient$CABARESC.OS <- data.Patient$CABARESC.clinical %>% 
  dplyr::left_join(data.Patient$CABARESC.Overview) %>% 
  dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
  dplyr::mutate(
    `Subject Number` = as.character(`Subject Number`),
    daysFromPreScreeningToEnd = `Date: Last contact (FU)` - `Date: Registration`, 
    monthsFromPreScreeningToEnd = daysFromPreScreeningToEnd / (365.25 / 12),
    `Dichotomized CTC count` = ifelse(`CTC Count` >= 5, 'CTC Count ≥5', 'CTC Count <5'),
    `WHO status (Pooled)` = ifelse(`WHO/ECOG PS at registration` %in% c('1', '2'), '1-2', '0')
  ) %>% 
  dplyr::select(`Genome-wide status`, `Dichotomized CTC count`, `WHO status (Pooled)`, monthsFromPreScreeningToEnd, Survival) %>% 
  dplyr::mutate(
    `Genome-wide status` = gsub('Genome-wide Z-', 'Aneuploidy ', `Genome-wide status`),
    `Combined Scores` = paste(`Genome-wide status`, `Dichotomized CTC count`, sep = ' & ')
  ) %>% 
  dplyr::mutate(subset = ifelse(!is.na(`Genome-wide status`), T, F))


# Multivariate analysis. ----

data.Patient$CABAV7.OS %>% 
  dplyr::filter(!is.na(Survival)) %>% 
  survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ `Genome-wide status (Baseline)` + `Dichotomized CTC count (Baseline)` + `WHO status (Pooled)`, data = ., ties = 'breslow') %>% 
  plotHR(.)


# Survival Analysis (Cox regression) ----

## CABA-V7

survData <- data.frame(data.Patient$CABAV7.OS)
plotFits <- list()

## Dichotomized mFAST-SeqS. ----
fit.mFASTSeqs <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status..Baseline., data = survData)
names(fit.mFASTSeqs$strata) <-  base::gsub('.*=', '', names(fit.mFASTSeqs$strata))
plotFits$fit.CABAV7.mFASTSeqs <- plotSurvival(fit.mFASTSeqs, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status..Baseline., data = survData), data = survData, ylim = 45, palette = c('#648FFF', '#FE6100'))


## Dichotomized CTC Counts. ----
fit.CTC <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count..Baseline., data = survData)
names(fit.CTC$strata) <-  base::gsub('.*=', '', names(fit.CTC$strata))
plotFits$fit.CABAV7.CTC <- plotSurvival(fit.CTC, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count..Baseline., data = survData), data = survData, ylim = 45, palette = c('#f23005', '#ffbe73'))

## Dichotomized WHO. ----
fit.WHO <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.status..Pooled., data = survData)
names(fit.WHO$strata) <-  base::gsub('.*=', 'WHO: ', names(fit.WHO$strata))
plotFits$fit.CABAV7.WHO <- plotSurvival(fit.WHO, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.status..Pooled., data = survData), data = survData, ylim = 45, palette = c('#2C3D4F', '#1ABB9A'))

## CABARESC (Subset)

survData <- data.frame(data.Patient$CABARESC.OS %>% dplyr::filter(subset))

## Dichotomized mFAST-SeqS. ----
fit.mFASTSeqs <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status, data = survData)
names(fit.mFASTSeqs$strata) <-  base::gsub('.*=', '', names(fit.mFASTSeqs$strata))
plotFits$fit.CABARESC.mFASTSeqs <- plotSurvival(fit.mFASTSeqs, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Genome.wide.status, data = survData), data = survData, ylim = 45, palette = c('#648FFF', '#FE6100'))


## Dichotomized CTC Counts. ----
fit.CTC <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count, data = survData)
names(fit.CTC$strata) <-  base::gsub('.*=', '', names(fit.CTC$strata))
plotFits$fit.CABARESC.CTC <- plotSurvival(fit.CTC, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count, data = survData), data = survData, ylim = 45, palette = c('#f23005', '#ffbe73'))

## Dichotomized WHO. ----
fit.WHO <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.status..Pooled., data = survData)
names(fit.WHO$strata) <-  base::gsub('.*=', 'WHO: ', names(fit.WHO$strata))
plotFits$fit.CABARESC.WHO <- plotSurvival(fit.WHO, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.status..Pooled., data = survData), data = survData, ylim = 45, palette = c('#2C3D4F', '#1ABB9A'))


## Combine plots. ----

layout <- 'ABC
DEF
GHI
JKL'

plotFits$fit.CABAV7.mFASTSeqs$plot +
  plotFits$fit.CABAV7.CTC$plot +
  plotFits$fit.CABAV7.WHO$plot +
  plotFits$fit.CABAV7.mFASTSeqs$table +
  plotFits$fit.CABAV7.CTC$table +
  plotFits$fit.CABAV7.WHO$table +
  plotFits$fit.CABARESC.mFASTSeqs$plot +
  plotFits$fit.CABARESC.CTC$plot +
  plotFits$fit.CABARESC.WHO$plot +
  plotFits$fit.CABARESC.mFASTSeqs$table +
  plotFits$fit.CABARESC.CTC$table +
  plotFits$fit.CABARESC.WHO$table +
  patchwork::plot_layout(design = layout, heights = c(1, .2, 1, .2), guides = 'auto') +
  patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Roboto'))


## Compare OS: CABAV7 vs. CABARESC. ----

survData <- data.frame(
  dplyr::bind_rows(
    data.Patient$CABAV7.OS %>% dplyr::mutate(strata = 'CABA-V7<br><i>n</i>=131'),
    data.Patient$CABARESC.OS %>% dplyr::mutate(strata = 'CABARESC<br><i>n</i>=224'),
    data.Patient$CABARESC.OS %>% dplyr::filter(subset) %>% dplyr::mutate(strata = 'CABARESC (Validation)<br><i>n</i>=50')
  )
)


fit.between <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ strata, data = survData)
names(fit.between$strata) <-  base::gsub('strata=', '', names(fit.between$strata))
plotSurvival(fit.between, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ strata, data = survData), data = survData, ylim = 60, palette = c('#FF8E2B', '#2980B9', '#50AC5D'))


# Combine CTC and mFAST-SeqS ----

plotCombined <- function(data){
  data =  data.frame(data)
  fit.Combined <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Combined.Scores, data = data)
  names(fit.Combined$strata) <-  base::gsub('.*=', '', names(fit.Combined$strata))
  return(plotSurvival(fit.Combined, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Combined.Scores, data = data), data = data, ylim = 45, palette = c('#0073C2', '#EFC000', '#CD534C', '#00A100')))
  
}

plotFits$fitCombined.CABAV7 <- plotCombined(data.Patient$CABAV7.OS)
plotFits$fitCombined.CABARESC <- plotCombined(data.Patient$CABARESC.OS %>% dplyr::filter(subset))

layout <- 'AB
DE'

plotFits$fitCombined.CABAV7$plot +
  plotFits$fitCombined.CABARESC$plot +
  plotFits$fitCombined.CABAV7$table +
  plotFits$fitCombined.CABARESC$table +
  patchwork::plot_layout(design = layout, heights = c(1, .2), guides = 'auto') +
  patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Roboto'))


# Add HR. ----

data.Patient$CABARESC.OS %>% 
  dplyr::filter(subset) %>% 
  dplyr::filter(!is.na(Survival)) %>% 
  survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ `Combined Scores`, data = ., ties = 'breslow') %>% 
  gtsummary::tbl_regression(
    exponentiate = T, 
    add_estimate_to_reference_rows = T,
  ) %>%
  gtsummary::add_n() %>% 
  gtsummary::add_nevent() %>% 
  gtsummary::bold_p() %>% 
  gtsummary::bold_labels() %>% 
  gtsummary::italicize_levels() %>% 
  gtsummary::sort_p() %>% 
  bstfun::add_inline_forest_plot(header = '', spec_pointrange.args = list(lim = c(-3, 3), width = 550, cex = 1, col = 'black', pch = 1))
