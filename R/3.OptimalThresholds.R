# Author:                      Job van Riet
# Date:                        04-04-2023
# Function:                    Determine optimal mFAST-SeqS thresholds for stratifying Survival status

# Import libraries ----

library(dplyr)
library(cutpointr)
library(ggplot2)
library(survminer)
library(gtsummary)
library(extrafont)
library(patchwork)

# Helper themes and functions.
source('R/misc_Themes.R')
source('R/misc_Functions.R')

plotFits <- list()

#  Determine cutoff based on CABA-V7. ----

data.Cutoff <- readxl::read_xlsx('Misc./Suppl. Table 1 - Overview of Data.xlsx', sheet = 'Clinical (CABAV7)', trim_ws = T) %>% 
  dplyr::inner_join(readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Overview (CABAV7)')) %>% 
  # Convert dates.
  dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    dateCensor = ifelse(!is.na(`Date: Death`), `Date: Death`, na.omit(c(`Date: Last follow-up`, `Date: End of study`, `Date: Pre-screening`))),
    dateCensor = as.Date(dateCensor, origin = '1970-01-01'),
    daysFromPreScreeningToEnd = dateCensor - `Date: Pre-screening`,
    monthsFromPreScreeningToEnd = daysFromPreScreeningToEnd / (365.25 / 12)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::select(
    `Genome-Wide Z Score (Baseline)`,
    `CTC Count (Baseline – 7.5mL)`,
    Survival, monthsFromPreScreeningToEnd
  ) %>% 
  dplyr::filter(complete.cases(.))

# Determine optimal cutt-off by max. optimizing sens./spec. to stratify survival status (1 / 0)
mFastSeqS.Threshold <- cutpointr::cutpointr(data.Cutoff, x = `Genome-Wide Z Score (Baseline)`, class = Survival, method = maximize_metric, metric = sum_sens_spec, pos_class	= 1, direction = '>=', boot_runs = 10000, boot_stratify = T)

## Optimal mFAST-Seq classes - CABA-V7. ----
survData <- data.Cutoff %>% 
  dplyr::mutate(
    Z = ifelse(`Genome-Wide Z Score (Baseline)` >= mFastSeqS.Threshold$optimal_cutpoint, sprintf('Aneuploidy score ≥%s', mFastSeqS.Threshold$optimal_cutpoint), sprintf('Aneuploidy score <%s', mFastSeqS.Threshold$optimal_cutpoint)),
  )

fit <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Z, data = survData)
names(fit$strata) <-  base::gsub('.*=', '', names(fit$strata))
plotFits$fit.CABAV7.mFASTSeqs.Optimized <- plotSurvival(fit, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Z, data = survData), data = survData, ylim = 45, palette = c('#648FFF', '#FE6100'))

## Optimal mFAST-Seq classes - CABARESC. ----

survData <-  readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Overview (CABARESC)') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`)) %>% 
  dplyr::left_join(readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Clinical (CABARESC)') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))) %>% 
  dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
  dplyr::mutate(
    `Subject Number` = as.character(`Subject Number`),
    daysFromPreScreeningToEnd = `Date: Last contact (FU)` - `Date: Registration`, 
    monthsFromPreScreeningToEnd = daysFromPreScreeningToEnd / (365.25 / 12),
    Z = ifelse(`Genome-Wide Z Score` >= mFastSeqS.Threshold$optimal_cutpoint, sprintf('Aneuploidy score ≥%s', mFastSeqS.Threshold$optimal_cutpoint), sprintf('Aneuploidy score <%s', mFastSeqS.Threshold$optimal_cutpoint)),
  ) %>% 
  dplyr::select(Z, monthsFromPreScreeningToEnd, Survival) 

fit <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Z, data = survData)
names(fit$strata) <-  base::gsub('.*=', '', names(fit$strata))
plotFits$fit.CABARESC.mFASTSeqs.Optimized <- plotSurvival(fit, hr = survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Z, data = survData), data = survData, ylim = 45, palette = c('#648FFF', '#FE6100'))


## Plot kaplap of new mFAST-SeqS cutoff
plotFits$fit.CABAV7.mFASTSeqs.Optimized$plot +
  plotFits$fit.CABARESC.mFASTSeqs.Optimized$plot +
  plotFits$fit.CABAV7.mFASTSeqs.Optimized$table +
  plotFits$fit.CABARESC.mFASTSeqs.Optimized$table +
  patchwork::plot_layout(ncol = 2, heights = c(1, .2), guides = 'auto') +
  patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Roboto'))
