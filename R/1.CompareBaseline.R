# Author:                      Job van Riet
# Date:                        09-08-2022
# Function:                    Compare relevant baseline clinical characteristics of CABA-V7 and CABARESC.

# Import libraries ----

library(plyr)
library(dplyr)
library(patchwork)


# Import data. ----

data.CABAV7 <- readxl::read_xlsx('Misc./Suppl. Table 1 - Overview of Data.xlsx', sheet = 'Clinical (CABAV7)', trim_ws = T)
data.CABARESC <- readxl::read_xlsx('Misc./Suppl. Table 1 - Overview of Data.xlsx', sheet = 'Clinical (CABARESC)', trim_ws = T)
data.CABARESC_Mini <- readxl::read_xlsx('Misc./Suppl. Table 1 - Overview of Data.xlsx', sheet = 'Overview (CABARESC)', trim_ws = T, skip = 1)


data.Table1 <- dplyr::bind_rows(
    data.CABARESC %>% dplyr::inner_join(data.CABARESC_Mini) %>% dplyr::select(`Age at registration`, `WHO/ECOG PS at registration`, PSA = `PSA at baseline`, ANC = `ANC at baseline`, HBG = `hbm at baseline`, ALP = `alkphos at baseline`, LDH = `ldh at baseline`) %>% dplyr::mutate(cohort = 'CABARESC', `WHO/ECOG PS at registration` = as.character(`WHO/ECOG PS at registration`)),
    data.CABAV7 %>% dplyr::select(`Age at registration`, `WHO/ECOG PS at registration`, PSA = `PSA at primary diagnosis [ug/L]`, ANC, HBG, ALP, LDH) %>% dplyr::mutate(cohort = 'CABA-V7')
) %>% 
    dplyr::mutate(
        `WHO status (Pooled)` = ifelse(`WHO/ECOG PS at registration` %in% c(1,2), '1 - 2', `WHO/ECOG PS at registration`)
    )

data.Table1_CABARESC <- dplyr::bind_rows(
    data.CABARESC %>% dplyr::inner_join(data.CABARESC_Mini) %>% dplyr::select(`Age at registration`, `WHO/ECOG PS at registration`, PSA = `PSA at baseline`, ANC = `ANC at baseline`, HBG = `hbm at baseline`, ALP = `alkphos at baseline`, LDH = `ldh at baseline`) %>% dplyr::mutate(cohort = 'CABARESC (mini)', `WHO/ECOG PS at registration` = as.character(`WHO/ECOG PS at registration`)),
    data.CABARESC %>% dplyr::select(`Age at registration`, `WHO/ECOG PS at registration`, PSA = `PSA at baseline`, ANC = `ANC at baseline`, HBG = `hbm at baseline`, ALP = `alkphos at baseline`, LDH = `ldh at baseline`) %>% dplyr::mutate(cohort = 'CABARESC', `WHO/ECOG PS at registration` = as.character(`WHO/ECOG PS at registration`)),
) %>% 
    dplyr::mutate(
        `WHO status (Pooled)` = ifelse(`WHO/ECOG PS at registration` %in% c(1,2), '1 - 2', `WHO/ECOG PS at registration`)
    )


# Perform statistics on baseline characteristics. ----

tableone::CreateTableOne(
    vars = c('Age at registration', 'PSA', 'ANC', 'HBG', 'ALP', 'LDH', 'WHO status (Pooled)'),
    strata = 'cohort', factorVars = 'WHO status (Pooled)', includeNA = T, 
    data = data.Table1_CABARESC, 
    smd = F,
    test = T
) %>% 
    print(.,
          exact = c('n'),
          smd = F,
          minMax = T,
          missing = T,
          explain = T,
          showAllLevels = T,
          quote = F,
          noSpaces = T
    )


# Output table. ----

data.Table1 %>% 
    reshape2::melt(id.vars = 'cohort') %>%
    dplyr::filter(!grepl('WHO', variable)) %>% 
    dplyr::group_by(cohort, variable) %>%
    dplyr::mutate(value = as.numeric(value)) %>% 
    dplyr::summarise(
        mean = mean(na.omit(value)),
        sd = sd(na.omit(value)),
        median = median(na.omit(value)),
        min = min(na.omit(value)),
        max = max(na.omit(value)),
        outMean = sprintf('%s±%s', round(mean, 1), round(sd, 1)),
        outMedian = sprintf('%s (%s - %s)', round(median, 1), round(min, 1), round(max, 1)),
        missing = sum(is.na(value))
    ) %>% 
    dplyr::select(cohort, variable, outMean, outMedian)


# Selection of relevant clinical characteristics. ----

data.AIC <- readxl::read_xlsx('Misc./Suppl. Table 1 - Overview of Data.xlsx', sheet = 'Clinical (CABAV7)', trim_ws = T) %>% 
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
    dplyr::mutate(
        `WHO status (Pooled)` = ifelse(`WHO/ECOG PS at registration` %in% c(1,2), '1 - 2', `WHO/ECOG PS at registration`),
        `Dichotomized CTC count (Baseline)` = ifelse(`CTC Count (Baseline – 7.5mL)` >= 5, 'CTC Count ≥5', 'CTC Count <5'),
    ) %>% 
    dplyr::select(
        `Age at registration`,
        `WHO status (Pooled)`,
        `PSA at primary diagnosis [ug/L]`,
        HBG, ALB, LDH, ALP, ANC,
        `Genome-Wide Z Score (Baseline)`,
        `Dichotomized CTC count (Baseline)`,
        Survival, monthsFromPreScreeningToEnd
    ) %>% 
    dplyr::filter(complete.cases(.))

MASS::stepAIC(survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ `Age at registration` + `WHO status (Pooled)` + `PSA at primary diagnosis [ug/L]` + HBG + ALB + 
                                  LDH + ALP + ANC + `Genome-Wide Z Score (Baseline)` + `Dichotomized CTC count (Baseline)`, data = data.AIC, ties = 'breslow') , direction = 'backward')


# Outcomes:
#                                                   coef exp(coef)  se(coef)      z       p
# `WHO status (Pooled)`1 - 2                       0.729202  2.073426  0.294403  2.477 0.01325
# HBG                                             -0.507536  0.601977  0.190796 -2.660 0.00781
# `Genome-Wide Z Score (Baseline)`                 0.012943  1.013027  0.004631  2.795 0.00519
# `Dichotomized CTC count (Baseline)`CTC Count ≥5  1.082219  2.951223  0.342824  3.157 0.00160