# Author:                      Job van Riet
# Date:                        07-07-2022
# Function:                    Compare baseline characteristics of CABA-V7 and CABARESC.

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
)

data.Table1_CABARESC <- dplyr::bind_rows(
    data.CABARESC %>% dplyr::inner_join(data.CABARESC_Mini) %>% dplyr::select(`Age at registration`, `WHO/ECOG PS at registration`, PSA = `PSA at baseline`, ANC = `ANC at baseline`, HBG = `hbm at baseline`, ALP = `alkphos at baseline`, LDH = `ldh at baseline`) %>% dplyr::mutate(cohort = 'CABARESC (mini)', `WHO/ECOG PS at registration` = as.character(`WHO/ECOG PS at registration`)),
    data.CABARESC %>% dplyr::select(`Age at registration`, `WHO/ECOG PS at registration`, PSA = `PSA at baseline`, ANC = `ANC at baseline`, HBG = `hbm at baseline`, ALP = `alkphos at baseline`, LDH = `ldh at baseline`) %>% dplyr::mutate(cohort = 'CABARESC', `WHO/ECOG PS at registration` = as.character(`WHO/ECOG PS at registration`)),
)

# Perform statistics. ----

test.Table1 <- data.Table1 %>% 
    reshape2::melt(id.vars = 'cohort') %>%
    dplyr::filter(variable != 'WHO/ECOG PS at registration') %>% 
    dplyr::mutate(value = as.numeric(value)) %>% 
    dplyr::group_by(variable) %>% 
    rstatix::t_test(paired = F, detailed = T, formula = value ~ cohort, alternative = 'two', p.adjust.method = 'none') 


test.Table1_CABARESC <- data.Table1_CABARESC %>% 
    reshape2::melt(id.vars = 'cohort') %>%
    dplyr::filter(variable != 'WHO/ECOG PS at registration') %>% 
    dplyr::mutate(value = as.numeric(value)) %>% 
    dplyr::group_by(variable) %>% 
    rstatix::t_test(paired = F, detailed = T, formula = value ~ cohort, alternative = 'two', p.adjust.method = 'none') 


# Output table. ----

data.Table1_CABARESC %>% 
reshape2::melt(id.vars = 'cohort') %>%
    dplyr::filter(variable != 'WHO/ECOG PS at registration') %>% 
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
    dplyr::inner_join(test.Table1_CABARESC %>% dplyr::select(variable, p, method)) %>% 
    dplyr::select(cohort, variable, outMean, outMedian, p, method)


print(tableone::CreateTableOne(
    vars = 'WHO/ECOG PS at registration',
    strata = c("cohort"),
    data = data.Table1_CABARESC,
    smd = F,
    test = T
), showAllLevels = T)

                  