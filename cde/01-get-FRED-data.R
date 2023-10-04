# Program: 01-get-FRED-data
# Purpose: Retrieve data from FRED.
# 
# Because the data are revised, the saved dataset includes an optional date
# stamp. This feature allows the project to be replicated or updated using the
# latest data. This feature of the code can be adjusted by changing 
# date_stamp_true:
#   1 : include date stamp to generated .csv file
#   0 : do not include date stamp to generated .csv file
# 
# Date Started: 2023-08-21
# Date Revised: 2023-08-21
library(tidyverse)
library(lubridate)
library(tidyquant) # To retrieve data from FRED
library(here)

file_prg <- "01-get-FRED-data"
date_stamp_true <- 1

# Retrieve data from FRED -------------------------------------------------

# Hires: Total Nonfarm (JTSHIL)	DOWNLOAD 
# Observation: Apr 2020: 3,524 (+ more) 
# Updated: Jun 9, 2020
# Units: Level in Thousands, Seasonally Adjusted
# Frequency: Monthly
# 
# Hires: Total Private (JTS1000HIL)	DOWNLOAD 
# Observation: Apr 2020: 3,305 (+ more) 
# Updated: Jun 9, 2020
# Units: Level in Thousands, Seasonally Adjusted
# Frequency: Monthly
# 
# Job Openings: Total Nonfarm (JTSJOL)	DOWNLOAD 
# Observation: Apr 2020: 5,046 (+ more) 
# Updated: Jun 9, 2020
# Units: Level in Thousands, Seasonally Adjusted
# Frequency: Monthly
# 
# Unemployment Level (UNEMPLOY)	DOWNLOAD 
# Observation: May 2020: 20,985 (+ more) 
# Updated: Jun 5, 2020
# Units: Thousands of Persons, Seasonally Adjusted
# Frequency: Monthly
# 
# All Employees, Total Nonfarm (PAYEMS)
# Observation: Jul 2023: 156,342 (+ more)   Updated: Aug 4, 2023
# Units: Thousands of Persons, Seasonally Adjusted
# Frequency: Monthly 
# 
# Labor Force Flows Employed to Employed (LNS17000000)
# Observation: Jul 2023: 154,589 (+ more)   Updated: Aug 4, 2023
# Units: Thousands of Persons, Seasonally Adjusted
# Frequency: Monthly 

fred_series <- c("JTSHIL", "JTS1000HIL", "UNEMPLOY", "JTSJOL", "PAYEMS", "UNRATE", "LNS17000000")
dat <- tq_get(fred_series,
              get = "economic.data",
              from = "1800-01-01")

dat <- dat %>% 
  mutate(series = case_when(
    symbol == "JTSHIL"~ "hires",
    symbol == "JTS1000HIL" ~ "hires_private",
    symbol == "UNEMPLOY" ~ "unemp",
    symbol == "JTSJOL" ~ "vacancies",
    symbol == "PAYEMS"~ "empl",
    symbol == "UNRATE" ~ "ur",
    symbol == "LNS17000000" ~ "eeflow"
  )) %>% 
  rename(FRED_symbol = symbol,
         lvl = price)

if (date_stamp_true == 1) {
  fout <- paste0("dat_", file_prg, "_", today(), ".csv")
} else {
  fout <- paste0("dat_", file_prg, ".csv")
}

write_csv(dat, here("out", fout))
