library(tidyverse)
library(ggrepel)
library(assertr)
library(here)

file_prg <- "05-plot-replacement-rates"

mywidth <- 6
golden <- 0.5*(1 + sqrt(5))
myheight <- mywidth / golden

csub_blue <- rgb(0, 53, 148, maxColorValue = 255)
csub_gray <- rgb(112, 115, 114, maxColorValue = 255)


# Read in data ------------------------------------------------------------

# Replacement-rate files
rr_files <- c("2001-1997_UI_Replacement_Rates-April-19-2024.csv",
              "2008-2002_UI_Replacement_Rates-April-19-2024.csv",
              "2015-2009_UI_Replacement_Rates-April-19-2024.csv",
              "2023-2016_UI_Replacement_Rates-April-19-2024.csv")

read_in_dol_data <- function(fname) {
  ii_fin <- here("dta", fname)
  ii_dat <- read_csv(ii_fin)
  ii_dat |> filter(str_detect(Year, "^\\d\\d+"))
}

dat <- read_in_dol_data(rr_files[[1]])

for (ii in 2:length(rr_files)) {
  ii_dat <- read_in_dol_data(rr_files[[ii]])
  dat <- bind_rows(dat, ii_dat)
}

# Check
dat_plt <- dat |> 
  mutate(`...9` = str_replace_na(`...9`, " "), 
         `Average Weekly Wage` = str_c(`Average Weekly Wage`, `...9`, sep = ""),
         average_wba = as.numeric(str_replace(`Average WBA`, "[$]", "")),
         average_weekly_wage = as.numeric(str_replace(`Average Weekly Wage`, "[$]", "")),
         year = as.numeric(Year),
         check_rr2 = average_wba / average_weekly_wage)

ggplot(data = dat_plt) +
  geom_line(mapping = aes(x = year, y = `Replacement Ratio 2`), color = "black") +
  geom_line(mapping = aes(x = year, y = check_rr2), color = "red")


# Plot: Weekly Benefit Amount and Weekly Wage -----------------------------

ggplot(data = dat_plt) +
  geom_line(mapping = aes(x = year, y = average_wba), color = "red") +
  geom_line(mapping = aes(x = year, y = average_weekly_wage)) +
  labs(x = "", y = "Amount in dollars")

# Plot replacement rates --------------------------------------------------

# dat_plt <- dat_plt |> 
#   mutate(my_series_label = case_when(
#     year == 2003 ~ ""
#   ))

mean_rr1 <- mean(dat_plt$`Replacement Ratio 1`)
mean_rr2 <- mean(dat_plt$`Replacement Ratio 2`)

dat_mylabel <- tribble(
  ~x,   ~y,       ~mylabel,
  2021, mean_rr1, round(mean_rr1, 2),
  2007, mean_rr2, round(mean_rr2, 2)
)

myltype = "dotted"
ggplot(data = dat_plt) +
  # Replacement ratio 1
  geom_line(mapping = aes(x = year, y = `Replacement Ratio 1`), linewidth = 1.0, color = "black") +
  geom_hline(yintercept = mean_rr1, linetype = myltype, linewidth = 1.0) +
  geom_text_repel(data = filter(dat_plt, year == 2021),
                  mapping = aes(x = year, y = `Replacement Ratio 1`, label = "Replacement ratio 1"), color = "black") +
  # Replacement ratio 2
  geom_line(mapping = aes(x = year, y = `Replacement Ratio 2`), color = csub_blue, linewidth = 1.0) +
  geom_hline(yintercept = mean_rr2, color = csub_blue, linetype = myltype, linewidth = 1.0) +
  geom_text_repel(data = filter(dat_plt, year == 2005),
                  mapping = aes(x = year, y = `Replacement Ratio 2`, label = "Replacement ratio 2"), color = csub_blue, box.padding = 1.1) +  
  # Label averages
  geom_text_repel(data = dat_mylabel,
                  mapping = aes(x = x, y = y, label = mylabel)) +
  labs(x = "", y = "Replacement ratios") +
  theme_minimal() 

fout <- paste0("fig_", file_prg, ".pdf")
ggsave(here("out", fout), heigh = myheight, width = mywidth)

