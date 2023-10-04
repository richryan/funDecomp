# Program: 02-adjust-transition-rates.R
# Purpose: Adjust job-finding and job-separation rates. The adjustment uses the
# fact that data are available at discrete intervals but transitions happen
# continuously.
# 
# Date Started: 2023-08-22
# Date Revised: 2023-08-22
library(tidyverse)
library(tidyquant)
library(lubridate)
library(ggrepel)
library(here)

file_prg <- "02-adjust-transition-rates"

# Golden ratio plotting parameters
mywidth <- 6
golden <- 0.5*(1 + sqrt(5))
myheight <- mywidth / golden

csub_blue <- rgb(0, 53, 148, maxColorValue = 255)

# Read in the data --------------------------------------------------------

dat <- read_csv(here("out", "dat_01-get-FRED-data_2023-08-22.csv"))

dat_wide <- dat %>% 
  pivot_wider(id_cols = date, names_from = series, values_from = lvl) %>% 
  mutate(find_crude = hires / unemp,
         sep_crude = find_crude * (ur / 100) / (1 - (ur / 100)),
         lf = empl + unemp) %>% 
  drop_na(hires) %>% 
  arrange(date)

ggplot(data = drop_na(dat_wide, find_crude)) +
  geom_line(mapping = aes(x = date, y = find_crude)) +
  geom_line(mapping = aes(x = date, y = sep_crude), color = "red") +
  geom_hline(yintercept = 1, color = "blue")

# Adjustment of transition probabilities ----------------------------------

make_get_s <- function(e_t, eh_t, e_t1) {
  function(s) {
    sep <- 1 - exp(-s)
    T1 <- e_t * (1 - sep) + eh_t 
    T2 <- eh_t * (1 - sep / s)
    ret <- e_t1 - T1 + T2
  }
}

compute_s <- function(sep_crude, e_t, eh_t, e_t1) {
  if (is.na(sep_crude) | is.na(e_t) | is.na(eh_t) | is.na(e_t1)) {
    NA
  } else {
    s0 <- -log(1 - sep_crude)
    get_s <- make_get_s(e_t, eh_t, e_t1)
    sol <- uniroot(
      get_s,
      lower = 0.0001, 
      upper = 0.2, 
      extendInt = "yes",
      check.conv = TRUE,
      tol = 1e-8
    ) 
    return(sol$root)
  }
}

make_get_f <- function(e_t, e_t1, lf, s) {
  function(f) {
    e_t1 - f * lf * (1 - exp(-s -f)) / (s + f) - e_t * exp(-s - f)
  }
}

compute_f <- function(find_crude, e_t, e_t1, lf, s) {
  if (is.na(find_crude) | is.na(e_t) | is.na(e_t1) | is.na(lf) | is.na(s)) {
    NA
  } else {
    f0 <- -log(1 - find_crude)
    get_f <- make_get_f(e_t, e_t1, lf, s)
    sol <- uniroot(get_f,
                   lower = 0.1 * find_crude,
                   upper = 2.0 * find_crude,
                   extendInt = "yes",
                   check.conv = TRUE,
                   tol = 1e-8)
    return(sol$root)
  }
}

dat_wide <- dat_wide %>% 
  mutate(sep_approx = 1 - (lead(empl) - hires) / empl) %>% 
  mutate(s = pmap_dbl(list(sep_crude, empl, hires, lead(empl)), compute_s),
         sep = 1 - exp(-s),
         f = pmap_dbl(list(find_crude, empl, lead(empl), lf, s), compute_f),
         find = 1 - exp(-f))

plot(dat_wide$date, dat_wide$sep, type = "l", col = "black", main = "Separation Rates")
points(dat_wide$date, dat_wide$sep_approx, type = "l", col = "blue")
points(dat_wide$date, dat_wide$sep_crude, type = "l", col = "red")

plot(dat_wide$date, dat_wide$find_crude, type = "l", col = "blue", main = "Finding Rates")
points(dat_wide$date, dat_wide$find, type = "l", col = "black")
abline(a = 1, b = 0)

# Recession data for plotting ----------------------------------------------------------------

# Data for recession shading
recess <- tq_get("USRECM", get = "economic.data", from = "1800-01-01")
recess_dat <- recess %>% 
  arrange(date) %>% 
  mutate(same = 1 - (price == lag(price))) %>% 
  # Remove first row, an NA, for cumulative sum
  filter(date > min(recess$date)) %>% 
  mutate(era = cumsum(same)) %>% 
  # Filter only recessions
  filter(price == 1)

recess_dat <- recess_dat %>% 
  group_by(era) %>% 
  # Unncessary, but to be sure...
  arrange(date) %>% 
  filter(row_number() == 1 | row_number() == n())

# Now reshape the data wide.
# Each row will contain the start and end dates of a recession.
recess_dat <- recess_dat %>% 
  mutate(junk = row_number()) %>% 
  mutate(begin_end = case_when(
    junk == 1 ~ "begin",
    junk == 2 ~ "end"
  ))

recess_wide <- recess_dat %>%
  ungroup() %>% 
  select(symbol, price, date, era, begin_end) %>% 
  pivot_wider(names_from = begin_end, values_from = date)

# Plots of job-finding and job-separation probabilities -------------------

dat_sep <- dat_wide %>% 
  select(date, starts_with("sep")) %>% 
  pivot_longer(cols = starts_with("sep"), names_to = "series", values_to = "lvl") %>% 
  mutate(my_label = case_when(
    series == "sep_crude" & date == ymd("2009-01-01") ~ "Uncorrected",
    series == "sep_approx" & date == ymd("2021-01-01") ~ "Approximated",
    series == "sep" & date == ymd("2017-03-01") ~ "Corrected",
    TRUE ~ ""
  ))


chart_begin <- min(dat_sep$date)
chart_end <- max(dat_sep$date)
psep <- ggplot(data = dat_sep) +
  geom_rect(
    data = filter(recess_wide,
                  begin >= chart_begin,
                  begin <= chart_end),
    mapping = aes(
      xmin = begin,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = csub_blue,
    alpha = 0.2
  ) +      
  geom_line(mapping = aes(x = date, y = lvl, color = series, linetype = series)) +
  geom_text_repel(mapping = aes(x = date, y = lvl, label = my_label, color = series),
                  box.padding = 0.75,
                  max.overlaps = Inf) +
  scale_color_viridis_d(end = 0.6) +  
  theme_minimal() +
  labs(x = "", y = "Probability of separating") +  
  theme(legend.position = "none")

fname <- paste0("fig_", file_prg, "_sep.pdf")
ggsave(here("out", fname), heigh = myheight, width = mywidth)

psep + coord_cartesian(ylim=c(0.02, 0.06))
fname <- paste0("fig_", file_prg, "_sep-truncated.pdf")
ggsave(here("out", fname), heigh = myheight, width = mywidth)

dat_find <- dat_wide %>% 
  select(date, find, find_crude) %>% 
  pivot_longer(cols = starts_with("find"), names_to = "series", values_to = "lvl") %>% 
  mutate(my_label = case_when(
    series == "find" & date == ymd("2004-01-01") ~ "Corrected",
    series == "find_crude" & date == ymd("2017-01-01") ~ "Uncorrected",
    TRUE ~ ""
  ))

chart_begin <- min(dat_find$date)
chart_end <- max(dat_find$date)
ggplot(data = dat_find) +
  geom_rect(
    data = filter(recess_wide,
                  begin >= chart_begin,
                  begin <= chart_end),
    mapping = aes(
      xmin = begin,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = csub_blue,
    alpha = 0.2
  ) +    
  geom_hline(yintercept = 1.0, color = "black") +
  geom_line(mapping = aes(x = date, y = lvl, color = series, linetype = series),
            linewidth = 1.0) +
  geom_text_repel(mapping = aes(x = date, y = lvl, label = my_label, color = series),
                  max.overlaps = Inf, 
                  box.padding = 0.5) +
  labs(x = "", y = "Probability of finding a job") +
  scale_color_viridis_d(end = 0.4) +
  theme_minimal() +
  theme(legend.position = "none")

fname <- paste0("fig_", file_prg, "_find.pdf")
ggsave(here("out", fname), heigh = myheight, width = mywidth)

# Save dataset ------------------------------------------------------------

dat_wide <- dat_wide %>% 
  select(-s, -f, lf)

fout <- paste0("dat_", file_prg, ".csv")
write_csv(dat_wide, here("out", fout))

