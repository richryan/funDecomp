library(tidyverse)
library(docstring)
library(ggrepel)
library(assertr)
library(here)

file_prg <- "04-plot-elasticity"

mywidth <- 6
golden <- 0.5*(1 + sqrt(5))
myheight <- mywidth / golden

csub_blue <- rgb(0, 53, 148, maxColorValue = 255)
csub_gray <- rgb(112, 115, 114, maxColorValue = 255)

# Functions ---------------------------------------------------------------

compute_find <- function(ttheta, m0, ggamma) {
  #' Compute probability of finding a job
  #'
  #' @param ttheta Number that represents market tightness, 
  #' ratio of vacancies to number of unemployed
  #' @param m0 Number that parameterizes matching efficiency.
  #' @param ggamma Number that represents curvature of matching technology
  #'
  #' @return Returns a number, the probability of finding a job.
  #'
  #' @examples 
  #' probability_find <- compute_find(0.5, 0.04, 0.1)
  m0 * ttheta / ((1 + ttheta^ggamma)^(1 / ggamma))
}

compute_find_prime <- function(ttheta, m0, ggamma) {
  T1 <- 1 - ((1 + ttheta^ggamma)^(-1)) * ttheta^ggamma
  T2 <- (1 + ttheta^ggamma)^(1 / ggamma)
  m0 * T1 / T2
}

# Check that the numerical derivative is close to the coded derivative
x_check <- seq(0.1, 1.5, length.out = 100)
my_test_ggamma <- 0.5
stopifnot(all(near(numDeriv::grad(func = function(x) compute_find(x, m0 = 1.0, ggamma = my_test_ggamma), x = x_check) - compute_find_prime(x_check, m0 = 1.0, ggamma = my_test_ggamma), 0)))

compute_fill <- function(ttheta, m0, ggamma) {
  #' Compute probability of filling a vacancy.
  #'
  #' @param ttheta Number that represents market tightness, 
  #' ratio of vacancies to number of unemployed
  #' @param m0 Number that parameterizes matching efficiency.
  #' @param ggamma Number that represents curvature of matching technology
  #'
  #' @return Returns a number, the probability of filling a vacancy.
  #'
  #' @examples 
  #' probability_fill <- compute_find(ttheta = 0.5, m0 = 0.04, ggamma = 0.1)
  m0 / ((1 + ttheta^ggamma)^(1 / ggamma))
}

# Parameters --------------------------------------------------------------

compute_monthly2daily <- function(x) {
  #' Convert a monthly rate to a daily rate
  #'
  #' @param x A number or vector of numbers that are monthly rates.
  #'
  #' @return A number or vector that corresponds to the implied daily rates, given the inputted monthly rates.
  #' @export
  #'
  #' @examples
  #' r_monthly <- 0.004
  #' r <- compute_monthly2daily(r_monthly)
  1 - (1 - x)^(1 / 30)
}

compute_daily2monthly <- function(x) {
  #' Convert a daily rate to a monthly rate
  #'
  #' @param x A number or vector of numbers that are daily rates.
  #'
  #' @return A number or vector that corresponds to the implied monthly rates, given the inputted daily rates.
  #' @export
  #'
  #' @examples
  #' r_monthly <- 0.004
  #' r <- compute_monthly2daily(r_monthly)
  #' r_monthly_check <- compute_monthly2daily(r)  
  1 - (1 - x)^30
}

target_tight <- 0.72
target_find_monthly <- 0.594
target_find <- compute_monthly2daily(target_find_monthly)
r_monthly <- 0.004
r <- compute_monthly2daily(r_monthly)
bbeta <- 1 / (1 + r)
s_monthly <- 0.036
s <- compute_monthly2daily(s_monthly)

target_elasticity <- 7.5636

ggamma <- readRDS(here("out", "dat_03-est-matching.rds"))

y <- 1 
# Pissarides (2009, Econometrica) calibration
z <- 0.71 
# Shimer (2005, AER) calibration
z_shimer2005 <- 0.4
pphi <- 0.50
TOL <- 1e-8

# Get matching efficiency to target job-finding rate
npts <- 100
match_efficiency_lo <- 1e-3
match_efficiency_hi <- 100
v_match_efficiency <- seq(match_efficiency_lo, match_efficiency_hi, length.out = npts)
v_find <- compute_find(ttheta = target_tight, m0 = v_match_efficiency, ggamma = ggamma)
plot(v_match_efficiency, v_find, type = "l")  
abline(h = target_find)

get_m0 <- function(m0, ttheta, ggamma, target_find) {
  #' Helper function to compute matching efficiency associated a target unemployment rate.
  #' 
  #' The function `get_m0()` can be used to searh for matching efficiency associated with targeted tightness.
  #'
  #' @param m0 A number for matching efficiency.
  #' @param ttheta A number for tightness, the ratio of vacancies to unemployment.
  #' @param ggamma A number for matching curvature.
  #' @param target_find A number for targeted  tightness.
  #'
  #' @return A number equal to the difference between the implied job-finding probability and the targeted job-finding probability.
  #' @export
  #'
  #' @examples
  #' sol_match_efficiency <- uniroot(get_m0, ttheta = target_tight, ggamma = ggamma, target_find = target_find, lower = match_efficiency_lo, upper = match_efficiency_hi)
  #' match_efficiency <- sol_match_efficiency$root
  find <- compute_find(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  find - target_find
}

sol_match_efficiency <- uniroot(get_m0, ttheta = target_tight, ggamma = ggamma, target_find = target_find,
                                lower = match_efficiency_lo, 
                                upper = match_efficiency_hi)

match_efficiency <- sol_match_efficiency$root

find_monthly <- compute_daily2monthly(compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma))
cat('Implied monthly unemployment rate =', s_monthly / (s_monthly + find_monthly))


# Table: Model results at different combinations of job-creation c --------

compute_c <- function(ttheta, m0, y, z, bbeta, s, ct, r, ch, cl, pphi, ggamma) {
  #' Compute the cost of posting a vacancy, given labor-market tightness.
  #'
  #' @param ttheta A number for labor-market tightness, ratio of vacancies to unemployed.
  #' @param m0 A number for matching efficiency.
  #' @param y A number for worker productivity.
  #' @param z A number for the flow value of nonwork.
  #' @param bbeta A number for the discount factor.
  #' @param s A number for the probability of separation.
  #' @param ct A number for the layoff tax.
  #' @param r A number for the interest rate.
  #' @param ch A number for the fixed cost of hiring paid by a firm.
  #' @param cl A number for the fixed cost of hiring paid by a worker.
  #' @param pphi A number for a worker's relative bargaining strength.
  #' @param ggamma A number for the curvature of the matching technology.
  #'
  #' @return A number for the cost of posting a vacancy associated with a given level of labor-market tightness.
  #' @export
  #'
  #' @examples
  #' c <- compute_c(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct, r = r, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma)
  if (!near(bbeta, 1 / (1+r))) {
    stop("!! Check discount factor and the interest rate")
  }
  find <- compute_find(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  fill <- compute_fill(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  c_T1 <- y - z - bbeta * s * ct - bbeta * (r + s) * ch / (1 - pphi) - bbeta * find * (pphi * ch / (1 - pphi) - cl)
  c_T2 <- (1 - pphi) * fill / (r + s + pphi * find)
  return(c_T2 * c_T1)
}

compute_w <- function(ttheta, m0, y, z, bbeta, s, r, c, ct, ch, cl, pphi, ggamma) {
  #' Compute the equilibrium wage.
  #'
  #' @param ttheta A number for labor-market tightness, ratio of vacancies to unemployed.
  #' @param m0 A number for matching efficiency.
  #' @param y A number for worker productivity.
  #' @param z A number for the flow value of nonwork.
  #' @param bbeta A number for the discount factor.
  #' @param s A number for the probability of separation.
  #' @param r A number for the interest rate.
  #' @param c A number for the flow cost of posting a vacancy.
  #' @param ct A number for the layoff tax.  
  #' @param ch A number for the fixed cost of hiring paid by a firm.
  #' @param cl A number for the fixed cost of hiring paid by a worker.
  #' @param pphi A number for a worker's relative bargaining strength.
  #' @param ggamma A number for the curvature of the matching technology.
  #'
  #' @return A number for the equilibrium wage.
  #' @export
  #'
  #' @examples
  #' c <- compute_wage(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct, r = r, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma)
  if (!near(bbeta, 1 / (1+r))) {
    stop("!! Check discount factor and the interest rate")
  }
  find <- compute_find(ttheta = ttheta, m0 = m0, ggamma = ggamma)

  T1 <- z + pphi * (y - z - bbeta * s * ct + ttheta * c)
  T2 <- bbeta * find * (pphi * ch - (1 - pphi) * cl)
  return(T1 + T2)
}

compute_ttheta_hi <- function(pphi, y, z, bbeta, r, s, c, cl, ch, ct) {
  #' Compute upper end of range for equilibrium tightness (ratio of vacancies to unemployment)
  #' 
  #' @description The steady-state level of equilibrium tightness, the ratio of vacancies to
  #' unemployment, is guaranteed to fall in a range. The upper end of the range is
  #' computed by this function. The computed number is useful because it can be
  #' passed to functions that compute equilibrium tightness.
  #'
  #' @param pphi A number in [0, 1] that represents a worker's bargaining strength.
  #' @param y A number for output produced by a worker.
  #' @param z A number for the flow value of nonwork.
  #' @param bbeta A number for the discount factor.
  #' @param r A number for the rate of interest (bbeta = (1 + r)^-1).
  #' @param s A number for the separation probability.
  #' @param c A number for the flow cost of posting a vacancy.
  #' @param cl A number for the fixed cost paid by a working of finding a job.
  #' @param ch A number for the fixed cost paid by a firm of hiring a worker.
  #' @param ct A number for the fixed cost paid by a worker when hired.
  #'
  #' @return A number equal to the end point of a range of values in which equilibrium tightness must fall.
  #' @export
  #'
  #' @examples
  #' ttheta_hi <- compute_ttheta_hi(pphi = pphi, y = y, z = z, bbeta = bbeta, r = r, s = s, c = c, ch = ch, ct = ct)
  if (!near(bbeta, 1 / (1+r))) {
    stop("!! Check discount factor and the interest rate")
  }
  T_top <- (1 - pphi) * (y - z - bbeta * s * ct + bbeta * cl) - bbeta * (r + s) * ch
  T_bottom <- c * pphi
  return(T_top / T_bottom)
}

compute_eqm_tight0 <- function(ttheta, m0, y, z, bbeta, s, ct, r, c, ch, cl, pphi, ggamma) {
  #' Function that takes deep parameters of the DMP model and returns a value that is zero in equilibrium.
  #' 
  #' The function can be used to take different guesses of labor-market tightness
  #' and search for equilibrium labor-market tightness, taking as given other deep
  #' parameters.
  #'
  #' @param ttheta A number for labor-market tightness, ratio of vacancies to unemployed.
  #' @param m0 A number for matching efficiency.
  #' @param y A number for worker productivity.
  #' @param z A number for the flow value of nonwork.
  #' @param bbeta A number for the discount factor.
  #' @param s A number for the probability of separation.
  #' @param ct A number for the layoff tax.
  #' @param r A number for the interest rate.
  #' @param c A number for the flow cost of posting a vacancy.
  #' @param ch A number for the fixed cost of hiring paid by a firm.
  #' @param cl A number for the fixed cost of hiring paid by a worker.
  #' @param pphi A number for a worker's relative bargaining strength.
  #' @param ggamma A number for the curvature of the matching technology.
  #'
  #' @return A value that is zero in equilibrium.
  #' @export
  #'
  #' @examples  
  if (!near(bbeta, 1 / (1+r))) {
    stop("!! Check discount factor and the interest rate")
  }
  find <- compute_find(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  fill <- compute_fill(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  T1 <- ((1 - pphi) / c) * (y - z - bbeta * s * ct - bbeta * (r + s) * ch / (1 - pphi)) 
  T2 <- (r + s) / fill + pphi * ttheta + bbeta * find * (pphi * ch - (1 - pphi) * cl) / c
  return(T1 - T2)
}

compute_eqm_tight <- function(m0, y, z, bbeta, s, ct, r, c, ch, cl, pphi, ggamma, mytol) {
  #' Compute equilibrium tightness for deep parameter values of DMP model with fixed costs.
  #' 
  #' The function searches for equilibrium labor-market tightness. The search
  #' relies on `compute_eqm_tight0()` and `uniroot()`.
  #'
  #' @param m0 A number for matching efficiency.
  #' @param y A number for worker productivity.
  #' @param z A number for the flow value of nonwork.
  #' @param bbeta A number for the discount factor.
  #' @param s A number for the probability of separation.
  #' @param ct A number for the layoff tax.
  #' @param r A number for the interest rate. 
  #' @param c A number for the flow cost of posting a vacancy.
  #' @param ch A number for the fixed cost of hiring paid by a firm.
  #' @param cl A number for the fixed cost of hiring paid by a worker.
  #' @param pphi A number for a worker's relative bargaining strength.
  #' @param ggamma A number for the curvature of the matching technology.
  #' @param mytol A number for tolerance passed to `uniroot()`.
  #'
  #' @return A number for equilibrium tightness given deep parameters of the DMP model with fixed costs.
  #' @export
  #'
  #' @examples
  if (!near(bbeta, 1 / (1+r))) {
    stop("!! Check discount factor and the interest rate")
  }  
  sol_tight <- uniroot(compute_eqm_tight0, 
                       y = y, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct, r = r, 
                       c = c, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma,
                       lower = 1e-5,
                       upper = compute_ttheta_hi(pphi = pphi, y = y, z = z, bbeta = bbeta, r = r, s = s, c = c, cl = cl, ch = ch, ct = ct),
                       extendInt = "no",
                       tol = mytol,
                       maxiter = 10000,
                       check.conv = TRUE)
  sol_tight$root
  # if (sol_tight$f.root < 2 * mytol) {
  #   return(sol_tight$root)
  # } else  {
  #   stop("!!! Check convergence !!!")
  # }
}

compute_elasticity_m_u <- function(ttheta, ggamma) {
  #' Compute elasticity of matching with respect to unemployment.
  #' 
  #' The function `compute_elasticity_m_u()` computes the elasticity of matching
  #' with respect to unemployment for a particular functional form. The functional
  #' form is:
  #'    q(ttheta) = m0 * (1 + ttheta^ggamma)^(-1 / ggamma)
  #'
  #' @param ttheta A number for labor-market tightness, ratio of vacancies to
  #'   unemployed.
  #' @param ggamma A number for the curvature of the matching technology.
  #'
  #' @return A number equal to the elasticity of matching with respect to
  #'   unemployment.
  #' @export
  #'
  #' @examples  
  (ttheta^ggamma) / (1 + ttheta^ggamma)
}

compute_T1 <- function(ttheta, m0, y, z, bbeta, s, ct, r, c, ch, cl, pphi, ggamma) {
  #' Compute first term in the two-factor decomposition of of the elasticity of
  #' labor-market tightness with respect to productivity.
  #'
  #' @param ttheta A number for labor-market tightness.
  #' @param m0 A number for matching efficiency.
  #' @param y A number for worker productivity.
  #' @param z A number for the flow value of nonwork.
  #' @param bbeta A number for the discount factor.
  #' @param s A number for the probability of separation.
  #' @param ct A number for the layoff tax.
  #' @param r A number for the interest rate. 
  #' @param c A number for the flow cost of posting a vacancy.
  #' @param ch A number for the fixed cost of hiring paid by a firm.
  #' @param cl A number for the fixed cost of hiring paid by a worker.
  #' @param pphi A number for a worker's relative bargaining strength.
  #' @param ggamma A number for the curvature of the matching technology.
  #' 
  #' @return A number equal to Upsilon, the first term of the two-factor
  #'   decomposition of the elasticity of labor-market tightness with respect to
  #'   productivity.
  #' @export
  #'
  #' @examples  
  if (!near(bbeta, 1 / (1+r))) {
    stop("!! Check discount factor and the interest rate")
  }  
  find <- compute_find(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  fill <- compute_fill(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  eta_m_u <- compute_elasticity_m_u(ttheta = ttheta, ggamma = ggamma)
  
  cost_term <- (pphi * ch - (1 - pphi) * cl) / c
  
  Upsilon_T1 <- r + s + find * (pphi + bbeta * fill * cost_term)
  Upsilon_T2 <- (r + s) * eta_m_u + find * (pphi + bbeta * (1 - eta_m_u) * fill * cost_term)
  Upsilon <- Upsilon_T1 / Upsilon_T2
  
  return(Upsilon)
}

compute_elasticity <- function(ttheta, m0, y, z, bbeta, s, ct, r, c, ch, cl, pphi, ggamma) {
  #' Compute elasticity of labor-market tightness with respect to productivity.
  #'
  #' @param ttheta A number for labor-market tightness.
  #' @param m0 A number for matching efficiency.
  #' @param y A number for worker productivity.
  #' @param z A number for the flow value of nonwork.
  #' @param bbeta A number for the discount factor.
  #' @param s A number for the probability of separation.
  #' @param ct A number for the layoff tax.
  #' @param r A number for the interest rate. 
  #' @param c A number for the flow cost of posting a vacancy.
  #' @param ch A number for the fixed cost of hiring paid by a firm.
  #' @param cl A number for the fixed cost of hiring paid by a worker.
  #' @param pphi A number for a worker's relative bargaining strength.
  #' @param ggamma A number for the curvature of the matching technology.
  #' 
  #' @return A number equal to the elasticity of labor-market tightness with
  #'   respect to productivity.
  #' @export
  #'
  #' @examples  
  if (!near(bbeta, 1 / (1+r))) {
    stop("!! Check discount factor and the interest rate")
  }  

  Upsilon <- compute_T1(ttheta, m0, y, z, bbeta, s, ct, r, c, ch, cl, pphi, ggamma)  
  eta_ttheta_y_T1 <- Upsilon 
  eta_ttheta_y_T2 <- y / (y - z - bbeta * s * ct - bbeta * (r + s) * ch / (1 - pphi))
  eta_ttheta_y <- eta_ttheta_y_T1 * eta_ttheta_y_T2
  return(eta_ttheta_y)
}

compute_dwdy <- function(ttheta, m0, y, z, bbeta, s, ct, r, c, ch, cl, pphi, ggamma) {
  find <- compute_find(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  fill <- compute_fill(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  find_prime <- compute_find_prime(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  eta_m_u <- compute_elasticity_m_u(ttheta = ttheta, ggamma = ggamma)
  
  cost_term <- (pphi * ch * (1 - pphi) * cl) / c
  T1 <- pphi * (1 - pphi) * find + bbeta * (1 - pphi) * find * cost_term * find_prime
  T2 <- (r + s) * eta_m_u + find * (pphi + bbeta * (1 - eta_m_u) * fill * cost_term)
  dwdy <- pphi + T1 / T2
}

compute_wage_elasticity <- function(ttheta, m0, y, z, bbeta, s, ct, r, c, ch, cl, pphi, ggamma) {
  find <- compute_find(ttheta = ttheta, m0 = m0, ggamma = ggamma)
  dwdy <- compute_dwdy(ttheta = ttheta, m0 = m0, y = y, z = z, bbeta = bbeta, s = s, ct = ct, r = r, c = c, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma)
  
  wage <- z + pphi * (y - z - bbeta * s * ct + ttheta * c) + bbeta * find * (pphi * ch + (1 - pphi) * cl)
  
  ret <- dwdy * (y / wage)
  return(ret)
}

# For all the economies, ct = 0
ct <- 0.0

# BASELINE Flow-cost of maintaining a vacancy for baseline with no fixed costs
cl_baseline <- 0.0
ch_baseline <- 0.0

# HIGH H: Maximum fixed costs paid by firms
ct_hhi <- ct
ch_hhi_T1 <- (1 - pphi) * (y - z - bbeta * s * ct_hhi + bbeta * compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma) * cl_baseline)
ch_hhi_T2 <- bbeta * (r + s + pphi * compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma))
ch_hhi <- ch_hhi_T1 / ch_hhi_T2 - 0.001
cl_hhi <- cl_baseline

# MIDDLE H: Half the fixed costs paid by firms
# ct_hhi2 <- 0.0
ch_hhi2 <- ch_hhi / 2.0
cl_hhi2 <- cl_baseline

# SPLIT: Split fixed cost of matching between worker and firm
cH <- ch_hhi
ppsi <- 0.5
# ct_ppsi <- 0.0
ch_ppsi <- (1 - ppsi) * cH 
cl_ppsi <- ppsi * cH

dat_tab <- tribble(
  ~economy,    ~case,       ~ch,         ~cl,         
  "Baseline",   "baseline",  ch_baseline, cl_baseline,
  "Middle $h$", "middle_h",  ch_hhi2,     cl_hhi2,
  "High $h$",   "high_h",    ch_hhi,      cl_hhi,   
  "Split",      "split",     ch_ppsi,     cl_ppsi    
)

# Compute c
dat_tab <- dat_tab |> 
  mutate(c = pmap_dbl(list(cl = dat_tab$cl, dat_tab$ch), 
                            # function
                            compute_c,
                            # parameters that do not vary
                            ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma, ct = ct)) 


# Compute elasticities and equilibrium wage
dat_tab <- dat_tab |> 
  mutate(wage = pmap_dbl(list(c = dat_tab$c, ch = dat_tab$ch, cl = dat_tab$cl),
                         compute_w,
                         # parameters that do not vary
                         ttheta = target_tight, m0 = match_efficiency, 
                         y = y, z = z, bbeta = bbeta, s = s, r = r, 
                         pphi = pphi, ggamma = ggamma,
                         ct = ct),
         elasticity = pmap_dbl(list(c = dat_tab$c, ch = dat_tab$ch, cl = dat_tab$cl),
                               # function
                               compute_elasticity,
                               # parameters that do not vary
                               ttheta = target_tight, m0 = match_efficiency, 
                               y = y, z = z, bbeta = bbeta, s = s, r = r, 
                               pphi = pphi, ggamma = ggamma,
                               ct = ct),
         elasticity_w = pmap_dbl(list(c = dat_tab$c, ch = dat_tab$ch, cl = dat_tab$cl),
                                 # function
                                 compute_wage_elasticity,
                                 # parameters that do not vary
                                 ttheta = target_tight, m0 = match_efficiency, 
                                 y = y, z = z, bbeta = bbeta, s = s, r = r, 
                                 pphi = pphi, ggamma = ggamma,
                                 ct = ct))


check_dwdy0 <- pmap_dbl(list(c = dat_tab$c, ch = dat_tab$ch, cl = dat_tab$cl),
  compute_dwdy,
  ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, ggamma = ggamma, ct = ct, 
  pphi = 0)
stopifnot(near(check_dwdy0, 0))

check_dwdy1 <- pmap_dbl(list(c = dat_tab$c, ch = dat_tab$ch, cl = dat_tab$cl),
  compute_dwdy,
  ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, ggamma = ggamma, ct = ct, 
                     pphi = 1)
stopifnot(near(check_dwdy1, 1))

check_target_tight <- pmap_dbl(
  # Arguments that vary go first
  list(c = dat_tab$c,
       ch = dat_tab$ch,
       cl = dat_tab$cl),
  # function
  compute_eqm_tight,
  # non-varying arguments
  m0 = match_efficiency,
  y = y,
  z = z,
  bbeta = bbeta,
  s = s,
  ct = 0,
  r = r,
  pphi = pphi,
  ggamma = ggamma,
  mytol = TOL
)

stopifnot(near(check_target_tight, target_tight))

make_pretty_num <- function(x, ndigits, smalln) {
  fpretty <- function(z, smalln, ndigits) {
    candidaten <- as.numeric(formatC(z, digits = ndigits, format = "f"))
    if (near(z, 0)) {
      formatC(z, format = "d")
    } else if (candidaten < smalln) {
      # Scientific notation for small numbers
      formatC(z, digits = ndigits, format = "e")
    } else {
      candidaten
    }
  }
  sapply(x, fpretty, smalln = smalln, ndigits = ndigits)
}
  
# Save table as LaTeX output
ndigits <- 3
smalln <- 0.01
dat_tbl_output <- dat_tab %>% 
  select(-case, -wage) %>% 
  mutate(across(-economy, ~ make_pretty_num(.x, ndigits = ndigits, smalln = smalln))) %>%
  rename(`{$c$}` = c,
         `{$h$}` = ch,
         `{$\\ell$}` = cl,
         `{$\\eta_{\\theta,y}$}` = elasticity,
         `{$\\eta_{w,y}$}` = elasticity_w,
         Economy = economy)

# Create .tex table 
tbl_fout <- here("out", paste0("tbl_", file_prg, "-experiment", ".tex"))
kableExtra::kable(dat_tbl_output, format = "latex", booktabs = TRUE, align = c("l", rep("S", 5)),
                  caption = "\\label{tab:model-results} Model results at different combinations of job-creation costs.",
                  escape = FALSE,
                  table.envir = "table") %>% 
  kableExtra::footnote(general  = "The column $\\\\eta_{\\\\theta,y}$ is the elasticity of market tightness with respect to productivity. 
                       The column $\\\\eta_{w,y}$ is the elasticity of wages with respect to productivity.",
                       general_title = "Notes: ",
                       footnote_as_chunk = TRUE,
                       title_format = c("italic"),
                       escape = FALSE,
                       threeparttable = TRUE) %>%
    kableExtra::save_kable(file = tbl_fout)  


# Table: Model results with lower z ---------------------------------------

# BASELINE Flow-cost of maintaining a vacancy for baseline with no fixed costs
cl_baseline_shimer2005 <- 0.0
ch_baseline_shimer2005 <- 0.0

# HIGH H: Maximum fixed costs paid by firms
ct_hhi_shimer2005 <- ct
ch_hhi_T1_shimer2005 <- (1 - pphi) * (y - z_shimer2005 - bbeta * s * ct_hhi_shimer2005 + bbeta * compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma) * cl_baseline)
ch_hhi_T2_shimer2005 <- bbeta * (r + s + pphi * compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma))
ch_hhi_shimer2005 <- ch_hhi_T1_shimer2005 / ch_hhi_T2_shimer2005 - 0.001
cl_hhi_shimer2005 <- cl_baseline

# MIDDLE H: Half the fixed costs paid by firms
ch_hhi2_shimer2005 <- ch_hhi_shimer2005 / 2.0
cl_hhi2_shimer2005 <- cl_baseline_shimer2005

# SPLIT: Split fixed cost of matching between worker and firm
cH_shimer2005 <- ch_hhi_shimer2005
ppsi <- 0.5
ch_ppsi_shimer2005 <- (1 - ppsi) * cH_shimer2005
cl_ppsi_shimer2005 <- ppsi * cH_shimer2005

dat_tab_shimer2005 <- tribble(
  ~economy,    ~case,       ~ch,                     ~cl,         
  "Baseline",   "baseline",  ch_baseline_shimer2005, cl_baseline_shimer2005,
  "Middle $h$", "middle_h",  ch_hhi2_shimer2005,     cl_hhi2_shimer2005,
  "High $h$",   "high_h",    ch_hhi_shimer2005,      cl_hhi_shimer2005,   
  "Split",      "split",     ch_ppsi_shimer2005,     cl_ppsi_shimer2005    
)

# Compute c
dat_tab_shimer2005 <- dat_tab_shimer2005 |> 
  mutate(c = pmap_dbl(list(cl = dat_tab_shimer2005$cl, dat_tab_shimer2005$ch), 
                      # function
                      compute_c,
                      # parameters that do not vary
                      ttheta = target_tight, m0 = match_efficiency, y = y, z = z_shimer2005, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma, ct = ct)) 


# Compute elasticities and equilibrium wage
dat_tab_shimer2005 <- dat_tab_shimer2005 |> 
  mutate(wage = pmap_dbl(list(c = dat_tab_shimer2005$c, ch = dat_tab_shimer2005$ch, cl = dat_tab_shimer2005$cl),
                         compute_w,
                         # parameters that do not vary
                         ttheta = target_tight, m0 = match_efficiency, 
                         y = y, z = z_shimer2005, bbeta = bbeta, s = s, r = r, 
                         pphi = pphi, ggamma = ggamma,
                         ct = ct),
         elasticity = pmap_dbl(list(c = dat_tab_shimer2005$c, ch = dat_tab_shimer2005$ch, cl = dat_tab_shimer2005$cl),
                               # function
                               compute_elasticity,
                               # parameters that do not vary
                               ttheta = target_tight, m0 = match_efficiency, 
                               y = y, z = z_shimer2005, bbeta = bbeta, s = s, r = r, 
                               pphi = pphi, ggamma = ggamma,
                               ct = ct),
         elasticity_w = pmap_dbl(list(c = dat_tab_shimer2005$c, ch = dat_tab_shimer2005$ch, cl = dat_tab_shimer2005$cl),
                                 # function
                                 compute_wage_elasticity,
                                 # parameters that do not vary
                                 ttheta = target_tight, m0 = match_efficiency, 
                                 y = y, z = z_shimer2005, bbeta = bbeta, s = s, r = r, 
                                 pphi = pphi, ggamma = ggamma,
                                 ct = ct))


check_dwdy0 <- pmap_dbl(list(c = dat_tab_shimer2005$c, ch = dat_tab_shimer2005$ch, cl = dat_tab_shimer2005$cl),
                        compute_dwdy,
                        ttheta = target_tight, m0 = match_efficiency, y = y, z = z_shimer2005, bbeta = bbeta, s = s, r = r, ggamma = ggamma, ct = ct, 
                        pphi = 0)
stopifnot(near(check_dwdy0, 0))

check_dwdy1 <- pmap_dbl(list(c = dat_tab_shimer2005$c, ch = dat_tab_shimer2005$ch, cl = dat_tab_shimer2005$cl),
                        compute_dwdy,
                        ttheta = target_tight, m0 = match_efficiency, y = y, z = z_shimer2005, bbeta = bbeta, s = s, r = r, ggamma = ggamma, ct = ct, 
                        pphi = 1)
stopifnot(near(check_dwdy1, 1))

check_target_tight <- pmap_dbl(
  # Arguments that vary go first
  list(c = dat_tab_shimer2005$c,
       ch = dat_tab_shimer2005$ch,
       cl = dat_tab_shimer2005$cl),
  # function
  compute_eqm_tight,
  # non-varying arguments
  m0 = match_efficiency,
  y = y,
  z = z_shimer2005,
  bbeta = bbeta,
  s = s,
  ct = 0,
  r = r,
  pphi = pphi,
  ggamma = ggamma,
  mytol = TOL
)

stopifnot(near(check_target_tight, target_tight))

# Save table as LaTeX output
dat_tbl_output_shimer2005 <- dat_tab_shimer2005 %>% 
  select(-case, -wage) %>% 
  mutate(across(-economy, ~ make_pretty_num(.x, ndigits = ndigits, smalln = smalln))) %>%
  rename(`{$c$}` = c,
         `{$h$}` = ch,
         `{$\\ell$}` = cl,
         `{$\\eta_{\\theta,y}$}` = elasticity,
         `{$\\eta_{w,y}$}` = elasticity_w,
         Economy = economy)

# Create .tex table 
tbl_fout_shimer2005 <- here("out", paste0("tbl_", file_prg, "-experiment-shimer2005", ".tex"))
my_general_footnote = paste0("The column $\\\\eta_{\\\\theta,y}$ is the elasticity of market tightness with respect to productivity. The column $\\\\eta_{w,y}$ is the elasticity of wages with respect to productivity. The value of nonwork is $", z_shimer2005, "$.")
kableExtra::kable(dat_tbl_output_shimer2005, format = "latex", booktabs = TRUE, align = c("l", rep("S", 5)),
                  caption = "\\label{tab:model-results-shimer2005} Model results at different combinations of job-creation costs.",
                  escape = FALSE,
                  table.envir = "table") %>% 
  kableExtra::footnote(general  = ,
                       general_title = "Notes: ",
                       footnote_as_chunk = TRUE,
                       title_format = c("italic"),
                       escape = FALSE,
                       threeparttable = TRUE) %>%
  kableExtra::save_kable(file = tbl_fout_shimer2005)  

# Steady-state dynamics --------------------------------------------------

# Productivity
sizey <- 100
v_prod <- seq(0.8 * y, 1.2 * y, length.out = sizey)

TOL <- 1e-8

compute_dynamics <- function(v_prod, m0, z, bbeta, s, c, ch, cl, ct, r, pphi, ggamma, mytol) {
  df <- tibble(ttheta = map_vec(v_prod, compute_eqm_tight, c = c, ch = ch, cl = cl, m0 = m0, z = z, bbeta = bbeta, s = s, ct= ct, r = r, pphi = pphi, ggamma = ggamma, mytol = mytol),
               y = v_prod)
  return(df)
}

dat_dynamics <- dat_tab %>% 
  mutate(dynamics = pmap(list(c = c, ch = ch, cl = cl), 
                         compute_dynamics, 
                         # Arguments that do not change
                         v_prod = v_prod, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct, r = r, pphi = pphi, ggamma = ggamma, mytol = TOL)) 

dat_long <- dat_dynamics %>% 
  unnest(dynamics) %>% 
  mutate(find = compute_find(ttheta = ttheta, m0 = match_efficiency, ggamma = ggamma),
         ur = s_monthly / (s_monthly + compute_daily2monthly(find)),
         Upsilon = compute_T1(ttheta = ttheta, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct_baseline, r = r, c = c, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma),
         elasticity_m_u = compute_elasticity_m_u(ttheta = ttheta, ggamma = ggamma),
         elasticity_m_u_hi = compute_elasticity_m_u(ttheta = ttheta, ggamma = 1.27),
         Upsilon_bound = pmax(1 / elasticity_m_u, 1 / (1 - elasticity_m_u)),
         Upsilon_bound_hi = pmax(1 / elasticity_m_u_hi, 1 / (1 - elasticity_m_u_hi)),) %>% 
  verify(Upsilon <= Upsilon_bound) 

dat_long_label <- dat_long %>% 
  group_by(case) %>% 
  mutate(id = row_number()) %>% 
  mutate(my_label = case_when(
    case == "baseline" & id == 10 ~ economy,
    case == "middle_h" & id == 25 ~ "plain(Middle)~italic('h')",
    case == "high_h" & id == 10 ~ "plain(High)~italic('h')",
    case == "split" & id == 75 ~ economy,
    TRUE ~ ""
  )) %>% 
  filter(my_label != "")

ggplot(data = dat_long) +
  geom_hline(yintercept = s_monthly / (s_monthly + target_find_monthly))  +
  geom_line(mapping = aes(x = y, y = ur, color = case, linetype = case), linewidth = 1.2) +
  geom_text(mapping = aes(x = 0.85, y = s_monthly / (s_monthly + target_find_monthly)),
            label = paste0(round(s_monthly / (s_monthly + target_find_monthly),3)),
            nudge_y = -0.005) +
  geom_text_repel(data = dat_long_label,
                  mapping = aes(x = y, y = ur, label = my_label, color = case),
                  max.overlaps = Inf,
                  nudge_y = 0.05,
                  box.padding = 3.0,
                  parse = TRUE) +
  labs(x = expression(paste(plain(Productivity), ", ",  italic('y'))), 
       y = expression(paste(plain(Unemployment), " ", plain(rate), ", ",  italic('u')))) +
  scale_color_viridis_d(begin = 0.0, end = 0.85) +
  guides(color = "none", linetype = "none") +
  theme_minimal()

fout <- paste0("fig_", file_prg, ".pdf")
ggsave(here("out", fout), heigh = myheight, width = mywidth)

ggplot(data = dat_long) +
  geom_line(mapping = aes(x = y, y = Upsilon, color = case)) +
  geom_line(mapping = aes(x = y, y = Upsilon_bound, color = case, linetype = case))


# Figure: Bound for different levels of productivity ------------------------------

ggplot(data = dat_long) +
  geom_hline(yintercept = 2)  +
  geom_line(mapping = aes(x = y, y = Upsilon_bound_hi, color = case, linetype = case), linewidth = 1.2) +
  geom_text_repel(data = dat_long_label,
                  mapping = aes(x = y, y = Upsilon_bound_hi, label = my_label, color = case),
                  max.overlaps = Inf,
                  nudge_y = 0.05,
                  box.padding = 3.0,
                  parse = TRUE) +
  labs(x = expression(paste(plain(Productivity), ", ",  italic('y'))),
       y = "Bound") +
  scale_color_viridis_d(begin = 0.0, end = 0.85) +
  guides(color = "none", linetype = "none") +
  theme_minimal()

fout <- paste0("fig_", file_prg, "_bound-y.pdf")
ggsave(here("out", fout), heigh = myheight, width = mywidth)

# Wage Curve --------------------------------------------------------------

npts_wage <- 10000

ct_baseline <- ct
c_baseline <- dat_tab |> filter(economy == "Baseline") |> pull(c)

dat_eqm <- tibble(
  ttheta = seq(0.1 * target_tight, 1.9 * target_tight, length.out = npts_wage),
  job_creation = y - bbeta * s * ct_baseline - (r + s) * c_baseline / compute_fill(ttheta = ttheta, m0 = match_efficiency, ggamma = ggamma) - (r + s)/(1 + r) * ch_baseline,
  wage_curve_baseline = z + pphi * (y - z - bbeta * s * ct_baseline + ttheta * c_baseline) + bbeta * compute_find(ttheta = ttheta, m0 = match_efficiency, ggamma = ggamma) * (pphi * ch_baseline - (1 - pphi) * cl_baseline),
  wage_curve_cl = z + pphi * (y - z - bbeta * s * ct_baseline + ttheta * c_baseline) + bbeta * compute_find(ttheta = ttheta, m0 = match_efficiency, ggamma = ggamma) * (pphi * ch_baseline - (1 - pphi) * cl_ppsi)
) %>% 
  mutate(id = row_number(),
         mylabel_job_creation = case_when(
           id == 25 ~ "Job creation by firms",
           TRUE ~ ""
           ),
         mylabel_wage_curve = case_when(
           id == 8000 ~ "Job creation by workers",
           TRUE ~ ""
         ),
         mylabel_wage_hi_cl = case_when(
           id == 6000 ~ "A higher job-creation cost paid by workers\nlowers wages and\nincreases labor-market tightness",
           TRUE ~ ""
         ))

eqm_wages_size <- 1.1
w_eqm <- y - bbeta * s * ct_baseline - (r + s) * c_baseline / compute_fill(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma) - bbeta * (r + s) * ch_baseline
ggplot(data = dat_eqm) +
  geom_line(mapping = aes(x = ttheta, y = job_creation), color = "black", linewidth = eqm_wages_size) +
  geom_line(mapping = aes(x = ttheta, y = wage_curve_cl), color = csub_gray, linewidth = eqm_wages_size, linetype = "dashed") +
  geom_line(mapping = aes(x = ttheta, y = wage_curve_baseline), color = csub_blue, linewidth = eqm_wages_size) +
  geom_text_repel(mapping = aes(x = ttheta, y = job_creation, label = mylabel_job_creation),
                  max.overlaps = Inf) +
  geom_text_repel(mapping = aes(x = ttheta, y = wage_curve_baseline, label = mylabel_wage_curve),
                  max.overlaps = Inf) +  
  geom_text_repel(mapping = aes(x = ttheta, y = wage_curve_cl, label = mylabel_wage_hi_cl),
                  max.overlaps = Inf) +
  geom_vline(xintercept = target_tight, 
             linetype = "dotted") +
  geom_hline(yintercept = w_eqm,
             linetype = "dotted") +
  geom_text(mapping = aes(x = target_tight, y = min(wage_curve_baseline), label = paste0(target_tight)), nudge_y = -0.05) + 
  geom_text(mapping = aes(x = min(ttheta), y = w_eqm, label = paste0(round(w_eqm, 3))), nudge_y = 0) +   
  theme_minimal() +
  labs(x = expression(paste(plain('Labor-market tightness'), ", ",  theta)), 
       y = expression(paste(plain('Equilibrium wages,'), " ", italic('w')))) 

fout_eqm_w_tightness <- paste0("fig_", file_prg, "_eqm-wages-mkt-tightness.pdf")
ggsave(here("out", fout_eqm_w_tightness), heigh = 1.4 * myheight, width = mywidth)


# Write parameters for LaTeX ----------------------------------------------

# File for LaTeX parameters
file_latex <- here("out", paste0("tex_", file_prg, ".tex"))
writeLines(paste0("% Numbers generated by ", file_prg), file_latex)
# Open connection to append
CON <- file(file_latex, "a") 

writeLines(paste0("\\newcommand{\\y}{", y, "}"), con = CON)
writeLines(paste0("\\newcommand{\\z}{", z, "}"), con = CON)
writeLines(paste0("\\newcommand{\\zShimer2005}{", z_shimer2005, "}"), con = CON)
writeLines(paste0("\\newcommand{\\pphi}{", pphi, "}"), con = CON)
writeLines(paste0("\\newcommand{\\rmonthly}{", r_monthly, "}"), con = CON)
writeLines(paste0("\\newcommand{\\smonthly}{", s_monthly, "}"), con = CON)
writeLines(paste0("\\newcommand{\\targetTight}{", target_tight, "}"), con = CON)
writeLines(paste0("\\newcommand{\\findMonthly}{", find_monthly, "}"), con = CON)
writeLines(paste0("\\newcommand{\\uMonthly}{", round(100 * s_monthly / (s_monthly + find_monthly), digits = 1), "}"), con = CON)
# Baseline parameters
economy_baseline <- dat_tab |> 
  filter(economy == "Baseline")
writeLines(paste0("\\newcommand{\\ttau}{", ct, "}"), con = CON)
writeLines(paste0("\\newcommand{\\chBaseline}{", economy_baseline$ch, "}"), con = CON)
writeLines(paste0("\\newcommand{\\clBaseline}{", economy_baseline$cl, "}"), con = CON)
writeLines(paste0("\\newcommand{\\wBaseline}{", round(economy_baseline$wage, 3), "}"), con = CON)
writeLines(paste0("\\newcommand{\\flowGain}{", round(100 * (economy_baseline$wage / z - 1), 0), "}"), con = CON)

close(CON)