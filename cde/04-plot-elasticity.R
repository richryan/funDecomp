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

y <- 1 
z <- 0.71 
pphi <- 0.50

target_elasticity <- 7.5636

ggamma <- readRDS(here("out", "dat_03-est-matching.rds"))

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
stopifnot(all(near(numDeriv::grad(func = function(x) compute_find(x, m0 = 1.0, ggamma = ggamma), x = x_check) - compute_find_prime(x_check, m0 = 1.0, ggamma = ggamma), 0)))

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

# Steady-state dynamics --------------------------------------------------

compute_ttheta_hi <- function(pphi, y, z, bbeta, r, s, c, cl, ch, ct) {
  #' Compute upper end of range for equilibrium tightness (ratio of vacancies to unemployment)
  #' 
  #' The steady-state level of equilibrium tightness, the ratio of vacancies to
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


# Check that compute_eqm_tight() returns target_tight
# TOL <- 1e-8
# check_tight <- compute_eqm_tight(m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct, r = r, c = c, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma, mytol = TOL)
# sol_tight <- uniroot(compute_eqm_tight0,
#   y = y, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct, r = r, 
#   c = c, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma,
#   lower = 1e-5,
#   upper = compute_ttheta_hi(pphi = pphi, y = y, z = z, bbeta = bbeta, r = r, s = s, c = c, ch = ch, ct = ct),
#   tol = TOL
#   )
# cat("Difference between numerical and target tight:", sol_tight$root - target_tight)
# cat("Difference between numerical and target tight:", check_tight - target_tight)
# cat("Difference between numerical and target tight:", check_tight - sol_tight$root)

# Productivity
sizey <- 100
v_prod <- seq(0.8 * y, 1.2 * y, length.out = sizey)

TOL <- 1e-8

# BASELINE Flow-cost of maintaining a vacancy for baseline with no fixed costs
cl_baseline <- 0.0
ch_baseline <- 0.0
ct_baseline <- 0.0
c_baseline <- compute_c(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct_baseline, r = r, ch = ch_baseline, cl = cl_baseline, pphi = pphi, ggamma = ggamma)
ttheta_baseline <- compute_eqm_tight(m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct_baseline, r = r, c = c_baseline, ch = ch_baseline, cl = cl_baseline, pphi = pphi, ggamma = ggamma, mytol = TOL)
stopifnot(near(ttheta_baseline, target_tight))
elasticity_baseline <- compute_elasticity(ttheta = ttheta_baseline, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma,
                                          c = c_baseline, ct = ct_baseline, ch = ch_baseline, cl = cl_baseline)
elasticity_w_baseline <- compute_wage_elasticity(ttheta = ttheta_baseline, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct_baseline, r = r, c = c_baseline, ch = ch_baseline, cl = cl_baseline, pphi = pphi, ggamma = ggamma)

# Check to see uniqueness condition holds
# ttheta_bar <- (1 - pphi) * (y - z - bbeta * s * ct_baseline - bbeta * (r + s) * ch_baseline / (1 - pphi) + bbeta * cl_baseline) / pphi / c_baseline
# find_bar <- compute_find(ttheta = target_bar, m0 = match_efficiency, ggamma = ggamma)
# fill_bar <- compute_fill(ttheta = target_bar, m0 = match_efficiency, ggamma = ggamma)
# find_prime_bar <- compute_find_prime(ttheta = ttheta_bar, m0 = match_efficiency, ggamma = ggamma)
# check_unique <- pphi * find_bar / (r + s + pphi * find_bar) - find_prime_bar / fill_bar
# stopifnot(check_unique < 0)

# Confirm that dwdy evaluates to 1 when pphi = 1 and 0 when pphi = 0
dwdy0 <- compute_dwdy(ttheta = ttheta_baseline, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct_baseline, r = r, c = c_baseline, ch = ch_baseline, cl = cl_baseline, ggamma = ggamma,
                     pphi = 0)
stopifnot(near(dwdy0, 0))
dwdy1 <- compute_dwdy(ttheta = ttheta_baseline, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct_baseline, r = r, c = c_baseline, ch = ch_baseline, cl = cl_baseline, ggamma = ggamma,
                     pphi = 1)
stopifnot(near(dwdy1, 1))

# HIGH H: Maximum fixed costs paid by firms
ct_hhi <- 0.0
ch_hhi_T1 <- (1 - pphi) * (y - z - bbeta * s * ct_hhi + bbeta * compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma) * cl_baseline)
ch_hhi_T2 <- bbeta * (r + s + pphi * compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma))
ch_hhi <- ch_hhi_T1 / ch_hhi_T2 - 0.001
cl_hhi <- cl_baseline
c_hhi <- compute_c(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, ch = ch_hhi, ct = ct_hhi, cl = cl_hhi, pphi = pphi, ggamma = ggamma) 
ttheta_hhi <- compute_eqm_tight(m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct_hhi, r = r, c = c_hhi, ch = ch_hhi, cl = cl_hhi, pphi = pphi, ggamma = ggamma, mytol = TOL)
stopifnot(near(ttheta_hhi, target_tight))
elasticity_hhi <- compute_elasticity(ttheta = ttheta_baseline, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma,
                                     c = c_hhi, ct = ct_hhi, ch = ch_hhi, cl = cl_hhi)
elasticity_w_hhi <- compute_wage_elasticity(ttheta = ttheta_hhi, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma,
                                            c = c_hhi, ct = ct_hhi, ch = ch_hhi, cl = cl_hhi)

# MIDDLE H: Half the fixed costs paid by firms
ct_hhi2 <- 0.0
ch_hhi2 <- ch_hhi / 2.0
cl_hhi2 <- cl_baseline
c_hhi2 <- compute_c(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma,
                    ch = ch_hhi2, ct = ct_hhi2, cl = cl_hhi2) 
ttheta_hhi2 <- compute_eqm_tight(m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma, mytol = TOL,
                                c = c_hhi2, ct = ct_hhi2, ch = ch_hhi2, cl = cl_baseline)
stopifnot(near(ttheta_hhi2, target_tight))
elasticity_hhi2 <- compute_elasticity(ttheta = ttheta_baseline, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma,
                                      c = c_hhi2, ct = ct_hhi2, ch = ch_hhi2, cl = cl_hhi2)
elasticity_w_hhi2 <- compute_wage_elasticity(ttheta = ttheta_hhi, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma,
                                             c = c_baseline, ct = ct_hhi, ch = ch_hhi2, cl = cl_hhi2)

# SPLIT: Split fixed cost of matching between worker and firm
cH <- ch_hhi
ppsi <- 0.5
ct_ppsi <- 0.0
ch_ppsi <- (1 - ppsi) * cH 
cl_ppsi <- ppsi * cH
c_ppsi <- compute_c(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, ch = ch_ppsi, ct = ct_ppsi, cl = cl_ppsi, pphi = pphi, ggamma = ggamma) 
ttheta_ppsi <- compute_eqm_tight(m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct_ppsi, r = r, c = c_ppsi, ch = ch_ppsi, cl = cl_ppsi, pphi = pphi, ggamma = ggamma, mytol = TOL)
stopifnot(near(ttheta_ppsi, target_tight))
elasticity_ppsi <- compute_elasticity(ttheta = ttheta_baseline, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma,
                                      c = c_ppsi, ct = ct_ppsi, ch = ch_ppsi, cl = cl_ppsi)
elasticity_w_ppsi <- compute_wage_elasticity(ttheta = ttheta_hhi, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, pphi = pphi, ggamma = ggamma,
                                             c = c_ppsi, ct = ct_ppsi, ch = ch_ppsi, cl = cl_ppsi)

# Confirm that dwdy evaluates to 1 when pphi = 1 and 0 when pphi = 0
dwdy0 <- compute_dwdy(ttheta = ttheta_ppsi, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, ggamma = ggamma,
                      c = c_ppsi, ct = ct_ppsi, ch = ch_ppsi, cl = cl_ppsi, 
                     pphi = 0)
stopifnot(near(dwdy0, 0))
dwdy1 <- compute_dwdy(ttheta = ttheta_baseline, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, ggamma = ggamma,
                      c = c_ppsi, ct = ct_ppsi, ch = ch_ppsi, cl = cl_ppsi, 
                     pphi = 1)
stopifnot(near(dwdy1, 1))

dat <- tribble(
  ~economy,    ~case,      ~c,         ~ch,         ~cl,         ~elasticity,          ~elasticity_w,
  "Baseline",   "baseline", c_baseline, ch_baseline, cl_baseline, elasticity_baseline,  elasticity_w_baseline,
  "Middle $h$", "middle_h", c_hhi2,     ch_hhi2,     cl_hhi2,     elasticity_hhi2,      elasticity_w_hhi2,
  "High $h$",   "high_h",   c_hhi,      ch_hhi,      cl_hhi,      elasticity_hhi,       elasticity_w_hhi,
  "Split",      "split",    c_ppsi,     ch_ppsi,     cl_ppsi,     elasticity_ppsi,      elasticity_w_ppsi
)

compute_dynamics <- function(v_prod, m0, z, bbeta, s, c, ch, cl, ct, r, pphi, ggamma, mytol) {
  df <- tibble(ttheta = map_vec(v_prod, compute_eqm_tight, c = c, ch = ch, cl = cl, m0 = m0, z = z, bbeta = bbeta, s = s, ct= ct, r = r, pphi = pphi, ggamma = ggamma, mytol = mytol),
               y = v_prod)
  return(df)
}

dat <- dat %>% 
  mutate(dynamics = pmap(list(c = c, ch = ch, cl = cl), compute_dynamics, v_prod = v_prod, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct_baseline, r = r, pphi = pphi, ggamma = ggamma, mytol = TOL)) 

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
dat_tbl <- dat %>% 
  select(-dynamics, -case) %>% 
  mutate(across(-economy, ~ make_pretty_num(.x, ndigits = ndigits, smalln = smalln))) %>%
  rename(`{$c$}` = c,
         `{$h$}` = ch,
         `{$\\ell$}` = cl,
         `{$\\eta_{\\theta,y}$}` = elasticity,
         `{$\\eta_{w,y}$}` = elasticity_w,
         Economy = economy)

# Create .tex table 
tbl_out <- here("out", paste0("tbl_", file_prg, "-experiment", ".tex"))
kableExtra::kable(dat_tbl, format = "latex", booktabs = TRUE, align = c("l", rep("S", 5)),
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
    kableExtra::save_kable(file = tbl_out)  
  
dat_long <- dat %>% 
  unnest(dynamics) %>% 
  mutate(find = compute_find(ttheta = ttheta, m0 = match_efficiency, ggamma = ggamma),
         ur = s_monthly / (s_monthly + compute_daily2monthly(find)),
         Upsilon = compute_T1(ttheta = ttheta, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = ct_baseline, r = r, c = c, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma),
         elasticity_m_u = compute_elasticity_m_u(ttheta = ttheta, ggamma = ggamma),
         Upsilon_bound = pmax(1 / elasticity_m_u, 1 / (1 - elasticity_m_u))) %>% 
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


# Wage Curve --------------------------------------------------------------

npts_wage <- 10000

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

fout_eqm_w_tightness <- paste0("fig_", file_prg, "-eqm-wages-mkt-tightness.pdf")
ggsave(here("out", fout_eqm_w_tightness), heigh = 1.4 * myheight, width = mywidth)

stopifnot(33==12)

v_ttheta <- seq(0.01, 1.5, length.out = 100)
v_find <- compute_find(v_ttheta, m0 = match_efficiency, ggamma = ggamma)
v_fill <- compute_fill(v_ttheta, m0 = match_efficiency, ggamma = ggamma)
v_find_prime <- compute_find_prime(v_ttheta, m0 = match_efficiency, ggamma = ggamma)

plot(v_ttheta, v_find_prime)

plot(v_ttheta, v_fill)

plot(v_ttheta, v_find, type = "b")
points(v_ttheta, v_find_prime, type = "l", col = "blue")
points(v_ttheta, v_find - v_find_prime)
abline(v = target_tight)
abline(h = 0.0)

plot(v_theta, compute_fill(v_ttheta, m0 = match_efficiency, ggamma = ggamma))

test <- compute_dynamics(vprod = v_prod, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, c = c_baseline, ch = ch_baseline, cl = cl_baseline, r = r, pphi = pphi, ggamma = ggamma, ct = ct_baseline, mytol = TOL)

junk <- tribble(
  ~x, ~y,
  1, 2,
  3, 4
)

fun_junk <- function(x, y, nsim) {
  df <- tibble(junk = rnorm(nsim, mean = x, sd = y))
  return(df)
}

junk <- junk %>% 
  mutate(result_list = list(x, y)) %>% 
  mutate(result = pmap(list(x, y), fun_junk, nsim = 100)) %>% 
  unnest(result)

stopifnot(33==12)


dat <- dat %>% 
  mutate(dynamics = map_vec(v_prod, compute_eqm_tight, c = c, ch = ch, cl = cl, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct_baseline, r = r, pphi = pphi, ggamma = ggamma, mytol = TOL))

# compute_dynamics <- function(df, bbeta, m0) {
#   c <- df$c
#   ch <- df$ch
#   cl <- df$
# }

dat <- dat %>% 
  mutate(ttheta = pmap_dbl(list(), compute_eqm_tight))

ttheta_baseline <- map_dbl(v_prod, compute_eqm_tight, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct, r = r, c = c_baseline, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma, mytol = TOL)
find_baseline <- compute_find(ttheta_baseline, m0 = match_efficiency, ggamma = ggamma)

plot(v_prod, s_monthly / (s_monthly + compute_daily2monthly(find_baseline)), type = "l", xlab = "productivity", ylab = "unemployment rate", main = "Response of UR to prod shock")
abline(h = s_monthly / (s_monthly + target_find_monthly))

elasticity_baseline <- compute_elasticity()

stopifnot(33==12)

plot(ttheta_baseline, compute_find_prime(ttheta_baseline, m0 = match_efficiency, ggamma = ggamma))


cl_hi_T1 <- y - z - bbeta * s * ct - bbeta * target_find * pphi * ch / (1 - pphi) - bbeta * (r + s) * ch / (1 - pphi)
cl_hi_T2 <- (1 + r) / target_find
cl_hi_1 <- cl_hi_T1 * cl_hi_T2 - 0.01

# cl_hi_T3 <- c_baseline * pphi + bbeta * compute_find_prime(target_tight, m0 = match_efficiency, ggamma = ggamma) * pphi * ch 
# cl_hi_T4 <- bbeta * compute_find_prime(target_tight, m0 = match_efficiency, ggamma = ggamma) * (1 - pphi)
# cl_hi_2 <- cl_hi_T3 / cl_hi_T4 - 0.01

# check <- c_baseline * pphi + bbeta * compute_find_prime(target_tight, m0 = match_efficiency, ggamma = ggamma) * (pphi * ch - (1 - pphi) * cl_hi_2)

# cl_hi <- min(cl_hi_1, cl_hi_2) - 6.3
cl_hi <- cl_hi_1

ttheta_cl <- map_dbl(v_prod, compute_eqm_tight, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct, r = r, c = c_baseline, ch = ch, cl = cl_hi, pphi = pphi, ggamma = ggamma, mytol = TOL)
find_cl <- compute_find(ttheta_cl, m0 = match_efficiency, ggamma = ggamma)

check <- (c_baseline * pphi + bbeta * compute_find_prime(ttheta = ttheta_cl, m0 = match_efficiency, ggamma = ggamma) * pphi * ch) / (bbeta * compute_find_prime(ttheta = ttheta_cl, m0 = match_efficiency, ggamma = ggamma) * (1 - pphi))

lm(log(ttheta_cl) ~ log(v_prod))
lm(log(ttheta_baseline) ~ log(v_prod))

plot(v_prod, s_monthly / (s_monthly + compute_daily2monthly(find_baseline)), type = "l")
points(v_prod, s_monthly / (s_monthly + compute_daily2monthly(find_cl)), type = "l", col = "blue")
abline(h = s_monthly / (s_monthly + target_find_monthly))

stopifnot(33==12)

# High fixed cost paid by firm
ch_hi <- (1 - pphi) * (y - z - bbeta * s * ct + bbeta * compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma) * cl) / bbeta / (r + s + pphi * compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma))
ch_hi <- ch_hi - 0.01
c_hi <- compute_c(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, ch = ch_hi, ct = ct, cl = cl, pphi = pphi, ggamma = ggamma) 

cH <- ch_hi
ppsi <- 0.5
cl_cl <- ppsi * cH
ch_cl <- (1 - ppsi) * cH
stopifnot(pphi * ch_cl - (1 - pphi) * cl_cl >= 0)
c_cl <- compute_c(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, r = r, ch = ch, ct = ct, cl = cl, pphi = pphi, ggamma = ggamma) 

dat <- tibble(y = v_prod, 
              ttheta = map_dbl(v_prod, compute_eqm_tight, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct, r = r, c = c, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma, mytol = TOL),
              ttheta_hi = map_dbl(v_prod, compute_eqm_tight, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct, r = r, c = c_hi, ch = ch_hi, cl = cl, pphi = pphi, ggamma = ggamma, mytol = TOL),
              ttheta_cl = map_dbl(v_prod, compute_eqm_tight, m0 = match_efficiency, z = z, bbeta = bbeta, s = s, ct = ct, r = r, c = c_cl, ch = ch_cl, cl = cl_cl, pphi = pphi, ggamma = ggamma, mytol = TOL),
              find = compute_find(ttheta, m0 = match_efficiency, ggamma = ggamma),
              find_hi = compute_find(ttheta_hi, m0 = match_efficiency, ggamma = ggamma),
              find_cl = compute_find(ttheta_cl, m0 = match_efficiency, ggamma = ggamma),
              ur_monthly = s_monthly / (s_monthly + compute_daily2monthly(find)),
              ur_monthly_hi = s_monthly / (s_monthly + compute_daily2monthly(find_hi)),
              ur_monthly_cl = s_monthly / (s_monthly + compute_daily2monthly(find_cl)))

ggplot(data = dat) +
  geom_line(mapping = aes(x = y, y = ur_monthly), color = "black") +
  geom_line(mapping = aes(x = y, y = ur_monthly_hi), color = "blue") +
  geom_line(mapping = aes(x = y, y = ur_monthly_cl), color = "maroon") +
  geom_hline(yintercept = s_monthly / (s_monthly + target_find_monthly), color = "gray")

ggplot(data = dat) +
  geom_line(mapping = aes(x = y, y = log(ttheta))) +
  geom_line(mapping = aes(x = y, y = log(ttheta_hi)))

lm(log(ttheta) ~ log(y), data = dat)
lm(log(ttheta_hi) ~ log(y), data = dat)

stopifnot(33==12)

check <- y - z - bbeta * s * ct - bbeta * (r + s) * ch_hi / (1 - pphi) - bbeta * compute_find(ttheta = target_tight, m0 = match_efficiency, ggamma = ggamma) * (pphi * ch_hi / (1 - pphi) - cl)



v_find <- compute_find(v_ttheta, m0 = match_efficiency, ggamma = ggamma)
v_find_hi <- compute_find(v_ttheta_hi, m0 = match_efficiency, ggamma = ggamma)
plot(v_prod, s_monthly / (s_monthly + compute_daily2month(v_find)), type = "l")
points(v_prod, s_monthly / (s_monthly + compute_daily2month(v_find_hi)), type = "l", col = "red")
abline(h = s_monthly / (s_monthly + target_find_monthly), lty = 3)
abline(v = y, lty = 3)

lm(log(v_ttheta) ~ log(v_prod))
# plot(log(v_prod), log(v_ttheta))

stopifnot(33==12)

for (i in seq_along(v_ttheta)) {
  print(i)
  i_y <- v_prod[[i]]
  i_ttheta_lo <- 1e-5
  i_ttheta_hi <- get_ttheta_hi(pphi = pphi, y = i_y, z = z, bbeta = bbeta, r = r_monthly, s = s_monthly, c = c, ch = ch, ct = ct)
  stopifnot(i_ttheta_lo < i_ttheta_hi)
  i_sol <- uniroot(compute_eqm_tight,
                   y = i_y, 
                   m0 = match_efficiency, 
                   z = z, bbeta = bbeta, s = s, ct = ct, r = r, 
                   c = c, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma,
                   lower = i_ttheta_lo,
                   upper = i_ttheta_hi)
  i_ttheta <- i_sol$root
  v_ttheta[i] <- i_ttheta
}

stopifnot(33==12)



# Change fixed cost of hiring ---------------------------------------------
# c_pissarides <- 0.1
# ch <- (1 - pphi) * c_pissarides / (bbeta * (r + s))
ch <- 0.1
c_hi <- compute_c(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = 0, r = r, ch = ch, cl = 0.0, pphi = pphi, ggamma = ggamma)

v_ttheta_ch <- vector(mode = "double", length = npts)
for (i in seq_along(v_ttheta)) {
  print(i)
  i_y <- v_prod[[i]]
  i_ttheta_lo <- 1e-5
  i_ttheta_hi <- get_ttheta_hi(pphi = pphi, y = i_y, z = z, bbeta = bbeta, r = r, s = s, c = c_hi, ch = ch, ct = ct)
  stopifnot(i_ttheta_lo < i_ttheta_hi)
  i_sol <- uniroot(compute_eqm_tight,
                   y = i_y, 
                   m0 = match_efficiency, 
                   z = z, bbeta = bbeta, s = s, ct = ct, r = r, 
                   c = c_hi, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma,
                   lower = i_ttheta_lo,
                   upper = i_ttheta_hi)
  i_ttheta <- i_sol$root
  v_ttheta_ch[i] <- i_ttheta
}

v_find_ch <- compute_find(v_ttheta_ch, m0 = match_efficiency, ggamma = ggamma)
plot(v_prod, s / (s + v_find), type = "l", col = "black")
points(v_prod, s / (s + v_find_ch), type = "l", col = "red")
abline(h = s / (s + target_find), lty = 3)

lm(log(v_ttheta_ch) ~ log(v_prod))

ch <- 7.0
c_hi <- compute_c(ttheta = target_tight, m0 = match_efficiency, y = y, z = z, bbeta = bbeta, s = s, ct = 0, r = r, ch = ch, cl = 0.0, pphi = pphi, ggamma = ggamma)

v_ttheta_chhi <- vector(mode = "double", length = npts)
for (i in seq_along(v_ttheta)) {
  print(i)
  i_y <- v_prod[[i]]
  i_ttheta_lo <- 1e-5
  i_ttheta_hi <- get_ttheta_hi(pphi = pphi, y = i_y, z = z, bbeta = bbeta, r = r, s = s, c = c_hi, ch = ch, ct = ct)
  # stopifnot(i_ttheta_lo < i_ttheta_hi)
  i_sol <- uniroot(compute_eqm_tight,
                   y = i_y, 
                   m0 = match_efficiency, 
                   z = z, bbeta = bbeta, s = s, ct = ct, r = r, 
                   c = c_hi, ch = ch, cl = cl, pphi = pphi, ggamma = ggamma,
                   lower = i_ttheta_lo,
                   upper = i_ttheta_hi)
  i_ttheta <- i_sol$root
  v_ttheta_chhi[i] <- i_ttheta
}

v_find_chhi <- compute_find(v_ttheta_chhi, m0 = match_efficiency, ggamma = ggamma)
plot(v_prod, s_monthly / (s_monthly + compute_daily2month(v_find)), type = "l", col = "black")
points(v_prod, s_monthly / (s_monthly + compute_daily2month(v_find_ch)), type = "l", col = "red")
points(v_prod, s_monthly / (s_monthly + compute_daily2month(v_find_chhi)), type = "l", col = "blue")
abline(h = s_monthly / (s_monthly + compute_daily2month(target_find)), lty = 3)

lm(log(v_ttheta_chhi) ~ log(v_prod))
