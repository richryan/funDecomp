# funDecomp #

Code for the paper 
*Unemployment Volatility: When Workers Pay Costs upon Accepting Jobs*.
The paper is included in the top of the directory in 
`unemployment_volatility_when_workers_pay_costs.pdf`.

Here's the abstract:

> When a firm hires a worker,
>   adding the new hire to payroll is costly.
>   These costs
>   reduce the amount of resources that can go to recruiting workers
>   and % This % affects
>   amplify how unemployment responds to changes in productivity.
>   Workers also incur up-front costs upon accepting jobs.
>   Examples include moving expenses and regulatory fees.
>   I establish that workers' costs
>   lessen the response of unemployment to productivity changes and
>   do not subtract from resources available for recruitment.
>   The influence of workers' costs is bounded by properties of a matching function,
>   which describes how
>   job openings and unemployment produce hires.
>   Using data on job finding that are adjusted
>   for workers' transitions between employment and unemployment and 
>   for how the Job Openings and Labor Turnover Survey records hires,
>   I estimate a bound that ascribes limited influence to workers' costs. 
>   The results demonstrate that costs paid by workers upon accepting jobs
>   affect outcomes in the labor market
>   (firms threaten workers with paying the up-front costs again if wage
>   negotiations fail),
>   but their influence on volatility is less important than firms' costs.

Code that replicates that project is contained in `cde/`.
The scripts are

  * `01-get-FRED-data.R`: Retrieves data from [FRED](https://fred.stlouisfed.org/).
  * `02-adjust-transition-rates.R`: Uses the FRED data to adjust rates of 
  job finding and job separation to account for 
      * workers' transitions between employment and unemployment and
      * how the Job Openings and Labor Turnover Survey (JOLTS) records hires---all
	  hires for a month are recorded, 
	  even those who separate before the month ends.
  * `03-est-matching.R`: Estimates the nonlinear mathching technology,
	$M \left( u, v \right) = \mu \frac{uv}{\left( u^{\gamma} + v^{\gamma} \right)^{1/\gamma}}$.   Various plots are produced showing how 
  the Beveridge curve shifted since JOLTS data on vacancies 
  became available and the fit produced by the estimated model.
  * `04-plot-elasticity.R`: Plots the elasticity of matching
  with respect to unemployment, $\eta_{M,u}$, and
  the upper bound for $\Upsilon$,
  $\max \left( \frac{1}{\eta_{M,u}}, \frac{1}{1-\eta_{M,u}} \right)$.
  The factor $\Upsilon$ is one of the factors 
  in [Ljungqvist and Sargent's (2017)](https://www.aeaweb.org/articles?id=10.1257/aer.20150233) fundamental decomposition of the 
  elasticity of market tightness with respect to productivity, $\eta_{\theta,y}$.

