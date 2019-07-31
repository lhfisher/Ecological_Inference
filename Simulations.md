Simulations
================
Leigh Fisher
July 31 2019

  - [Introduction and notes](#introduction-and-notes)
  - [Simulations to assess simplifying assumptions in the absence of
    vaccination](#simulations-to-assess-simplifying-assumptions-in-the-absence-of-vaccination)
  - [Simulations in partially vaccinated
    populations](#simulations-in-partially-vaccinated-populations)
  - [Simulations to assess asymptotic behavior of ecological vaccine
    model](#simulations-to-assess-asymptotic-behavior-of-ecological-vaccine-model)

# Introduction and notes

  - This file contains code and analysis for the simulations presented
    in \`\`Ecological inference for infectious disease data, with
    application to vaccination strategies".
  - All functions are found in
ecofall\_sim\_functions.R

# Simulations to assess simplifying assumptions in the absence of vaccination

For all simulations in this section, we summarize results over 250
simulated populations, each with a population of 100,000 and we assume
the endemic component, \(\alpha_{EN}\), is -10.

## \(R_0\) = 2.5

``` r
alpha_ar = log(2.5)
# 
set.seed(47)
sim_r0big=replicate(Nsims, sim.simple(N=N, Nwks=3*52, alpha_ar=alpha_ar, alpha_en=alpha_en))

vals=c(alpha_ar, alpha_en)

res_r0big = plyr::ldply(1:8, function(i){
        m = paste0('mod', i)
        
        est=sapply(sim_r0big[m, ], function(x){x$par})
        ests = rowMeans(est)
        se.ests = apply(est, 1, sd)
        percentiles = apply(est, 1, function(x){quantile(x, c(0.025, 0.975))})
        
        bias = rowMeans(sapply(sim_r0big[m, ], function(x){(x$par - vals)}))
        mse = rowMeans(sapply(sim_r0big[m, ], function(x){(x$par - vals)^2}))
        c(mod = i, Est.ar = ests[1], pct.025.ar = as.numeric(percentiles[1, 1]), pct.975.ar = as.numeric(percentiles[2, 1]), Bias.ar = bias[1], MSE.ar = mse[1], 
          Est.en = ests[2], pct.025.en = as.numeric(percentiles[1, 2]), pct.975.en = as.numeric(percentiles[2, 2]), Bias.en = bias[2], MSE.en = mse[2])

})
```

Examples of simulated epidemic curves when \(R_0=2.5\)
<img src="Simulations_files/figure-gfm/eg_r0_big-1.png" width="55%" style="display: block; margin: auto;" />

<img src="Simulations_files/figure-gfm/bias_mse_r0_big-1.png" title="Red lines indicate true parameter value; blue line is at zero." alt="Red lines indicate true parameter value; blue line is at zero." width="33%" /><img src="Simulations_files/figure-gfm/bias_mse_r0_big-2.png" title="Red lines indicate true parameter value; blue line is at zero." alt="Red lines indicate true parameter value; blue line is at zero." width="33%" /><img src="Simulations_files/figure-gfm/bias_mse_r0_big-3.png" title="Red lines indicate true parameter value; blue line is at zero." alt="Red lines indicate true parameter value; blue line is at zero." width="33%" />

| mod |  Est.ar | pct.025.ar | pct.975.ar | Bias.ar | MSE.ar |   Est.en | pct.025.en | pct.975.en | Bias.en | MSE.en |
| --: | ------: | ---------: | ---------: | ------: | -----: | -------: | ---------: | ---------: | ------: | -----: |
|   1 |   0.916 |      0.910 |      0.923 |   0.000 |  0.000 | \-10.003 |   \-10.244 |    \-9.790 | \-0.003 |  0.014 |
|   2 |   0.073 |      0.072 |      0.074 | \-0.843 |  0.711 | \-12.606 |   \-12.839 |   \-12.416 | \-2.606 |  6.802 |
|   3 |   0.726 |      0.722 |      0.730 | \-0.191 |  0.036 |  \-9.874 |   \-10.115 |    \-9.660 |   0.126 |  0.031 |
|   4 | \-0.020 |    \-0.021 |    \-0.019 | \-0.936 |  0.877 | \-12.573 |   \-12.806 |   \-12.374 | \-2.573 |  6.632 |
|   5 |   0.916 |      0.909 |      0.922 |   0.000 |  0.000 | \-10.003 |   \-10.244 |    \-9.785 | \-0.003 |  0.014 |
|   6 |   0.089 |      0.087 |      0.091 | \-0.827 |  0.685 | \-12.612 |   \-12.845 |   \-12.417 | \-2.612 |  6.834 |
|   7 |   0.743 |      0.739 |      0.747 | \-0.173 |  0.030 |  \-9.890 |   \-10.133 |    \-9.678 |   0.110 |  0.027 |
|   8 | \-0.002 |    \-0.003 |    \-0.001 | \-0.919 |  0.844 | \-12.576 |   \-12.805 |   \-12.346 | \-2.576 |  6.649 |

Summary of parameter estimates from 250 simulations, where alpha\_ar =
0.916 and alpha\_en = -10

## \(R_0 = 1\)

``` r
alpha_ar = log(1)

set.seed(47)
sim_r0mod=replicate(Nsims, sim.simple(N=N, Nwks=3*52, alpha_ar=alpha_ar, alpha_en=alpha_en))

vals=c(alpha_ar, alpha_en)

res_r0mod = plyr::ldply(1:8, function(i){
        m = paste0('mod', i)
        
        est=sapply(sim_r0mod[m, ], function(x){x$par})
        ests = rowMeans(est)
        se.ests = apply(est, 1, sd)
        percentiles = apply(est, 1, function(x){quantile(x, c(0.025, 0.975))})
        
        bias = rowMeans(sapply(sim_r0mod[m, ], function(x){(x$par - vals)}))
        mse = rowMeans(sapply(sim_r0mod[m, ], function(x){(x$par - vals)^2}))
        c(mod = i, Est.ar = ests[1], pct.025.ar = as.numeric(percentiles[1, 1]), pct.975.ar = as.numeric(percentiles[2, 1]), Bias.ar = bias[1], MSE.ar = mse[1], 
          Est.en = ests[2], pct.025.en = as.numeric(percentiles[1, 2]), pct.975.en = as.numeric(percentiles[2, 2]), Bias.en = bias[2], MSE.en = mse[2])

})
```

Examples of simulated epidemic curves when \(R_0=1\)
<img src="Simulations_files/figure-gfm/eg_r0_mod-1.png" width="55%" style="display: block; margin: auto;" />

<img src="Simulations_files/figure-gfm/bias_mse_r0_mod-1.png" title="Red lines indicate true parameter value; blue line is at zero." alt="Red lines indicate true parameter value; blue line is at zero." width="33%" /><img src="Simulations_files/figure-gfm/bias_mse_r0_mod-2.png" title="Red lines indicate true parameter value; blue line is at zero." alt="Red lines indicate true parameter value; blue line is at zero." width="33%" /><img src="Simulations_files/figure-gfm/bias_mse_r0_mod-3.png" title="Red lines indicate true parameter value; blue line is at zero." alt="Red lines indicate true parameter value; blue line is at zero." width="33%" />

| mod |   Est.ar | pct.025.ar | pct.975.ar |  Bias.ar | MSE.ar |    Est.en | pct.025.en | pct.975.en |  Bias.en | MSE.en |
| --: | -------: | ---------: | ---------: | -------: | -----: | --------: | ---------: | ---------: | -------: | -----: |
|   1 | \-0.0097 |   \-0.0758 |     0.0259 | \-0.0097 | 0.0007 |  \-9.9181 |  \-10.4353 |   \-9.3373 |   0.0819 | 0.0851 |
|   2 | \-0.0416 |   \-0.1019 |   \-0.0150 | \-0.0416 | 0.0022 | \-10.5054 |  \-11.3577 |   \-9.6039 | \-0.5054 | 0.4303 |
|   3 | \-0.0104 |   \-0.0784 |     0.0258 | \-0.0104 | 0.0007 |  \-9.9150 |  \-10.4178 |   \-9.3384 |   0.0850 | 0.0851 |
|   4 | \-0.0416 |   \-0.0956 |   \-0.0162 | \-0.0416 | 0.0021 | \-10.5051 |  \-11.3486 |   \-9.6409 | \-0.5051 | 0.4158 |
|   5 | \-0.0097 |   \-0.0758 |     0.0262 | \-0.0097 | 0.0007 |  \-9.9181 |  \-10.4353 |   \-9.3373 |   0.0819 | 0.0851 |
|   6 | \-0.0416 |   \-0.1019 |   \-0.0150 | \-0.0416 | 0.0022 | \-10.5048 |  \-11.3577 |   \-9.6039 | \-0.5048 | 0.4296 |
|   7 | \-0.0104 |   \-0.0784 |     0.0258 | \-0.0104 | 0.0007 |  \-9.9149 |  \-10.4178 |   \-9.3384 |   0.0851 | 0.0849 |
|   8 | \-0.0416 |   \-0.0956 |   \-0.0162 | \-0.0416 | 0.0021 | \-10.5050 |  \-11.3486 |   \-9.6409 | \-0.5050 | 0.4155 |

Summary of parameter estimates from 250 simulations, where alpha\_ar =
0.000 and alpha\_en = -10

## \(R_0 = 0.85\)

``` r
alpha_ar = log(0.85)

set.seed(47)
sim_r0lo=replicate(Nsims, sim.simple(N=N, Nwks=3*52, alpha_ar=alpha_ar, alpha_en=alpha_en))

vals=c(alpha_ar, alpha_en)

res_r0lo = plyr::ldply(1:8, function(i){
        m = paste0('mod', i)
        
        est=sapply(sim_r0lo[m, ], function(x){x$par})
        ests = rowMeans(est)
        se.ests = apply(est, 1, sd)
        percentiles = apply(est, 1, function(x){quantile(x, c(0.025, 0.975))})
        
        bias = rowMeans(sapply(sim_r0lo[m, ], function(x){(x$par - vals)}))
        mse = rowMeans(sapply(sim_r0lo[m, ], function(x){(x$par - vals)^2}))
        c(mod = i, Est.ar = ests[1], pct.025.ar = as.numeric(percentiles[1, 1]), pct.975.ar = as.numeric(percentiles[2, 1]), Bias.ar = bias[1], MSE.ar = mse[1], 
          Est.en = ests[2], pct.025.en = as.numeric(percentiles[1, 2]), pct.975.en = as.numeric(percentiles[2, 2]), Bias.en = bias[2], MSE.en = mse[2])

})
```

Examples of simulated epidemic curves when \(R_0=0.85\)
<img src="Simulations_files/figure-gfm/eg_r0_lo-1.png" width="55%" style="display: block; margin: auto;" />

<img src="Simulations_files/figure-gfm/bias_mse_r0_lo-1.png" title="Red lines indicate true parameter value; blue line is at zero." alt="Red lines indicate true parameter value; blue line is at zero." width="33%" /><img src="Simulations_files/figure-gfm/bias_mse_r0_lo-2.png" title="Red lines indicate true parameter value; blue line is at zero." alt="Red lines indicate true parameter value; blue line is at zero." width="33%" /><img src="Simulations_files/figure-gfm/bias_mse_r0_lo-3.png" title="Red lines indicate true parameter value; blue line is at zero." alt="Red lines indicate true parameter value; blue line is at zero." width="33%" />

| mod |   Est.ar | pct.025.ar | pct.975.ar |  Bias.ar | MSE.ar |   Est.en | pct.025.en | pct.975.en | Bias.en | MSE.en |
| --: | -------: | ---------: | ---------: | -------: | -----: | -------: | ---------: | ---------: | ------: | -----: |
|   1 | \-0.3137 |   \-0.4674 |   \-0.1993 | \-0.0261 | 0.0065 | \-9.9595 |  \-10.3147 |   \-9.6225 |  0.0405 | 0.0372 |
|   2 | \-0.3274 |   \-0.4820 |   \-0.2076 | \-0.0397 | 0.0075 | \-9.9741 |  \-10.3469 |   \-9.6231 |  0.0259 | 0.0368 |
|   3 | \-0.3139 |   \-0.4674 |   \-0.1993 | \-0.0262 | 0.0065 | \-9.9592 |  \-10.3147 |   \-9.6232 |  0.0408 | 0.0371 |
|   4 | \-0.3275 |   \-0.4843 |   \-0.2071 | \-0.0399 | 0.0075 | \-9.9741 |  \-10.3493 |   \-9.6229 |  0.0259 | 0.0369 |
|   5 | \-0.3138 |   \-0.4674 |   \-0.1993 | \-0.0261 | 0.0065 | \-9.9594 |  \-10.3147 |   \-9.6225 |  0.0406 | 0.0372 |
|   6 | \-0.3274 |   \-0.4820 |   \-0.2076 | \-0.0397 | 0.0075 | \-9.9741 |  \-10.3469 |   \-9.6231 |  0.0259 | 0.0368 |
|   7 | \-0.3139 |   \-0.4674 |   \-0.1993 | \-0.0262 | 0.0065 | \-9.9592 |  \-10.3147 |   \-9.6232 |  0.0408 | 0.0371 |
|   8 | \-0.3276 |   \-0.4843 |   \-0.2071 | \-0.0399 | 0.0075 | \-9.9740 |  \-10.3493 |   \-9.6229 |  0.0260 | 0.0370 |

Summary of parameter estimates from 250 simulations, where alpha\_ar =
-0.288 and alpha\_en =
-10

## Summary

<img src="Simulations_files/figure-gfm/sum_sim_novacc-1.png" width="50%" /><img src="Simulations_files/figure-gfm/sum_sim_novacc-2.png" width="50%" /><img src="Simulations_files/figure-gfm/sum_sim_novacc-3.png" width="50%" /><img src="Simulations_files/figure-gfm/sum_sim_novacc-4.png" width="50%" />

# Simulations in partially vaccinated populations

``` r
source('ecofall_sim_functions.R')

N=rep(1e5, 5)
x=c(0.65, 0.7, 0.75, 0.8, 0.85) 
phi0=0
phi=0.8
alpha_ar = log(2)
alpha_en =   1
exp(alpha_ar)*(1-phi*x)
Nwks=52*3 # Number of weeks
# Nsims=1e2
sts=c(0.5, 0, 1)
```

Below, we plot examples of a population consisting of 5
areas

<img src="Simulations_files/figure-gfm/sim_partial_vacc-1.png" width="33%" /><img src="Simulations_files/figure-gfm/sim_partial_vacc-2.png" width="33%" /><img src="Simulations_files/figure-gfm/sim_partial_vacc-3.png" width="33%" />

## All-or-none vaccine

| Model            | Est AR | 2.5% AR | 97.5% AR | Bias AR | MSE AR | Est EN | 2.5% EN | 97.5% EN | Bias EN | MSE EN | Est Phi | 2.5% Phi | 97.5% Phi | Bias Phi | MSE Phi |
| :--------------- | -----: | ------: | -------: | ------: | -----: | -----: | ------: | -------: | ------: | -----: | ------: | -------: | --------: | -------: | ------: |
| All-or-none      |  0.680 |   0.588 |    0.742 | \-0.013 |  0.002 |  1.009 |   0.769 |    1.216 |   0.009 |  0.013 |   0.798 |    0.775 |     0.822 |  \-0.002 |   0.000 |
| Leaky            |  0.680 |   0.588 |    0.741 | \-0.013 |  0.002 |  1.010 |   0.771 |    1.218 |   0.010 |  0.013 |   0.799 |    0.776 |     0.822 |  \-0.001 |   0.000 |
| Ecological       |  0.661 |   0.353 |    0.894 | \-0.032 |  0.019 |  1.004 |   0.566 |    1.428 |   0.004 |  0.051 |   0.787 |    0.571 |     0.923 |  \-0.013 |   0.009 |
| Epidemic-Endemic |  0.513 |   0.165 |    0.861 | \-0.181 |  0.062 |  0.056 | \-0.174 |    0.266 | \-0.944 |  0.905 |   0.631 |    0.564 |     0.701 |  \-0.169 |   0.030 |

Summary of parameter estimates from 250 simulations, where alpha\_ar =
0.693 alpha\_en = 1 and an all-or-none vaccine effect of
0.8

<img src="Simulations_files/figure-gfm/summary_aon_sims-1.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_aon_sims-2.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_aon_sims-3.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_aon_sims-4.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_aon_sims-5.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_aon_sims-6.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_aon_sims-7.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_aon_sims-8.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_aon_sims-9.png" width="33%" />

## Leaky vaccine

| Model            | Est AR | 2.5% AR | 97.5% AR | Bias AR | MSE AR | Est EN | 2.5% EN | 97.5% EN | Bias EN | MSE EN | Est Phi | 2.5% Phi | 97.5% Phi | Bias Phi | MSE Phi |
| :--------------- | -----: | ------: | -------: | ------: | -----: | -----: | ------: | -------: | ------: | -----: | ------: | -------: | --------: | -------: | ------: |
| All-or-none      |  0.681 |   0.582 |    0.744 | \-0.012 |  0.002 |  0.999 |   0.786 |    1.221 | \-0.001 |  0.013 |   0.798 |    0.772 |     0.819 |  \-0.002 |   0.000 |
| Leaky            |  0.681 |   0.582 |    0.743 | \-0.013 |  0.002 |  1.001 |   0.789 |    1.223 |   0.001 |  0.013 |   0.799 |    0.773 |     0.820 |  \-0.001 |   0.000 |
| Ecological       |  0.657 |   0.365 |    0.912 | \-0.036 |  0.021 |  0.990 |   0.555 |    1.441 | \-0.010 |  0.050 |   0.784 |    0.600 |     0.951 |  \-0.016 |   0.009 |
| Epidemic-Endemic |  0.508 |   0.191 |    0.872 | \-0.186 |  0.070 |  0.045 | \-0.158 |    0.284 | \-0.955 |  0.924 |   0.630 |    0.570 |     0.704 |  \-0.170 |   0.030 |

Summary of parameter estimates from 250 simulations, where alpha\_ar =
0.693 alpha\_en = 1 and a leaky vaccine effect of
0.8

<img src="Simulations_files/figure-gfm/summary_lk_sims-1.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_lk_sims-2.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_lk_sims-3.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_lk_sims-4.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_lk_sims-5.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_lk_sims-6.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_lk_sims-7.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_lk_sims-8.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_lk_sims-9.png" width="33%" />

# Simulations to assess asymptotic behavior of ecological vaccine model

``` r
Nsims=100
pmt=proc.time()
set.seed(47)
asym_sims=replicate(Nsims, sim.one(vacc='AoN', N=N, Nwks=52*20, alpha_ar=alpha_ar, alpha_en=alpha_en, phi=phi, x=x))
proc.time()-pmt
```

user system elapsed 2720.43 32.19 2848.21

``` r
# 
```

Example of epidemic curve
![](Simulations_files/figure-gfm/asymptotic_eg-1.png)<!-- -->

Simulation
results

| Model            | Est AR | 2.5% AR | 97.5% AR | Bias AR | MSE AR | Est EN | 2.5% EN | 97.5% EN | Bias EN | MSE EN | Est Phi | 2.5% Phi | 97.5% Phi | Bias Phi | MSE Phi |
| :--------------- | -----: | ------: | -------: | ------: | -----: | -----: | ------: | -------: | ------: | -----: | ------: | -------: | --------: | -------: | ------: |
| All-or-none      |  0.692 |   0.653 |    0.717 | \-0.001 |  0.000 |  1.007 |   0.922 |    1.088 |   0.007 |  0.002 |   0.800 |    0.791 |     0.810 |    0.000 |   0.000 |
| Leaky            |  0.692 |   0.653 |    0.716 | \-0.001 |  0.000 |  1.011 |   0.927 |    1.093 |   0.011 |  0.002 |   0.804 |    0.794 |     0.814 |    0.004 |   0.000 |
| Ecological       |  0.614 |   0.523 |    0.718 | \-0.079 |  0.009 |  0.927 |   0.769 |    1.089 | \-0.073 |  0.012 |   0.767 |    0.698 |     0.835 |  \-0.033 |   0.002 |
| Epidemic-Endemic |  0.449 |   0.337 |    0.597 | \-0.244 |  0.064 |  0.041 | \-0.033 |    0.123 | \-0.959 |  0.922 |   0.621 |    0.598 |     0.650 |  \-0.179 |   0.032 |

Summary of parameter estimates from 100 simulations, where alpha\_ar =
0.693 alpha\_en = 1 and a leaky vaccine effect of
0.8

<img src="Simulations_files/figure-gfm/summary_asymp_sims-1.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_asymp_sims-2.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_asymp_sims-3.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_asymp_sims-4.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_asymp_sims-5.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_asymp_sims-6.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_asymp_sims-7.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_asymp_sims-8.png" width="33%" /><img src="Simulations_files/figure-gfm/summary_asymp_sims-9.png" width="33%" />
