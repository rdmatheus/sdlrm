# `sdlrm`: Sdl Regression for Integer-Valued and Paired Count Data

[![Build Status](https://travis-ci.com/travis-ci/travis-web.svg?branch=master)](https://travis-ci.com/travis-ci/travis-web)

> Rodrigo M. R. Medeiros
> <rodrigo.matheus@live.com>, IME-USP

Implementation of the SDL regression model proposed by Medeiros and Bourguignon (2020). Provide a set of functions for a complete analysis of integer-valued data, in which it is assumed that the dependent variable follows a Skew Discrete Laplace (SDL) distribution. This regression model is also useful for the analysis of experimental studies in which a paired count data is observed.

## Installation

You can install the current development version of `sdlrm` from [GitHub](https://github.com/rdmatheus/sdlrm) with:

``` r
devtools::install_github("rdmatheus/sdlrm")
```
To run the above command, it is necessary that the `devtools` package is previously installed on R. If not, install it using the following command:

``` r
install.packages("devtools")
```
After installing the devtools package, if you are using Windows, install the most current [RTools](https://cran.r-project.org/bin/windows/Rtools/) program. Finally, run the command `devtools::install_github("rdmatheus/sdlrm")`, and then the package will be installed on your computer.

## Example

This package provide complete estimation and inference for the parameters as well as simulation envelope plots, useful for assessing the goodness-of-fit of the model. The implementation is straightforward and similar to other popular packages, like `betareg` and `glm`, where the main function is `sdlrm()`. Below is an example of some functions usage and available methods.

``` r
library(sdlrm)

# Data visualization (Description: ?pss)
```

![alt text](figures/fig1a.png)
![alt text](figures/fig1b.png)

``` r
# Fit a double model (mean and dispersion) with a test for constant dispersion

fit <- sdlrm(Diff ~ Group | Group, data = pss, disp_test = TRUE)

# Print

fit
#> 
#> Call:
#> sdlrm(formula = Diff ~ Group | Group, data = pss, disp_test = TRUE)
#> 
#> mean Coefficients:
#> [1]   7.363636 -11.296970
#> 
#> dispersion Coefficients:
#> [1]  2.6555938 -0.5875912
#> # Summary

summary(fit)
#>
#> Call:
#> sdlrm(formula = Diff ~ Group | Group, data = pss, disp_test = TRUE)
#> 
#> Summary for residuals:
#>       Mean       Sd  Skewness Kurtosis
#>   0.012713 0.966192 -0.177681 1.800122
#> 
#> ----------------------------------------------------------------
#> Mean:
#> Coefficients:
#>             Estimate Std. Error t value  Pr(>|t|)   
#> (Intercept)   7.3636     3.6010  2.0449  0.040867 * 
#> GroupSport  -11.2970     4.0119 -2.8159  0.004865 **
#> ---
#> Signif. codes:  
#>   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ----------------------------------------------------------------
#> Dispersion:
#> 
#> Link function: log 
#> Coefficients:
#>             Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)  2.65559    0.32089  8.2758    <2e-16 ***
#> GroupSport  -0.58759    0.43006 -1.3663    0.1718    
#> ---
#> Signif. codes:  
#>   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ----------------------------------------------------------------
#> 
#> Test for constant dispersion:
#>                 S         W        LR         G
#> Value   1.9292900 1.8667800 1.8920900 1.8964500
#> P value 0.1648361 0.1718443 0.1689664 0.1684762
#>                 
#> In addition, Log-lik value: -88.56555 
#> AIC: 185.1311 and BIC: 190.1635

# Plot

plot(fit)
```
![alt text](figures/fig2a.png)

``` r
# Envelope

envel_sdl(fit)
```

![alt text](figures/fig3a.png)
![alt text](figures/fig3b.png)
