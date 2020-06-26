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

