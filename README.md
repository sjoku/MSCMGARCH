# MSCMGARCH

## Disclaimer

This repository is a fork of the original repo ``MJFullness/MSCMGARCH``. The intention of this repo is to expand on the documentation and working reproducible examples, aligned to the original paper. 

## Introduction

The MSCMGARCH package provides a Markov-switching multivariate GARCH model with copula-distributed innovations. It allows for the dynamic assessment of financial market interdependencies, capturing potential regime shifts in asymmetric tail dependence structures among volatility innovations in speculative return series.

## Installation

To install the MSCMGARCH package, you can use the following commands in R:

```R
# Install the devtools package if you don't have it already
install.packages("devtools")

# Install the MSCMGARCH package from GitHub
devtools::install_github("nutle/MSCMGARCH")
```

## Quickstart Guide

Here is a quickstart guide with examples of the main functions in the MSCMGARCH package.

### Simulate Data

You can simulate data from the Markov-switching multivariate GARCH model using the `Sim_MSCMGARCH` function:

```R
# Simulate data
sim_data <- Sim_MSCMGARCH(type = "Normal Gumbel", n = 1000, amount = 10, true_par = c(0.9, 0.2, 12), lower_bound = c(0, 0, 1), upper_bound = c(1, 1, 17), nc = 1, seed = 123)
```

### Estimate Model Parameters

You can estimate the model parameters using the `MLE_MSCMGARCH` function:

```R
# Estimate model parameters
estimation <- MLE_MSCMGARCH(r = sim_data, type = "Normal Gumbel", nc = 1)
```

### Perform Rolling Window Analysis

You can perform rolling window analysis for forecasting using the `rolling_window` function:

```R
# Perform rolling window analysis
forecast <- rolling_window(series = sim_data, type = "MS_CMGARCH", copula_type = c(4, 4, 4), asymmetric = FALSE, window_length = 250, portfolio_weights = c(0.5, 0.5), signs = NULL, nc = 1)
```

## Contributing

If you would like to contribute to the MSCMGARCH package, please follow these steps:

1. Fork the repository on GitHub.
2. Create a new branch with a descriptive name.
3. Make your changes and commit them with clear and concise messages.
4. Push your changes to your forked repository.
5. Create a pull request to the main repository.

## Vignettes

You can find detailed vignettes for the MSCMGARCH package at the following links:

- [Quickstart Guide](vignettes/quickstart.html)
- [Experiment Reproduction](vignettes/experiment.html)
