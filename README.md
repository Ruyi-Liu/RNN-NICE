# RNN-NICE
Recurrent Neural Network for Non-Iterative Conditional Expectation for Time-Varying Treatment Effect

<img src="Stabilized Weighted RNN.png" alt="Model Architecture"  width="800" height="250"/>


## Overview

This repository provides implementations of four models for time-varying treatment effect estimation. Each model is organized in its own directory, and all necessary dependencies are listed within each Jupyter notebook.

Included contents:
- Code for model development
- Data loading
- Reproducible steps for large-scale data
- Final causal effect estimation and performance metrics calculation

## Directory Structure

- **Example Data Files**: Example datasets for demonstration and reproduction.
- **RNN_original_model**: Implementation of the original RNN model.
- **Sequential model**: Sequential modeling approach for time-varying covariates and outcomes.
- **Stabilized Weighted RNN**: Model integrating stabilized weights to improve robustness against positivity violations.
- **gfoRmula**: Baseline method using the R package `gfoRmula` for comparison.

## Notes

- Large-scale data files are not included in this repository due to size limitations.
- The baseline `gfoRmula` method is included as a reference for comparison.




