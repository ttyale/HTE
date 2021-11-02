## Title: HTE
## Author: Guangyu Tong
## Description: This folder includes the simulation code and examples for the paper entitled "Accounting for unequal cluster sizes in designing cluster randomized trials to detect treatment effect heterogeneity". The folder contains three files and the details of each are explained below.

## Files:
  1. Simulation_Functions_HTE.R: This file includes functions estimating sample sizes and power for HTE in cluster randomized trials without considering variable cluster sizes [sample_size_cal(), power_cal()]; considering variable cluster sizes [sample_size_cal_new(), power_cal_new()]; generating simulation data [datagen()]; performing simulations[sim()]; summarizing sample size requirement [sample_size()]. 
  2. Simulation_Functions_Mar.R: This file includes functions estimating sample sizes and power for conditional average treatment effect induced from HTE in cluster randomized trials without considering variable cluster sizes [sample_size_cal_mar(), power_cal_mar()]; considering variable cluster sizes [sample_size_cal_new_mar(), power_cal_new_mar()]; estimating sample sizes and power for marginal average treatment effect induced from HTE in cluster randomized trials without considering variable cluster sizes [sample_size_cal_Y(), power_cal_Y()]; considering variable cluster sizes [sample_size_cal_new_Y(), power_cal_new_Y()], generating simulation data [datagen()]; performing simulations[sim()]; summarizing sample size requirement [sample_size_mar()].
  3. Simulation_Examples.R: This file includes four simulation examples each considering mean cluster sizes of {20, 50, 100}, outcome intraclass correlation coefficients of {0.01,0.05,0.10}, covariate intraclass correlation coefficients of {0.10, 0.25, 0.50}, and 5000 simulations. 
      Example 1 considers powering the heterogeneous treatment effect for a continuous endpoint with delta=0.10 and cv=0.0. 
      Example 2 considers powering the heterogeneous treatment effect for a binary endpoint with delta=0.45 and cv=0.6. 
      Example 3 considers powering the conditional overall effect for a continuous endpoint with delta=0.25 and cv=0.3.
      Example 4 considers powering the marginal overall effect for a continuous endpoint with delta=0.15 and cv=0.0.
