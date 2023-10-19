# mpox-engage-sex-networks
Code to analyze data on sexual networks and its impacts on mpox transmission during the 2022-2023 outbreak.

The results are presented in the pre-print titled _[Characteristics of the sexual networks of gay, bisexual, and other men who have sex with men in Montréal, Toronto, and Vancouver: implications for the transmission and control of mpox in Canada](https://www.medrxiv.org/content/10.1101/2023.08.31.23294912v1)_.

## code structure
Briefly, the first number in the file prefix means the following
- `1`: descriptive analyses
- `2`: main analyses (degree distribution, mpox _R<sub>0</sub>_ estimation)
- `3`: sensitivity and supplementary analyses (alternative outcome: anal partnerships, restriction, covariate standardization), mpox _R<sub>t</sub>_ from mpox incidence data.

The second number is simply to maintain files in order.

NB: _R<sub>0</sub>_ is the basic reproduction number and _R<sub>t</sub>_ is the effective or instantaneous reproduction number.

### descriptive analyses
- [`11_descr_table1.R`](11_descr_table1.R) produces __Table 1__ (correlates of & summary statistics of number of partners).
- [`12_descr_retention_ipcw.R`](12_descr_retention_ipcw.R) computes the IPCWs (details in Supplementary Materials section _RDS-II weights and inverse probability of censoring weights_) and produces __Table S2__ (retention of participants).

### main analyses
- Scripts `21` through `23` fit the negative binomial regression models ([`21_distr_fit_negbin_coef.R`](21_distr_fit_negbin_coef.R)), perform the post-stratification ([`src-stan/regression_negbin_aggregate.stan`](src-stan/regression_negbin_aggregate.stan)), and compute the fitted population-wide distribution of the number of sexual partners based on the post-stratified samples (script [`22_distr_fit_negbin_predictive.R`](22_distr_fit_negbin_predictive.R)).
	- Script [`21`](21_distr_fit_negbin_coef.R) also produces __Table S3__ (regression coefficients).
	- The post-stratification procedure is explained in detail in Supplementary Materials section _Distribution of sexual partner numbers (Bayesian regression and post-stratification)_.
	- [`23_distr_fit_negbin_figs.R`](23_distr_fit_negbin_figs.R) generates [__Figure 1__](fig/fig_1_cdf_main_model.png) (CDF of the number of sexual partners in P6M in each city) and [__Figure S1__](fig/fig_S1_cdf_model_fit_to_data.png) (model fit).
- [`24_distr_r0_ngm_case.R`](24_distr_r0_ngm_case.R) computes the mpox _R<sub>0</sub>_ using two separate methods: the next-generation matrix (NGM) method based on the fitted distribution of sexual partner numbers estimated above, and the growth-rate method using the reported mpox cases over time. These two _R<sub>0</sub>_ estimation methods are then combined to estimate the mpox per-partnership secondary attack rate (SAR).
	- This script produces [__Figure S2__](fig/fig_2_r0_main.png) (cumulative mpox incidence in each province).
	- This script also projects the mpox _R<sub>0</sub>_ based no pre-COVID-19 sexual activity levels based on the NGM method and the average SAR across cities.
	- [`25_distr_r0_figs.R`](25_distr_r0_figs.R) generates [__Figure 2__](fig/fig_2_r0_main.png) (mpox _R<sub>0</sub>_ results and 'projected' _R<sub>0</sub>_).

### sensitivity and supplementary analyses
- [`31_supp_restriction.R`](31_supp_restriction.R) performs the restriction sensitivity analysis, which fits the regression models restricting the analytical sample to participants who had visits at all three timepoints. Results are in [__Figure S3__](fig/fig_S3_cdf_main_vs_restriction.png).
- [`32_supp_standardization.R`](32_supp_standardization.R) performs the regression-based standardization. Results are in [__Figure S4__](fig/fig_S4_cdf_main_vs_standardization.png).
- Scripts `21` through `23` described in the previous section can be modified to use the alternative anal partners in P6M outcome variable. Results are shown in [__Figure S5__](fig/fig_S5_cdf_main_all_vs_anal.png).
- [`33_supp_sensitivity_figs.R`](33_supp_sensitivity_figs.R) generates __Figures S3--S5__ from the analyses described above.
- [`34_supp_rt_case.R`](34_supp_rt_case.R) estimates the _R<sub>t</sub>_ using the _EpiEstim_ package, and results of which are shown in [__Figure S6__](fig/fig_S6_rt.png).
