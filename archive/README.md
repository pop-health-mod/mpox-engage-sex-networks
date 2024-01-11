## archive
Code from previous versions of the analyses. 

- [`24_distr_r0_ngm_case.R`](./archive/24_distr_r0_ngm_case.R) computes the mpox _R<sub>0</sub>_ using two separate methods: the next-generation matrix (NGM) method based on the fitted distribution of sexual partner numbers estimated above, and the growth-rate method using the reported mpox cases over time. These two _R<sub>0</sub>_ estimation methods are then combined to estimate the mpox per-partnership secondary attack rate (SAR).
	- This script also projects the mpox _R<sub>0</sub>_ based no pre-COVID-19 sexual activity levels based on the NGM method and the average SAR across cities.
	- [`25_distr_r0_figs.R`](./archive/25_distr_r0_figs.R) generates [__Figure 2__] (mpox _R<sub>0</sub>_ results and 'projected' _R<sub>0</sub>_).

NB: _R<sub>0</sub>_ is the basic reproduction number and _R<sub>e</sub>_ is the effective or instantaneous reproduction number.
