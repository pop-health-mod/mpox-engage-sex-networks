empiric_nb_partn_p6m
	time_pt:	time period (pre-pandemic, pandemic, or post-restrictions);
	x:			number of partners in P6M;
	n:			number of individuals reporting x partners;
	p:			proportion of individuals reporting x partners;
				.rds_ipw suffix indicate that n and p are weighted with the RDS-II weights (pre-pandemic period) or RDS-IPC weights (pandemic and post-restrictions periods), hence those distributions should be representative of the whole population.

fitted_pmf_nb_partn_p6m
	time_pt:	time period (pre-pandemic, pandemic, or post-restrictions);
	y_pred:		number of partners in P6M;
	mean:		proportion of individuals reporting 'y_pred' partners (RDS-II or RDS-IPC weighted, depending on the time period);
	cr.i_low:	lower 95% credible interval around 'mean';
	cr.i_upp:	upper 95% credible interval around 'mean'.

fitted_cdf_nb_partn_p6m
	same as the pmf file above except that 'mean' is the proportion of individuals reporting >=y_pred partners.
