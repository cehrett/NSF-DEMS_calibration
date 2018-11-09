To reproduce figure 2 from the paper, run the script "reproduce_figure.m". 

This will call "MCMC_settings.m" to package settings for the calibration.

"MCMC_calib", also called by "reproduce_figure.m", performs the actual calibration. 

"Ex_sim.m" is the example computer model in this case, which is used directly (rather than a code surrogate). 

"gp_cov.m" builds the covariance matrix for the Gaussian process discrepancy function. 

"direct_fn_output.mat" stores the output of "Ex_sim.m" over a 1000x1000 grid on the parameter space [0,3]x[0,6].