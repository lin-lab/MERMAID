
require(data.table)

ndays <- 365

## These seeds are used to ensure covariates are identical 
## across simulation replicates.

seeds <- c(422200,911591,37997,223277,876269,230747,651954,609235)

runif(1);
print(head(.Random.seed))

rpois_const_seed <- function(n, lam, seed_val) {
	old <- .Random.seed
	set.seed(seed_val)
	res <- rpois(n, lam)
	.Random.seed <<- old
	res
}

rnorm_const_seed <- function(n, seed_val) {
  old <- .Random.seed
  set.seed(seed_val)
  res <- rnorm(n)
  .Random.seed <<- old
  res
}

max_nr <- 100

z1_list <- split(rnorm_const_seed(max_nr * ndays, seed_val = seeds[1]), rep(1:max_nr, each = ndays))
z2_list <- split(rnorm_const_seed(max_nr * ndays, seed_val = seeds[2]), rep(1:max_nr, each = ndays))

zr_list <- split(rnorm_const_seed(max_nr * 8, seed_val = seeds[3]), rep(1:max_nr, each = 8))
zi_list <- split(0.05 * rnorm_const_seed(max_nr * 1, seed_val = seeds[4]), rep(1:max_nr, each = 1))

library(splines)

options(stringsAsFactors=FALSE)

args <- commandArgs(TRUE)
num_regions <- as.integer(args[1])
sim_rep <- as.integer(args[2])


# Parameters are correctly specified
lag_p_fit_size <- 5
mu_fit <- 4.7
sigma_fit <- 2.9

P_sero <- 0.01 

out_dir <- "simulation_output_custom_PC1"

system(paste0("mkdir -p ", out_dir))

out_file <- paste0(out_dir, "/main_regions_", num_regions, "_rep_", sim_rep, "_out.rds")

init_Lam <- 50

n_prev <- 6

prev_w <- max(1.0, 0.25/P_sero)
reweight_prev <- FALSE

source("src/helper_functions.R")
source("src/plotting_functions.R")
source("src/cumulative_incidence.R")
source("src/fit_MERMAID.R")

source("simulation_scripts/simulation_functions_simple_custom.R")

rt_spline_df <- 8

date <- as.Date("2020-01-01") + seq(ndays) - 1

X_rt <- model.matrix(~ cut(date, breaks = 12 ))[,-1]

# -------------------------------------

# Lag probabilities (now reverse fitted vs simulated)

lag_p_fit <- dnbinom(0:20, mu = 5, size= 5)
lag_p_fit <- lag_p_fit/sum(lag_p_fit)

lag_p <- dnbinom(0:20, mu = 5, size= lag_p_fit_size)
lag_p <- lag_p/sum(lag_p)

# Serial interval weights (now reversed)

si_weights_fit <- discr_si(1:30, mu = 4.7, sigma =  2.9)
si_weights <- discr_si(1:30, mu = mu_fit, sigma = sigma_fit)

# ------------------------------------

lag_p_weekday <- NULL

use_brownian <- FALSE

intervention_eff <- 0.00

intervention_df_list <- lapply(1:max_nr, function(i){
	out <- data.frame(Intervention = rep(zi_list[[i]], ndays))
	out
})

ppos_vec_list <- lapply(z1_list, function(z1){
	out <- 0.25 * (2 + sin( (36/ndays) * pi * (1:ndays) ) + 0.25 * z1 ) 
	out - mean(out)
})
ntest_vec_list <- lapply(z2_list, function(z2){
	out <- 5 * (1 + cos( (19/ndays) * pi * (1:ndays) ) + 0.25 * z2 ) + sqrt(1:ndays)
	out - mean(out)
})

mean_ascertain <- 0.5
beta_ascertain <- c(0.05, 0.00)

rt_beta_vec_list <- list(
c( -0.3859,-0.4079,-0.4810,-0.5631,-0.3805,-0.3212,-0.2718,-0.0988,-0.5504,-0.5593,-0.9889 ),
c( -0.275,-0.600,-0.614,-0.738,-0.649,-0.466,-0.376,-0.288,-0.569,-0.696,-1.135 ),
c( -0.192,-0.489,-0.545,-0.527,-0.555,-0.372,-0.411,-0.200,-0.404,-0.578,-0.927 ),
c( -0.132,-0.493,-0.375,-0.553,-0.519,-0.190,-0.228,-0.235,-0.172,-0.640,-0.870 )
)

pop_size <- 8000000

N_sero <- round(P_sero * pop_size)

rt_global_mean_val <-  1.00

n_sim <- 0
sim_list <- list()
sero_list <- list()

print(paste("Num regions =", num_regions ))

while( n_sim < num_regions ){

	repeat{
		options(na.action = 'na.omit')

		rt_mean_val <- rt_global_mean_val

		ntest_vec <- ntest_vec_list[[ n_sim + 1 ]]
		ppos_vec <- ppos_vec_list[[ n_sim + 1 ]]

		rt_beta_vec <- rt_beta_vec_list[[ n_sim + 1 ]]
		intervention_df <- intervention_df_list[[ n_sim + 1 ]]	

		sim <- simulate_data(time_points = ndays, X_rt, Lambda_0 = init_Lam, mean_rt = rt_mean_val, rt_beta = rt_beta_vec, mean_pi = mean_ascertain, ntest = ntest_vec, ppos = ppos_vec, alpha_tests = beta_ascertain[1], alpha_ppos = beta_ascertain[2], lag_probs = lag_p, weekday_lag_probs = lag_p_weekday, XC = matrix(intervention_df$Intervention, ncol = 1), Beta_XC = intervention_eff, si_weights = si_weights )
	
		sim$data$Intervention <- intervention_df$Intervention

		print( quantile(sim$data$C_lag, 0.05))
		print( sum(sim$data$Y) )
		if(all(!is.na(sim$data$C + sim$data$Y)) ){
			if( sum(sim$data$Y) < 0.5*pop_size & sum(sim$data$Y) >= 0.03*pop_size & quantile(sim$data$C_lag, 0.05) >= 5 ){
				break;
			}
		}else{
			print("NAs")
		}
	}
	
	n_sim <- n_sim + 1
	
	sim$data$pop_size <- pop_size 
	
	select_dates <-  function(x,n) x[round(quantile(seq_along(x), c(1:n)/(n+1)))]
	
	prev_dates <- select_dates( sim$data$date, n_prev )
	
	sim$data[,csY := cumsum(Y),]
	sim$data[,prev := cumsum(Y)/pop_size,]
	print( summary(sim$data$prev))
	
	sim$data[,rep_prev := cumsum(C_lag)/pop_size,]
	
	sero_dat <- subset(sim$data, date %in% prev_dates)[,list(ntested = N_sero, csY, prev),by=date][,list(
		npos = rhyper(1, k = ntested, m = csY, n = pop_size - csY ), region = n_sim
	),by=list(date, ntested)]
	
	sero_list[[n_sim]] <- sero_dat
	
	sim$data$region <- n_sim
	
	sim_list[[n_sim]] <- sim
}

sero_data <- rbindlist(sero_list)

sim_data <- rbindlist(lapply(sim_list, function(x) x$data))

print(head(sim_data))

if( num_regions == 1 ){

	rt_eqn <- ~ cut(date, breaks = 12)
	pi_eqn <- ~ ntest

}else{

	rt_eqn <- ~ Intervention + factor(region)*cut(date, breaks = 12)
	pi_eqn <- ~ factor(region)*ntest

}


sim_data <- as.data.frame(copy(sim_data))
sero_data <- as.data.frame(copy(sero_data))

sim_par_list <- list(
	rt_formula = rt_eqn, 
	pi_formula = pi_eqn, 
	pi_intercept = 0, 
	pi_init_mean = mean(sim_data$pi_t_true),
	data = as.data.frame(sim_data), 
	prev_data= as.data.frame(sero_data),
	prev_weight = prev_w,
	lag_probs = lag_p_fit,
	fix_lag_probs = TRUE,
	estimate_weekday_lag = FALSE,
	si_weights = si_weights_fit,
	si_nts = length(si_weights_fit),
	confirmed_cases_var = "C_lag", 
	date_var = "date", 
	init_Lambda = init_Lam,
	adj_Lambda = 1,
	subset_var = NULL, 
	region_var = "region", 
	pop_size_var = NULL, 
	max_it_E_step = 5,
	tol = 1e-3,
	max_it = 1000,
	plot_every = NA,
	recovered_offset = FALSE,
	recycle_start = TRUE,
	return_vcov = TRUE,
	accelerated_em = TRUE,
	reweight_prev = reweight_prev,
	prev_method = 'hyper',
	print_mem = TRUE
)

sim_fit <- do.call(fitMERMAID, sim_par_list)

sim_fit$data[,`:=`(
	'num_regions'=num_regions,
	'sim_rep'=sim_rep,
	'lag_p_fit_size'=lag_p_fit_size,
	'mu_fit'=mu_fit,
	'sigma_fit'=sigma_fit
),]

saveRDS(
	object = list(sim_list=sim_list, sim_part_list = sim_par_list, sim_data=sim_data, sero_data = sero_data, sim_fit = sim_fit),
	file = out_file
)



