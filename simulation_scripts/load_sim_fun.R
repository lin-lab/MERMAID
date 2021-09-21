
source("src/helper_functions.R")

source("src/cumulative_incidence.R")

init_Lam <- 50
si_weights <- discr_si(1:30, mu = 4.7, sigma =  2.9)

simple_deconv <- function(x, probs, lags = seq_along(probs) - 1){
        out <- rep(0, length(x))

        for(i in 1:length(probs)){
                out <- out + data.table::shift(x, -lags[i], fill =  tail(x,1) )*probs[i]
        }
        out
}

calc_Lambda <- function(yy, si_w = si_weights, Lambda_0 = init_Lam){
	n <- length(yy)
	Lambda <- rep(NA, n)

	for(i in 1:n){
	  Lambda[i] <- calc_Lam(yy, si_w, i, fill = Lambda_0)
	}
	Lambda
}

require(data.table)

process_simulation_output <- function(rds_file, rep_name = gsub(".rds", "", rds_file), orthog_beta = TRUE, glm_beta = FALSE){
	require(Matrix)
 
	logit <- function(x) log(x) - log(1-x)
 
	if(!file.exists(rds_file)) return(NULL)

	d <- readRDS(rds_file)

	nr <- length(unique(d$sim_data$region))

	rt_beta_list <- lapply(d$sim_list, function(x){
		x$par$rt_beta
	} ) 


	lag_p <- dnbinom(0:20, mu = 5, size= 5)
	lag_p <- lag_p/sum(lag_p)
	

	d$sim_data <- as.data.table(d$sim_data)[order(region,time)]

	d$sim_data[,Y_0 := C_lag,]
	d$sim_data[,Y_1 := simple_deconv(C_lag, probs = lag_p, lags = 0:20),]
	
	d$sim_data[,Lambda_0 := calc_Lambda(Y_0, Lambda_0 = 50),]
	d$sim_data[,Lambda_1 := calc_Lambda(Y_1, Lambda_0 = 50),]
	
	rt_spline_df <- length(rt_beta_list[[1]]) - 2
	
	if( any(grepl("^cut", d$sim_fit$Rt_coef$Variable)) ){
		
		X_rt <- model.matrix(~ cut(date, breaks = 12), d$sim_data)
		
		rt_formula <- Y ~ offset(log(Lambda)) + cut(date, breaks = 12)
		
		rt_formula_0 <- Y_0 ~ offset(log(Lambda_0)) + cut(date, breaks = 12)
		rt_formula_1 <- Y_1 ~ offset(log(Lambda_1)) + cut(date, breaks = 12)
		
	}else{
		
		if( grepl("BS1", rds_file) ){
			print("Reverse!")
			
			X_rt <- model.matrix(~ splines::bs(max(date) - date, degree = 3, df = rt_spline_df), d$sim_data)
			
			rt_formula <- Y ~ offset(log(Lambda)) + splines::bs(max(date) - date, degree = 3, df = rt_spline_df)
				
			rt_formula_0 <- Y_0 ~ offset(log(Lambda_0)) + splines::bs(max(date) - date, degree = 3, df = rt_spline_df)
			rt_formula_1 <- Y_1 ~ offset(log(Lambda_1)) + splines::bs(max(date) - date, degree = 3, df = rt_spline_df)
			
		}else{
					
			X_rt <- model.matrix(~ splines::bs(date, degree = 3, df = rt_spline_df), d$sim_data)
			
			rt_formula <- Y ~ offset(log(Lambda)) + splines::bs(date, degree = 3, df = rt_spline_df)
						
			rt_formula_0 <- Y_0 ~ offset(log(Lambda_0)) + splines::bs(date, degree = 3, df = rt_spline_df)
			rt_formula_1 <- Y_1 ~ offset(log(Lambda_1)) + splines::bs(date, degree = 3, df = rt_spline_df)
		}
		
	}
	
	
	rt_fit <- stats::glm(rt_formula, family = poisson(), d$sim_data)
		
	rt_fit_0 <- stats::glm(rt_formula_0, family = poisson(), d$sim_data)
		
	rt_fit_1 <- stats::glm(rt_formula_1, family = poisson(), d$sim_data)
		
	
	rt_y_hat <- predict(rt_fit, type = 'response')/d$sim_data$Lambda
	rt_y_hat_se <- predict(rt_fit, type = 'response', se.fit=TRUE)$se.fit/d$sim_data$Lambda
	
	rt_0_hat <- predict(rt_fit_0, type = 'response')/d$sim_data$Lambda_0
	rt_0_hat_se <- predict(rt_fit_0, type = 'response', se.fit=TRUE)$se.fit/d$sim_data$Lambda_0
	
	rt_1_hat <- predict(rt_fit_1, type = 'response')/d$sim_data$Lambda_1
	rt_1_hat_se <- predict(rt_fit_1, type = 'response', se.fit=TRUE)$se.fit/d$sim_data$Lambda_1
	
	y_0_hat <- d$sim_data$Y_0
	y_0_hat_se <- 0 + NA
	
	y_1_hat <- d$sim_data$Y_1 
	y_1_hat_se <- 0 + NA
	
	NR <- nrow(X_rt) - 1
	
	X_svd <- svd(X_rt)
	
	A <- diag(X_svd$d) %*% t(X_svd$v) #/sqrt(NR)
	
	Vb <- d$sim_fit$vcov[1:ncol(X_rt),1:ncol(X_rt)]
	b <- d$sim_fit$Rt_coef$Estimate
	
	Beta <- head(d$sim_list[[1]]$par$rt_beta,ncol(X_rt))
	
	Vg <- A %*% Vb %*% t(A)
	g <- (A %*% b)[,1]
	g_SE <- sqrt(diag(Vg))
	
	Gamma <- (A %*% Beta)[,1]
	
	
	if( nr > 1 ){
		
		Intercepts <- c("(Intercept)", paste0("factor(region)", 2:nr))
		prefixes <- c("", paste0("factor(region)", 2:nr, ":"))
		Slopes <- paste0("bs(date, degree = 3, df = rt_spline_df)", 1:8)
		
		rt_beta <- c()
		intervention_beta <- 0
		rt_terms <- c("Intervention")
		for(i in 1:nr){
			rt_terms <- c(rt_terms, Intercepts[i], paste0(prefixes[i], Slopes))
			rt_beta <- c(rt_beta, head(rt_beta_list[[i]], -1))
			intervention_beta <- intervention_beta + tail(rt_beta_list[[i]], 1)/nr
		}
		rt_beta_true <- c(intervention_beta, rt_beta)

		rt_k <- match(rt_terms, d$sim_fit$Rt_coef$Variable)
		
	}else{
			
		Intercepts <- c("(Intercept)")
		prefixes <- c("")
		Slopes <- paste0("bs(date, degree = 3, df = rt_spline_df)", 1:8)
		
		rt_beta <- c()
		intervention_beta <- 0
		rt_terms <- c()
		for(i in 1:nr){
			rt_terms <- c(rt_terms, Intercepts[i], paste0(prefixes[i], Slopes))
			rt_beta <- c(rt_beta, head(rt_beta_list[[i]], -1))
			intervention_beta <- intervention_beta + tail(rt_beta_list[[i]], 1)/nr
		}
		rt_beta_true <- c(rt_beta)

		rt_k <- match(rt_terms, d$sim_fit$Rt_coef$Variable)
		

	}

	pi_beta_true <- d$sim_list[[1]]$par$pi_beta

	p_0 <- as.data.table(d$sero_data)[,list(date, region, prev = npos/ntested, n = ntested),][order(region,date)]

	p_1 <- merge(p_0[,list(date,region),], as.data.table(d$sim_fit$data)[,list(date,region,prev_1=Prevalence,prev_1_se = Prevalence_SE),], by = c("date", "region"))[order(region,date)]

	delta <- max(abs(p_0$prev - p_1$prev_1))
	print(delta)
	
	Vcov <- solve(d$sim_fit$J_w)
		
		
	J_1 <- as.matrix(d$sim_fit$J_1)
	set_0 <- which(diag(J_1) <= 0)
	J_1[set_0,] <- 0
	J_1[,set_0] <- 0
	Vcov_1 <- MASS::ginv(J_1)
	
	N <- nrow(Vcov)
	
	np_c <- 2
	
	i_a <- (N - np_c + 1):N
	i_r <- 1:(N - np_c)
	
	Vcov_pi <- Vcov[i_a, i_a]
	
	Vcov_rt <- Vcov[i_r, i_r]
	Vcov_rt_1 <- Vcov_1[i_r, i_r]
	
	if( any(is.na(diag(Vcov_rt_1)) ) ){
		Vcov_rt_1 <- solve(J_1[i_r, i_r])
	}else if( any(diag(Vcov_rt_1) <= 0) ){
		Vcov_rt_1 <- solve(J_1[i_r, i_r])
	}
	
	pi_t <- d$sim_fit$data$pi_t
	
	X <- cbind(1,d$sim_fit$data$ntest, d$sim_fit$data$ppos )[,1:np_c]

	
	SE_pi <- pi_t*(1 - pi_t) * sqrt(rowSums(tcrossprod(X, Vcov_pi) * X))
	
	SE_Rt_ww <- d$sim_fit$data$R_t * sqrt(rowSums(tcrossprod(X_rt, Vcov_rt) * X_rt))
	
	SE_Rt_w1 <- d$sim_fit$data$R_t * sqrt(rowSums(tcrossprod(X_rt, Vcov_rt_1) * X_rt))
	
	if( orthog_beta ){
		dt_beta <- data.table(
			param = "Rt_beta",
			variable = paste0("Gamma_", 1:length(g)), #[rt_k],
			theta_true =  Gamma, 
			theta_hat = g, #[rt_k],
			theta_hat_SE = g_SE #[rt_k]
		)
	}else{
		dt_beta <- data.table(
			param = "Rt_beta",
			variable = d$sim_fit$Rt_coef$Variable, #[rt_k],
			theta_true =  head(d$sim_list[[1]]$par$rt_beta,-1), 
			theta_hat = d$sim_fit$Rt_coef$Estimate, #[rt_k],
			theta_hat_SE = d$sim_fit$Rt_coef$SE #[rt_k]
		)
	}
	
	
	SE_weighted_info <- sqrt(diag(Vcov))
	
	pi_beta_SE <- tail(SE_weighted_info, np_c)
	Rt_beta_SE <- head(SE_weighted_info, -np_c)

	out <- rbindlist(list(
		dt_beta,
		data.table(
			param = "pi_beta",
			variable = d$sim_fit$pi_coef$Variable,
			theta_true = d$sim_list[[1]]$par$pi_beta[1:np_c],
			theta_hat = d$sim_fit$pi_coef$Estimate,
			theta_hat_SE = pi_beta_SE #d$sim_fit$pi_coef$SE
		),
		d$sim_fit$data[,list(
			param = "Rt", 
			variable = paste(region, date, sep = ":"),
			theta_true = R_t_true, 
			theta_hat = R_t, 
			theta_hat_SE = SE_Rt_w1
		),],
		d$sim_fit$data[,list(
			param = "Rt_ww", 
			variable = paste(region, date, sep = ":"),
			theta_true = R_t_true, 
			theta_hat = R_t, 
			theta_hat_SE = SE_Rt_ww
		),],
		d$sim_fit$data[,list(
			param = "Rt_fit", 
			variable = paste(region, date, sep = ":"),
			theta_true = R_t_true, 
			theta_hat = rt_y_hat, 
			theta_hat_SE = rt_y_hat_se
		),],
		d$sim_fit$data[,list(
			param = "Rt_fit_0", 
			variable = paste(region, date, sep = ":"),
			theta_true = R_t_true, 
			theta_hat = rt_0_hat, 
			theta_hat_SE = rt_0_hat_se
		),],
		d$sim_fit$data[,list(
			param = "Rt_fit_1", 
			variable = paste(region, date, sep = ":"),
			theta_true = R_t_true, 
			theta_hat = rt_1_hat, 
			theta_hat_SE = rt_1_hat_se
		),],
		d$sim_fit$data[,list(
			param = "pi_t", 
			variable = paste(region, date, sep = ":"),
			theta_true = pi_t_true, 
			theta_hat = pi_t, 
			theta_hat_SE = SE_pi # pi_t_SE
		),],
		d$sim_fit$data[,list(
			param = "pi_t_w1", 
			variable = paste(region, date, sep = ":"),
			theta_true = pi_t_true, 
			theta_hat = pi_t, 
			theta_hat_SE = pi_t_SE
		),],
		d$sim_fit$data[,list(
			param = "Infections", 
			variable = paste(region, date, sep = ":"),
			theta_true = Y, 
			theta_hat = Infections, 
			theta_hat_SE = Infections_SE
		),],
		d$sim_fit$data[,list(
			param = "Infections_0", 
			variable = paste(region, date, sep = ":"),
			theta_true = Y, 
			theta_hat = y_0_hat, 
			theta_hat_SE = y_0_hat_se
		),],
		d$sim_fit$data[,list(
			param = "Infections_1", 
			variable = paste(region, date, sep = ":"),
			theta_true = Y, 
			theta_hat = y_1_hat, 
			theta_hat_SE = y_1_hat_se
		),],
		d$sim_fit$data[,list(
			param = "Prevalence", 
			variable = paste(region, date, sep = ":"),
			theta_true = cumsum(Y)/pop_size, 
			theta_hat = Prevalence, 
			theta_hat_SE = Prevalence_SE
		),],
		d$sim_fit$data[,list(
			param = "Prevalence_0", 
			variable = paste(region, date, sep = ":"),
			theta_true = cumsum(Y)/pop_size, 
			theta_hat = cumsum(y_0_hat)/pop_size, 
			theta_hat_SE = as.numeric(NA)
		),],
		d$sim_fit$data[,list(
			param = "Prevalence_1", 
			variable = paste(region, date, sep = ":"),
			theta_true = cumsum(Y)/pop_size, 
			theta_hat = cumsum(y_1_hat)/pop_size, 
			theta_hat_SE = as.numeric(NA)
		),]
	))

	out[,`:=`(
		'theta_hat_z' = (theta_hat - theta_true)/theta_hat_SE,'num_regions'=d$sim_fit$data$num_regions[1],'sim_rep'=d$sim_fit$data$sim_rep[1],'lag_p_fit_size'=d$sim_fit$data$lag_p_fit_size[1],'mu_fit'=d$sim_fit$data$mu_fit[1],'sigma_fit'=d$sim_fit$data$sigma_fit[1]
	),]

	out

}

