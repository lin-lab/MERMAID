


expit <- function(x) 1/(1 + exp(-x))
logit <- function(x) log(x) - log(1-x)

# simple "deconvolution" (technically, convolution)
simple_deconv <- function(x, probs, lags = seq_along(probs) - 1){
	out <- rep(0, length(x))

	for(i in 1:length(probs)){
		out <- out + data.table::shift(x, -lags[i], fill = 0)*probs[i]
	}
	out
}

logit <- function(x) log(x) - log(1 - x)

calc_shift <- function(vec, shifts, skip_cols = NULL){
	nk <- length(shifts)
	mat <- matrix(c(vec), ncol = nk)
	out <- rep(0, nrow(mat))
	for(i in 1:nk){
		if( i %in% skip_cols ) next;
		out <- out + data.table::shift(mat[,i],shifts[i],fill=0)
	}
	out
}

lag_mat <- function(i,n, reverse = FALSE){
	require(Matrix)
	if( reverse ){
		Matrix::Diagonal(n = n) - lag_mat(i=i, n=n, reverse = FALSE)
	}else{
		Matrix::sparseMatrix(
			i = (1:(n-i))+i, j = (1:(n-i)),
			dims = c(n,n)
		)
	}
}

# Calculate the marginal moments of Y = C + U
calc_Y_moments <- function(R_t, w_t, Lambda_0 = 1, max_window = Inf, calc_var = TRUE){

	require(Matrix)

	nr <- length(R_t)
	nb <- length(w_t)
	
	if( length(Lambda_0) == 1 ){
		Lambda_0 <- c(1,rep(0, nr-1))
	}

	Sigma <- NULL
	EY <- rep(NA, nr)
	
	for(i in 1:nr){
		
		if( i == 1 ){ 
			EY[i] <- R_t[i] * Lambda_0[i]
		}else{
			EY[i] <- R_t[i] * (sum(  w_t[min(nb, max(i-1,1)):1] * EY[max(i - nb, 1):(i-1)] ) + Lambda_0[i])
		}
	}
	
	if( calc_var ){
		Sigma <- matrix(NA, nr, nr)
		for(i in 1:nr){
		
			if( i == 1 ){ 
				Sigma[i,i] <- EY[i] 
			}else{
				ww <- w_t[min(nb, max(i-1,1)):1]
				ii <- max(i - nb, 1):(i-1)
				Sigma[i,i] <- EY[i] + (R_t[i]^2) * sum(crossprod(Sigma[ii,ii], ww) * ww)
			}
		
		
			if( i == nr ) break;
			for( j in (i+1):nr ){
			
				if( j - i > max_window ){ 
					Sigma[i,j:nr] <- Sigma[j:nr,i] <- 0
					break;
				}
			
				Sigma[i,j] <- Sigma[j,i] <- R_t[j] * sum( w_t[min(nb, max(j-1,1)):1] * Sigma[i,max(j - nb, 1):(j-1)] )
			
			}
		
		}
		list( mu = EY, Sigma = as(Sigma, "RsparseMatrix"))
	}else{
	
		list( mu = EY)
	}

}

# Calculate marginal moments of Y = C + U across regions
calc_Y_moments_by_region <- function(R_t, w_t, region, Lambda_0 = 1, max_window = Inf, calc_var = TRUE){

	require(Matrix)

	Lambda_0_list <- NULL
	if( length(Lambda_0) != 1 ){
		Lambda_0_list <- base::split(x = Lambda_0, f = region)
	}

	ii <- 0
	out_list <- lapply(base::split(x = R_t, f = region), function(x){
		ii <<- ii + 1
		if( is.null(Lambda_0_list) ){
			Lambda_0_region <- Lambda_0
		}else{
			Lambda_0_region <- Lambda_0_list[[ii]]
		}
		calc_Y_moments(R_t = x, w_t = w_t, Lambda_0 = Lambda_0_region, max_window = max_window, calc_var = calc_var)
	})

	list(
		mu = Reduce(c, lapply(out_list, function(x) x$mu)),
		Sigma = if(calc_var){
			Reduce(Matrix::bdiag, lapply(out_list, function(x) x$Sigma))
		}else{
			NULL
		}
	)
}

# Calculate Lambda recursively
calc_Lam <- function(y, w, t, fill){
	if( t <= 1 ){
		return(fill)
	}
	np <- min(length(w), t - 1)	
	yy <- y[(t - 1):(t - np)]
	if(length(yy) < length(w)) yy <- c(yy, rep(fill, length(w) - length(yy)))
	sum(yy * w)
}

# Make weight matrix for calculating Lambda
makeLTW <- function(vals, nr, npt = NA){
	require(Matrix)
	
	npt <- min(npt, length(vals), na.rm=TRUE)
	vals <- vals[1:npt]

	dd <- do.call(rbind.data.frame, lapply( 1:nr, function(ii){
		data.frame( ii = ii, jj =  ii - (1:npt), xx = vals)
	} ))

	with(
			subset(dd, jj > 0 & jj < ii),
			Matrix::sparseMatrix( i = ii, j = jj, x = xx, dims = c(nr, nr))
	)
}



calc_log_lik <- function(C, U, pi_t, R_t, W){
	Y <- U + C
	Lambda <- (W %*% Y)[,1] + 1
	RL <- R_t * Lambda
	eta_pi <- log(pi_t) - log(1- pi_t)
	
	logLik_Y <- sum(Y * log( RL )) - sum(RL) - sum(lgamma(Y + 1))
	logLik_C_Y <- sum( C * eta_pi ) + sum( Y * log(1 - pi_t) ) + sum(lgamma(Y + 1)) - sum(lgamma(Y-C + 1)) - sum(lgamma(C + 1))
	
	logLik_Y + logLik_C_Y
}

get_log_lik_fun <- function(pi_t, R_t, W){

	function(x){
		if(min(x) < 0) return(-Inf)
	
		U <- x[1:length(pi_t)]
		C <- x[(1:length(pi_t)) + length(pi_t)]
		Y <- U + C
		Lambda <- (W %*% Y)[,1] + 1
		RL <- R_t * Lambda
		eta_pi <- log(pi_t) - log(1- pi_t)
		
		logLik_Y <- sum(Y * log( RL )) - sum(RL) - sum(lgamma(Y + 1))
		logLik_C_Y <- sum( C * eta_pi ) + sum( Y * log(1 - pi_t) ) + sum(lgamma(Y + 1)) - sum(lgamma(Y-C + 1)) - sum(lgamma(C + 1))
		
		logLik_Y + logLik_C_Y
	}
}

get_grad_fun <- function(pi_t, R_t, W, const = 1){

	function(x){
		U <- x[1:length(pi_t)]
		C <- x[(1:length(pi_t)) + length(pi_t)]
		Y <- U + C
		Lambda <- (W %*% Y)[,1] + const
		
		eta0 <- log(1- pi_t)
		eta1 <- log(pi_t) - eta0
		dg_Y_C <- digamma(Y - C + 1)
		
		vec <- (t(W) %*% diag(1/Lambda)  %*% Y )[,1]
		
		grad_Y <- log(Lambda) + log(R_t) + vec - dg_Y_C  + eta0 - (t(W) %*% R_t)[,1]
		
		grad_C <- eta1 + dg_Y_C - digamma(C + 1)
		
		c(
			grad_Y,
			grad_C
		)
	}
}

get_hessian_fun <- function(pi_t, R_t, W, const = 1){

	function(x){
		U <- x[1:length(pi_t)]
		C <- x[(1:length(pi_t)) + length(pi_t)]
		#U <- Y - C
		Y <- U + C
		Lambda <- (W %*% Y)[,1] + const
		
		eta0 <- log(1- pi_t)
		eta1 <- log(pi_t) - eta0
		dg_Y_C <- digamma(Y - C + 1)
		
		vec <- (t(W) %*% diag(1/Lambda)  %*% Y )[,1]
		
		grad_Y <- log(Lambda) + log(R_t) + vec - dg_Y_C  + eta0 - (t(W) %*% R_t)[,1]
		
		grad_C <- eta1 + dg_Y_C - digamma(C + 1)
		
		Wtfp <- crossprod(W, diag(1/Lambda))
		
		hess_YC <- Matrix::Diagonal( x = trigamma(Y - C + 1))
		hess_C <- (-1)*hess_YC - Matrix::Diagonal( x = trigamma(C + 1))
		hess_Y <- crossprod(W, diag(-Y/(Lambda^2) ) %*% W ) + Wtfp + t(Wtfp) - hess_YC
		
		rbind(
			cbind(hess_Y, hess_YC),
			cbind(hess_YC, hess_C)
		)
	}
}
