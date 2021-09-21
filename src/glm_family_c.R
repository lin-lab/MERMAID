

# Poisson and binomial families that do not 
# return errors for non-integer response. 

dpois_c <- function(x, lambda, log = FALSE){
	out <- x*log(lambda) - lambda - lgamma(x+1)
	if( !log ){
		exp(out)
	}else{
		out
	}
}

lchoose_c <- function(n, k) lgamma(n+1) - lgamma(n-k+1) - lgamma(k+1)

dbinom_c <- function(x, size, prob, log = FALSE){
	out <- lgamma(size+1) - lgamma(size-x+1) - lgamma(x+1) + x * log(prob) + (size-x) * log(1-prob)
	if( !log ){
		exp(out)
	}else{
		out
	}
}

# The following functions are copied 
# from stats::poisson and modified to
# omit error messages for non-integer
# values 

poisson_c <- function (link = "log")
{
    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    okLinks <- c("log", "identity", "sqrt")
    if (linktemp %in% okLinks)
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    }
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for poisson family; available links are %s",
                linktemp, paste(sQuote(okLinks), collapse = ", ")),
                domain = NA)
        }
    }
    variance <- function(mu) mu
    validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)
    dev.resids <- function(y, mu, wt) {
        r <- mu * wt
        p <- which(y > 0)
        r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
        2 * r
    }
    aic <- function(y, n, mu, wt, dev) -2 * sum(dpois_c(y, mu,
        log = TRUE) * wt)
    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the 'Poisson' family")
        n <- rep.int(1, nobs)
        mustart <- y + 0.1
    })
    simfun <- function(object, nsim) {
        wts <- object$prior.weights
        if (any(wts != 1))
            warning("ignoring prior weights")
        ftd <- fitted(object)
        rpois(nsim * length(ftd), ftd)
    }
    structure(list(family = "poisson", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
        validmu = validmu, valideta = stats$valideta, simulate = simfun),
        class = "family")
}

binomial_c <- function (link = "logit")
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) linktemp <- deparse(linktemp)
    okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
    if (linktemp %in% okLinks)
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    } else {
        ## what else shall we allow?  At least objects of class link-glm.
        if(inherits(link, "link-glm")) {
            stats <- link
            if(!is.null(stats$name)) linktemp <- stats$name
        } else {
	    stop(gettextf('link "%s" not available for binomial family; available links are %s',
			  linktemp, paste(sQuote(okLinks), collapse =", ")),
	     domain = NA)
        }
    }
    variance <- function(mu) mu * (1 - mu)
    validmu <- function(mu) all(is.finite(mu)) && all(mu>0 &mu<1)
    dev.resids <- function(y, mu, wt) .Call(stats:::C_binomial_dev_resids, y, mu, wt)
    aic <- function(y, n, mu, wt, dev) {
        m <- if(any(n > 1)) n else wt
	-2*sum(ifelse(m > 0, (wt/m), 0)*
               dbinom_c(m*y, m, mu, log=TRUE))
    }
    initialize <- expression({
	if (NCOL(y) == 1) {
	    ## allow factors as responses
	    ## added BDR 29/5/98
	    if (is.factor(y)) y <- y != levels(y)[1L]
	    n <- rep.int(1, nobs)
            ## anything, e.g. NA/NaN, for cases with zero weight is OK.
            y[weights == 0] <- 0
	    if (any(y < 0 | y > 1))
		stop("y values must be 0 <= y <= 1")
            mustart <- (weights * y + 0.5)/(weights + 1)
            m <- weights * y
            # if(any(abs(m - round(m)) > 1e-3))
                # warning("non-integer #successes in a binomial glm!")
	}
	else if (NCOL(y) == 2) {
            # if(any(abs(y - round(y)) > 1e-3))
                # warning("non-integer counts in a binomial glm!")
	    n <- y[, 1] + y[, 2]
	    y <- ifelse(n == 0, 0, y[, 1]/n)
	    weights <- weights * n
            mustart <- (n * y + 0.5)/(n + 1)
	}
	else stop("for the 'binomial' family, y must be a vector of 0 and 1\'s\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
    })
    simfun <- function(object, nsim) {
        ftd <- fitted(object)
        n <- length(ftd)
        ntot <- n*nsim
        wts <- object$prior.weights
        if (any(wts %% 1 != 0))
            stop("cannot simulate from non-integer prior.weights")
        ## Try to fathom out if the original data were
        ## proportions, a factor or a two-column matrix
        if (!is.null(m <- object$model)) {
            y <- model.response(m)
            if(is.factor(y)) {
                ## ignote weights
                yy <- factor(1+rbinom(ntot, size = 1, prob = ftd),
                             labels = levels(y))
                split(yy, rep(seq_len(nsim), each = n))
            } else if(is.matrix(y) && ncol(y) == 2) {
                yy <- vector("list", nsim)
                for (i in seq_len(nsim)) {
                    Y <- rbinom(n, size = wts, prob = ftd)
                    YY <- cbind(Y, wts - Y)
                    colnames(YY) <- colnames(y)
                    yy[[i]] <- YY
                }
                yy
            } else
            rbinom(ntot, size = wts, prob = ftd)/wts
        } else rbinom(ntot, size = wts, prob = ftd)/wts
    }
    structure(list(family = "binomial",
		   link = linktemp,
		   linkfun = stats$linkfun,
		   linkinv = stats$linkinv,
		   variance = variance,
		   dev.resids = dev.resids,
		   aic = aic,
		   mu.eta = stats$mu.eta,
		   initialize = initialize,
		   validmu = validmu,
		   valideta = stats$valideta,
                   simulate = simfun),
	      class = "family")
}


