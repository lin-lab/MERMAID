
# Download raw data 

setwd("data")

# get CDC seroprevalence survey data
sero_url <- system('curl -k https://healthdata.gov/dataset/covid-19-diagnostic-laboratory-testing-pcr-testing-time-series | tr " " "\\n" | grep "https://healthdata.gov/sites/default/files/covid-19_diagnostic_lab_testing"', intern = TRUE)
sero_url <- gsub('"', '', sero_url)
sero_url <- gsub('.*=', '', sero_url)

system(paste("wget --no-check-certificate", sero_url))

# get CDC state-level positive count data
system("curl -k https://data.cdc.gov/api/views/9mfq-cb36/rows.csv?accessType=DOWNLOAD > cdc_positive_by_state.csv")

# get USAFacts data
system("rm -f covid_confirmed_usafacts.csv")
system("wget https://static.usafacts.org/public/data/covid-19/covid_confirmed_usafacts.csv")

# get COVIDTracking data
system("rm -f all-states-history.csv")
system("wget https://covidtracking.com/data/download/all-states-history.csv")

## get CDC HealthData.gov COVID-19 PCR time series data
system("curl https://healthdata.gov/api/views/j8mb-icvb/rows.csv?accessType=DOWNLOAD > data/covid-19_diagnostic_lab_testing.csv")

## get CDC seroprevalence survey data
system("curl https://data.cdc.gov/api/views/d2tw-32xv/rows.csv?accessType=DOWNLOAD  > data/CDC_seroprevalence.csv")



# ------------------------------------------

library(data.table)
library(ggplot2)

du <- fread("data/covid_confirmed_usafacts.csv")

du <- melt(du, id.vars = names(du)[1:4], variable.name = "date", value.name = "cases")

du[,date := as.Date(date, format = "%Y-%m-%d"),]

setorder(du, State, countyFIPS, date)

du[,new_c := pmax(0,cases - shift(cases, 1, fill = 0)),by=list(State,countyFIPS)]

du <- du[,list(incid_US = sum(new_c)),by=list(state=State,date)][order(state,date)]

# ------------------------------------------

# https://covidtracking.com/about-data/data-definitions

dt_ct <- fread("data/all-states-history.csv")

dt_ct[,date := as.Date(date),]

getIncrements <- function(x, fill = 0, nneg = TRUE){
	wna <- which(is.na(x))
	x[wna] <- fill
	if(nneg) x <- pmax(0, x)
	out <- x - data.table::shift(x, 1, fill = fill)
	if(nneg) out <- pmax(0, out)
	out
}

dt_ct <- dt_ct[order(state,date)]

plot(subset(dt_ct, state=="ME")$totalTestsPeopleViralIncrease)

dt_ct_ss <- dt_ct[,list(
	date=date,
	pos_viral_CT = getIncrements(positiveTestsViral), 
	neg_viral_CT = getIncrements(negativeTestsViral),
	pos_1 = positiveIncrease,
	neg_1 = negativeIncrease
),by=list(state)]

merge_states <- c("CT",  "NJ", "TX", "IA", "IN", "MN", "MO", "SD", "WY")

dt_ct <- dt_ct[,list(
	state, date,
	incid_CT = (positiveIncrease), 
	npos_CT = getIncrements(positiveTestsViral),
	ntest_CT = (totalTestResultsIncrease),
	ntest_indiv_CT = (totalTestsPeopleViralIncrease)
),]

subset(dt_ct, npos_CT>ntest_CT)

states_no_CT_pos <- subset(dt_ct[,max(npos_CT),by=state], V1 == 0)$state

setorder(dt_ct, state, date)

# ------------------------------------------

dt <- fread("data/covid-19_diagnostic_lab_testing.csv")
dt[,date := as.Date(as.character(date), format = "%Y/%m/%d"),]

dt <- dcast(dt, state + date ~ tolower(overall_outcome), value.var = "new_results_reported")
dt[is.na(dt)] <- 0

dt[,ntest_spec_CDC := negative + positive,]

setnames(dt, c("inconclusive", "negative", "positive"), paste0( c("inconclusive", "negative", "positive"), "_CDC"))

d_sero <- fread("data/CDC_seroprevalence.csv")

sero_cn <- c(
	"Site" = "state", 
	"Date Range of Specimen Collection" = "date_range", 
	"Catchment population" = "pop_size", 
	"Rate (%) [Cumulative Prevalence]" = "pct_pos",
	"Lower CI [Cumulative Prevalence]" = "pct_pos_lb",
	"Upper CI [Cumulative Prevalence]" = "pct_pos_ub",
	"n [Cumulative Prevalence]" = "sample_size"
)

setnames(d_sero, 
	names(sero_cn), sero_cn
)

d_sero <- d_sero[,list(
	state,
	date_s  = as.Date(paste0(trimws(substr(date_range,1,6)),",2020"),format="%b %d, %y"),
	date_e  = as.Date(trimws(substr(date_range,9,nchar(date_range))),format="%b %d, %y"),
	pop_size,
	sample_size,
	pos_count = (pct_pos/100) * sample_size,
	prev_est = (pct_pos/100),
	prev_lb = (pct_pos_lb/100),
	prev_ub = (pct_pos_ub/100)
),]

d_sero[,date := date_s + round((date_e - date_s)/2),]

dc <- fread("data/cdc_positive_by_state.csv")

dc[,date := as.Date(as.character(submission_date), format = "%m/%d/%Y"),]

setorder(dc, state, date)
dc[is.na(dc)] <- 0

dc[,incid_t_CDC := pmax(0, tot_cases - shift(tot_cases, 1, fill = 0)),by=state]
dc[,incid_p_CDC := pmax(0, prob_cases - shift(prob_cases, 1, fill = 0)),by=state]
dc[,incid_c_CDC := pmax(0,incid_t_CDC - incid_p_CDC),by=state]

dc[,deaths_t_CDC := pmax(0, tot_death - shift(tot_death, 1, fill = 0)),by=state]
dc[,deaths_p_CDC := pmax(0, prob_death - shift(prob_death, 1, fill = 0)),by=state]
dc[,deaths_c_CDC :=  pmax(0,deaths_t_CDC - deaths_p_CDC),by=state]

dm <- merge(
	dc[,list(state, date, incid_c_CDC, incid_t_CDC, deaths_c_CDC, deaths_t_CDC),],
	merge(dt, du, by = c("state", "date"), all = TRUE),
	by = c("state", "date"), all = TRUE
)[order(state,date)]

dm <- merge(dm, dt_ct, by = c("state", "date"), all = TRUE)

dm <- subset(dm, date < "2021-01-01" & date >= "2020-02-01")

dm <- merge(dm, d_sero[,list(pop_size=tail(pop_size,1)),by=state], by = "state", all = FALSE)

dm <- merge(dm, dt_ct_ss, by = c("state", "date"))

M_smooth <- 3
K_smooth <- 2

cbn_series <- function(p1, p2, n1, n2){
	p1 <- pmax(0, p1, na.rm=TRUE)
	p2 <- pmax(0, p2, na.rm=TRUE)
	n1 <- pmax(0, n1, na.rm=TRUE)
	n2 <- pmax(0, n2, na.rm=TRUE)
	
	s1 <- kza::kz(p1 + n1, m = M_smooth, k = K_smooth)
	s2 <- kza::kz(p2 + n2, m = M_smooth, k = K_smooth)

	w1 <- mean(s1)/(max(s1)+1)
	w2 <- mean(s2)/(max(s2)+1)
	list(
		(p1 * w1 + p2 * w2)/(w1 + w2),
		(n1 * w1 + n2 * w2)/(w1 + w2)
	)
}

dm <- dm[order(state,date)]

dm[,
	c("npos_2S", "nneg_2S") := cbn_series(pos_viral_CT, positive_CDC, neg_viral_CT, negative_CDC)
,]


diffMat <- function(X){
	n <- ncol(X)
	out <- matrix(NA, n, n)
	X <- as.matrix(X)
	for(i in 1:n){
		for(j in 1:n){
			out[i,j] <- mean( abs(X[,i] - X[,j]) )
		}
	}
	out
}

merge_counts <- function(..., name = ""){
	Xdt <- data.table(...)
	X <- do.call(cbind, list(...))
	X[is.na(X)] <- 0
	colnames(X) <- names(Xdt)
	
	X <- apply(X, 2, function(x)  cumsum(kza::kz(x, m = M_smooth, k = K_smooth)))
	
	R <- 1/(1 + diffMat(X)) 
	diag(R) <- 0
	R[is.na(R)] <- 0
	max_R <- apply(R, 1, max, na.rm = TRUE)
	use_vars <- which(max_R > min(max_R, na.rm =TRUE))
	if( length(use_vars) >0 ){
		print(paste0(name, ": ", paste( colnames(X)[use_vars], collapse = ', ' )))
		X <- X[,use_vars]
	}
	
	xv <- apply(X, 1, mean, na.rm = TRUE)
	xv - data.table::shift(xv, n = 1, fill = 0)
}

setorder(dm, state, date)
dm[,incid_M := merge_counts(incid_t_CDC, incid_CT, incid_US, name = state),by=state]


# --------------------------------------------------------

dm[,npos_S2S := kza::kz(pmax(npos_2S, 0, na.rm=TRUE), m = M_smooth, k = K_smooth),by=state]

dm[,nneg_S2S := kza::kz(pmax(nneg_2S, 0, na.rm=TRUE), m = M_smooth, k = K_smooth),by=state]

dm[,ntest_S2S := npos_S2S + nneg_S2S,by=state]

dm[,ppos_S2S := npos_S2S/(1.00 + ntest_S2S),by=state]

dm[,incid_M_1S := 1.00,by=state]
dm$incid_M_1S[dm$incid_M <= 0 & dm$date <= "2020-03-25"] <- 0.00

dm <- dm[order(state,date)]

dm[,incid_S2S := pmax(0, predict(lm( incid_M ~ 0 + splines::bs(date, degree = 3, df = 3 + 12):npos_S2S:incid_M_1S ), type = "response")),by=state]


# --------------------------------------------------------
#  Plot merged case counts
# --------------------------------------------------------

tmp <- data.table(
"state" = state.abb,
"state_division"= state.division,
"stateName"= state.name,
"state_region" = state.region
)

dd <- melt(data=dm, id.vars = c("state", "date"), measure.vars = c("incid_c_CDC", "incid_t_CDC", "positive_CDC", "incid_US", "incid_CT", "incid_M2"))


dd <- merge(dd, tmp, by = "state")

labs <- c("incid_c_CDC" = "CDC confirmed", "incid_t_CDC" = "CDC probable", "positive_CDC" = "CDC positive tests", "incid_US" = "USAFacts cases", "incid_CT" = "COVIDTracking cases", "incid_M2" = "Merged")

dd[,type := labs[variable],]

dd[,med := median(value,na.rm = TRUE),by=list(state,date)]

dd[,max_med := max(med),by=list(state)]

dd$val2 <- dd$value
dd$val2[dd$val2 > 4*dd$max_med] <- NA
dd$val2[dd$val2 < 0 ] <- NA

require(ggplot2)

pdf("division_counts_v5.pdf", 8.5, 11)

plot_list <- list()

divi <- sort(unique(dd$state_division))
for(x in divi){

tmp <- data.table(val = c("CDC positive tests", "Merged series"), val2 = NA, date = head(dd$date,2), type = c("CDC confirmed"), stateName = subset(dd, state_division == x)$stateName[1])

print(
plot_list[[x]] <- ggplot(subset(dd, state_division == x & type != "Merged" & type != "CDC positive tests"), aes(x = date, y = val2, colour = type)) + geom_line(size = 0.2) + theme_bw() + xlab(NULL) + ylab("Daily cases") + guides(colour = guide_legend(title=NULL)) + ggtitle(x) + 
geom_line(data = subset(dd, state_division == x & type == "Merged"),  aes(x = date, y = val2), colour = "black", size = 0.2, linetype = 'dashed') + 
geom_line(data = subset(dd, state_division == x & type == "CDC positive tests"),  aes(x = date, y = val2), colour = "red", size = 0.2, linetype = 'dotted') +
geom_line(data = tmp, aes(y = as.numeric(val2), x = date, linetype = val)) +  
facet_wrap(~stateName, scales = 'free', ncol = 1) + scale_linetype_manual(name=NULL,values= c("CDC positive tests" = "dotted", "Merged series" = "dashed"), guide = guide_legend(drop=FALSE), drop = FALSE) + theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.1))
)
}


dev.off()


# --------------------------------------------------------

tmp <- data.table(
"state" = state.abb,
"state_division"= state.division,
"stateName"= state.name,
"state_region" = state.region
)

dd <- melt(data=dm, id.vars = c("state", "date"), measure.vars = c("incid_c_CDC", "incid_t_CDC", "positive_CDC", "incid_US", "incid_CT",  "incid_M", "incid_S2S"))

dd <- merge(dd, tmp, by = "state")

labs <- c("incid_c_CDC" = "CDC confirmed", "incid_t_CDC" = "CDC probable", "positive_CDC" = "CELR positive tests", "incid_US" = "USAFacts cases", "incid_CT" = "COVIDTracking cases", "incid_M" = "Consensus cases", "incid_S2S" = "Final cases (scaled CELR)" )

dd[,type := labs[variable],]

dd <- dd[order(state,type,date)]

dd[,value := (function(x,v){
	if(v[1] == "Merged" ){x}else{ kza::kz(x, m = M_smooth, k = K_smooth) }
})(value,type),by=list(state,type)]

dd[,med := median(value,na.rm = TRUE),by=list(state,date)]

dd[,max_med := max(med),by=list(state)]

dd$val2 <- dd$value
dd$val2[dd$val2 > 4*dd$max_med] <- NA
dd$val2[dd$val2 < 0 ] <- NA

require(ggplot2)

plot_list <- list()

consensus_names <- names(c("CELR positive tests" = "dotted", "Consensus cases" = "dashed", "Final cases (scaled CELR)" = "dotdash"))

x_brks <- as.Date(c(paste0("2020-0",2:9,"-01"), paste0("2020-",10:12,"-01"), "2021-01-01"))
x_labs <- c("Feb", "", "Apr", "", "Jun", "", "Aug", "", "Oct", "", "Dec", "")


# --------------------------------------------------------
# CONFIRMED CASE PLOTS
# --------------------------------------------------------

pdf("division_counts_smooth_final.pdf", 8.5, 11)

divi <- sort(unique(dd$state_division))
for(x in divi){

	tmp <- data.table(val = c("CELR positive tests", "Consensus cases", "Final cases (scaled CELR)"), val2 = NA, date = head(dd$date,3), type = c("CDC confirmed"), stateName = subset(dd, state_division == x)$stateName[1])

	print(
		plot_list[[x]] <- (
			ggplot(subset(dd, state_division == x & !(type %in% consensus_names)), aes(x = date, y = val2, colour = type)) + geom_line(size = 0.2) + theme_bw() + xlab(NULL) + ylab("Daily cases") + guides(colour = guide_legend(title=NULL)) + ggtitle(x) + 
			scale_x_date(breaks = x_brks, labels = x_labs, limits = c(min(brks), max(brks))) + 
			facet_wrap(~stateName, scales = 'free', ncol = 1) + theme(panel.grid.minor = element_line(size = 0.1), panel.grid.major = element_line(size = 0.1))  + 
			geom_line(data = subset(dd, state_division == x & type == "Consensus cases"),  aes(x = date, y = val2), colour = "black", size = 0.2, linetype = 'dashed') + 
			geom_line(data = subset(dd, state_division == x & type == "CELR positive tests"),  aes(x = date, y = val2), colour = "black", size = 0.2, linetype = 'dotted') +
			geom_line(data = subset(dd, state_division == x & type == "Final cases (scaled CELR)"),  aes(x = date, y = val2), colour = "black", size = 0.2, linetype = 'solid') + 
			geom_line(data = tmp, aes(y = as.numeric(val2), x = date, linetype = val)) + 
			scale_linetype_manual(name=NULL,values= c("CELR positive tests" = "dotted", "Consensus cases" = "dashed", "Final cases (scaled CELR)" = "solid"), guide = guide_legend(drop=FALSE), drop = FALSE) + 
			theme(strip.background = element_rect(colour = NA), strip.text.x = element_text(margin = margin(0,0,0,0, "lines")))
		)
	)
}

dev.off()

