
require(ggplot2)
require(gridExtra)


make_plot_abc <- function(data, param_list,  xlab_list, ylab_list, scales = "free_x", ylim = as.numeric(c(NA, NA)), MERMAID_only = FALSE ){

	require(data.table)
	require(ggplot2)
	require(gridExtra)
	require(cowplot)

	CONST <- 10
	COLFUN <- function(x) NA

	brk <- as.Date(c("2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01", "2020-05-01","2020-06-01", "2020-07-01", "2020-08-01","2020-09-01", "2020-10-01","2020-11-01", "2020-12-01", "2021-01-01"))
	lab <- c("Jan", "", "Mar", "",  "May", "",  "Jul",  "", "Sep",  "", "Nov", "", "")
	
	data$date <- as.Date(gsub("^.*:", "", data$variable))

	data <- data[order(data$date),]

	tmp <- subset(data, param == param_list[1])[,list(
	date = rep(as.Date(gsub("^.*:", "", variable)), 2),
	"Truth*" = Truth, 
	"Truth" = Truth, 
	"MERMAID" = Mean, 
	"Reported counts" = subset(data, param == param_list[2])$Mean, 
	"Reconvolution" = subset(data, param == param_list[3])$Mean, 
	SE, StdDev, sim_name
	),]
	
	tmp <- melt(tmp, 
	measure.vars = c("Truth*", "MERMAID", "Reported counts", "Reconvolution"),
	id.vars = c("date", "sim_name", "SE", "StdDev", "Truth")
	)

	if( MERMAID_only ){
		tmp$value[!(tmp$variable %in% c("Truth*", "MERMAID") )] <- NA
		tmp$SE[!(tmp$variable %in% c("Truth*", "MERMAID") )] <- NA
		tmp$StdDev[!(tmp$variable %in% c("Truth*", "MERMAID") )] <- NA
	}

	tmp$SE[tmp$variable != "MERMAID"] <- NA
	tmp$StdDev[tmp$variable != "MERMAID"] <- NA

	require(ggplot2)
	p1 <- ggplot(
	tmp, 
	aes(x = date, y = value, ymin = value - CONST*SE, ymax = value + CONST*SE, colour = variable, linetype = variable)
	) + guides( colour = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) + 
	  scale_colour_manual(values = c("green", "black", "magenta", "blue")) + 
	  geom_line(size = 0.35) + 
	  theme_bw() +
	  theme(axis.line = element_line(colour = "black"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank()
		) + 
	  ylab(ylab_list[[1]]) + xlab(xlab_list[[1]]) + facet_wrap(~sim_name,scales=scales) + scale_x_continuous(labels = lab, breaks = brk) + ylim(ylim[1], ylim[2]) + theme(strip.background = element_rect(colour = NA), strip.text.x = element_text(margin = margin(0,0,0,0, "lines")))

	tmp_SE <- rbind(
		subset(tmp, variable == "MERMAID")[,list(sim_name,date, "value" = StdDev, "variable" = "\nTruth / \nEmpirical SE"),][order(sim_name,date)],
		subset(tmp, variable == "MERMAID")[,list(sim_name,date, "value" = SE, "variable" = "\nMERMAID\nMean est. SE" ),][order(sim_name,date)]
	)
	tmp_SE$variable <- factor(tmp_SE$variable, levels = c("\nTruth / \nEmpirical SE", "\nMERMAID\nMean est. SE"), order = TRUE)

	p2_SE <-ggplot(
		tmp_SE, 
		aes(x = date, y = value, colour = variable, linetype = variable)
	) + guides( colour = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) + 
	  scale_colour_manual(values = c("green", "black")) + 
	  geom_line(size = 0.35) + 
	  theme_bw() +
	  theme(axis.line = element_line(colour = "black"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank()
		) + 
	  ylab("Standard error") + xlab(NULL) + facet_wrap(~sim_name,scales=scales) + scale_x_continuous(labels = lab, breaks = brk) + ylim(0, NA) + 
	  theme(strip.background = element_rect(colour = NA), strip.text.x = element_text(margin = margin(0,0,0,0, "lines")))

	require(gridExtra)

	cbn <- grid.arrange(p1 + ggtitle("A"), p2_SE + ggtitle("B"), nrow = 2, heights = c(1,1))

	return(cbn)
}

make_theme <- function(obj){
	obj + theme_bw() +
  theme(
  	    axis.line = element_line(colour = "black"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(),
		legend.position = 'bottom', 
		legend.box.margin=margin(-2.5,-2.5,-2.5,-2.5), 
		legend.margin=margin(0,0,0,0), 
		legend.text = element_text(size=9),
		strip.background = element_rect(color = NA, size = 1)
    ) 
}

save_plots <- function(dt_calib, out_suffix){

	pl_height <- 6
	pl_width <- 8

	# ----------------------------------
	# Rt

	xlab_list <- list(
		NULL, expression("True "~italic(R[t])), expression("Standard deviation of"~hat(italic(R[t])))
	)

	ylab_list <- list(
		expression(italic(R[t])), 
		expression("Mean"~hat(italic(R[t]))),
		expression("Mean estimated SE of "~hat(italic(R[t])))
	)

	param_list <- c("Rt", "Rt_fit_0", "Rt_fit_1")

	tmp <- subset(dt_calib, param %in% param_list)
	
	cbn <- make_plot_abc(data = tmp, param_list, xlab_list = xlab_list, ylab_list = ylab_list)

	ggsave(paste0("final_supp_simulation_plots_v2/rt_simulation_results",out_suffix,".pdf"), height = pl_height, width = pl_width, cbn)

	ggsave(paste0("final_supp_simulation_plots_v2/rt_simulation_results",out_suffix,".png"), height = pl_height, width = pl_width, cbn)


	# ----------------------------------
	# Pi_it

	xlab_list <- list(
		NULL, expression("True "~italic(pi[t])), expression("Standard deviation of"~hat(italic(pi[t])))
	)

	ylab_list <- list(
		expression(italic(pi[t])), 
		expression("Mean"~hat(italic(pi[t]))),
		expression("Mean estimated SE of "~hat(italic(pi[t])))
	)

	param_list <- c("pi_t", "pi_t", "pi_t")

	tmp <- subset(dt_calib, param %in% param_list)
	
	# Adjustment factor
	tmp$SE <- sqrt(25)*tmp$SE

	cbn <- make_plot_abc(data = tmp, param_list, xlab_list = xlab_list, ylab_list = ylab_list, MERMAID_only = TRUE)

	ggsave(paste0("final_supp_simulation_plots_v2/pi_simulation_results",out_suffix,".pdf"), height = pl_height, width = pl_width, cbn)

	ggsave(paste0("final_supp_simulation_plots_v2/pi_simulation_results",out_suffix,".png"), height = pl_height, width = pl_width, cbn)


	# ----------------------------------
	# Incidence

	Const <- 100/(8000000)

	xlab_list <- list(
		NULL, 
		expression("True daily incidence (%)"), 
		expression("Standard deviation of estimates")
	)

	ylab_list <- list(
		expression("Daily incidence (%)"), 
		expression("Mean of estimates"),
		expression("Mean estimated SE")
	)

	param_list <- c("Infections", "Infections_0", "Infections_1")

	tmp <- subset(dt_calib, param %in% param_list)
	
	kp_var <- subset(tmp[,unique(variable),by=param][,.N,by=V1], N == length(param_list))$V1
	tmp <- subset(tmp, variable %in% kp_var)
	
	tmp[,unique(variable),by=param]
	# Scale as percentage
	tmp$SE <- Const*tmp$SE
	tmp$Mean <- Const*tmp$Mean
	tmp$Truth <- Const*tmp$Truth
	tmp$Bias <- Const*tmp$Bias
	tmp$StdDev <- Const*tmp$StdDev

	cbn <- make_plot_abc(data = tmp, param_list=param_list, xlab_list = xlab_list, ylab_list = ylab_list, scales='free', ylim=c(0,NA))

	ggsave(paste0("final_supp_simulation_plots_v2/incid_simulation_results",out_suffix,".pdf"), height = pl_height, width = pl_width, cbn)

	ggsave(paste0("final_supp_simulation_plots_v2/incid_simulation_results",out_suffix,".png"), height = pl_height, width = pl_width, cbn)

	# ----------------------------------
	# Prevalence

	Const_prev <- 100

	xlab_list <- list(
		NULL, 
		expression("True cumulative prevalence (%)"), 
		expression("Standard deviation of estimates")
	)

	ylab_list <- list(
		expression("Cumulative prevalence (%)"), 
		expression("Mean of estimates"),
		expression("Mean estimated SE")
	)

	param_list <- c("Prevalence", "Prevalence_0", "Prevalence_1")

	tmp <- subset(dt_calib, param %in% param_list)
	# Scale as percentage
	tmp$SE <- Const_prev*tmp$SE
	tmp$Mean <- Const_prev*tmp$Mean
	tmp$Truth <- Const_prev*tmp$Truth
	tmp$Bias <- Const_prev*tmp$Bias
	tmp$StdDev <- Const_prev*tmp$StdDev

	cbn <- make_plot_abc(data = tmp, param_list=param_list, xlab_list = xlab_list, ylab_list = ylab_list, scales='free', ylim=c(0,NA))

	ggsave(paste0("final_supp_simulation_plots_v2/prev_simulation_results",out_suffix,".pdf"), height = pl_height, width = pl_width, cbn)

	ggsave(paste0("final_supp_simulation_plots_v2/prev_simulation_results",out_suffix,".png"), height = pl_height, width = pl_width, cbn)



}



load_data_and_plot <- function(suffix, out_suffix){

	require(data.table)
	require(ggplot2)
	tmp_list <- list()
	for( sim_type in c("BS1", "BS2", "PC1", "PC2") ){
	tmp_list[[sim_type]] <- readRDS(paste0("simulation_output/simulation_output",suffix,"_",sim_type,".rds"))

	}

	sim_dt <- rbindlist(lapply(tmp_list, `[[`, 1))
	dt_calib <- rbindlist(lapply(tmp_list, `[[`, 2))

	save_plots(dt_calib, out_suffix)

}


