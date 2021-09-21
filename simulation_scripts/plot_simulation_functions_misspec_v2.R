
require(data.table)
require(ggplot2)
require(gridExtra)
require(cowplot)

make_plot_abc <- function(data,  xlab_list, ylab_list, scales = "free_x", ylim = as.numeric(c(NA,NA)) ){

	require(data.table)
	require(ggplot2)
	require(gridExtra)
	require(cowplot)

	CONST <- 10
	COLFUN <- function(x) NA

	brk <- as.Date(c("2020-01-01", "2020-02-01", "2020-03-01", "2020-04-01", "2020-05-01","2020-06-01", "2020-07-01", "2020-08-01","2020-09-01", "2020-10-01","2020-11-01", "2020-12-01", "2021-01-01"))
	lab <- c("Jan", "", "Mar", "",  "May", "",  "Jul",  "", "Sep",  "", "Nov", "", "")
	
	tmp <- subset(data, lag_p_fit_size == 1)[,list(
	date = rep(as.Date(gsub("^.*:", "", variable)), 2),
	"Truth*" = Truth, 
	"Truth" = Truth, 
	"MERMAID\nTrue lag r = 1" = Mean,
	"MERMAID\nTrue lag r = 999" = subset(data, lag_p_fit_size != 1)$Mean ,
	SE, StdDev, sim_name =  si_config #gsub("SI Mean", "SI\nMu", si_config)
	),]

	tmp <- melt(tmp, 
		measure.vars = c("Truth*", "MERMAID\nTrue lag r = 1", "MERMAID\nTrue lag r = 999"),
		id.vars = c("date", "sim_name", "SE", "StdDev", "Truth")
	)

	require(ggplot2)
	p1 <- ggplot(
	tmp, 
	aes(x = date, y = value, ymin = value - CONST*SE, ymax = value + CONST*SE, colour = variable, linetype = variable)
	) + guides( colour = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) + 
	  scale_colour_manual(values = c("green", "red", "blue")) + 
	  geom_line(size = 0.35) + 
	  theme_bw() +
	  theme(axis.line = element_line(colour = "black"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank()
		) + theme(strip.background = element_rect(colour = NA), strip.text.x = element_text(margin = margin(0,0,0,0, "lines"))) + 
	  ylab(ylab_list[[1]]) + xlab(xlab_list[[1]]) + facet_wrap(~sim_name,scales=scales) + scale_x_continuous(labels = lab, breaks = brk) + ylim(ylim[1], ylim[2])


	tmp_SE <- rbind(
		subset(tmp, variable == "MERMAID\nTrue lag r = 1")[,list(sim_name,date, "value" = StdDev, type = "Truth /\nEmpirical SE", "variable" = "True lag r = 1"),],
		subset(tmp, variable == "MERMAID\nTrue lag r = 999")[,list(sim_name,date, "value" = StdDev, type = "Truth /\nEmpirical SE","variable" = "True lag r = 999" ),],
		subset(tmp, variable == "MERMAID\nTrue lag r = 1")[,list(sim_name,date, "value" = SE, type = "Mean est. SE", "variable" = "True lag r = 1"),],
		subset(tmp, variable == "MERMAID\nTrue lag r = 999")[,list(sim_name,date, "value" = SE, type = "Mean est. SE","variable" = "True lag r = 999" ),]
	)

	p2_SE <-ggplot(
		tmp_SE, 
		aes(x = date, y = value, colour = variable, linetype = type)
	) + guides( colour = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) + 
	  scale_colour_manual(values = c("red", "blue")) + 
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




