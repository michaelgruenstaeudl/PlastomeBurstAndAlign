#!/usr/bin/env Rscript

library(plyr) # For rounding up
library(ggplot2)
library(svglite)

inFn = file.choose()
ouDr = dirname(inFn)
plotData = read.csv(inFn, header=T)

max_time = round_any(max(plotData$TIME), 10, f=ceiling)

base_plot = ggplot(data=plotData, aes(x=CPUS, y=TIME, group=MARKER))

myPlot = base_plot + 
    geom_point(data=plotData, aes(color=DATASET)) +
    geom_line(data=plotData, aes(color=DATASET)) +
    geom_text(data=plotData, aes(label=paste(TIME, " sec"), color=DATASET), size=4, vjust=-0.35, hjust=-0.15) +
    scale_x_continuous(breaks=c(1,5,10), 
                       minor_breaks=c(2,3,4,6,7,8,9,11,12), 
                       labels=c(1,5,10),
                       expand=expansion(0.01),
                       limits = c(1,12)) +
    scale_y_log10(breaks=c(100, seq(500, max_time, 500)),
                  minor_breaks=NULL,
                  labels=c(100, seq(500, max_time, 500)),
                  limits = c(10,max_time)) +
    #scale_linetype_manual(values = c("dotted", "solid")) +
    scale_color_manual(labels = c("Asteraceae (n=155)", "monocots (n=733)"), values = c("blue", "maroon")) +
    xlab("\nNumber of CPUs") + 
    ylab("Computation Time (seconds) - log-scale\n") +
    ggtitle("Speed improvements through MPI implementation\n") +
    theme_minimal() + 
    theme(plot.title = element_text(size=18, hjust=0.5),
        plot.subtitle = element_text(size=16, face="italic"),
        axis.text=element_text(size=12, color="gray"),
        axis.title=element_text(size=14),
        #plot.margin=unit(c(0.1,1.0,0.1,0.1),"cm"),  # Note: margin(t=0, r=0, b=0, l=0, unit="pt")
        legend.key.width=unit(1,"cm"),
        panel.grid.major=element_line(color="gray"), 
        panel.grid.minor=element_line(color="gray")
  )


out_fn=paste(ouDr, "/PlastomeBurst_SpeedImprovements.svg", sep='')
svglite(out_fn, width=12, height=8) #height=29.7)
myPlot
dev.off()

