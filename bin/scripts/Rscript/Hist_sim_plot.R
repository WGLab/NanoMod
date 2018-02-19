Hist_sim_plot <- function(mdfperc, spvector, rankstrvector){
   histplot <- ggplot(mdfperc, aes(x=MixedPerc, y=Fraction, fill=Percentile)) + 
              theme_bw() +
	      geom_bar(stat = 'identity', position = 'stack') +
	      scale_fill_manual(values=c("green", "green3", "blue","blue4","grey","yellow","red","magenta")) + 
	      theme(axis.text.x = element_text(angle = 90, hjust = 1))

   grid.arrange(histplot)
}


