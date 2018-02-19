Hist_sim_plot9 <- function(mdfperc, spvector, rankstrvector, xlabs){
   methods <- unique(mdfperc$Method)

   lapply(1:length(methods), function(m_ind){
      curm <- methods[m_ind]
      cplotDat <- mdfperc[mdfperc$Method==curm,]

      mp <- ggplot(cplotDat, aes(x=MixedPerc, y=Fraction, fill=Percentile)) +
            theme_bw() +
	    geom_bar(stat = 'identity', position = 'stack') +
	    scale_fill_manual(values=c("green", "green3", "blue","blue4","grey","yellow","red","magenta")) +
	    facet_wrap(~factor(CaseControl), ncol = 3) +
	    xlab(xlabs) +
	    #ggtitle(ccp) +
            theme(panel.margin=unit(0, "lines")) +
	    theme(axis.title.x=element_blank())
	    #theme(axis.text.x=element_blank())

      #return (mp)
      grid.arrange(mp, ncol=1)
   })
}


