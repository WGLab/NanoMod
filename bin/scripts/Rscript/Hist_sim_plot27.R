Hist_sim_plot27 <- function(mdfperc, spvector, rankstrvector, xlabs){
   methods <- unique(mdfperc$Method)
   
   casecontrolpairs <- unique(mdfperc$CaseControl)

   #ps <- lapply(1:length(casecontrolpairs), function(ccp_ind){
   lapply(1:length(casecontrolpairs), function(ccp_ind){
      ccp <- casecontrolpairs[ccp_ind]
      cplotDat <- mdfperc[mdfperc$CaseControl==ccp,]

      mp <- ggplot(cplotDat, aes(x=MixedPerc, y=Fraction, fill=Percentile)) +
            theme_bw() +
	    geom_bar(stat = 'identity', position = 'stack') +
	    scale_fill_manual(values=c("green", "green3", "blue","blue4","grey","yellow","red","magenta")) +
            facet_grid(~factor(Method)) +
	    xlab(xlabs) +
	    #ggtitle(ccp) +
            theme(panel.margin=unit(0, "lines")) +
	    theme(axis.title.x=element_blank())
	    #theme(axis.text.x=element_blank())

      #return (mp)
      grid.arrange(mp, ncol=1)
   })
   
   #grid.arrange(ps, ncol=1);
}


