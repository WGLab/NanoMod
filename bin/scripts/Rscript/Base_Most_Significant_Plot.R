mscaley <- function(my) sprintf("%6.2f", my)

Base_Most_Significant_Plot <- function(plotDat, stru, strt, strks, strcb, mrtitle, mhasbox) {
   dftmargin = 0.1

   volinplot <- ggplot(plotDat) +
                facet_grid(~factor(Position)) +
		scale_fill_manual(values=c('1'='blue', '2'='green'))+
		ylab('Mean raw signal') + theme_bw() +
		theme(panel.margin=unit(0, "lines")) +
		geom_violin(aes(x=DS, y=Signal, fill=DS)) +
		#geom_violin(aes(x=factor(DS), y=Signal, fill=factor(DS))) +
		theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) +
		theme(plot.margin = unit(c(0.4+dftmargin,dftmargin,dftmargin,dftmargin+0.2),"cm")) + 
		labs(title = mrtitle) + 
		theme(legend.margin=margin(t=0, r=0, b=0, l=-0.2, unit="cm")) +
		#scale_y_continuous(labels=mscaley, oob = rescale_none) + 
		theme(plot.title = element_text(lineheight=.7, face="bold")) 

   mdflist <- list(strks, stru);
   if (nrow(strcb)>0){ 
      mdflist <- list(strks, strcb);
   }
   mmin_pv <- c(0, 0);
   for (dfi in 1:length(mdflist)){
      if (mmin_pv[dfi] > min(mdflist[[dfi]]$Pvalue)){
         mmin_pv[dfi] <- min(mdflist[[dfi]]$Pvalue)
      }
   }

   if (nrow(strcb)>0){
      if (mhasbox==1){ yl4 <- c("log(KS P-value)", "log(Comb P-value)") }
      else {yl4 <- c("log(KS Pv)", "log(Comb Pv)") }
   }else{
      if (mhasbox==1){ yl4 <- c("log(KS P-value)") }#, "log(MWU P-value)") }
      else {yl4 <- c("log(KS Pv)") }#, "log(MWU Pv)") }
   }
   mcolor <- c("cyan", "darkgoldenrod2")

   pvs <- lapply(1:length(mdflist), function(pt){
      cpst <- mdflist[[pt]]

      mp <- ggplot(cpst, aes(x=factor(Position), y=Pvalue)) + 
            ylab(yl4[pt]) + theme_bw() + theme(axis.title.y = element_text(size = rel(0.9))) +
	    theme(panel.margin=unit(0, "lines")) +
	    theme(axis.ticks.x=element_blank(),axis.title.x=element_blank()) +
	    theme(plot.margin = unit(c(0.1,1.0,0.1,0.4), "cm")) +
	    theme(axis.title.y=element_text(size=7)) +
	    geom_bar(stat="identity", width=0.33, fill=mcolor[pt], color=mcolor[pt]) +
	    scale_y_continuous(oob = rescale_none, limits=c(mmin_pv[pt],0)) +
	    #scale_y_continuous(labels=mscaley, oob = rescale_none) + 
	    theme(axis.text.x=element_blank())


      return(mp)
   })
   
   if (mhasbox==1){
      boxplot <- ggplot(plotDat) +
                 geom_boxplot(aes(x=DS, y=Signal, fill=DS), outlier.colour="red",outlier.shape=16, outlier.size=0.3, notch=TRUE) +
                 #geom_boxplot(aes(x=factor(DS), y=Signal, fill=factor(DS)), outlier.colour="red",outlier.shape=16, outlier.size=0.3, notch=TRUE) +
                 facet_grid(~factor(Position)) +
		 scale_fill_manual(values=c('1'='blue', '2'='green'))+
		 ylab('Mean Raw signal') + theme_bw() +
		 theme(axis.title.y=element_text(size=9)) +
		 theme(panel.margin=unit(0, "lines")) +
		 theme(legend.margin=margin(t=0, r=0, b=0, l=-0.2, unit="cm")) +
		 theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) +
		 theme(plot.margin = unit(c(dftmargin,dftmargin,dftmargin,dftmargin+0.2),"cm")) +
		 #scale_y_continuous(labels=mscaley, oob = rescale_none) +
		 theme(strip.background = element_blank(),strip.text.x = element_blank())
      if (nrow(strcb)>0){
         grid.arrange(volinplot, boxplot, pvs[[1]], pvs[[2]], ncol = 1, heights=c(2.5,1.25,0.75,0.75));
      }else{
         grid.arrange(volinplot, boxplot, pvs[[1]], ncol = 1, heights=c(2.5,1.25,0.75));
      }
   }
   else{
     if (nrow(strcb)>0){
        grid.arrange(volinplot, pvs[[1]], pvs[[2]], ncol = 1, heights=c(2,0.7,0.7));
     }else{
        grid.arrange(volinplot, pvs[[1]], ncol = 1, heights=c(2,1));
     }
   }

}

