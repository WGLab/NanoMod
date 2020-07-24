Base_Most_Significant_Plot <- function(plotDat, stru, strt, strks, strcb, mrtitle, mhasbox, mplotType) {
dftmargin = 0.1
  print( mplotType)
  print( mplotType=="Density")
  if(mplotType=="Violin"){

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
  } else{
  print("In density Plot") #Density Plot

    plotDat$Base <- str_extract(as.character(plotDat$Position),"[ATCG]" )

    plotDat$Location <- strsplit(as.character(plotDat$Position), "/")
    plotDat$Location <- do.call("rbind", plotDat$Location)
    plotDat$Location <- plotDat$Location[,1] %>% as.numeric()

    kspv <- strsplit(as.character(plotDat$Position), "\n")
    kspv <- do.call("rbind", kspv)
    kspv <- kspv[,2]
    kspv <- as.numeric(kspv)
    cpv <-  strsplit(as.character(plotDat$Position), "\n")
    cpv <- do.call("rbind", cpv)
    cpv <- cpv[,3]
    cpv <- as.numeric(cpv)
    plotDat$kspv <- log10(kspv)
    plotDat$cpv <- log10(cpv)
    pvpltDat <- plotDat[!duplicated(plotDat[,'Location']),]

    densDat <- lapply(split(
      plotDat, paste0(plotDat$Position,
                      plotDat$DS)),
      function(pgDat){
        pgDens <- density(pgDat$Signal)
        nDens <- length(pgDens$x)
        data.frame(Position=pgDens$y,
                   Signal=pgDens$x,
                   gPos=rep(pgDat$Location[1], nDens),
                   Group=rep(pgDat$DS[1], nDens))
      })

    maxDens <- max(unlist(lapply(densDat, function(x) x$Position)))
    normDensDat <- do.call(
      rbind.data.frame,
      lapply(densDat, function(posDens){
        pos <- ifelse(posDens$Group == 1, posDens$Position * -1,
                      posDens$Position) / (maxDens * 2)
        data.frame(Position=pos + posDens$gPos[1],
                   Signal=posDens$Signal,
                   gPos=posDens$gPos,
                   Group=posDens$Group)
      }))

    dftmargin = 0.1
    base_pos <- min(plotDat$Signal)-1.5
    p <- ggplot(normDensDat) +
      geom_polygon(aes(x=Position + 0.5, y=Signal, fill=as.factor(Group),
                       group=paste0(Group, gPos)),
                   size=0, show.legend=FALSE) +
      geom_text(
        aes(x=Location+0.5, y=base_pos, label=Base, color=Base),
        data=plotDat, hjust=0.5, size=2.5, show.legend=FALSE,
      ) +
      labs(title = mrtitle) +
      scale_fill_manual(values=c('1'='#e34a33', '2'='#404040')) +
      geom_vline(
        xintercept=
          min(plotDat$Location):
          (max(plotDat$Location) + 1),
        size=0.01) +theme_bw() + theme(
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin=unit(c(0.1,1.0,0.1,0.4), "cm"),
          axis.title.y=element_text(size=7)) +
      theme(plot.title = element_text(lineheight=.7, face="bold"))


    ks_p <- ggplot(pvpltDat) + geom_bar(stat="identity",fill = "#ef8a62",aes(x=Location + 0.5, y=kspv)) +
      geom_vline(
        xintercept=
          min(plotDat$Location):
          (max(plotDat$Location) + 1),
        size=0.01) +theme_bw() + theme(panel.grid.major.x=element_blank(),
                                       panel.grid.minor.x=element_blank(),
                                       panel.grid.minor.y=element_blank(),
                                       axis.title.x=element_blank(),axis.text.x=element_blank(),
                                       axis.ticks.x=element_blank(),
                                       plot.margin=unit(c(0.1,1.0,0.1,0.4), "cm"),
                                       axis.title.y=element_text(size=7)) +
      ylab("log10(KS Pv)")

    cs_p <- ggplot(pvpltDat) + geom_bar(stat="identity",fill = "#91bfdb",aes(x=Location + 0.5, y=cpv)) +
      geom_vline(
        xintercept=
          min(plotDat$Location):
          (max(plotDat$Location) + 1),
        size=0.01) +theme_bw() + theme(axis.text.x=element_text(hjust=0),
                                       panel.grid.major.x=element_blank(),
                                       panel.grid.minor.x=element_blank(),
                                       panel.grid.minor.y=element_blank(),
                                       plot.margin=unit(c(0.1,1.0,0.1,0.4), "cm"),
                                       axis.title.y=element_text(size=7)) +
      ylab("log10(Comb Pv)") + xlab("Position")
    print(head(normDensDat))
    grid.arrange(p, ks_p, cs_p, ncol = 1)

  }
}
