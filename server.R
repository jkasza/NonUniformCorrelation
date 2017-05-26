###################################################################
# R code to accompany "Impact of non-uniform correlation structure
#    on sample size and power in multiple-period cluster randomised
#    trials", by J Kasza et al.
# 
# This file contains the server for the Shiny app 
###################################################################


library(shiny)
library(ggplot2)
library(reshape2)
library(plyr)
library(rootSolve)
library(plotly)

source("NonUnifCorrFunctions.R", local=TRUE)

shinyServer(function(input, output) {
  
#Output design matrices and explanations:  
  
  output$xmat <- renderTable({
    if(input$design == 1) Xdes <- SWdesmat2(input$T)
    if(input$design == 2) Xdes <- plleldesmat2(input$T)
    if(input$design == 3) Xdes <- pllelbasedesmat2(input$T)
    if(input$design == 4) Xdes <- crxodesmat2(input$T)
    head(Xdes, n=input$T)
    
  },digits=0)
  
  
  output$text1 <- renderText({ 
    "Design matrix:"
  })
  
  
  output$text1 <- renderText({ 
    "Stepped wedge design matrix:"
  })
  output$text2 <- renderText({ 
    "Parallel design matrix:"
  })
  output$text3 <- renderText({ 
    "Parallel w/ baseline design matrix:"
  })
  output$text4 <- renderText({ 
    "CRXO design matrix:"
  })
  
  
  
  output$SWxmat <- renderTable({
    head(SWdesmat2(input$T), n=input$T)
  },digits=0)
  output$pllelxmat <- renderTable({
    head(plleldesmat2(input$T), n=input$T)
  },digits=0)
  output$pllelbxmat <- renderTable({
    head(pllelbasedesmat2(input$T), n=input$T)
  },digits=0)
  output$crxoxmat <- renderTable({
    head(crxodesmat2(input$T), n=input$T)
  },digits=0)
  
  
  
  output$HHrelvarexplan <- renderText({
    "Relative variances of the treatment effect estimator for each design assuming an exponential decay
    within-cluster correlation structure for a range of decay parameters, relative to the variance 
    calculated assuming the Hussey and Hughes model."
  })
  
  output$HGrelvarexplan <- renderText({
    "Relative variances of the treatment effect estimator for selected design assuming an exponential decay
    within-cluster correlation structure for a range of decay parameters, relative to the variance 
    calculated assuming the Hooper/Girling model for a range of alpha parameter values."
  })
  
  output$ExpDecayExplan <- renderText({
    if(input$select==1){
      "Variance of the treatment coefficient assuming an exponential decay 
      between-period correlation structure for a range of decay parameters, 
      for each of the stepped wedge, parallel with baseline,
      parallel and cross-over designs, given the user-specified design parameters."
    }
    
    else if(input$select==2){
      "Relative variances of the treatment coefficient for each design assuming an exponential decay
      between-period correlation structure for a range of decay parameters, relative to the variance 
      calculated assuming an exchangeable correlation structure (the Hussey and Hughes model)."
    }
    else if(input$select==3){
      "Design effects for each design, assuming an expoenential decay between-period correlation structure for 
      a range of decay parameters, given the user-specified design parameters."
    }
    else if(input$select==4){
      "Power to detect the user-specified effect size for each design, assuming an expoenential decay between-period correlation structure for 
      a range of decay parameters, given the user-specified design parameters."
    }
    })  
  
  output$Hooperincltext<- renderText({
    if(input$HooperYN==2){
      "The quantity corresponding to a model with constant between-period correlation structure
      is also indicated, for a user-specified value of the constant between-period correlation."
    }
    })  

#Create plots   
 
  # Exponential decay versus Hooper/Girling model:
   output$Contourplot <- renderPlot({
    
    if(input$design == 1) Xdes <- SWdesmat2(input$T)
    if(input$design == 2) Xdes <- plleldesmat2(input$T)
    if(input$design == 3) Xdes <- pllelbasedesmat2(input$T)
    if(input$design == 4) Xdes <- crxodesmat2(input$T)
    
    Xdes <- matrix(data=as.vector(t(Xdes)), nrow=input$nclust*nrow(Xdes), ncol=ncol(Xdes), byrow = TRUE)
    
    decay<- seq(0,1,0.0025)
    r <- 1-decay
    
    decay1vec <- as.vector(matrix(data= decay, nrow=length(r), ncol=length(r), byrow=FALSE))
    decay2vec <- as.vector(matrix(data= decay, nrow=length(r), ncol=length(r), byrow=TRUE))
    
    r1vec <- as.vector(matrix(data= r, nrow=length(r), ncol=length(r), byrow=FALSE))
    r2vec <- as.vector(matrix(data= r, nrow=length(r), ncol=length(r), byrow=TRUE))
    
    color_paletteBW <-colorRampPalette(c( "white", "black"))(8)
    expdecayvars <- sapply(r, ExpDecayVar,  Xmat=Xdes,  m=input$m, rho0=input$rho0, simplify = "array")
    hgvars <- sapply(r, HooperVar, Xmat=Xdes,  m=input$m, rho0=input$rho0, simplify = "array")
    
    #want to divide each value of expdecayvars by each value of hgvars
    relexphg <- matrix(data= expdecayvars, nrow= length(expdecayvars), ncol=length(expdecayvars), byrow = FALSE)
    relexphg<- sweep(relexphg, 2, hgvars, "/")
    relexphgvec<- as.vector(relexphg)
    relexphgdf <- data.frame(r1 = r1vec, r2= r2vec, decay1= decay1vec, decay2= decay2vec, relvar = relexphgvec)
    
    brks<-cut(relexphgdf$relvar, breaks=c(0, 0.5, 0.75, 1, 2, 3, 4, 5, 20))
    brks <- gsub(",", "-",brks,fixed=TRUE)
    relexphgdf$brks<-gsub("\\(|\\]","",brks)
    
    contourplot<- ggplot(relexphgdf[relexphgdf$decay1<=input$maxr & relexphgdf$decay2<=input$maxr,], aes(decay1,decay2)) +
      geom_raster(aes(fill=brks)) +
      geom_abline(slope=1, intercept=0) +
      scale_fill_manual("Relative var." , values = color_paletteBW) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme(legend.position ="bottom", legend.key.size = unit(1, "cm"), legend.text=element_text(size=12), legend.background = element_rect(fill="grey95") )  + 
      guides(fill=guide_legend(nrow=1, keywidth=2, unit="cm"))  +
      ylab("Hooper/Girling alpha") +  xlab("Exponential decay (1-r)") #+
  #  coord_fixed() + 

    contourplot
    
  })

  
 #Variance, power, design effect for the exponential decay model
  output$varplotlyexp <- renderPlotly({
   
    Xsw <- SWdesmat2(input$T)
    Xpllel <- plleldesmat2(input$T)
    Xpllelbase <- pllelbasedesmat2(input$T)
    Xcrxo <- crxodesmat2(input$T)
    
    Xsw <- matrix(data=as.vector(t(Xsw)), nrow=input$nclust*nrow(Xsw), ncol=ncol(Xsw), byrow = TRUE)
    Xpllel <- matrix(data=as.vector(t(Xpllel)), nrow=input$nclust*nrow(Xpllel), ncol=ncol(Xpllel), byrow = TRUE)
    Xpllelbase <- matrix(data=as.vector(t(Xpllelbase)), nrow=input$nclust*nrow(Xpllelbase), ncol=ncol(Xpllelbase), byrow = TRUE)
    Xcrxo <- matrix(data=as.vector(t(Xcrxo)), nrow=input$nclust*nrow(Xcrxo), ncol=ncol(Xcrxo), byrow = TRUE)
    
    if((input$T-1)%%2 == 1) { correction <- (input$T)/(input$T-1) }
    if((input$T-1)%%2 == 0) { correction <- 1}
    
    
    myvars<- data.frame(decay = seq(1,0,-0.01), 
                        SteppedWedge= simplify2array(lapply(seq(0, 1, 0.01), ExpDecayVar, Xmat = Xsw, m=input$m, rho0=input$rho0)),
                        Parallel = correction*simplify2array(lapply(seq(0, 1, 0.01), ExpDecayVar, Xmat = Xpllel, m=input$m, rho0=input$rho0)),
                        ParallelBaseline = correction*simplify2array(lapply(seq(0, 1, 0.01), ExpDecayVar, Xmat = Xpllelbase, m=input$m, rho0=input$rho0)),
                        CRXO = correction*simplify2array(lapply(seq(0, 1, 0.01), ExpDecayVar,Xmat = Xcrxo, m=input$m, rho0=input$rho0)))
  
      constantdecayvar <- c(HooperVar(1-input$rCons, Xsw, input$m, input$rho0),
                          correction*HooperVar(1-input$rCons, Xpllel, input$m, input$rho0),
                          correction*HooperVar(1-input$rCons, Xpllelbase, input$m, input$rho0),
                          correction*HooperVar(1-input$rCons, Xcrxo, input$m, input$rho0))
    
    if(input$select == 3) {
      #Design effects: compare to an individually randomised trial of same size
      #total variance here = 1
      # 4*sigma2/(NmT) 
      #Need different iRCT variances for the different designs
      #For SW: K = T-1
      variRCT <- 4*1/((input$T-1)*input$m*input$T)   
      myvars$SteppedWedge <-  myvars$SteppedWedge/variRCT
      myvars$Parallel <-  myvars$Parallel/variRCT
      myvars$ParallelBaseline <-  myvars$ParallelBaseline/variRCT
      myvars$CRXO <-  myvars$CRXO/variRCT
      
      #For other designs, K depends on T in a more complex way, but that doesn't matter
      #since the variances have all already been re-scaled
      
      constantdecayvar <- constantdecayvar/variRCT
    
    }
    
    if(input$select == 4) {
      
      myvars$SteppedWedge <-  pnorm( -1.96 + sqrt(1/myvars$SteppedWedge)*input$effsize )
      myvars$Parallel <-  pnorm( -1.96 + sqrt(1/myvars$Parallel)*input$effsize )
      myvars$ParallelBaseline <-  pnorm( -1.96 + sqrt(1/myvars$ParallelBaseline)*input$effsize )
      myvars$CRXO <- pnorm( -1.96 + sqrt(1/myvars$CRXO)*input$effsize )
      
      constantdecayvar <- pnorm( -1.96 + sqrt(1/constantdecayvar)*input$effsize )
      
      
    }
    
      myvars$SteppedWedge <- round(myvars$SteppedWedge, 4)
      myvars$Parallel <- round(myvars$Parallel, 4)
      myvars$ParallelBaseline <- round(myvars$ParallelBaseline, 4)
      myvars$CRXO <- round(myvars$CRXO, 4)  
      constantdecayvar <- round(constantdecayvar, 4)
      
    myplot<- plot_ly(myvars, x = ~decay, y = ~SteppedWedge, name = 'Stepped wedge', type = 'scatter', mode = 'lines',
                          line = list(color = "black", width = 4)) %>%
        add_trace(y = ~Parallel, name = 'Parallel', line = list(color = "darkred", width = 4, dash = 'dash')) %>%
        add_trace(y = ~ParallelBaseline, name = 'Pllel + BL', line = list(color = "darkblue", width = 4, dash = 'dot')) %>%
        add_trace(y = ~CRXO, name = 'CRXO', line = list(color = "darkgreen", width = 4, dash = 'dashdot')) %>%
        layout(xaxis = list(title = "Decay (1-r)"), yaxis = list (title = ""),
               legend=list(orientation="h")) 
    
    if(input$select==4) myplot <- myplot %>% add_trace(x = c(0, 1), y = c(1, 1), line=list(width = 0.5))
    
    
    if(input$HooperYN==2) {
      myplot <- myplot  %>%
        add_trace(x = c(0, 1), y = c(constantdecayvar[1],constantdecayvar[1]), name = 'Stepped wedge', line=list(color="grey",width = 2, dash="solid"), showlegend=FALSE) %>%
        add_trace(x=input$rCons, y=constantdecayvar[1], marker=list(color="black"), name = 'Stepped wedge', showlegend=FALSE) %>%
        add_trace(x = c(0, 1), y = c(constantdecayvar[2],constantdecayvar[2]), name="Parallel", line=list(color="grey",width = 2, dash="solid"), showlegend=FALSE) %>%
        add_trace(x=input$rCons, y=constantdecayvar[2], marker=list(color="darkred"), name="Parallel", showlegend=FALSE) %>%
        add_trace(x = c(0, 1), y = c(constantdecayvar[3],constantdecayvar[3]), name = 'Pllel + BL', line=list(color="grey",width = 2, dash="solid"), showlegend=FALSE) %>%
        add_trace(x=input$rCons, y=constantdecayvar[3], marker=list(color="darkblue"), name = 'Pllel + BL', showlegend=FALSE) %>%
        add_trace(x = c(0, 1), y = c(constantdecayvar[4],constantdecayvar[4]), name="CRXO", line=list(color="grey",width = 2, dash="solid"), showlegend=FALSE) %>%
        add_trace(x=input$rCons, y=constantdecayvar[4], marker=list(color="darkgreen"), name="CRXO", showlegend=FALSE) 
    }
    
    
    myplot
    
  })
  
  #Comparing the exponential decay model to the Hussey and Hughes model
  output$HHrelvarplot <- renderPlotly({
    #Generate the design matrices:
    Xsw <- SWdesmat2(input$T)
    Xpllel <- plleldesmat2(input$T)
    Xpllelbase <- pllelbasedesmat2(input$T)
    Xcrxo <- crxodesmat2(input$T)
    
    Xsw <- matrix(data=as.vector(t(Xsw)), nrow=input$nclust*nrow(Xsw), ncol=ncol(Xsw), byrow = TRUE)
    Xpllel <- matrix(data=as.vector(t(Xpllel)), nrow=input$nclust*nrow(Xpllel), ncol=ncol(Xpllel), byrow = TRUE)
    Xpllelbase <- matrix(data=as.vector(t(Xpllelbase)), nrow=input$nclust*nrow(Xpllelbase), ncol=ncol(Xpllelbase), byrow = TRUE)
    Xcrxo <- matrix(data=as.vector(t(Xcrxo)), nrow=input$nclust*nrow(Xcrxo), ncol=ncol(Xcrxo), byrow = TRUE)
    
    if((input$T-1)%%2 == 1) { correction <- (input$T)/(input$T-1) }
    if((input$T-1)%%2 == 0) { correction <- 1}
    
    
    myvars<- data.frame(decay = seq(1,0,-0.01), 
                        SteppedWedge= simplify2array(lapply(seq(0, 1, 0.01), ExpDecayVar, Xmat = Xsw, m=input$m, rho0=input$rho0)),
                        Parallel = correction*simplify2array(lapply(seq(0, 1, 0.01), ExpDecayVar, Xmat = Xpllel, m=input$m, rho0=input$rho0)),
                        ParallelBaseline = correction*simplify2array(lapply(seq(0, 1, 0.01), ExpDecayVar, Xmat = Xpllelbase, m=input$m, rho0=input$rho0)),
                        CRXO = correction*simplify2array(lapply(seq(0, 1, 0.01), ExpDecayVar,Xmat = Xcrxo, m=input$m, rho0=input$rho0)))
    

    
    
    constantdecayvar <- c(HooperVar(1-input$rCons, Xsw, input$m, input$rho0),
                          correction*HooperVar(1-input$rCons, Xpllel, input$m, input$rho0),
                          correction*HooperVar(1-input$rCons, Xpllelbase, input$m, input$rho0),
                          correction*HooperVar(1-input$rCons, Xcrxo, input$m, input$rho0))
    
    
    #Compare to a standard H+H design!
    myvarsHH<- data.frame(sw= simplify2array(lapply(rep(1,length(seq(0, 1, 0.01))), ExpDecayVar, Xmat = Xsw, m=input$m, rho0=input$rho0)),
                          pllel = correction*simplify2array(lapply(rep(1,length(seq(0, 1, 0.01))), ExpDecayVar, Xmat = Xpllel, m=input$m, rho0=input$rho0)),
                          pllelbase = correction*simplify2array(lapply(rep(1,length(seq(0, 1, 0.01))), ExpDecayVar, Xmat = Xpllelbase, m=input$m, rho0=input$rho0)),
                          crxo = correction*simplify2array(lapply(rep(1,length(seq(0, 1, 0.01))), ExpDecayVar,Xmat = Xcrxo, m=input$m, rho0=input$rho0)))
    
    myvars$SteppedWedge <-  myvars$SteppedWedge/myvarsHH$sw
    myvars$Parallel <-  myvars$Parallel/myvarsHH$pllel
    myvars$ParallelBaseline <-  myvars$ParallelBaseline/myvarsHH$pllelbase
    myvars$CRXO <-  myvars$CRXO/myvarsHH$crxo
    
    constantdecayvar[1] <- constantdecayvar[1]/myvarsHH$sw[1]
    constantdecayvar[2] <- constantdecayvar[2]/myvarsHH$pllel[1]
    constantdecayvar[3] <- constantdecayvar[3]/myvarsHH$pllelbase[1]
    constantdecayvar[4] <- constantdecayvar[4]/myvarsHH$crxo[1]
    
    myvars$SteppedWedge <- round(myvars$SteppedWedge, 3)
    myvars$Parallel <- round(myvars$Parallel, 3)
    myvars$ParallelBaseline <- round(myvars$ParallelBaseline, 3)
    myvars$CRXO <- round(myvars$CRXO, 3)  
    constantdecayvar <- round(constantdecayvar, 3)
    
    myplot<- plot_ly(myvars, x = ~decay, y = ~SteppedWedge, name = 'Stepped wedge', type = 'scatter', mode = 'lines',
                     line = list(color = "black", width = 4)) %>%
      add_trace(y = ~Parallel, name = 'Parallel', line = list(color = "darkred", width = 4, dash = 'dash')) %>%
      add_trace(y = ~ParallelBaseline, name = 'Pllel + BL', line = list(color = "darkblue", width = 4, dash = 'dot')) %>%
      add_trace(y = ~CRXO, name = 'CRXO', line = list(color = "darkgreen", width = 4, dash = 'dashdot')) %>%
      layout(xaxis = list(title = "Decay (1-r)"),
             yaxis = list (title = ""), legend=list(orientation="h"))
    
    
    if(input$HooperYN==2) {
      myplot <- myplot  %>%
        add_trace(x = c(0, 1), y = c(constantdecayvar[1],constantdecayvar[1]), name = 'Stepped wedge', line=list(color="grey",width = 2, dash="solid"), showlegend=FALSE) %>%
        add_trace(x=input$rCons, y=constantdecayvar[1], marker=list(color="black"), name = 'Stepped wedge', showlegend=FALSE) %>%
        add_trace(x = c(0, 1), y = c(constantdecayvar[2],constantdecayvar[2]), name="Parallel", line=list(color="grey",width = 2, dash="solid"), showlegend=FALSE) %>%
        add_trace(x=input$rCons, y=constantdecayvar[2], marker=list(color="darkred"), name="Parallel", showlegend=FALSE) %>%
        add_trace(x = c(0, 1), y = c(constantdecayvar[3],constantdecayvar[3]), name = 'Pllel + BL', line=list(color="grey",width = 2, dash="solid"), showlegend=FALSE) %>%
        add_trace(x=input$rCons, y=constantdecayvar[3], marker=list(color="darkblue"), name = 'Pllel + BL', showlegend=FALSE) %>%
        add_trace(x = c(0, 1), y = c(constantdecayvar[4],constantdecayvar[4]), name="CRXO", line=list(color="grey",width = 2, dash="solid"), showlegend=FALSE) %>%
        add_trace(x=input$rCons, y=constantdecayvar[4], marker=list(color="darkgreen"), name="CRXO", showlegend=FALSE) 
    }
    
    
    myplot 
    
  })
  
  
  })


