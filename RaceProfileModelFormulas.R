dataToPCA <- function(data, basisNo)
{
  dataTable = data[,7:(7+basisNo)]
  numberEvents = nrow(data)
  length = as.numeric(colnames(data[1,][7+basisNo]))
  
  sortData <- function(n_curves = numberEvents){
    ID <- 1:n_curves
    x <- vector(mode = "list", length = n_curves)
    t <- vector(mode = "list", length = n_curves)
    
    for (i in 1:n_curves){
      t[i] <- list(seq(0, length, length / basisNo))
      x[i] <- list(unlist(unname(data[i,7:(7 + basisNo)])))
    }
    df <- tibble(ID,t,x)
    names(df) <- c("ID", "Time", "Curve")
    return(df)
  }
  
  df = sortData(numberEvents)
  
  df_1 <- df %>% select(!c(ID,Curve)) %>% unnest_longer(Time) 
  df_2 <- df %>% select(!c(ID,Time)) %>% unnest_longer(Curve)
  ID <- sort(rep(1:numberEvents,(basisNo + 1)))
  df_l <- cbind(ID,df_1,df_2)
  
  
  knots    = c(seq(0,length, length/basisNo)) #Location of knots
  n_knots   = length(knots) #Number of knots
  n_order   = 4 # order of basis functions: for cubic b-splines: order = 3 + 1
  n_basis   = length(knots) + n_order - 2;
  basis = create.bspline.basis(rangeval = c(0,length), n_basis)
  
  argvals <- matrix(df_l$Time, nrow = (basisNo + 1), ncol = numberEvents)
  y_mat <- matrix(as.numeric(df_l$Curve), nrow = basisNo + 1, ncol = numberEvents)
  
  W.obj <- Data2fd(argvals = argvals, y = y_mat, basisobj = basis, lambda = 0.5)
  fun_pca <- pca.fd(W.obj, nharm = 4)
  return(fun_pca)
}

dataTransform <- function(data, pcaObject, all, basisNo)
{
  pcData = as.data.frame(cbind(rownames(data), data[,2], data[,1], data[,3], as.numeric(pcaObject$scores[,1]), as.numeric(pcaObject$scores[,2]), as.numeric(pcaObject$scores[,3]), as.numeric(pcaObject$scores[,4]), data[,5], data[,4], data[,basisNo + 8]))
  colnames(pcData) <- c("ID", "Level", "Athlete", "Type", "PC1", "PC2", "PC3", "PC4", "Date", "Year", "Age")
  pcData$PC1 = as.numeric(pcData$PC1)
  pcData$PC2 = as.numeric(pcData$PC2)
  pcData$PC3 = as.numeric(pcData$PC3)
  pcData$PC4 = as.numeric(pcData$PC4)
  
  pcData$Level <- factor(pcData$Level , levels=c("Domestic", "World Cup / Juniors", "WC / Olympics"))
  
  if(all == FALSE)
  {
    pcData = pcData[!(pcData$Level == "Domestic" & pcData$Type == "Heat") & !(pcData$Level == "Domestic" & pcData$Type == "Semi"),]
  }
  
  pcDataSort = pcData
  pcDataSort$Date = as.Date(pcDataSort$Date, "%d/%m/%Y")
  pcDataSort = pcDataSort[order(pcDataSort$Date),]
  pcDataSort = pcDataSort[order(pcDataSort$Athlete),]
  
  for(i in 2:nrow(pcDataSort))
  {
    if(pcDataSort$Date[i] <= pcDataSort$Date[i-1] & pcDataSort$Athlete[i] == pcDataSort$Athlete[i-1])
    {
      pcDataSort$Date[i] = pcDataSort$Date[i-1] + 10
    }
  }
  return(pcDataSort)
}


PCAtoHMM <- function(pcData, iterations)
{
  pcData$Athlete = as.numeric(pcData$Athlete)
  pcData$Age = as.factor(pcData$Age)
  pcDataCount = aggregate(pcData$ID, by = list(pcData$Athlete), FUN = "length")
  
  set.seed(100)
  model <-
    depmix(
      list(PC1 ~ Age + Level, 
           PC2 ~ Age + Level, 
           PC3 ~ Age + Level,
           PC4 ~ Age + Level),
      data = pcData,
      nstates = 4,
      ntimes = pcDataCount[,2],
      family =  list(gaussian(),gaussian(),gaussian(),gaussian()),
      transition = ~ 1,
      instart = runif(4)
    )
  fit = multistart(model, nstart = iterations, initIters = 200, emc = em.control(rand = TRUE))
  
  return(fit)
}

exportCareerStatePlot <- function(data, HMM, minDate, maxDate, minimumRaces, gender)
{
  pcDataAdd = cbind(data,data.frame(t=1:nrow(data),variable='state.multi',value=HMM@posterior$state), posterior(HMM, type = "viterbi"))
  pcDataAdd$response = HMM@response[[1]][[1]]@y
  data$ID = as.numeric(data$ID)
  data$Age = as.factor(data$Age)
  dataCount = aggregate(data$ID, by = list(data$Athlete), FUN = "length")
  dataCount = dataCount[dataCount$x > minimumRaces - 1,]
  colnames(dataCount) = c("ID", "Number")
  
  for(i in dataCount$ID)
  {
    ggplot(pcDataAdd[pcDataAdd$Athlete %in% i,],aes(x=Date,y=value, colour = Age, shape = Age)) +
      geom_jitter(size = 4, height = 0, width = 1) + scale_shape_manual(values= c(16, 17, 18, 20)) + scale_x_date(date_breaks = "1 years", date_labels = ("%Y"), 
                                                                   limits = as.Date(c(paste(minDate,'/01/01',sep=""),paste(maxDate + 1,'/01/01',sep="")))) + ylab("Predicted State") + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=13), legend.text=element_text(size=13),title = element_text(size = 12), 
            strip.text = element_text(size=12)) + ggtitle(paste("Athlete", i, "Career State Tracking")) + ylim(1, 4)
    if(gender == "Men")
    {
      ggsave(file = file.path("Plot Outputs", "Men", "Career State Tracking", paste("Athlete ",  i, " - ", " Career State Tracking.png", sep = "")), height = 10, width = 17, bg = "white")
    }
    else
    {
      ggsave(file = file.path("Plot Outputs", "Women", "Career State Tracking", paste("Athlete ",  i, " - ", " Career State Tracking.png", sep = "")), height = 10, width = 17, bg = "white")
    }
  }
}

exportPCPlot <- function(data, HMM, minDate, maxDate, minimumRaces, gender)
{
  pcDataAdd = cbind(data,data.frame(t=1:nrow(data),variable='state.multi',value=HMM@posterior$state), posterior(HMM, type = "viterbi"))
  pcDataAdd$response = HMM@response[[1]][[1]]@y
  data$ID = as.numeric(data$ID)
  data$Age = as.factor(data$Age)
  dataCount = aggregate(data$ID, by = list(data$Athlete), FUN = "length")
  dataCount = dataCount[dataCount$x > minimumRaces - 1,]
  colnames(dataCount) = c("ID", "Number")
  
  for(i in dataCount$ID)
  {
    a <- ggplot(pcDataAdd[pcDataAdd$Athlete %in% i,],aes(x=Date,y=PC1)) +
    geom_jitter(size = 4, height = 0, width = 1, colour = "blue") + scale_shape_manual(values= c(16, 17, 18, 20)) + scale_x_date(date_breaks = "1 years", date_labels = ("%Y"), 
                                                                                                                limits = as.Date(c(paste(minDate,'/01/01',sep=""),paste(maxDate + 1,'/01/01',sep="")))) + ylab("Predicted State") + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=13), legend.text=element_text(size=13),title = element_text(size = 12), 
          strip.text = element_text(size=12)) + ggtitle("PC1") + ylim(-1, 1) +geom_smooth(method = "lm", formula = y~x)
    b <- ggplot(pcDataAdd[pcDataAdd$Athlete %in% i,],aes(x=Date,y=PC2)) +
      geom_jitter(size = 4, height = 0, width = 1, colour = "blue") + scale_shape_manual(values= c(16, 17, 18, 20)) + scale_x_date(date_breaks = "1 years", date_labels = ("%Y"), 
                                                                                                                                   limits = as.Date(c(paste(minDate,'/01/01',sep=""),paste(maxDate + 1,'/01/01',sep="")))) + ylab("Predicted State") + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=13), legend.text=element_text(size=13),title = element_text(size = 12), 
            strip.text = element_text(size=12)) + ggtitle("PC2") + ylim(-1, 1) +geom_smooth(method = "lm", formula = y~x)
    c <- ggplot(pcDataAdd[pcDataAdd$Athlete %in% i,],aes(x=Date,y=PC3)) +
      geom_jitter(size = 4, height = 0, width = 1, colour = "blue") + scale_shape_manual(values= c(16, 17, 18, 20)) + scale_x_date(date_breaks = "1 years", date_labels = ("%Y"), 
                                                                                                                                   limits = as.Date(c(paste(minDate,'/01/01',sep=""),paste(maxDate + 1,'/01/01',sep="")))) + ylab("Predicted State") + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=13), legend.text=element_text(size=13),title = element_text(size = 12), 
            strip.text = element_text(size=12)) + ggtitle("PC3") + ylim(-1, 1) +geom_smooth(method = "lm", formula = y~x)
    d <- ggplot(pcDataAdd[pcDataAdd$Athlete %in% i,],aes(x=Date,y=PC4)) +
      geom_jitter(size = 4, height = 0, width = 1, colour = "blue") + scale_shape_manual(values= c(16, 17, 18, 20)) + scale_x_date(date_breaks = "1 years", date_labels = ("%Y"), 
                                                                                                                                   limits = as.Date(c(paste(minDate,'/01/01',sep=""),paste(maxDate + 1,'/01/01',sep="")))) + ylab("Predicted State") + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=13), legend.text=element_text(size=13),title = element_text(size = 12), 
            strip.text = element_text(size=12)) + ggtitle("PC4") + ylim(-1, 1) +geom_smooth(method = "lm", formula = y~x)
    
    e <- ggarrange(a, b, c, d, ncol = 2, nrow = 2)
    annotate_figure(e, top = text_grob(paste("Athlete", i, "Principal Component Tracking"), 
                                          color = "black", face = "bold", size = 18))
    
    if(gender == "Men")
    {
      ggsave(file = file.path("Plot Outputs", "Men", "Principal Component Tracking", paste("Athlete ",  i, " - ", " Principal Component Tracking.png", sep = "")), height = 10, width = 17, bg = "white")
    }
    else
    {
      ggsave(file = file.path("Plot Outputs", "Women", "Principal Component Tracking", paste("Athlete ",  i, " - ", " Principal Component Tracking.png", sep = "")), height = 10, width = 17, bg = "white")
    }
  }
  
}


