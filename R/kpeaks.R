genpolygon <- function(x, binrule, nbins, disp=FALSE){
  if(missing(x)) 
    stop("Missing data")
  if(!is.numeric(x)) 
    stop("Non-numeric data")
  if(!is.vector(x)) 
     stop("Data x must be a vector")
  n <- length(x)
  if(n<2) 
    stop("Data vector should contain at least two observations")
  if(missing(nbins)) 
    nbins <- 30
  if(!is.numeric(nbins))
    stop("Number of classes should be a positive integer")
  if(!missing(nbins) & missing(binrule))
    binrule="usr"
  if(missing(binrule))
    binrule="sturges"
  xsd <- sqrt(var(x))
  xrange <- diff(range(x))
  xskew <- n*(sum((x-mean(x))^3)/(xsd^3))/((n-1)*(n-2))
  switch(binrule,
    sqr = {nbins <- floor(sqrt(n))},
    sturges = {nbins <- ceiling(1+log(n, 2))},
    huntsberger = {nbins <- ceiling(1+3.332*log(n, 10))},
    bc = {nbins <- ceiling(5*log(n, 10))},
    cencov = {nbins <- ceiling(n^(1/3))},
    rice = {nbins <- ceiling(2*n^(1/3))},
    ts = {nbins <- ceiling((2*n)^(1/3))},
    scott = {nbins <- ceiling(xrange/(3.5*xsd*n^(-1/3)))},
    fd = {nbins <- ceiling(xrange/(2*IQR(x)*n^(-1/3)))},
    doane = {nbins <- ceiling(1+log(n,2)+log(1+abs(xskew)/(6*(n-2)/((n+1)*(n+3))^(1/2)),2))},
    cebeci = {nbins <- ceiling(log(n)/2*pi)},
    usr = {nbins <- nbins},
    stop("Enter a valid binning rules: 'sqr', 'sturges', 'scott', 'huntsberger', 'bc', 'cencov', 'rice', 'ts', 'fd', 'doane', 'cebeci' and 'usr'.")) 
  hdata <- hist(x, breaks=nbins, plot=disp)
  frepoly <- list()
    frepoly$freqs <- hdata$counts
    frepoly$mids <- hdata$mids
    frepoly$nbins <- length(hdata$mids)
  return(frepoly)
}

plotpolygon <- function(x, nbins, ptype, bcol="gray", pcol="blue"){
  if(missing(x)) 
    stop("Missing data vector")
  if(!is.numeric(x)) 
    stop("Non-numeric data")
  if(!is.vector(x)) 
     stop("Data should be a vector")
  n <- length(x)
  if(n<2) 
    stop("Data vector should contain at least two observations")
  if(missing(nbins)) 
    nbins <- 30
  if(!is.numeric(nbins))
    stop("Number of bins should be a positive integer")
  if(missing(ptype)) 
    ptype <- "sp"
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5))
   switch(ptype,
    p  = {hx <- hist(x, breaks=nbins, border=FALSE, main=paste("Frequency Polygon"))},
    ph = {hx <- hist(x, breaks=nbins, col=bcol, border=TRUE, main="Freq. Polygon & Histogram")},
    sp = {par(new=TRUE); hx <- hist(x, breaks=nbins, main="", xlab="", ylab="", axes = FALSE)},
    {stop("Valid plot types are 'p', 'ph' and 'sp'.")})
  m <- hx$breaks[2]-hx$breaks[1]
  xm <- c(min(hx$mids)-m, hx$mids, max(hx$mids)+m)
  xc <- c(0, hx$counts, 0)
  lines(xm, xc, col=pcol, lwd=2)
  abline(h=0)
}

findpolypeaks <- function(xm, xc, tcmethod, tc){
  if(missing(xm))
    stop("Missing data vector for the middles of classes")
  if(!is.numeric(xm)) 
    stop("The middles of classes should be numeric")
  if(!is.vector(xm)) 
    stop("Data type of the middles of classes should be a vector")
  if(length(xm)<1) 
    stop("Data vector for the middles of classes should contain at least one value")
  if(missing(xc)) 
    stop("Missing data for the frequencies of classes")
  if(!is.numeric(xc)) 
    stop("The frequencies of classes should be numeric")
  if(!is.vector(xc)) 
     stop("Data type of the frequencies of classes should be a vector")
  if(length(xc)<1) 
    stop("Data vector for the frequencies of classes should contain at least one value")
  if(length(xm)!=length(xc)) 
    stop("Data vectors for the frequencies and middles of classes should be equal length")
  if(missing(tcmethod)) 
    tcmethod <- "usr"
  if(missing(tc)) 
    tc <- 1
  if(!is.numeric(tc)) 
    stop("Threshold frequency value should be numeric")
  switch(tcmethod,
    sd1 = {tc <- sqrt(var(xc))},
    sd2 = {tc <- 0.25*sqrt(var(xc))},
    q1 = {tc <- unname(quantile(xc, 0.25))},
    iqr = {tc <- unname(0.25*(quantile(xc, 0.75)-quantile(xc, 0.25)))},
    min = {tc <- min(xc[which(xc > 0)])},
    min2 = {tc <- 2*min(xc[which(xc > 0)])},
    log2 = {tc <- log(length(xc)/10,2)},
    avg = {tc <- 0.10*mean(xc)},
    usr = {tc <- tc},
    stop("The valid methods are 'sd1', 'sd2', 'q1', 'iqr', 'min', 'min2', 'log2', 'avg' and 'usr' for threshold value calculation."))
  xm <- xm[xc >= tc]
  xc <- xc[xc >= tc]
  nc <- length(xc)
  pfreqs <- c()
  pvalues <- c()
  pidx <- 1
  if(nc > 1){
    if(xc[1] > xc[2]){
          pvalues[1] <- xm[1]
          pfreqs[1] <- xc[1]
          pidx <- 2
    }
    for(i in 2:(nc-1)){
      if(xc[i] != xc[i-1]){
      	if((xc[i] > xc[i-1]) && (xc[i] >= xc[i+1])){
          pvalues[pidx] <- xm[i]
          pfreqs[pidx] <- xc[i]
          pidx <- pidx + 1
	}
      }
    }
    if(xc[nc] > xc[nc-1]){
          pvalues[pidx] <- xm[nc]
          pfreqs[pidx] <- xc[nc]
     }
  }
  else{
    pvalues <- xm
    pfreqs <- xc
  }
  peaks <- list()
    peaks$pm <- cbind(pvalues, pfreqs)
    peaks$np <- length(pvalues)
  return(peaks)
}

rmshoulders <- function(xm, xc, trmethod, tv){
  if(missing(xm)) 
    stop("Data for the middles of classes is missing")
  if(!is.numeric(xm)) 
    stop("The middles of classes should be numeric")
  if(!is.vector(xm)) 
     stop("Data type of the middles of classes should be a vector")
  if(length(xm)<1) 
    stop("Data vector for the middles of classes should contain at least one value")
  if(missing(xc)) 
    stop("Missing data for the frequencies of classes")
  if(!is.numeric(xc)) 
    stop("The frequencies of classes should be numeric")
  if(!is.vector(xc)) 
     stop("Data type of the frequencies of classes should be a vector")
  if(length(xc)<1) 
    stop("Data vector for the frequencies of classes should contain at least one value")
  if(length(xm)!=length(xc)) 
    stop("Data vectors should be equal length for the frequencies and middles of classes")
  if(missing(tv) & missing(trmethod))
    trmethod <- "avg"
  if(missing(tv))
    tv <- 1
  if(!missing(tv) & missing(trmethod))
    trmethod <- "usr"
  if(!is.numeric(tv)) 
    stop("Threshold distance value should be numeric")
  npc <- length(xc)
  pfreqs <- c(); pvalues <- c()
  deleted <- c(); delidx <- 0	
  switch(trmethod,
    sd = {tv <- sqrt(var(abs(diff(xm))))},
    q1 = {tv <- unname(quantile(abs(diff(xm)), 0.25))},
    q3 = {tv <- unname(quantile(abs(diff(xm)), 0.75))},
    iqr = {tv <- 0.25*IQR(abs(diff(xm)))},
    avg = {tv <- mean(abs(diff(xm)))},
    med = {tv <- median(abs(diff(xm)))},
    usr = {tv <- tv},
    {stop("The valid options are 'sd', 'q1', 'q3', 'iqr', 'avg', 'med' and 'usr' for calculation of the threshold distance value.")})
  if(npc > 2){
    for(i in 1:(npc-1)){
       if(abs(xm[i] - xm[i+1]) <= tv){
        delidx <- delidx + 1
	if(xc[i] < xc[i+1]){
	   deleted[delidx] <- i
        }
        else{
           deleted[delidx] <- i+1
        }
      }
    }
    if(delidx >= 1){
      xm <- xm[-deleted]
      xc <- xc[-deleted]
    }
  }
  pvalues <- xm ; pfreqs <- xc
  peaks <- list()
   peaks$pm <- cbind(pvalues, pfreqs)
   peaks$np <- length(pvalues)
  return(peaks)
}

.estimatek <- function(pc, rcs=TRUE, tpc=1){
  pc <- as.integer(pc)
  dst <- "Full"
  if(rcs){
    pc <- subset(pc, pc > tpc)
    dst <- "Reduced"
  }
  am <- med <- mod <- cr <- ciqr <- mppc <- mq3m <- mtl <- ptk <- am <- avgk <- modk <- mtlk <- NA
  p <- length(pc)
  if(p>1){
    pcm <- as.matrix(table(pc)) 
    pcmval <- as.integer(rownames(pcm))
    pcmfq <- pcm[,1]
    sumpc <- 0; npc <- 0 
    for(v1 in 1:(p-1)){
      for(v2 in (v1+1):p){	
        sumpc <- sumpc + (pc[v1]+pc[v2])/2
        npc <- npc + 1
      }
    }
    sortedpc <- sort(pc, decreasing=TRUE)
    quantpc <- unname(quantile(pc))
    am <- ceiling(mean(pc))
    med <- ceiling(median(pc))
    if(length(unique(pcmfq)) != 1){
    	mod <- ceiling((pcmval[which.max(pcmfq)]))
    }
    else{
	mod <- NA
    }
    cr <- round((max(pc) + min(pc))/2)
    ciqr <- round((quantpc[4] + quantpc[2])/2)
    mppc <- round(sumpc/npc)
    mq3m <- round((quantpc[4] + quantpc[5])/2)
    mtl <- round((sortedpc[1] + sortedpc[2])/2)
    estimates <- sort(c(am,med,mod,cr,ciqr,mppc,mq3m,mtl), decreasing=TRUE) 
    estm <- as.matrix(table(estimates))
    estmval <- as.integer(rownames(estm))
    estmfq <- estm[,1]
    avgk <- round(mean(estimates))
    modk <- round((estmval[which.max(estmfq)]))
    mtlk <- round((estimates[1]+estimates[2])/2)
  } else{
    am <- pc
    med <- pc
  }
  k <- list()
    k$am <- am
    k$med <- med
    k$mod <- mod
    k$mppc <- mppc
    k$cr <- cr
    k$ciqr <- ciqr
    k$mq3m <- mq3m
    k$mtl <- mtl
    k$avgk <- avgk
    k$modk <- modk
    k$mtlk <- mtlk
    k$dst <- dst
    k$pcounts <- pc
  return(k)
}

findk <- function(x, binrule, nbins, tcmethod, tc, trmethod, tv, rms=FALSE, rcs=FALSE, tpc=1){
  if(missing(x)) 
    stop("Missing dataset")
  if (is(x, "data.frame"))
    x <- as.matrix(x)
  if (!is(x, "matrix"))
    stop("Dataset must be a numeric data frame or a matrix")
  if(missing(binrule))
    binrule <- "sturges"
  if(missing(nbins))
    nbins <- 20
  if(missing(tc))
    tc <- 0
  if(missing(tcmethod)){
    tcmethod <- "usr"
    tc <- 1
  }
  if(missing(tv))
    tv <- 0
  if(missing(trmethod)){
    trmethod <- "min"
  }
  pcounts <- c()
  for(i in 1:ncol(x)){
    hres <- genpolygon(x[,i], binrule=binrule, nbins=nbins)
    fpres <- findpolypeaks(hres$mids, hres$freqs, tcmethod=tcmethod, tc=tc)
    npc <- fpres$np 
    if(rms){
      fpres <- rmshoulders(fpres$pm[,1], fpres$pm[,2], trmethod=trmethod, tv=tv)
      npc <- fpres$np 
    }
    pcounts <- c(pcounts, npc)
  }
  k <- list()
  k <- .estimatek(pc=pcounts, rcs=rcs, tpc=1)
  return(k)
}