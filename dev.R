#development script for popgen simulations (R translation of Joe Felsenstein's PopG, adding more visualization and summary stats)
library(plyr);library(reshape);library(ggplot2);library(magrittr)

#binomial draw function for new generations
binomialDraw <- function(n,p){ 
  draws <- runif(n)
  bnl <- length(draws[draws<p])
  return(bnl)
}

#main simulation function
runPopSim <- function(gen=100,p=0.5,Waa=1,Wab=1,Wbb=1,n=100,nPop=2,m=0,stats=c("p","Fst")){
  allele.freq <- data.frame(matrix(ncol=3*nPop))
  allele.freq[1,] <- rep(p,nPop)
  
  for(i in 1:gen){
    #initialize starting allele frequencies
    mean.p <- as.numeric(rowMeans(allele.freq[i,]))
    for(j in 1:nPop){
      p <- allele.freq[i,j]
      p <- p*(1-m)+m*mean.p
      if(p>0 && p<1){
        
        q <- 1-p
        w <- p*p*Waa+2*p*q*Wab+q*q*Wbb
        freq.aa <- (p*p*Waa)/w
        freq.ab <- (2*p*q*Wab)/w
        #draw gametes for the next generation
          Naa <- binomialDraw(n,freq.aa)
        if(freq.aa<1 & Naa<n){
          Nab <- binomialDraw((n-Naa),(freq.ab/(1-freq.aa)))
        } else {
          Nab <- 0
        }
        allele.freq[(i+1),j] <- ((2*Naa)+Nab)/(2*n)
        allele.freq[(i+1),(j+nPop)] <- Nab/n
        allele.freq[(i+1),(j+2*nPop)] <- freq.ab
        
      } else {
        
        if(p<=0){
          p <- 0
        } else {
          p <- 1
        }
        allele.freq[(i+1),j] <- p
        allele.freq[(i+1),(j+nPop)] <- 0
        allele.freq[(i+1),(j+2*nPop)] <- 0
        
      }
    }
  }
  names <- c()
  for(i in 1:nPop){names[i]<-paste0("p",i)}
  for(i in (nPop+1):(2*nPop)){names[i]<-paste0("Ho",i-nPop)}
  for(i in (nPop+nPop+1):(3*nPop)){names[i]<-paste0("He",i-2*nPop)}
  colnames(allele.freq) <- names
  #calculate Fis
  allele.freq$meanHo <- rowMeans(allele.freq[(nPop+1):(nPop*2)])
  allele.freq$meanHe <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
  allele.freq$Fis <- abs(1-(allele.freq$meanHo/allele.freq$meanHe))
  allele.freq$Hs <- rowMeans(allele.freq[(2*nPop+1):(3*nPop)])
  allele.freq$mean.p <- rowMeans(allele.freq[1:nPop])
  allele.freq$Ht <- 2*allele.freq$mean.p*(1 - allele.freq$mean.p)
  allele.freq$Fst <- 1-allele.freq$Hs/allele.freq$Ht
  allele.freq$Fst[allele.freq$Fst<0] <- 0
  allele.freq$gen <- 0:gen
  df <- melt(allele.freq,id.vars = "gen")
  df$dataType <- c(rep("p",(nPop*(gen+1))),rep("Ho",nPop*(gen+1)),rep("He",nPop*(gen+1)),rep("meanHo",(gen+1)),rep("meanHe",(gen+1)),
                   rep("Fis",(gen+1)),rep("Hs",(gen+1)),rep("mean.p",(gen+1)),rep("Ht",(gen+1)),rep("Fst",(gen+1)))
  df <- subset(df,dataType %in% stats)
  return(df)
}

runPopSim.noMelt <- function(gen=100,p=0.5,Waa=1,Wab=1,Wbb=1,n=100,nPop=2,m=0,stats=c("p","Fst")){
  allele.freq <- data.frame(matrix(ncol=3*nPop))
  allele.freq[1,] <- rep(p,nPop)
  
  for(i in 1:gen){
    #initialize starting allele frequencies
    mean.p <- as.numeric(rowMeans(allele.freq[i,]))
    for(j in 1:nPop){
      p <- allele.freq[i,j]
      p <- p*(1-m)+m*mean.p
      if(p>0 && p<1){
        
        q <- 1-p
        w <- p*p*Waa+2*p*q*Wab+q*q*Wbb
        freq.aa <- (p*p*Waa)/w
        freq.ab <- (2*p*q*Wab)/w
        #draw gametes for the next generation
        Naa <- binomialDraw(n,freq.aa)
        if(freq.aa<1 & Naa<n){
          Nab <- binomialDraw((n-Naa),(freq.ab/(1-freq.aa)))
        } else {
          Nab <- 0
        }
        allele.freq[(i+1),j] <- ((2*Naa)+Nab)/(2*n)
        allele.freq[(i+1),(j+nPop)] <- Nab/n
        allele.freq[(i+1),(j+2*nPop)] <- freq.ab
        
      } else {
        
        if(p<=0){
          p <- 0
        } else {
          p <- 1
        }
        allele.freq[(i+1),j] <- p
        allele.freq[(i+1),(j+nPop)] <- 0
        allele.freq[(i+1),(j+2*nPop)] <- 0
        
      }
    }
  }
  names <- c()
  for(i in 1:nPop){names[i]<-paste0("p",i)}
  for(i in (nPop+1):(2*nPop)){names[i]<-paste0("Ho",i-nPop)}
  for(i in (nPop+nPop+1):(3*nPop)){names[i]<-paste0("He",i-2*nPop)}
  colnames(allele.freq) <- names
  #calculate Fis
  allele.freq$meanHo <- rowMeans(allele.freq[(nPop+1):(nPop*2)])
  allele.freq$meanHe <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
  allele.freq$Fis <- abs(1-(allele.freq$meanHo/allele.freq$meanHe))
  allele.freq$Hs <- rowMeans(allele.freq[(2*nPop+1):(3*nPop)])
  allele.freq$mean.p <- rowMeans(allele.freq[1:nPop])
  allele.freq$Ht <- 2*allele.freq$mean.p*(1 - allele.freq$mean.p)
  allele.freq$Fst <- 1-allele.freq$Hs/allele.freq$Ht
  allele.freq$Fst[allele.freq$Fst<0] <- 0
  allele.freq$gen <- 0:gen
  return(allele.freq)
}

plotSingleRun <- function(df,nPop,gen){
  print(ggplot(df,aes(x=gen,y=value,col=variable))+ylim(0,1)+facet_wrap(~dataType)+geom_line())
}

samplePopSim <- function(nReps=100,gen=100,p=0.5,Waa=1,Wab=1,Wbb=1,n=100,nPop=2,m=0){
  sumTable <- data.frame(matrix(ncol=3*nPop+8))
  for(i in 1:nReps){
    df <- runPopSim(gen=gen,p=p,Waa=Waa,Wab=Wab,Wbb=Wbb,n=n,nPop=nPop,m=m)
    names(sumTable) <- names(df)
    sumTable[i,] <- df[nrow(df),]
  }
  tbl <- colMeans(sumTable,na.rm=T)
  tbl <- tbl[c("Fis","Hs","Ht","Fst")]
  return(tbl)
}


