#development script for popgen simulations (+/- R translation of PopG, adding more visualization and summary stats)
library(plyr);library(reshape);library(ggplot2);library(magrittr)

#binomial draw for new genotype frequencies
binomialDraw <- function(n,p){ 
  draws <- runif(n)
  bnl <- length(draws[draws<p])
  return(bnl)
}

#main simulation function
runPopSim <- function(gen=100,p=0.5,Waa=1,Wab=1,Wbb=1,n=100,nPop=2,m=0,stats=c("p","Fst"),drift=T){
  allele.freq <- data.frame(matrix(ncol=3*nPop))
  allele.freq[1,(1:nPop)] <- rep(p,nPop) #starting allele freqs
  for(i in 1:gen){ 
    mean.p <- as.numeric(rowMeans(allele.freq[i,(1:nPop)]))
    for(j in 1:nPop){
      p <- allele.freq[i,j]
      q <- 1-p
      pp <- p #store parental allele freq's
      p <- p*(1-m)+m*mean.p # migration
      if(p>0 && p<1){ #if alleles are not fixed
        w <- p*p*Waa+2*p*q*Wab+q*q*Wbb #population average fitness
        freq.aa <- (p*p*Waa)/w #get post-selection genotype frequencies (ie weighted by relative fitness)
        freq.ab <- (2*p*q*Wab)/w
        if(drift==T){ 
            Naa <- binomialDraw(n,freq.aa) # binomial draw for number of new genotype counts (ie drift)
          if(freq.aa<1){ 
            Nab <- binomialDraw((n-Naa),(freq.ab/(1-freq.aa)))
          }
          else {
            Nab <- 0
          }
          allele.freq[(i+1),j] <- ((2*Naa)+Nab)/(2*n)#new p after drift in columns 1:nPop
          allele.freq[(i+1),(j+nPop)] <- Nab/n #Ho in columns (nPop+1):(nPop*2)
          allele.freq[(i+1),(j+2*nPop)] <- 2*pp*q #He in columns (nPop*2+1):nPop*3
        } 
        else { #no drift (infinite population) conditions
          allele.freq[(i+1),j] <- freq.aa+(freq.ab/2)
          allele.freq[(i+1),(j+nPop)] <- freq.ab 
          allele.freq[(i+1),(j+2*nPop)] <- 2*pp*q
        }
      } else { #if alleles are fixed
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
  allele.freq$meanHo <- rowMeans(allele.freq[(nPop+1):(nPop*2)])
  allele.freq$meanHe <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
  allele.freq$Fis <- abs(1-(allele.freq$meanHo/allele.freq$meanHe))
  allele.freq$mean.p <- rowMeans(allele.freq[1:nPop]) 
  allele.freq$Hs <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
  allele.freq$Ht <- 2*allele.freq$mean.p*(1 - allele.freq$mean.p)
  allele.freq$Fst <- (allele.freq$Ht-allele.freq$Hs)/allele.freq$Ht
  allele.freq$Fst[allele.freq$Fst<0] <- 0
  allele.freq$gen <- 0:gen
  df <- melt(allele.freq,id.vars = "gen")
  df$dataType <- c(rep("p",(nPop*(gen+1))),rep("Ho",nPop*(gen+1)),rep("He",nPop*(gen+1)),rep("meanHo",(gen+1)),rep("meanHe",(gen+1)),
                   rep("Fis",(gen+1)),rep("mean.p",(gen+1)),rep("Hs",(gen+1)),rep("Ht",(gen+1)),rep("Fst",(gen+1)))
  df <- subset(df,dataType %in% stats)
  return(df)
}

#duplicate function with wide-format output for printing summary tables.
runPopSim.noMelt <- function(gen=100,p=0.5,Waa=1,Wab=1,Wbb=1,n=100,nPop=2,m=0,stats=c("p","Fst"),drift=T){
  allele.freq <- data.frame(matrix(ncol=3*nPop))
  allele.freq[1,(1:nPop)] <- rep(p,nPop) #starting allele freqs
  for(i in 1:gen){ 
    mean.p <- as.numeric(rowMeans(allele.freq[i,(1:nPop)]))
    for(j in 1:nPop){
      p <- allele.freq[i,j]
      pp <- p
      p <- p*(1-m)+m*mean.p # migration
      if(p>0 && p<1){ #if alleles are not fixed
        q <- 1-p
        w <- p*p*Waa+2*p*q*Wab+q*q*Wbb
        freq.aa <- (p*p*Waa)/w # selection
        freq.ab <- (2*p*q*Wab)/w
        if(drift==T){ 
          Naa <- binomialDraw(n,freq.aa) # drift (popG magic)
          if(freq.aa<1){ 
            Nab <- binomialDraw((n-Naa),(freq.ab/(1-freq.aa)))
          }
          else {
            Nab <- 0
          }
          allele.freq[(i+1),j] <- ((2*Naa)+Nab)/(2*n)#new p after drift
          allele.freq[(i+1),(j+nPop)] <- Nab/n #Ho
          allele.freq[(i+1),(j+2*nPop)] <- 2*pp*q #He
        } 
        else { #no drift (infinite population) conditions
          allele.freq[(i+1),j] <- freq.aa+(freq.ab/2)
          allele.freq[(i+1),(j+nPop)] <- freq.ab 
          allele.freq[(i+1),(j+2*nPop)] <- 2*pp*q
        }
      } else { #if alleles are fixed
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
  allele.freq$meanHo <- rowMeans(allele.freq[(nPop+1):(nPop*2)])
  allele.freq$meanHe <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
  allele.freq$Fis <- abs(1-(allele.freq$meanHo/allele.freq$meanHe))
  allele.freq$mean.p <- rowMeans(allele.freq[1:nPop]) 
  allele.freq$Hs <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
  allele.freq$Ht <- 2*allele.freq$mean.p*(1 - allele.freq$mean.p)
  allele.freq$Fst <- (allele.freq$Ht-allele.freq$Hs)/allele.freq$Ht
  allele.freq$Fst[allele.freq$Fst<0] <- 0
  allele.freq$gen <- 0:gen
  return(allele.freq)
}

#plotting function
plotSingleRun <- function(df,nPop,gen){
  print(ggplot(df,aes(x=gen,y=value,col=variable))+ylim(0,1)+facet_wrap(~dataType)+geom_line())
}


