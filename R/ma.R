ma <-   function(d, partition=list(c(length(d[1,])),c(1:length(d[1,])-1)), 
                 ht=43.6978644, hp=0.8120818, hs=6.0049711 ){
  
    # ht 43.6978644  hp 0.8120818 hs 6.0049711
    # (-1/(ht A + 1) + 1)*(n^hp)/hs
  
    # a helper function
    concat <- function(lis){
      ret <- lis[[1]]
      if(length(lis)>1){
        for(i in 2:length(lis)){
          ret <- append(ret,lis[[i]])
        }
      }
      return(ret)
    }
    
    # another one
    getLikeVec <- function(dataSet,kerWidth,loobParam){
      n <- length(dataSet[,1])
      kvec <- c(kerWidth)
      loob <- c(loobParam)
      liVec <- rep(0.0,times=n)
      six <- order(dataSet[,1])
      ms<-length(dataSet[1,])
      ds <- dataSet[six,]
      ret <- .C("likelihoodVec", np=as.integer(n), mp=as.integer(ms), 
                rdp=as.double(ds), kvecp=as.double(kvec), loobp=as.double(loob),
                liVecp=as.double(liVec))
      liVec <- ret$liVecp[order(six)]
      # print(paste("liVec[1] = ",toString(liVec[1])))
      return(liVec)
    }

  d <- d[complete.cases(d),]
  if(nrow(d)<2) stop("need at least 2 complete cases to compute A")
    
  all <- concat(partition)

  n<-length(d[,1])
  mp <- length(partition)
  m<-length(d[1,])

  rd <- rwt(d)
  
    nullLogLike <- 0
 
    # a function to be optimized over kernel width.
    # the fuction returns the likelihood of the null model.
    nullScore <- function(tkw){
      nliVec <- rep(1.0,times=n)
      for (j in 1:mp){
        rd1 <- cbind(rd[,partition[[j]]])
        nliVec <- nliVec * getLikeVec(rd1,tkw,1.0)
        # print(paste("j = ",toString(j)," nliVec[1] = ",toString(nliVec[1])))
      }
      nullLogLike <- sum(sapply(nliVec,log))
      # print(paste("tkw = ",toString(tkw)))
      # print(paste("nullLogLike = ",toString(-nullLogLike)))
      return(-nullLogLike)
    }
  
    # get the optimum kernel width for the null model
    okwl <- 0.01
    okwr <- n / 2
    opt <- optimize(f=nullScore,lower=okwl,upper=okwr,tol=0.01,maximum=FALSE)    
    nkw <- opt$minimum
    nullLogLike <- -opt$objective
    # print(paste("nkw = ",toString(nkw)))
    # print(paste("nullLogLike = ",toString(nullLogLike)))
    
    # make one call to get data likelihood vector at optimal kernel widfths.
    nliVec <- rep(1.0,times=n)
    for (j in 1:mp){
      rd1 <- cbind(rd[,partition[[j]]])
      nliVec <- nliVec * getLikeVec(rd1,nkw,1.0)
    }
    
    # a function to be optimized over kernel width and mixture weight.
    # the function returns a mixture likelihood for the alternate model
    wllScore <- function(v) {
      # nkw <- 1 + n * v[1] / (40 - 38 * v[1])
      akw <- 1 + n * v[1] / (40 - 38 * v[1])
      weight <- v[2]
      
      # make one call to getLikeVec to get likelihood vector for the alternate model
      rd1 <- cbind(rd[,all])
      liVec <- getLikeVec(rd1,akw,1.0)

      # now compute weighted log likelihood of null and alt models
      wScore <- sum(sapply(nliVec * weight + liVec * (1-weight),log))
      return(-wScore)
      
    }
  
    # options(warn=-1)
    # get the optimum kernel width and weight for the alternate model
    opt <- nmkb(c(0.5,0.5),wllScore,lower=c(0,0),upper=c(1,1),
                control = list(tol=0.01))   
    fv <- opt$par  
    # recover the optimal parameters
    okw <- 1 + n * fv[1] / (40 - 38 * fv[1])
    optWeight <- fv[2]
    wLogLike <- - opt$value
    
    # LRS is required for computing P values
    LRS <- -2*nullLogLike +2*wLogLike
    
    rawA <- 1 - exp(-(2/n)*(wLogLike-nullLogLike))
    if(rawA<0.0) rawA <- 0.0
    
    # hyperbolic correction
    os <- (1.0 - 1.0/(1.0 + ht * rawA) ) * (n^hp) / hs
    awLogLike <- wLogLike + os 
    Ascore <- 1 - exp(-(2/n)*(awLogLike-nullLogLike))
    if(Ascore<0.0) Ascore <- 0.0
  
    # OLD code for applying the loop (leave one out proportion)  
    # two more calls to compute wLogLike and nullLogLike 
    # with a proportion of loob (leave one out)
    # loob <- lp
    
      # the first call
      # nliVec <- rep(1.0,times=n)
      # for (j in 1:mp){
      #   rd1 <- cbind(rd[,partition[[j]]])
      #   nliVec <- nliVec * getLikeVec(rd1,nkw,loob)
      # }
      # nullLogLike <- sum(sapply(nliVec,log))
  
      # the second call
      # rd1 <- cbind(rd[,all])
      # liVec <- getLikeVec(rd1,okw,loob)
    
      # now compute weighted log likelihood of null and alt models
      # wLogLike <- sum(sapply(nliVec * optWeight + liVec * (1-optWeight),log))
 
      # estimate r^2 from wLogLike and nullLogLike
      # Ascore <- 1 - exp(-(2/n)*(wLogLike-nullLogLike))
    
    # negative A is not allowed
    # Ascore <- abs(Ascore)
    
    
    # print(paste("Ascore = ",toString(Ascore)))

    rdf = list(A=Ascore,rawA=rawA,jointKW=okw,altLL=wLogLike,nullLL=nullLogLike,
              marginalKW=nkw,weight=optWeight,LRstat=LRS,nRows=n,mCols=m,partition=partition)

  return(rdf)
  
}
