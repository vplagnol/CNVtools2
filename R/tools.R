
getQualityScore <- function(x, hyp = 'H0'){

  if (hyp == 'H0') posterior <- x@best.fit.H0
  if (hyp == 'H1') posterior <- x@best.fit.H1
  
  Qvec <- c()
  for (coh in levels (posterior$batch)) {
    posterior2 <- subset(posterior, posterior$batch == coh)
    sds <- tapply(posterior2$signal, FUN=sd, INDEX = posterior2$cn)
    means <- tapply(posterior2$signal, FUN=mean, INDEX = posterior2$cn)
    freq <- table(posterior2$cn)/length(posterior2$cn)

    ########removes some missing categories
    sds <- subset(sds, table(posterior2$cn) > 4)
    means <- subset(means, table(posterior2$cn) > 4)
    freq <- subset(freq, table(posterior2$cn) > 4)


    l <- length(means)

    if (l == 1) {return (NA)} else {   ############an exception for the one component situation
      dmeans <- abs(means[ 1: (l-1) ] - means[ 2:l ])    ##the differences of means for pairs of adjacent clusters
      av.sds <- (freq[1:(l-1)] * sds [ 1:(l-1) ]  + freq[ 2:l ] * sds [ 2:l ])/ ( freq[1:(l-1)] + freq[ 2:l ])   #the average std dev for pairs of adjacent clusters
      weights <- freq[1:(l-1)]*freq[ 2:l ]   ###weights for each pair of clusters
      Q <- sum(weights*dmeans/av.sds)/sum(weights)    ##the quality score
      Qvec <- append(Qvec, values = Q )
    }
  }
  return (Qvec)
}

