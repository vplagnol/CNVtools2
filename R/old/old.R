


######################################################################### the workhorse function
CNV.fitModel <- function(x,
                         nind,
                         hyp = "H0",
                         pi.model = 0,
			 mix.model = 10,
                         force.replace = FALSE,
                         start.mean = NULL,
                         start.var = NULL,
			 control=list(tol=1e-5, max.iter = 3000, min.freq=4)) {

  message('Fitting under ', hyp)
  
  pi.mod <- as.integer(pi.model)
  mix.mod <- as.integer(mix.model)
  logit.offset = rep(0, x@nsamples*x@ncomp)


  ################# starting points
  km.estimate <- FALSE
  if (is.null(start.mean)) {
    km.estimate <- TRUE
  }  

  start.mean <- quantile(x = x@signal, probes = (1:x@ncomp - 0.5)/x@ncomp)
  start.alpha <- rep(1/x@ncomp, x@ncomp)  
  start.var <- rep((sd(x@signal)/x@ncomp)^2, 3)



  if (km.estimate) {
    message('Using k-means estimates as starting point')
    km <- kmeans( test@signal, centers = x@ncomp, nstart = 3)
    my.order <- order(km$centers)
    start.mean <- km$centers [ order(km$centers) ]
    start.var <- (km$withinss / km$size)[ order(km$centers) ]
    start.alpha <- (km$size / sum(km$size))[ order(km$centers) ]
    
                                        #print(x@ncomp)
                                        #print(start.mean)
                                        #print(sqrt(start.var))
                                        #print(start.alpha)    
  }

  
  clean.data <- data.frame(subject = rep(x@covariates$sample.names, 3),
                           strata.association = x@strata.association,
                           strata.mean = x@strata.mean,
                           strata.var = x@strata.var,
                           trait = rep(x@trait, x@ncomp),
                           cn = x@expanded.covariates$cn,
                           signal = rep(x@signal, 3),
                           mean.start = start.mean[ x@expanded.covariates$cn ],
                           alpha.start = start.alpha[ x@expanded.covariates$cn ],
                           var.start = start.var[  x@expanded.covariates$cn ],
                           nu.start = 3,
                           batch = rep(x@covariates$batch) )

###### Here we make sure that the numbering of the strata variance sense
  initial <- sort(unique(clean.data$strata.association))
  clean.data$strata.association <- as.integer(rank(initial)[ match (x = clean.data$strata.association, table = initial) ])

  initial <- sort(unique(clean.data$strata.var))
  clean.data$strata.var <- as.integer(rank(initial)[ match (x = clean.data$strata.var, table = initial) ])

  initial <- sort(unique(clean.data$strata.mean))
  clean.data$strata.mean <- as.integer(rank(initial)[ match (x = clean.data$strata.mean, table = initial) ])

  if(mix.model < 10) {stop("Specification of mix.model incorrect\n")}
  res <- .Call("C_fitmodel", 
	       as.integer(x@ncomp), 
               as.integer(x@nsamples), 
	       hyp, 
	       clean.data, 
	       logit.offset,
               as.matrix(x@design.matrix.mean[, -1]), 
               as.matrix(x@design.matrix.var[,-1]), 
               as.matrix(x@design.matrix.association[, -1]), 
               control,
               mix.mod,	
               pi.mod)

  
  new.data <- res[[1]]
  clean.data$posterior <- new.data[,1]
  clean.data$mean <- new.data[,2]
  clean.data$var <- new.data[,3]
  clean.data$pr <- new.data[,4]
  clean.data$alpha <- new.data[,5]
  clean.data$pdc <- new.data[,6]
  clean.data$nu <- new.data[,7]	

  my.params <-  getparams(clean.data)
  my.params[['convergence']] <- res[[2]]
  replace <- force.replace
  
  if (hyp == 'H0') {
    if ( length(x@mle.H0) == 0 ) replace = TRUE else {
      if (my.params$lnL >  x@mle.H0$lnL || (my.params$convergence == 'C' &&  x@mle.H0$convergence != 'C')) replace = TRUE
    }
    
    if (replace) {
      x@best.fit.H0 <- compact.data.frame(clean.data)
      x@mle.H0 <- my.params
    }    
  }

   if (hyp == 'H1') {
    if ( length(x@mle.H1) == 0 ) replace = TRUE else {
      if (my.params$lnL >  x@mle.H1$lnL || (my.params$convergence == 'C' &&  x@mle.H1$convergence != 'C')) replace = TRUE
    }
    
    if (replace) {
      x@best.fit.H1 <- compact.data.frame(clean.data)
      x@mle.H1 <- my.params
    }    
  }
    
  return(x)
}
#####################################################################
