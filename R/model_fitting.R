setGeneric("expand", def = function(object)
           standardGeneric('expand'))


setGeneric("fit", def = function(object,
                    ncomp,
                    hyp = 'H0',
                    EM.starting.point = 'kmeans',
                    gaussian.or.T = 'gaussian',
                    model.association.H0 = NULL,
                    model.association.H1 = NULL,
                    model.mean = NULL,
                    model.var = NULL,
                    start.mean = NULL,
                    start.var = NULL,
                    start.alpha = NULL,
                    control = list(tol = 1e-5, max.iter = 500, min.freq = 4, logP.outliers = 8, nstart.kmeans = 40),
                    force.replace = FALSE)
           standardGeneric('fit'))



#############################################################################
compact.data.frame <- function (full.frame) 
{
  full.frame <- full.frame[order(full.frame$cn),]
  full.frame.mod <- do.call(rbind.data.frame, split(x = full.frame$posterior, f = full.frame$IID))
  names(full.frame.mod) <- paste('P', c(1:(dim(full.frame.mod)[2])), sep='')       
  full.frame.mod$cn <-apply(full.frame.mod, FUN=which.max, MAR=1)    
  full.frame.mod$IID <- row.names(full.frame.mod)                              
  full.frame <- subset( full.frame[ ,c('IID', 'batch', 'signal', 'trait', 'proba.not.outlier')], full.frame$cn == 1)
  full.frame <- merge(full.frame, full.frame.mod)        
  return (full.frame)
}


getparams <- function(d)
{
  p <- list()
  p[["ns"]] <- length( levels(d$batch) )
  p[["nc"]] <- range(d$cn)[2]
  p[["nind"]] <- dim(d)[1]/p[["nc"]]

  maxLike <- tapply(d$pr, d$IID, max)  ##takes the max likelihood
  p[["lnL"]] <- sum ( maxLike + log( tapply(exp(d$pr - maxLike[ d$IID]), FUN=sum, d$IID)) )
  
  p[["alpha"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  p[["mean"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  p[["var"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  p[["nu"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  p[["pdc"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  
  lev <- levels(d$batch)
  for(j in 1:p$ns) {
    for(i in 1:p$nc) {
      p$mean[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$mean)
      p$alpha[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$alpha)
      p$var[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$var)
      p$nu[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$nu)
      p$pdc[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$pdc)
    }	
  }
  return (p)
}

#############################################################################
setMethod("expand", "CNVtools", function(object) {
  
  object@expanded.covariates <- data.frame() 
  for (i in 1:object@ncomp) object@expanded.covariates <- rbind.data.frame ( object@expanded.covariates, object@covariates)
  object@expanded.covariates$cn <- sort(rep( 1:object@ncomp, object@nsamples))


  ########################### Now build the design matrices for the clustering model
  model.mean <- object@model.mean
  model.var <- object@model.var
  model.association.H0 <- object@model.association.H0
  model.association.H1 <- object@model.association.H1
  #object@design.matrix.association <- model.matrix(data = object@expanded.covariates, object = model.association)
  
  special <- c("strata")
    
  for (design in c('var', 'mean', 'association.H0', 'association.H1')) {
    
    my.formula <- as.formula(get(paste('model', design, sep = '.')))
    Terms <- terms(my.formula, special, data= object@expanded.covariates)
    strats <- attr(Terms, "specials")$strata
    
    if (!is.null(strats)) {
      if (design %in% c('association.H0', 'association.H1')) stop('No strata term allowed in the association test')
      m <- list()
      m[[1]] <- as.name("model.frame")
      m[[2]] <- Terms
      names(m)[2] <- "formula"
      m[['data']] <-  object@expanded.covariates
      m <- as.call(c(as.list(m), list(na.action=as.symbol("na.omit"))))
      m <- eval(m)
      
      temps <- untangle.specials(Terms, "strata", 1)
      assign(x = paste('strata', design, sep = '.'), as.integer(strata(m[, temps$vars], shortlabel=TRUE)))
      nstrata <- length(unique(get(paste('strata', design, sep = '.'))))
      
      ########## here I essentially make the design matrices equal to 0 if strata is included. I should do better
      if (nstrata == 1) {
        assign(x = paste('design.matrix', design, sep = '.'), value = model.matrix(data = object@expanded.covariates, object = as.formula(' ~ 1') ))
      } else assign(x = paste('design.matrix', design, sep = '.'), value = matrix(0))      
    } else {
      assign(x = paste('design.matrix', design, sep = '.'), value = model.matrix(data = object@expanded.covariates, object = my.formula ))
      assign(x = paste('strata', design, sep = '.'), value = as.integer(rep(1, object@ncomp*object@nsamples)))
    }
  }

  object@design.matrix.mean <- design.matrix.mean
  object@design.matrix.var <- design.matrix.var
  object@design.matrix.association.H0 <- design.matrix.association.H0
  object@design.matrix.association.H1 <- design.matrix.association.H1

  object@strata.mean <- strata.mean
  object@strata.var <- strata.var

  return (object)
})


#############################################################################
setMethod("fit", "CNVtools", function(object,
                                      ncomp,
                                      hyp,
                                      EM.starting.point,
                                      gaussian.or.T,
                                      model.association.H0,
                                      model.association.H1,
                                      model.mean,
                                      model.var,
                                      start.mean,
                                      start.var,
                                      start.alpha,
                                      control,
                                      force.replace) {
  message('Fitting under ', hyp)

  object@ncomp <- ncomp
  pi.model <- as.integer(0)
  mix.model <- as.integer(10)
  logit.offset = rep(0, object@nsamples*object@ncomp)

  if (!is.null(model.mean)) {
    if (object@model.mean != model.mean) force.replace <- TRUE
    object@model.mean <- model.mean    
  }
  
  if (!is.null(model.var)) {
    if (object@model.var != model.var) force.replace <- TRUE
    object@model.var <- model.var
  }
  
  if (!is.null(model.association.H0)) {
    if (object@model.association.H0 != model.association.H0) force.replace <- TRUE
    object@model.association.H0 <- model.association.H0
  }

  if (!is.null(model.association.H1)) {
    if (object@model.association.H1 != model.association.H1) force.replace <- TRUE
    object@model.association.H1 <- model.association.H1
  }
  
  object <- expand(object)
  
  
################# starting points
  km.estimate <- FALSE
  manual <- FALSE
  basic <- FALSE
  H1.from.H0 <- FALSE
  
  if (EM.starting.point == 'basic') basic <- TRUE
  if (EM.starting.point == 'kmeans') km.estimate <- TRUE
  if (EM.starting.point == 'manual') manual <- TRUE
  if (EM.starting.point == 'H1.from.H0') {
    if ((hyp == 'H1') && (length(object@mle.H0) > 0)) H1.from.H0 <- TRUE
    if ((hyp == 'H1') && (length(object@mle.H0) == 0)) {stop('You selected H1.from.H0 as a starting point but this CNVtools object does not seem to have a fitted H0')}
    if (hyp == 'H0') {stop('The EM starting point H1.from.H0 only makes sense if the hypothesis is H1')}
  }
    
#####
  if (basic) {
    start.mean <- quantile(x = object@signal, probes = (1:object@ncomp - 0.5)/object@ncomp)
    start.alpha <- rep(1/object@ncomp, object@ncomp)  
    start.var <- rep((sd(object@signal)/object@ncomp)^2, object@ncomp)
  }
  
  if (manual) {
    if (is.null(start.mean)) stop("Starting means must be provided")
    start.alpha <- rep(1/object@ncomp, object@ncomp)  
    start.var <- rep((sd(object@signal)/object@ncomp)^2, object@ncomp)
  }

  if (km.estimate) {
    nstart <- ifelse ( 'nstart.kmeans' %in% names(control), control$nstart.kmeans, 40)
    message('Using k-means estimates as starting point of the EM algorithm, using ', nstart, ' starting points')
    km <- kmeans( object@signal, centers = object@ncomp, nstart = nstart)
    my.order <- order(km$centers)
    start.mean <- km$centers [ order(km$centers) ]
    start.var <- (km$withinss / km$size)[ order(km$centers) ]
    start.alpha <- (km$size / sum(km$size))[ order(km$centers) ]
    message('Done with k-means')
  }

  if (H1.from.H0) {
    start.mean <- as.numeric(object@mle.H0$mean)
    start.var <- as.numeric(object@mle.H0$var)
    start.alpha <- as.numeric(object@mle.H0$alpha)
  }

  clean.data <- data.frame(IID = rep(object@IID, object@ncomp),
                           strata.mean = object@strata.mean,
                           strata.var = object@strata.var,
                           trait = as.numeric(rep(object@trait, object@ncomp)),
                           cn = object@expanded.covariates$cn,
                           signal = rep(object@signal, object@ncomp),
                           mean.start = start.mean[ object@expanded.covariates$cn ],
                           alpha.start = start.alpha[ object@expanded.covariates$cn ],
                           var.start = start.var[  object@expanded.covariates$cn ],
                           nu.start = 3,
                           batch = rep(object@covariates$batch) )

###### Here we make sure that the numbering of the strata variance sense
  initial <- sort(unique(clean.data$strata.var))
  clean.data$strata.var <- as.integer(rank(initial)[ match (x = clean.data$strata.var, table = initial) ])

  initial <- sort(unique(clean.data$strata.mean))
  clean.data$strata.mean <- as.integer(rank(initial)[ match (x = clean.data$strata.mean, table = initial) ])

  if(mix.model < 10) {stop("Specification of mix.model incorrect\n")}


  if (hyp == 'H0') {
    res <- .Call("C_fitmodel", 
                 as.integer(object@ncomp), 
                 as.integer(object@nsamples), 
                 hyp, 
                 clean.data, 
                 logit.offset,
                 as.matrix(object@design.matrix.mean[, -1]), 
                 as.matrix(object@design.matrix.var[,-1]), 
                 as.matrix(object@design.matrix.association.H0), 
                 control,
                 mix.model,	
                 pi.model)
  }

  
  if (hyp == 'H1') {
    res <- .Call("C_fitmodel", 
                 as.integer(object@ncomp), 
                 as.integer(object@nsamples), 
                 hyp, 
                 clean.data, 
                 logit.offset,
                 as.matrix(object@design.matrix.mean[, -1]), 
                 as.matrix(object@design.matrix.var[,-1]), 
                 as.matrix(object@design.matrix.association.H1), 
                 control,
                 mix.model,	
                 pi.model)
  }

  

  
  
  new.data <- res[[1]]
  clean.data$posterior <- new.data[,1]
  clean.data$mean <- new.data[,2]
  clean.data$var <- new.data[,3]
  clean.data$pr <- new.data[,4]
  clean.data$alpha <- new.data[,5]
  clean.data$pdc <- new.data[,6]
  clean.data$nu <- new.data[,7]
  clean.data$proba.not.outlier <- new.data[,8]	

  
  my.params <-  getparams(clean.data)
  my.params[['convergence']] <- res[[2]]
  replace <- force.replace
  
  if (hyp == 'H0') {
    if ( length(object@mle.H0) == 0 ) replace = TRUE else {
      if (my.params$lnL >  object@mle.H0$lnL || (my.params$convergence == 'C' &&  object@mle.H0$convergence != 'C')) replace = TRUE
    }
    
    if (replace) {
      object@best.fit.H0 <- compact.data.frame(clean.data)
      object@best.fit.H0 <- object@best.fit.H0 [ match( object@covariates$IID, table = object@best.fit.H0$IID), ]
      object@mle.H0 <- my.params
      object@cn.calls <- object@best.fit.H0$cn
    }    
  }

   if (hyp == 'H1') {
    if ( length(object@mle.H1) == 0 ) replace = TRUE else {
      if (my.params$lnL >  object@mle.H1$lnL || (my.params$convergence == 'C' &&  object@mle.H1$convergence != 'C')) replace = TRUE
    }
    
    if (replace) {
      object@best.fit.H1 <- compact.data.frame(clean.data)
      object@best.fit.H1 <- object@best.fit.H1 [ match( object@covariates$IID, table = object@best.fit.H1$IID), ]
      object@mle.H1 <- my.params
    }    
  }
    
  return(object)
  
})



