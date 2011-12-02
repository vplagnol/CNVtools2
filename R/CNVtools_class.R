library(survival)
source('R/model_fitting.R')


setClass("CNVtools",
         representation(ncomp = 'numeric',
                        nprobes = 'numeric',
                        nsamples = 'numeric',
                        signal = "numeric",
                        trait = 'numeric',
                        batch = 'factor',                      
                        IID = 'character',
                        FID = 'factor',
                        f.member = 'numeric',
                        covariates = 'data.frame',
                        expanded.covariates = 'data.frame',  ## the matrix expanded for each copy number
                        multiprobe.signal = 'matrix',
                        cn.calls = "numeric",
                        model.association = "formula",
                        model.mean = "formula",
                        model.var = "formula",
                        design.matrix.association = 'matrix',
                        design.matrix.mean = 'matrix',
                        design.matrix.var = 'matrix',
                        strata.mean = 'integer',
                        strata.var = 'integer',
                        binary.trait = 'logical',
                        best.fit.H0 = 'data.frame',
                        mle.H0 = 'list',
                        best.fit.H1 = 'data.frame',
                        mle.H1 = 'list'))



validCNVtools <- function( object ) {
  if (length(object@signal) != object@nsamples) {message("Lengths of signal and nsamples do not match"); stop()}
  if (length(object@trait) != object@nsamples) {message("Lengths of trait and nsamples do not match"); stop()}
  if (length(object@batch) != object@nsamples) {message("Lengths of batch and nsamples do not match"); stop()}
  
  if ((nrow(object@covariates) > 0) && (nrow(object@covariates) != object@nsamples)) {message("The covariates matrix does not have the right size"); stop()}
}

setValidity("CNVtools", validCNVtools)



#############################################################################
setMethod("initialize", "CNVtools", function(.Object,
                                             signal,
                                             trait,
                                             batch,
                                             model.association = formula(' ~ cn'),
                                             model.mean = formula ('~ strata(cn)'),
                                             model.var = formula ('~ strata(cn)'),
                                             IID = NULL,
                                             FID = NULL,
                                             f.member = NULL,  ##1 father, 2 mother, after that children
                                             covariates = data.frame()) {

  .Object@nsamples = length(signal)
  
  if (is.numeric(batch)) batch <- factor(as.character(batch))
  if (is.character(batch)) batch <- factor(batch)

  
  if (is.character(model.association)) model.association <- as.formula(model.association)
  if (is.character(model.mean)) model.mean <- as.formula(model.mean)
  if (is.character(model.var)) model.var <- as.formula(model.var)
  
  .Object@model.association <- model.association
  .Object@model.mean <- model.mean
  .Object@model.var <- model.var

  
  .Object@signal = as(signal, 'numeric')
  .Object@trait = as(trait, 'numeric')

  
  .Object@batch = factor(batch)


  if (is.null(IID)) .Object@IID <- paste('sample', 1:.Object@nsamples, sep = '.') else .Object@IID = as(IID, 'character')
  
  
  sd.signal <-  sd(.Object@signal)
  if ( sd.signal > 0.3 ||   sd.signal < 1.5 ) {
    message('The standard deviation of the signal is quite different from 1. For numeric stability purpose it will be rescaled to sd = 1.')
    .Object@signal <- .Object@signal / sd.signal
  }

  
  if ( sum(! .Object@trait %in% c(0,1)) == 0 ) .Object@binary.trait <- TRUE else .Object@binary.trait <- FALSE
  

  
  if (nrow(covariates) > 0) {
    .Object@covariates = covariates    
    .Object@covariates$Y <- .Object@trait
    .Object@covariates$batch <- .Object@batch
    .Object@covariates$IID  <- .Object@IID
  } else {
    .Object@covariates <- data.frame(IID = .Object@IID,
                                     Y = .Object@trait,
                                     batch = .Object@batch)
  }


  ############# Now the family aspect
  if (!is.null(FID)) {
    if (is.null(f.member)) stop('If the FID is included one must also add the f.member field')
    if (is.numeric(FID) || is.character(FID)) FID <- factor(FID)
    .Object@f.member <- f.member
    .Object@FID <- FID
    .Object@covariates$FID <- FID     
  }


  
  validCNVtools ( .Object )
  .Object
})
                                               
setMethod("show", "CNVtools", function(object) {
  cat('Number of samples:', object@nsamples, '\n\n')

  message('Model for the means of the intensity clusters')
  print(object@model.mean,  showEnv = FALSE)

  message('Model for the variances of the intensity clusters')
  print(object@model.var, showEnv = FALSE)

  message('Model for the association')  
  print(object@model.association, showEnv = FALSE)

  if (object@binary.trait) {
    message('Binary trait, here is the case control distribution per batch.')
    print(table(object@trait, object@batch))
  }
  
})




#############################################################################
setGeneric("TDT" ,function(.Object){standardGeneric("TDT")}) 

setMethod("TDT", "CNVtools", function(.Object) {  ##Intensity based TDT test

  small.cov <- data.frame ( signal = .Object@signal, FID = .Object@FID, f.member = .Object@f.member)

  
  fathers <- subset(small.cov[, c('FID', 'signal')], small.cov$f.member == 1)
  if (max(table(fathers$FID)) > 1) stop('There is a family with 2 or more fathers')
  names(fathers)[2] <- 'signal.father'

  mothers <- subset(small.cov[, c('FID', 'signal')], small.cov$f.member == 2)
  names(mothers)[2] <- 'signal.mother'    
  if (max(table(mothers$FID)) > 1) stop('There is a family with 2 or more mothers')
  
  small.cov <- merge( small.cov, fathers, by = 'FID', all.x = TRUE)
  small.cov <- merge( small.cov, mothers, by = 'FID', all.x = TRUE)
  small.cov <- subset(small.cov, f.member >= 3)
  
  my.tab <- table(0.5*(small.cov$signal.father + small.cov$signal.mother) > small.cov$signal)  
  return (chisq.test(my.tab))
})
 

#############################################################################
setGeneric("apply.pca",function(.Object, multiprobe.signal = NULL){standardGeneric("apply.pca")}) 

setMethod("apply.pca", "CNVtools", function(.Object, multiprobe.signal) {
 
  if (!is.null(multiprobe.signal)) {
    if (class(multiprobe.signal) != 'matrix') stop('Multi-probe signal must be a matrix')
    if (nrow(multiprobe.signal) != .Object@nsamples) stop('Number of rows of multiprobe signal must be equal to the number of samples')
    .Object@multiprobe.signal <- multiprobe.signal
  }
  
  pca <- prcomp(.Object@multiprobe.signal, scale = TRUE)$x[,1]
  pca <- pca/sd(pca)
  if (cor(pca, .Object@signal) < 0) .Object@signal <- - pca else .Object@signal <- pca

  ######### clear the previously fitted objects that are not relevant anymore
  .Object@mle.H0 <- list()
  .Object@mle.H1 <- list()
  .Object@best.fit.H0 <- data.frame()
  .Object@best.fit.H1 <- data.frame()
  
  return(.Object)
})

#############################################################################
setGeneric("plot.cnv",function(.Object, hist.or.dens = "histogram", hyp = 'H0', batch = NULL, subset = NULL, freq = NULL, plot.outliers = FALSE, ... ){standardGeneric("plot.cnv")}) 

setMethod("plot.cnv", "CNVtools", function(.Object, hist.or.dens, hyp, batch, subset, freq, ...) {

  if (!is.null(subset) && ( (length(subset) != .Object@nsamples) || class(subset) != 'logical')) stop('If specified, subset must be a logical vector of length equal to the number of samples')
  if (!hyp %in% c('H0', 'H1') || length(hyp) != 1) stop('The argument hyp must of of length 1 and equal to H0 or H1')
  if (!is.null(batch) & !is.null(subset))  stop()
  
  
  if (hyp == 'H0') posterior <- .Object@best.fit.H0
  if (hyp == 'H1') posterior <- .Object@best.fit.H1
  
  if (!is.null(batch))  posterior <- posterior[posterior$batch %in% batch, ]
  if (!is.null(subset))  posterior <- posterior[posterior$batch %in% batch, ]
  
  posterior <- .Object@best.fit.H0
  
  posterior <- posterior[order(posterior$signal), ]
  if (hist.or.dens == "density") {
    dens <- density(posterior$signal)
    plot(dens, ...)
    my.max <- max(dens$y)
  }
  if (hist.or.dens == "histogram") {
    my.hist <- hist(posterior$signal, freq = freq, ...)
    my.max <- max(my.hist$counts)
    if (!is.null(freq)) {if (freq == FALSE) my.max <- max(my.hist$density)}
  }
  
  col <- 1
  ncomp <- max(posterior$cn)
  for (i in c(1:.Object@ncomp)) {
    if (!plot.outliers) pr <- ifelse ( posterior$proba.not.outlier < 0.1, 0, posterior[, paste("P", i, sep = "")])
    if ( plot.outliers) pr <- posterior[, paste("P", i, sep = "")]

    lines(posterior$signal, my.max*pr, type = "l", col = col)
    col <- col + 2
  }
  
})





CNVtools.multivariate.signal <- function (multivariate.signal,
                                          trait,
                                          batch,
                                          model.mean = formula ('~ strata(cn)'),
                                          model.var = formula ('~  strata(cn)'),
                                          model.association = formula ('~  cn'),
                                          IID = NULL,
                                          FID = NULL,
                                          f.member = NULL,
                                          covariates = data.frame()) {

  if(class(multivariate.signal) != 'matrix') stop('Multivariate input signal must be a matrix')
  
  mean.signal <-  apply(multivariate.signal, MAR = 1, FUN = mean, na.rm = TRUE)
  res <- new('CNVtools',
             trait = trait,
             signal = mean.signal,
             batch = batch,
             IID = IID,
             FID = FID,
             f.member = f.member,
             model.mean = model.mean,
             model.var = model.var,
             model.association = model.association)
  
  res@multiprobe.signal <- multivariate.signal
  return(res)
}
  
