library(survival)




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
                        sex = 'numeric',
                        covariates = 'data.frame',
                        expanded.covariates = 'data.frame',  ## the matrix expanded for each copy number
                        family.frame = 'data.frame',
                        family.test = 'list',
                        multiprobe.signal = 'matrix',
                        probe.correlations = 'matrix',  ##stores the correlations between parents and offsprings: row 1 mid parents, row 2 fathers/sons, row 3 mothers/daughters
                        signal.correlations = 'numeric', ##stores the correlations parents/offsprings using the one dimensional summary data (3 numbers)
                        cn.calls = "numeric",
                        model.association.H0 = "formula",
                        model.association.H1 = "formula",
                        model.mean = "formula",
                        model.var = "formula",
                        design.matrix.association.H0 = 'matrix',
                        design.matrix.association.H1 = 'matrix',
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
                                             model.association.H0 = formula(' ~ 1'),
                                             model.association.H1 = formula(' ~ cn'),
                                             model.mean = formula ('~ strata(cn)'),
                                             model.var = formula ('~ strata(cn)'),
                                             IID = NULL,
                                             FID = NULL,
                                             f.member = NULL,  ##1 father, 2 mother, after that 3 for affected offspring
                                             sex = NULL,
                                             covariates = data.frame()) {

  .Object@nsamples = length(signal)

  if (missing(batch)) batch <- rep('A',  .Object@nsamples)
  if (is.numeric(batch)) batch <- factor(as.character(batch))
  if (is.character(batch)) batch <- factor(batch)

  
  if (is.character(model.association.H0)) model.association <- as.formula(model.association.H0)
  if (is.character(model.association.H1)) model.association <- as.formula(model.association.H1)
  if (is.character(model.mean)) model.mean <- as.formula(model.mean)
  if (is.character(model.var)) model.var <- as.formula(model.var)
  
  .Object@model.association.H0 <- model.association.H0
  .Object@model.association.H1 <- model.association.H1
  .Object@model.mean <- model.mean
  .Object@model.var <- model.var

  
  .Object@signal = as(signal, 'numeric')
  .Object@trait = as(trait, 'numeric')
  .Object@batch = factor(batch)

  
  if (is.null(IID)) .Object@IID <- paste('sample', 1:.Object@nsamples, sep = '.') else {
    if (length(IID) != .Object@nsamples) stop('Nb of samples seems inconsistent between the vectors signal and IID')
    .Object@IID = as(IID, 'character')
  }
  
  if (is.null(sex)) .Object@sex <- as.numeric(rep(NA, .Object@nsamples)) else {
    if (!class(sex) %in% c('integer', 'numeric')) stop('Sex must be a numeric variable')
    .Object@sex = sex
  }

  if (sum(! .Object@sex %in% c(NA, 1, 2))) stop('Sex must be either NA, 1 (male) or 2 (female)')
  
  sd.signal <-  sd(.Object@signal)
  if ( sd.signal < 0.3 ||   sd.signal > 1.5 ) {
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




CNVtools.multivariate.signal <- function (multivariate.signal,
                                          trait,
                                          batch = NULL,
                                          model.mean = formula ('~ strata(cn)'),
                                          model.var = formula ('~  strata(cn)'),
                                          model.association.H0 = formula ('~  1'),
                                          model.association.H1 = formula ('~  cn'),
                                          IID = NULL,
                                          FID = NULL,
                                          f.member = NULL,
                                          sex = NULL,
                                          covariates = data.frame()) {

  if(class(multivariate.signal) != 'matrix') stop('Multivariate input signal must be a matrix')
  if (is.null(batch)) batch <- rep('A', nrow(multivariate.signal))
  
  mean.signal <-  apply(multivariate.signal, MAR = 1, FUN = mean, na.rm = TRUE)

  res <- new('CNVtools',
             trait = trait,
             signal = mean.signal,
             batch = batch,
             IID = IID,
             FID = FID,
             f.member = f.member,
             sex = sex,
             model.mean = model.mean,
             model.var = model.var,
             model.association.H0 = model.association.H0,
             model.association.H1 = model.association.H1)


  res@nprobes <- ncol(multivariate.signal)
  res@multiprobe.signal <- multivariate.signal
  return(res)
}
  



##########################################################################
setGeneric("familial.scaling",function(.Object){standardGeneric("familial.scaling")})


setMethod("familial.scaling", "CNVtools", function(.Object) {

  if (is.null(.Object@multiprobe.signal)) stop('We need the probe level intensity for this function to make sense')
  probe.names <- dimnames(.Object@multiprobe.signal)[[2]]

  .Object@family.frame <- data.frame ( signal = .Object@signal, FID = .Object@FID, f.member = .Object@f.member, status = .Object@trait, sex = .Object@sex)
  .Object@family.frame <- cbind.data.frame( .Object@family.frame, .Object@multiprobe.signal)

  fathers <- subset(.Object@family.frame[, c('FID', probe.names)], .Object@family.frame$f.member == 1)
  if (max(table(fathers$FID)) > 1) stop('There is a family with 2 or more fathers')
  names(fathers)[2:(1+length(probe.names)) ] <- paste('father', probe.names, sep = '.')

  mothers <- subset(.Object@family.frame[, c('FID', probe.names)], .Object@family.frame$f.member == 2)
  if (max(table(mothers$FID)) > 1) stop('There is a family with 2 or more mothers')
  names(mothers)[2:(1+length(probe.names)) ] <- paste('mother', probe.names, sep = '.')

  .Object@family.frame <- merge( .Object@family.frame, fathers, by = 'FID', all.x = TRUE)
  .Object@family.frame <- merge( .Object@family.frame, mothers, by = 'FID', all.x = TRUE)
  .Object@family.frame <- subset( .Object@family.frame, f.member >= 3)

  .Object@probe.correlations <- matrix(nrow = 3, ncol = .Object@nprobes, dimnames = list(c('midparent', 'men', 'women'), probe.names))
  
  
  for (i in 1:.Object@nprobes) {
    mid.parent <- 0.5*( .Object@family.frame[, paste('father', probe.names[i], sep = '.')] + .Object@family.frame[, paste('mother', probe.names[i], sep = '.')])
    .Object@probe.correlations[1,i] <- cor(  .Object@family.frame[, probe.names[i]],  mid.parent, use = 'pairwise.complete.obs')

    men <- subset( .Object@family.frame, sex == 1)
    .Object@probe.correlations[2,i] <- cor(  men[, probe.names[i]], men[, paste('father', probe.names[i], sep = '.')] , use = 'pairwise.complete.obs')
    
    women <- subset( .Object@family.frame, sex == 2)
    .Object@probe.correlations[3,i] <- cor(  women[, probe.names[i]], women[, paste('father', probe.names[i], sep = '.')] , use = 'pairwise.complete.obs')                                           
  }

  if (sum(.Object@probe.correlations[1,] > 0) > 0) {
    scaling <- .Object@probe.correlations[1,] / max(.Object@probe.correlations[1,])
    scaling <- pmax (scaling, 0.01) ##sorts out negative values
  } else {
    scaling <- rep(1, test@nprobes)
  }
  message('Scaling factor for the probe data: ')
  print(scaling)

  message('Now renormalizing the intensity data: ')
  for (i in 1:.Object@nprobes) {.Object@multiprobe.signal[,i] <- scaling[ i ] * .Object@multiprobe.signal[, i]}

  message('Now computing the mean intensity data: ')
  .Object@signal <-  apply(.Object@multiprobe.signal, MAR = 1, FUN = mean, na.rm = TRUE)

  
  return(.Object)
})





setMethod("show", "CNVtools", function(object) {
  cat('Number of samples:', object@nsamples, '\n\n')

  if (sum(!is.na(object@sex)) > 0) {
    message('Sex distribution (rows, 1 is male and 2 is female) across family member status (columns, 1 father, 2 mother, 3 offspring)')
    print(table(object@sex, object@f.member))
  }
  
  message('Model for the means of the intensity clusters')
  print(object@model.mean,  showEnv = FALSE)

  message('Model for the variances of the intensity clusters')
  print(object@model.var, showEnv = FALSE)

  message('Model for the association')  
  print(object@model.association.H0, showEnv = FALSE)
  print(object@model.association.H1, showEnv = FALSE)

  if (object@binary.trait) {
    message('Binary trait, here is the case control distribution per batch.')
    print(table(object@trait, object@batch))
  }
  
})




#############################################################################
setGeneric("FamilyTest" ,function(.Object, type = 'FBAT', chromosome = 'autosome'){standardGeneric("FamilyTest")}) 

setMethod("FamilyTest", "CNVtools", function(.Object, type, chromosome ) {  ##Intensity based TDT test

  if (! type %in% c('TDT', 'FBAT', 'binomial')) stop('Type of TDT test must be either TDT, FBAT or binomial')
  
  .Object@family.frame <- data.frame ( signal = .Object@signal, FID = .Object@FID, f.member = .Object@f.member, status = .Object@trait, sex = .Object@sex)

  if (length(.Object@cn.calls) == nrow(.Object@family.frame)) { ##do we have genotype calls?
    .Object@family.frame$cn.calls <- .Object@cn.calls
  } else {
    .Object@family.frame$cn.calls <- rep(NA, .Object@n.samples)
    if (type == 'TDT') {
      message('Error: Cannot run a TDT test of the genotypes have not been called')
      stop()
    }
  }

  if (type == 'TDT' && .Object@ncomp != 3) {
    message('Error: TDT can only be ran on CNVs with 3 components (i.e. di-allelic)')
    stop()
  }
  
  my.sd <- sd(.Object@family.frame$signal)
  
  fathers <- subset(.Object@family.frame[, c('FID', 'signal', 'cn.calls')], .Object@family.frame$f.member == 1)
  if (max(table(fathers$FID)) > 1) stop('There is a family with 2 or more fathers')
  names(fathers)[2] <- 'signal.father'
  names(fathers)[3] <- 'cn.calls.father'

  mothers <- subset(.Object@family.frame[, c('FID', 'signal', 'cn.calls')], .Object@family.frame$f.member == 2)
  names(mothers)[2] <- 'signal.mother'
  names(mothers)[3] <- 'cn.calls.mother'
  if (max(table(mothers$FID)) > 1) stop('There is a family with 2 or more mothers')
  
  .Object@family.frame <- merge( .Object@family.frame, fathers, by = 'FID', all.x = TRUE)
  .Object@family.frame <- merge( .Object@family.frame, mothers, by = 'FID', all.x = TRUE)
  .Object@family.frame <- subset(.Object@family.frame, f.member >= 3 & (status == 1))
 

############ TDT test
  if (type == 'TDT') {
    .Object@family.frame$TDT.score <- 0
    
###if both parents are hets
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 2 & .Object@family.frame$cn.calls.mother == 2 & .Object@family.frame$cn.calls == 3, 2, .Object@family.frame$TDT.score)
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 2 & .Object@family.frame$cn.calls.mother == 2 & .Object@family.frame$cn.calls == 1, -2, .Object@family.frame$TDT.score)

###### if the mother is het
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 3 & .Object@family.frame$cn.calls.mother == 2 & .Object@family.frame$cn.calls == 2, -1, .Object@family.frame$TDT.score)
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 3 & .Object@family.frame$cn.calls.mother == 2 & .Object@family.frame$cn.calls == 3, 1, .Object@family.frame$TDT.score)
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 1 & .Object@family.frame$cn.calls.mother == 2 & .Object@family.frame$cn.calls == 1, -1, .Object@family.frame$TDT.score)
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 1 & .Object@family.frame$cn.calls.mother == 2 & .Object@family.frame$cn.calls == 2, 1, .Object@family.frame$TDT.score)
    
###### if the father is het
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 2 & .Object@family.frame$cn.calls.mother == 3 & .Object@family.frame$cn.calls == 2, -1, .Object@family.frame$TDT.score)
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 2 & .Object@family.frame$cn.calls.mother == 3 & .Object@family.frame$cn.calls == 3, 1, .Object@family.frame$TDT.score)
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 2 & .Object@family.frame$cn.calls.mother == 1 & .Object@family.frame$cn.calls == 1, -1, .Object@family.frame$TDT.score)
    .Object@family.frame$TDT.score <- ifelse ( .Object@family.frame$cn.calls.father == 2 & .Object@family.frame$cn.calls.mother == 1 & .Object@family.frame$cn.calls == 2, 1, .Object@family.frame$TDT.score)

    .Object@family.frame$incompatible <-  (.Object@family.frame$cn.calls.father == 3 &  .Object@family.frame$cn.calls.mother %in% c(2, 3) & .Object@family.frame$cn.calls == 1) |
    (.Object@family.frame$cn.calls.father == 1 &  .Object@family.frame$cn.calls.mother %in% c(1, 2) & .Object@family.frame$cn.calls == 3) | 
    (.Object@family.frame$cn.calls.mother == 1 &  .Object@family.frame$cn.calls.father %in% c(1, 2) & .Object@family.frame$cn.calls == 3) |
    (.Object@family.frame$cn.calls.mother == 3 &  .Object@family.frame$cn.calls.father %in% c(2, 3) & .Object@family.frame$cn.calls == 1)

    
    if (sum( .Object@family.frame$incompatible, na.rm = TRUE) > 0) {
      message('Warning: Some genotype are incompatible between parents and offsprings. See table below:')
      print(subset( .Object@family.frame,  .Object@family.frame$incompatible))      
    } else {message('Good news: all genotypes are compatible between parents and offsprings')}


    my.tab <- c( 2 *sum(.Object@family.frame$TDT.score == -2, na.rm = TRUE) + sum(.Object@family.frame$TDT.score == -1, na.rm = TRUE),
                 2 *sum(.Object@family.frame$TDT.score == 2, na.rm = TRUE) + sum(.Object@family.frame$TDT.score == 1, na.rm = TRUE))
    my.test <- chisq.test(my.tab)

    .Object@family.test <- list(transmissions = my.tab,
                             n.compatible.transmissions = sum( !.Object@family.frame$incompatible, na.rm = TRUE),
                             n.incompatible.transmissions = sum( .Object@family.frame$incompatible, na.rm = TRUE),
                             statistic = my.test$statistic,
                             p.value = my.test$p.value)
  }

############### binomial test  
  if (type == 'binomial') {

    .Object@family.frame$U <- .Object@family.frame$signal - 0.5*(.Object@family.frame$signal.father + .Object@family.frame$signal.mother)
    my.tab <- table( .Object@family.frame$U > 0 )
    
    my.test <- chisq.test(my.tab)
    .Object@family.test <- list(transmissions = my.tab,
                             statistic = my.test$statistic,
                             p.value = my.test$p.value)
  }

############### FBAT test
  if (type == 'FBAT') {

    if (chromosome != 'X') {
      .Object@family.frame$U <- .Object@family.frame$signal - 0.5*(.Object@family.frame$signal.father + .Object@family.frame$signal.mother)
      .Object@family.frame$V <- .Object@family.frame$U^2
      T <- sum(.Object@family.frame$U, na.rm = TRUE)^2/ sum(.Object@family.frame$V, na.rm = TRUE)

          
    .Object@family.test <- list(statistic = T,
                             p.value = pchisq(q = T, lower.tail = FALSE, df = 1))
    
    }

    if (chromosome == 'X') {
      male <- subset( .Object@family.frame, sex == 1)
      male$U <- male$signal - male$signal.father
      male$V <- male$U^2
      
      female <- subset( .Object@family.frame, sex == 2)
      female$U <- female$signal - female$signal.mother
      female$V <- female$U^2

      T.male <- sum(male$U, na.rm = TRUE)^2 / sum(male$V, na.rm = TRUE)
      T.female <- sum(female$U, na.rm = TRUE)^2 / sum(female$V, na.rm = TRUE)
      T <- (sum(male$U, na.rm = TRUE) + sum(female$U, na.rm = TRUE))^2/(sum(male$V, na.rm = TRUE) + sum(female$V, na.rm = TRUE))

          
    .Object@family.test <- list(statistic = T,
                             p.value = pchisq(q = T, lower.tail = FALSE, df = 1),
                             statistic.male = T.male,
                             statistic.female = T.female)
    
    }

  }
  
  return (.Object)
})
 

#############################################################################
setGeneric("apply.pca",function(.Object, multiprobe.signal = NULL, scale = FALSE){standardGeneric("apply.pca")})

setMethod("apply.pca", "CNVtools", function(.Object, multiprobe.signal, scale) {
 
  if (!is.null(multiprobe.signal)) {
    if (class(multiprobe.signal) != 'matrix') stop('Multi-probe signal must be a matrix')
    if (nrow(multiprobe.signal) != .Object@nsamples) stop('Number of rows of multiprobe signal must be equal to the number of samples')
    .Object@multiprobe.signal <- multiprobe.signal
  }
  
  pca <- prcomp(.Object@multiprobe.signal, scale = scale, center = TRUE)$x[,1]
  pca <- pca/sd(pca)
  if (cor(pca, .Object@signal) < 0) .Object@signal <- - pca else .Object@signal <- pca

  ######### clear the previously fitted objects that are not relevant anymore
  .Object@mle.H0 <- list()
  .Object@mle.H1 <- list()
  .Object@best.fit.H0 <- data.frame()
  .Object@best.fit.H1 <- data.frame()
  
  return(.Object)
})



###############################


setGeneric("plot.cnv",function(.Object, hist.or.dens = "histogram", hyp = 'H0', freq = TRUE, plot.outliers = FALSE, subset, ...){standardGeneric("plot.cnv")}) 

setMethod(f = "plot.cnv", signature ( .Object = "CNVtools", hist.or.dens = 'ANY', hyp = 'ANY',  freq = 'ANY', plot.outliers = 'ANY', subset = 'ANY'), function(.Object, hist.or.dens, hyp, freq, plot.outliers, subset, ...) { ####

  if (!missing(subset)) {
    if (! class(subset) %in% c('logical')) stop('The subset argument must be a logical vector')
    if (length(subset) != .Object@nsamples) stop('If specified, subset must be a logical vector of length equal to the number of samples')
  }

  if (!missing(plot.outliers)) {
    if (! class(plot.outliers) %in% c('logical')) stop('The plot.outliers argument must be a logical vector')
  }

    
  if (!hyp %in% c('H0', 'H1') || length(hyp) != 1) stop('The argument hyp must of of length 1 and equal to H0 or H1')
  #if (!is.null(batch) & !is.null(subset))  stop()
    
  if (hyp == 'H0') posterior <- .Object@best.fit.H0
  if (hyp == 'H1') posterior <- .Object@best.fit.H1

  if (!missing(subset))  posterior <- posterior[subset, ]
  #if (!is.null(batch))  posterior <- posterior[posterior$batch %in% batch, ]

  posterior <- posterior[order(posterior$signal), ]
  if (hist.or.dens == "density") {
    dens <- density(posterior$signal)
    plot(dens, ...)
    my.max <- max(dens$y)
  }
  
  if (hist.or.dens == "histogram") {
    my.hist <- hist(posterior$signal, freq = freq, ...)
    my.max <- max(my.hist$counts)
    if (freq == FALSE) my.max <- max(my.hist$density)
  }


  col <- 1
  ncomp <- max(posterior$cn)
  for (i in c(1:.Object@ncomp)) {
    if (!plot.outliers) pr <- ifelse ( posterior$proba.not.outlier < 0.1, 0, posterior[, paste("P", i, sep = "")])
    if ( plot.outliers) pr <- posterior[, paste("P", i, sep = "")]

    lines(x = posterior$signal, y = my.max*pr, type = "l", col = col)
    col <- col + 2
  }
  
})



