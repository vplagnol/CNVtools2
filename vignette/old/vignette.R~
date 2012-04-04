

lib <- TRUE
if (lib) {library(CNVtools)} else {
  library(survival)      
  source('R/CNVtools_class.R')
  source('R/model_fitting.R')
  dyn.load('src/CNVtools.so')
}


manu.CNVs <- function() {
  my.files <- list.files('examples/log_ratio', recursive = TRUE, pattern = 'Rdata', full.names = TRUE)
  ped.file <- read.table('examples/annotation_files/Vanguard.987Samples.PEDFILE.txt', header = TRUE, sep = '\t')
  
  pdf('test.pdf', width = 8, height = 4)
  par(mfrow = c(1, 2))
  for (file in my.files) {
    print(file)
    load(file)
    probe.data <- probe_dat_sub[, - (1:3)]
    my.match <-  match(table = ped.file$Sample.Name, names(probe.data))

    ############ get some info about the samples
    FID <- ped.file$FamID[ my.match ]
    mid <- ped.file$mid[ my.match ]
    fid <- ped.file$fid[ my.match ]
    sex <- ped.file$sex[ my.match ]
    trait <- ped.file$t1d[ my.match ] - 1
    trait <- ifelse(is.na(trait) | trait < 0, 0, trait)
    
    f.member <- NA
    f.member <- ifelse( fid == 0 & mid == 0 & sex == 1, 1, f.member)
    f.member <- ifelse( fid == 0 & mid == 0 & sex == 2, 2, f.member)
    f.member <- ifelse( fid != 0 & mid != 0, 3, f.member)


    ############### Now create the CNVtools object
    test <- CNVtools.multivariate.signal(IID = names(probe.data),
                                         multivariate.signal = t(probe.data),
                                         trait = trait,
                                         model.association = '~ as.factor(cn)',
                                         FID = FID,
                                         f.member = f.member,
                                         batch = rep('A', ncol(probe.data)))
    test <- apply.pca(test)
    
    test <- fit(test, ncomp = 3, hyp = 'H0', EM.starting.point = 'kmeans')
    test <- fit(test, ncomp = 3, hyp = 'H1', EM.starting.point = 'H1.from.H0')
    
    plot.cnv (test, col = 'red', breaks = 20, hyp = 'H0', main = paste(gsub(pattern = '.Rdata', replacement = '', basename(file)), 'H0'))
    plot.cnv (test, col = 'red', breaks = 20, hyp = 'H1', main = paste(gsub(pattern = '.Rdata', replacement = '', basename(file)), 'H1'))
    title (paste('TDT test: ', signif(TDT(test)$p.value, 3)), outer = TRUE, line = -1)
    print(TDT(test))
  }
  
  dev.off()

  
  T <- 2*(test@mle.H1$lnL - test@mle.H0$lnL)

  #print(test@mle.H0$mean)
  #print(test@mle.H1$mean)
  #print(T)
  #print(test@mle.H0$mean)
  return(test)
}

test <- manu.CNVs()






basic.check <- function(niter = 1) {
  #set.seed(seed = 0)
  if (file.exists('test.tab')) file.remove('test.tab')

  n <- 500
  odds.ratio <- 1.3

  
  for (i in 1:niter) {
    X <- rbinom(size = 2, n = n, prob = 0.3)
    nu <- log(odds.ratio)*X
    prob <- exp(nu)/(1 + exp(nu))
    trait <- rbinom(size = 1, n = n, prob = prob)

    mod1 <- glm( trait ~ X, family = binomial)
    mod0 <- glm( trait ~ 1, family = binomial)
    dl <-     anova(mod0, mod1, test = 'Chisq')$Dev[2]

    test <- new('CNVtools',
                trait = trait,
                model.association = '~ cn',
                signal = X + rnorm(sd = 0.1, n = n),
                batch = rep('A', n))
    
    test <- fit(test, ncomp = 3, hyp = 'H0', EM.starting.point = 'manual', start.mean = c(0, 1, 2), start.var = rep(0.1^2, 3))
    test <- fit(test, ncomp = 3, hyp = 'H1', EM.starting.point = 'H1.from.H0')
    T <- 2*(test@mle.H1$lnL - test@mle.H0$lnL)
    
    cat(T, dl, '\n')
    cat(T, dl, '\n', file = 'test.tab', append = TRUE )
  }
  print(test)
  print(test@mle.H0$mean)
  print(test@mle.H1$mean)
  print(T)
}
#basic.check(10)
