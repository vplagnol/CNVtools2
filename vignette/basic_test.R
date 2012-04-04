#library(CNVtools) 
source('R/CNVtools_class.R')
source('R/model_fitting.R')

dyn.load('src/CNVtools.so')

load('examples/log_ratio/chr1/CNVR397.1.Rdata')
ped.file <- read.table('examples/annotation_files/Vanguard.987Samples.PEDFILE.txt', header = TRUE, sep = '\t')

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
f.member <- ifelse( fid == 0 & mid == 0 & sex == 1, 1, f.member) ##father
f.member <- ifelse( fid == 0 & mid == 0 & sex == 2, 2, f.member) ##mother
f.member <- ifelse( fid != 0 & mid != 0, 3, f.member)


test <- CNVtools.multivariate.signal(IID = names(probe.data),
                                     multivariate.signal = t(probe.data),
                                     trait = trait,
                                     model.association = '~ as.factor(cn)',
                                     FID = FID,
                                     f.member = f.member,
                                     batch = rep('A', ncol(probe.data)))


test <- apply.pca(test)    
#test@signal[50] <- 10
test <- fit(test, ncomp = 3, hyp = 'H0', EM.starting.point = 'kmeans')
