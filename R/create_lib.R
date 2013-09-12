###create the CNVtools2 library

base  <- '/ugi/home/shared/vincent/libraries/R/working/CNVtools'


system('rm /ugi/home/shared/vincent/libraries/R/working/CNVtools/R/*')
package.skeleton(name="CNVtools",
                 code_files = c('R/CNVtools_class.R', 'R/model_fitting.R', 'R/tools.R'),
                 path='/ugi/home/shared/vincent/libraries/R/working',
                 force=TRUE)

for (folder in c('data', 'src', 'inst', 'inst/doc', 'inst/doc/fig')) {
  my.folder <- paste('/ugi/home/shared/vincent/libraries/R/working/CNVtools/', folder, sep = '')
  if (!file.exists(my.folder)) dir.create(path = my.folder)
}



############### data for case control
data <- read.table("data/A112.dat", header=TRUE)
A112 <- subset(data, data$cohort %in% c('58C', 'NBS'))
save(list = c('A112'), file="/ugi/home/shared/vincent/libraries/R/working/CNVtools/data/A112.RData")

############### prepare family object
pedfile <- read.table('examples/annotation_files/Vanguard.987Samples.PEDFILE.txt', header = TRUE, sep = '\t')
load('examples/log_ratio/chr1/CNVR397.1.Rdata')
probe.data <- probe_dat_sub[, - (1:3)]
my.match <-  match(table = pedfile$Sample.Name, names(probe.data))
FID <- pedfile$FamID[ my.match ]
mid <- pedfile$mid[ my.match ]
fid <- pedfile$fid[ my.match ]
sex <- pedfile$sex[ my.match ]
trait <- pedfile$t1d[ my.match ] - 1
trait <- ifelse(is.na(trait) | trait < 0, 0, trait)
f.member <- NA
f.member <- ifelse( fid == 0 & mid == 0 & sex == 1, 1, f.member) ##father
f.member <- ifelse( fid == 0 & mid == 0 & sex == 2, 2, f.member) ##mother
f.member <- ifelse( fid != 0 & mid != 0, 3, f.member)
save(list = c('probe.data', 'trait', 'FID', 'f.member'), file="/ugi/home/shared/vincent/libraries/R/working/CNVtools/data/family_CNV.RData")



file.copy(from = c('vignette/vignette-family.Rnw', 'vignette/vignette-case_control.Rnw'),  '/ugi/home/shared/vincent/libraries/R/working/CNVtools/inst/doc', overwrite = TRUE)

system(paste("cp src/*.cpp src/*.c src/*.h ", base, '/src', sep = ''))
file.copy(from = 'doc/NAMESPACE', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/NAMESPACE', overwrite = TRUE)
file.copy(from = 'doc/DESCRIPTION', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/DESCRIPTION', overwrite = TRUE)




file.copy(from = 'doc/CNVtools.multivariate.signal.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/CNVtools.multivariate.signal.Rd', overwrite = TRUE)
file.copy(from = 'doc/FamilyTest.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/FamilyTest.Rd', overwrite = TRUE)
file.copy(from = 'doc/familial.scaling.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/familial.scaling.Rd', overwrite = TRUE)
file.copy(from = 'doc/apply.pca.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/apply.pca.Rd', overwrite = TRUE)
file.copy(from = 'doc/plot.cnv.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/plot.cnv.Rd', overwrite = TRUE)
file.copy(from = 'doc/validCNVtools.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/validCNVtools.Rd', overwrite = TRUE)
file.copy(from = 'doc/compact.data.frame.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/compact.data.frame.Rd', overwrite = TRUE)
file.copy(from = 'doc/CNVtools.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/CNVtools.Rd', overwrite = TRUE)
file.copy(from = 'doc/getparams.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/getparams.Rd', overwrite = TRUE)
file.copy(from = 'doc/expand.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/expand.Rd', overwrite = TRUE)
file.copy(from = 'doc/fit.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/fit.Rd', overwrite = TRUE)
file.copy(from = 'doc/getQualityScore.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/getQualityScore.Rd', overwrite = TRUE)
file.copy(from = 'doc/CNVtools-package.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/CNVtools-package.Rd', overwrite = TRUE)

                                        #file.copy(from = 'doc/expand.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/expand.Rd', overwrite = TRUE)

system("R CMD build /ugi/home/shared/vincent/libraries/R/working/CNVtools")
system("R CMD INSTALL CNVtools_2.2.0.tar.gz")
