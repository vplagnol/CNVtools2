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

data <- read.table("data/A112.dat", header=TRUE)
A112 <- subset(data, data$cohort %in% c('58C', 'NBS'))
save(list = 'A112', file="/ugi/home/shared/vincent/libraries/R/working/CNVtools/data/A112.RData")


system(paste("cp src/*.cpp src/*.c src/*.h ", base, '/src', sep = ''))
file.copy(from = 'doc/NAMESPACE', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/NAMESPACE', overwrite = TRUE)
file.copy(from = 'doc/DESCRIPTION', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/DESCRIPTION', overwrite = TRUE)




file.copy(from = 'doc/CNVtools.multivariate.signal.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/CNVtools.multivariate.signal.Rd', overwrite = TRUE)
file.copy(from = 'doc/TDT.Rd', to = '/ugi/home/shared/vincent/libraries/R/working/CNVtools/man/TDT.Rd', overwrite = TRUE)
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
system("R CMD INSTALL CNVtools_2.1.0.tar.gz")
