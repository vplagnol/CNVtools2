###create the CNVtools2 library

base <- '/ugi/home/shared/vincent/libraries/R/working/CNVtools'

package.skeleton(name="CNVtools",
                 code_files = c('R/CNVtools_class.R', 'R/model_fitting.R'),
                 path='/ugi/home/shared/vincent/libraries/R/working',
                 force=TRUE)

for (folder in c('data', 'src', 'inst', 'inst/doc', 'inst/doc/fig')) {
  my.folder <- paste('/ugi/home/shared/vincent/libraries/R/working/CNVtools/', folder, sep = '')
  if (!file.exists(my.folder)) dir.create(path = my.folder)
}

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

system("R CMD build /ugi/home/shared/vincent/libraries/R/working/CNVtools")
system("R CMD INSTALL CNVtools_2.0.0.tar.gz")
