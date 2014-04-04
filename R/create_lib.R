###create the CNVtools2 library

base  <- '/cluster/project4/vyp/vincent/libraries/R/working/CNVtools'

system('rm /cluster/project4/vyp/vincent/libraries/R/working/CNVtools/R/*')
package.skeleton(name="CNVtools",
                 code_files = c('R/CNVtools_class.R', 'R/model_fitting.R', 'R/tools.R'),
                 path='/cluster/project4/vyp/vincent/libraries/R/working',
                 force=TRUE)

for (folder in c('data', 'src', 'inst', 'inst/doc', 'inst/doc/fig')) {
  my.folder <- paste(base, '/', folder, sep = '')
  print(my.folder)
  if (!file.exists(my.folder)) dir.create(path = my.folder)
}



############### data for case control
data <- read.table("data/A112.dat", header=TRUE)
A112 <- subset(data, data$cohort %in% c('58C', 'NBS'))
save(list = c('A112'), file=paste(base, "/data/A112.RData", sep = ''))

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
save(list = c('probe.data', 'trait', 'FID', 'f.member'), file=paste(base, "/data/family_CNV.RData", sep = ''))



file.copy(from = c('vignette/vignette-family.Rnw', 'vignette/vignette-case_control.Rnw'),  paste(base, '/inst/doc', sep = ''), overwrite = TRUE)

system(paste("cp src/*.cpp src/*.c src/*.h ", base, '/src', sep = ''))
file.copy(from = 'doc/NAMESPACE', to = paste(base, '/NAMESPACE', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/DESCRIPTION', to = paste(base, '/DESCRIPTION', sep = ''), overwrite = TRUE)




file.copy(from = 'doc/CNVtools.multivariate.signal.Rd', to = paste(base, '/man/CNVtools.multivariate.signal.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/FamilyTest.Rd', to = paste(base, '/man/FamilyTest.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/familial.scaling.Rd', to = paste(base, '/man/familial.scaling.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/apply.pca.Rd', to = paste(base, '/man/apply.pca.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/plot.cnv.Rd', to = paste(base, '/man/plot.cnv.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/validCNVtools.Rd', to = paste(base, '/man/validCNVtools.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/compact.data.frame.Rd', to = paste(base, '/man/compact.data.frame.Rd',sep = ''), overwrite = TRUE)
file.copy(from = 'doc/CNVtools.Rd', to = paste(base, '/man/CNVtools.Rd',sep = ''), overwrite = TRUE)
file.copy(from = 'doc/getparams.Rd', to = paste(base, '/man/getparams.Rd',sep = ''), overwrite = TRUE)
file.copy(from = 'doc/expand.Rd', to = paste(base, '/man/expand.Rd',sep = ''), overwrite = TRUE)
file.copy(from = 'doc/fit.Rd', to = paste(base, '/man/fit.Rd',sep = ''), overwrite = TRUE)
file.copy(from = 'doc/getQualityScore.Rd', to = paste(base, '/man/getQualityScore.Rd',sep = ''), overwrite = TRUE)
file.copy(from = 'doc/CNVtools-package.Rd', to = paste(base, '/man/CNVtools-package.Rd',sep = ''), overwrite = TRUE)



system(paste("R CMD build ", base, sep = ''))
system("R CMD INSTALL CNVtools_2.2.1.tar.gz")
