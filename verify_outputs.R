
arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
root <- dirname(normalizePath(sub('^--file=','',arg[1]),winslash='/'))
args <- commandArgs(TRUE)
if(length(args)!=2L || args[1]!='--output')stop('Usage: Rscript verify_outputs.R --output DIRECTORY')
out <- normalizePath(args[2],winslash='/',mustWork=TRUE)
tab <- function(path)read.csv(file.path(out,path),check.names=FALSE,stringsAsFactors=FALSE)
checks <- data.frame(check=character(),status=character())
check <- function(name,condition) {
  if(!isTRUE(condition))stop('FAIL: ',name)
  checks <<- rbind(checks,data.frame(check=name,status='PASS'))
}
r <- tab('run_summary.csv')
check('All six included modules exited zero',nrow(r)==6L && all(r$exit_code==0L) && identical(r$module,paste0('Fig',1:6)))

expected <- read.csv(file.path(root,'Fig1/expected/summary.csv'))
actual <- tab('Fig1/summary.csv')
check('Fig1 metadata summary matches independent reference',
  identical(actual$metric,expected$metric) && isTRUE(all.equal(actual$value,expected$value,tolerance=1e-12)))
check('Fig1 distinct technical branch generated',file.exists(file.path(out,'Fig1/e_technical_demo/pcoa_technical_demo.pdf')))

wide <- read.csv(file.path(root,'Fig2/input/subtype_abundance.csv'),check.names=FALSE)
s <- tab('Fig2/sample_summaries.csv')
check('Fig2 all 72 abundance totals reconciled',nrow(s)==72L && nrow(wide)==2646L &&
  isTRUE(all.equal(s$Total_abundance,colSums(wide[-1])[match(s$Sample,names(wide)[-1])],check.attributes=FALSE,tolerance=1e-12)))
check('Fig2 richness independently reconciled',all(s$Richness==colSums(wide[-1]>0)[match(s$Sample,names(wide)[-1])]))
coverage2 <- tab('Fig2/panel_coverage.csv')
check('Fig2 included a-h branches and composite reported as executed',
      nrow(coverage2)==9L && all(coverage2$Status=='RUN') &&
      setequal(coverage2$Panel,
               c('Fig2a','Fig2b','Fig2c','Fig2b/c plots','Fig2d','Fig2e','Fig2f','Fig2g/h','Fig2 composite')))
tier <- read.csv(file.path(root,'Fig5/input/abundance.csv'),check.names=FALSE)
ij <- cbind(match(tier$Subtype,wide$Subtype),match(tier$Sample,names(wide)[-1]))
check('Fig2/Fig5 sample aliases and retained subtype values agree',!anyNA(ij) &&
  isTRUE(all.equal(tier$Abundance,as.matrix(wide[-1])[ij],check.attributes=FALSE,tolerance=1e-12)))
check('Fig2/Fig5 metadata are identical',identical(
  read.csv(file.path(root,'Fig2/input/sample_metadata.csv')),
  read.csv(file.path(root,'Fig5/input/sample_metadata.csv'))))


rscript <- file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
Sys.setenv(LC_ALL='C')
commands <- list(
  Fig3=c(file.path(root,'Fig3/verify_outputs.R'),file.path(out,'Fig3')),
  Fig4=c(file.path(root,'Fig4/verify_saved_outputs.R'),'--output',file.path(out,'Fig4')),
  Fig5=c(file.path(root,'Fig5/test_demo.R'),file.path(out,'Fig5')),
  Fig6=c(file.path(root,'Fig6/tests/verify_outputs.R'),'--output',file.path(out,'Fig6')))
for(nm in names(commands)) {
  status <- system2(rscript,c('--vanilla',vapply(commands[[nm]],shQuote,character(1))))
  check(paste(nm,'saved numerical checks'),status==0L)
}
for(nm in paste0('Fig',1:6)) {
  pdf <- list.files(file.path(out,nm),pattern='[.]pdf$',full.names=TRUE,recursive=TRUE)
  check(paste(nm,'nonempty figure files'),length(pdf)>0L && all(file.info(pdf)$size>1000))
}
write.csv(checks,file.path(out,'verification_summary.csv'),row.names=FALSE)
cat('PASS: integrated output totals, module checks and included figures verified.\n')
