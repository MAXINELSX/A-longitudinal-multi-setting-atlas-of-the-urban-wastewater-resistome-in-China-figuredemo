


args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
root <- dirname(normalizePath(sub('^--file=','',script_arg[1]),winslash='/'))
if(length(args) && !(length(args)==2L && args[1]=='--output'))
  stop('Usage: Rscript --vanilla run_all.R [--output NEW_DIRECTORY]')
out <- if(length(args)) args[2] else file.path(root,'generated_output')
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop('Output directory is not empty. Choose a new directory to preserve prior results.')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,winslash='/',mustWork=TRUE)
log_dir <- file.path(out,'logs')
dir.create(log_dir,showWarnings=FALSE)
Sys.setenv(LC_ALL='C')
extra <- Sys.getenv('DEMO_R_LIB','')
if(nzchar(extra)) .libPaths(c(extra,.libPaths()))
source(file.path(root,'check_dependencies.R'),local=TRUE)
inventory <- check_dependencies()
write.csv(inventory,file.path(out,'package_versions.csv'),row.names=FALSE)
rscript <- file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
summary <- data.frame(module=character(),exit_code=integer(),elapsed_seconds=numeric())
for(module in paste0('Fig',1:6)) {
  message(module)
  started <- proc.time()[['elapsed']]
  log <- file.path(log_dir,paste0(module,'.log'))
  status <- system2(rscript,c('--vanilla',shQuote(file.path(root,module,'run_demo.R')),
                              '--output',shQuote(file.path(out,module))),stdout=log,stderr=log)

  log_text <- readLines(log,warn=FALSE)
  for(p in c(out,gsub('/','\\',out,fixed=TRUE)))
    for(prefix in c(p,encodeString(p)))
      log_text <- gsub(prefix,'OUTPUT',log_text,fixed=TRUE)
  for(p in c(root,gsub('/','\\',root,fixed=TRUE)))
    for(prefix in c(p,encodeString(p)))
      log_text <- gsub(prefix,'PACKAGE',log_text,fixed=TRUE)
  writeLines(log_text,log)
  summary <- rbind(summary,data.frame(module=module,exit_code=status,
                         elapsed_seconds=round(proc.time()[['elapsed']]-started,3)))
  write.csv(summary,file.path(out,'run_summary.csv'),row.names=FALSE)
  if(status!=0L) {
    cat(tail(readLines(log,warn=FALSE),30),sep='\n')
    stop(module,' failed; inspect its log. Remaining modules were not run.')
  }
}
session <- capture.output(sessionInfo())
session <- session[!grepl('Time zone|tzcode',session,ignore.case=TRUE)]
writeLines(session,file.path(out,'sessionInfo.txt'))
print(summary,row.names=FALSE)
