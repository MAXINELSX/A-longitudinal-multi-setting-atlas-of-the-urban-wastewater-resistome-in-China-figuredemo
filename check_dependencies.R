check_dependencies <- function() {
  extra <- Sys.getenv('DEMO_R_LIB','')
  if(nzchar(extra)) .libPaths(c(extra,.libPaths()))
  required <- c('cowplot','data.table','dplyr','forcats','ggalluvial','ggnewscale',
    'ggplot2','lme4','lmerTest','multcompView','openxlsx','patchwork','purrr',
    'randomForest','readr','readxl','rfPermute','rstatix','scales','stringr',
    'tibble','tidyr','vegan','VennDiagram')
  missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop('Install required packages first: ',paste(missing,collapse=', '))
  if(!capabilities('cairo'))stop('R Cairo graphics support is required.')
  optional <- 'performance'
  available <- c(required,optional[vapply(optional,requireNamespace,logical(1),quietly=TRUE)])
  data.frame(package=available,version=vapply(available,function(p)as.character(packageVersion(p)),character(1)),
             required=available %in% required,row.names=NULL)
}
if(sys.nframe()==0L)print(check_dependencies(),row.names=FALSE)
