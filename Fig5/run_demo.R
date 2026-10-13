
args <- commandArgs(trailingOnly=TRUE)
flag <- grep('^--file=', commandArgs(FALSE), value=TRUE)
root <- dirname(normalizePath(sub('^--file=', '', flag[1]), winslash='/'))
if(length(args) && !(length(args)==2L && args[1]=='--output')) stop('Usage: Rscript run_demo.R [--output DIR]')
out <- if(length(args)) args[2] else file.path(root,'demo_output')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
extra <- Sys.getenv('DEMO_R_LIB','')
if(nzchar(extra)) .libPaths(c(extra,.libPaths()))
packages <- c('ggplot2','dplyr','tidyr','tibble','patchwork','ggalluvial','VennDiagram','cowplot','vegan','lme4','multcompView')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install packages first: ',paste(missing,collapse=', '))
started <- proc.time()[['elapsed']]
source(file.path(root,'source/panel_helpers.R'), encoding='UTF-8')
long <- read.csv(file.path(root,'input/abundance.csv'),stringsAsFactors=FALSE,check.names=FALSE)
meta <- read.csv(file.path(root,'input/sample_metadata.csv'),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(long)==7344L,nrow(meta)==72L,all(meta$Metadata_status=='SYNTHETIC_DEMO_ONLY'),
          all(is.finite(long$Abundance)),all(long$Abundance>=0), !anyDuplicated(long[c('Sample','Subtype')]))
data <- prepare_fig5_data(long,meta,expected_samples=72L,expected_subtypes=102L)
write.csv(data$samples,file.path(out,'sample_summaries.csv'),row.names=FALSE)
message('Fig5 b/c')
bc <- run_fig5_bc(file.path(root,'input'),file.path(out,'bc'))
stopifnot(nrow(bc$risk)==32L,nrow(bc$tier_i)==8L)
message('Fig5 d')
d <- run_fig5d(data$long,data$meta,file.path(out,'d'),permutations=99L,seed=123L,dpi=150)
p_d <- d$plot
ggplot2::ggsave(file.path(out,'d/Fig5d_PCoA.pdf'),p_d,width=5.5,height=5)
ggplot2::ggsave(file.path(out,'d/Fig5d_PCoA.png'),p_d,width=5.5,height=5,dpi=150)
message('Fig5 e')
e <- run_fig5e(data$samples,file.path(out,'e'))
p_e <- e$plot
ggplot2::ggsave(file.path(out,'e/Fig5e.pdf'),p_e,width=5.3,height=4.8)
ggplot2::ggsave(file.path(out,'e/Fig5e.png'),p_e,width=5.3,height=4.8,dpi=150)
stopifnot(nrow(e$contrasts)==6L,all(is.finite(e$contrasts$estimate)),all(is.finite(e$contrasts$q_BH)))
message('Fig5 f/g')
fg <- run_fig5_fg(data$long,data$meta,file.path(out,'fg'))
pf <- fg$plots$f
pg <- fg$plots$g
for(ext in c('pdf','png')) {
  ggplot2::ggsave(file.path(out,'fg',paste0('Fig5f_TierI_city_composition.',ext)),pf,width=11,height=8,dpi=150)
  ggplot2::ggsave(file.path(out,'fg',paste0('Fig5g_TierI_monthly_bubble.',ext)),pg,width=6.5,height=9.5,dpi=150)
}
coverage <- data.frame(Panel=c('b','c','d','e','f','g'),
  Status=c('PASS_EXISTING_CLASSIFICATION_SUBSET','PASS_EXISTING_CLASSIFICATION_SUBSET',
           'PASS_SYNTHETIC_GROUPS',if(isTRUE(e$diagnostics$Singular)) 'PASS_WITH_SINGULAR_FIT' else 'PASS_SYNTHETIC_GROUPS',
           'PASS_SYNTHETIC_GROUPS','PASS_SYNTHETIC_GROUPS'))
write.csv(coverage,file.path(out,'panel_status.csv'),row.names=FALSE)
write.csv(data.frame(package=packages,version=vapply(packages,function(x)as.character(packageVersion(x)),character(1))),
          file.path(out,'package_versions.csv'),row.names=FALSE)
session <- capture.output(sessionInfo())
session <- session[!grepl('Time zone|tzcode',session,ignore.case=TRUE)]
writeLines(session,file.path(out,'sessionInfo.txt'))
write.csv(data.frame(seconds=round(proc.time()[['elapsed']]-started,3),samples=72L,subtypes=102L,
                     class_records=32L,synthetic_grouping=TRUE),file.path(out,'run_metadata.csv'),row.names=FALSE)
cat('PASS: Fig5 b-g\n')
