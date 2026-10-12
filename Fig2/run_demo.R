

started <- proc.time()[['elapsed']]
if (.Platform$OS.type == 'windows') invisible(suppressWarnings(Sys.setlocale('LC_ALL','English_United States.utf8')))
args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
base <- dirname(normalizePath(sub('^--file=', '', script_arg[1]), winslash='/'))
if(length(args) && !(length(args)==2L && args[1]=='--output')) stop('Usage: Rscript run_demo.R [--output OUTPUT_DIRECTORY]')
out <- if(length(args)) args[2] else file.path(base,'output')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,winslash='/',mustWork=TRUE)
extra_lib <- Sys.getenv('DEMO_R_LIB','')
if(nzchar(extra_lib)) .libPaths(c(extra_lib,.libPaths()))
required <- c('readxl','openxlsx','ggplot2','dplyr','tidyr','stringr','purrr','scales',
              'tibble','patchwork','lme4','multcompView','rstatix','readr','forcats')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install required packages: ',paste(missing,collapse=', '))
save_csv <- function(x,name) write.csv(x,file.path(out,name),row.names=FALSE,na='',fileEncoding='UTF-8')
meta <- read.csv(file.path(base,'input/sample_metadata.csv'),check.names=FALSE,stringsAsFactors=FALSE)
wide <- read.csv(file.path(base,'input/subtype_abundance.csv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(identical(names(wide)[-1],meta$Sample),nrow(meta)==72L,nrow(wide)==2646L,
          !anyDuplicated(meta$Sample),!anyDuplicated(wide$Subtype),!anyNA(wide))
values <- as.matrix(wide[-1])
stopifnot(is.numeric(values),all(is.finite(values)),all(values>=0),
          all(meta$Metadata_status=='SYNTHETIC_DEMO_ONLY'))
meta$Sample_Date <- as.Date(meta$Sample_Date)
long <- data.frame(Sample=rep(meta$Sample,each=nrow(wide)),Subtype=rep(wide$Subtype,nrow(meta)),
                   Abundance=as.vector(values),stringsAsFactors=FALSE) |>
  dplyr::left_join(meta,by='Sample')
samples <- long |>
  dplyr::group_by(Sample) |>
  dplyr::summarise(Richness=sum(Abundance>0),Total_abundance=sum(Abundance),.groups='drop') |>
  dplyr::left_join(meta,by='Sample') |>
  dplyr::mutate(Raw_blank_count=0L)
cohort <- list(long=long,meta=meta,samples=samples)
stopifnot(nrow(unique(meta[c('City','Setting','PhysicalSite')]))==24L,
          all(table(meta$Setting)==18L),all(table(meta$PhysicalSite)==3L),
          all(table(meta$City,meta$Setting,meta$SamplingMonth)==2L),
          qr(model.matrix(~Setting+City+SamplingMonth,meta))$rank==8L)
save_csv(samples,'sample_summaries.csv')
save_csv(meta,'demo_metadata.csv')
save_csv(data.frame(Check=c('72_unique_source_sample_columns','2646_subtypes_per_sample',
 'finite_explicit_nonnegative_abundance','no_imputation_or_duplicated_sample_columns',
 'balanced_four_assigned_groups','24_assigned_sites_three_repeats_each',
 'full_rank_fixed_effect_design','all_assigned_metadata_marked_synthetic'),Status='PASS'),'input_validation.csv')

fig2 <- new.env(parent=globalenv())
source(file.path(base,'upstream/analysis_functions.R'),local=fig2,encoding='UTF-8')
fig2$.canonical_long <- long
fig2$.canonical_meta <- meta
fig2$.figure_results_dir <- out
fig2$FIG2_DATA_ROOT <- file.path(base,'input')
fig2$audited_city_order <- sort(unique(meta$City))
fig2$.f_type_colours <- c(fig2$.f_type_colours, capO='#D0D0D0')
fig2$.external_path <- function(...) file.path(base,'input',...)


fig2$demo_ggsave <- function(filename,plot=ggplot2::last_plot(),...) {
  dots <- list(...)
  if(grepl('\\.png$',filename,ignore.case=TRUE)) dots$dpi <- 130
  do.call(ggplot2::ggsave,c(list(filename=filename,plot=plot),dots))
}
fig2$demo_write_table <- function(x,file,...) {
  csv_file <- sub('\\.xlsx$','.csv',file,ignore.case=TRUE)
  if(is.data.frame(x) || is.matrix(x)) write.csv(x,csv_file,row.names=FALSE,na='',fileEncoding='UTF-8')
  else if(is.list(x)) for(nm in names(x)) write.csv(x[[nm]],sub('\\.csv$',paste0('_',nm,'.csv'),csv_file),row.names=FALSE,na='')
  else stop('Unsupported upstream table output type')
  invisible(csv_file)
}
Sys.setenv(ARG_EXPORT_FULL_TABLES='0',ARG_GLOBAL_DUPLICATES='drop')
coverage <- list(); warnings_seen <- list()
run_panel <- function(panel,expression,method) {
  panel_warnings <- character()
  value <- tryCatch(withCallingHandlers(force(expression),warning=function(w) {
    panel_warnings <<- c(panel_warnings,conditionMessage(w));invokeRestart('muffleWarning')
  }),error=function(e)e)
  success <- !inherits(value,'error')
  coverage[[length(coverage)+1L]] <<- data.frame(Panel=panel,Status=if(success)'RUN' else 'FAILED',
    Method=method,Detail=if(success)'' else conditionMessage(value))
  if(length(panel_warnings)) warnings_seen[[length(warnings_seen)+1L]] <<- data.frame(Panel=panel,Warning=unique(panel_warnings))
  save_csv(dplyr::bind_rows(coverage),'panel_coverage.csv')
  if(!success) message(panel,': ',conditionMessage(value))
  invisible(value)
}
run_panel('Fig2a',fig2$run_panel_a(),'COMPOSITION')
b <- run_panel('Fig2b',fig2$fit_fig2_model(samples,'Richness','negative_binomial',file.path(out,'panel_b_NB')),
               'NB_LMM_BH')
c <- run_panel('Fig2c',fig2$fit_fig2_model(samples,'Total_abundance','gaussian_original',file.path(out,'panel_c_LMM')),
               'GAUSSIAN_LMM_BH')
if(!inherits(b,'error') && !inherits(c,'error')) run_panel('Fig2b/c plots',fig2$run_panels_bc(cohort),
               'VIOLIN')
run_panel('Fig2d',fig2$run_panel_d(),
          'KRUSKAL_DUNN')
run_panel('Fig2e',fig2$run_panel_e(),'CORE_PREVALENCE')
has_original_style <- requireNamespace('ggnewscale',quietly=TRUE)
f <- run_panel('Fig2f',fig2$run_panel_f_original(make_plots=has_original_style),
               if(has_original_style)'CIRCLE_UPSET' else 'TABLES_ONLY')
if(!has_original_style && !inherits(f,'error')) run_panel('Fig2f fallback plot',fig2$audited_panel_f(make_plots=TRUE),
          'CIRCLE_UPSET_ALTERNATIVE')
run_panel('Fig2g/h',fig2$run_panels_gh(),
          'DETECTION_MONTHLY')
save_csv(if(length(warnings_seen))dplyr::bind_rows(warnings_seen) else data.frame(Panel=character(),Warning=character()),'warnings.csv')
save_csv(data.frame(Package=c(required,if(has_original_style)'ggnewscale'),
 Version=vapply(c(required,if(has_original_style)'ggnewscale'),function(p)as.character(packageVersion(p)),character(1))),
 'package_versions.csv')
session_info <- capture.output(sessionInfo())
session_info <- session_info[!grepl('time zone:|tzcode source:',session_info,ignore.case=TRUE)]
writeLines(session_info,file.path(out,'R_session_info.txt'))
save_csv(data.frame(Elapsed_seconds=round(proc.time()[['elapsed']]-started,2),Samples=72L,
 Subtypes=2646L,Real_abundance_cells=190512L,Synthetic_metadata=TRUE,Original_f_rendering=has_original_style),'run_record.csv')
failed <- dplyr::bind_rows(coverage)$Status=='FAILED'
if(any(failed)) stop('Some panels failed. See panel_coverage.csv; successful outputs are retained.')
cat('PASS', '\n')
