compose_fig2 <- function(out) {
  needed <- c('png', 'grid')
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) stop('Install required packages: ', paste(missing, collapse=', '))
  panel <- function(...) file.path(out, ...)
  first_existing <- function(primary, fallback) if(file.exists(primary)) primary else fallback
  source_images <- list(
    a=panel('panel_a','panel_a_type_mean_abundance.png'),
    b=panel('panels_b_c','half_violin_site_richness_refstyle_narrow60.png'),
    c=panel('panels_b_c','half_violin_site_copiescell_refstyle_narrow60_fullrange.png'),
    d=panel('panel_d_global','copies_per_cell_by_region_with_thisstudy_region07.png'),
    e=panel('panel_e_core_scatter','core_subtype_distribution_by_site.png'),
    f=first_existing(panel('panel_f_original_style','ARG_type_3quarter_circle_groupedTree_publication_ready.png'),
      panel('panel_f_audited','Fig2f_circle_ALL_core_subtypes.png')),
    f_upset=first_existing(panel('panel_f_original_style','ARG_subtype_upset_ALL_core_subtypes_FIXED_label_order.png'),
      panel('panel_f_audited','Fig2f_upset_ALL_core_subtypes.png')),
    g=panel('panels_gh','figures','Fig2g_clinical9_detection_corrected.png'),
    h=panel('panels_gh','figures','Fig2h_caption8_monthly_abundance_corrected.png')
  )
  absent <- unlist(source_images)[!file.exists(unlist(source_images))]
  if(length(absent)) stop('Missing panel images: ', paste(absent, collapse=', '))
  read_trim <- function(path) {
    img <- png::readPNG(path)
    if(length(dim(img))==2L) img <- array(rep(img,3L), dim=c(dim(img),3L))
    mask <- apply(img[,,seq_len(min(dim(img)[3],3L)),drop=FALSE] < 0.985, c(1,2), any)
    keep_y <- which(rowSums(mask)>0L)
    keep_x <- which(colSums(mask)>0L)
    if(!length(keep_y) || !length(keep_x)) return(img)
    pad <- 5L
    yr <- max(1L,min(keep_y)-pad):min(dim(img)[1],max(keep_y)+pad)
    xr <- max(1L,min(keep_x)-pad):min(dim(img)[2],max(keep_x)+pad)
    img[yr,xr,,drop=FALSE]
  }
  images <- lapply(source_images,read_trim)
  draw <- function(name,x,y,w,h) {
    grid::pushViewport(grid::viewport(x=x,y=y,width=w,height=h,just=c('left','bottom')))
    grid::grid.raster(images[[name]],x=0.5,y=0.5,width=1,height=1,default.units='npc',interpolate=TRUE)
    grid::popViewport()
  }
  letter <- function(label,x,y) grid::grid.text(label,x=x,y=y,default.units='npc',just=c('left','top'),gp=grid::gpar(fontsize=16,fontface='bold'))
  render <- function() {
    grid::grid.newpage()
    grid::grid.rect(gp=grid::gpar(fill='white',col=NA))
    draw('a',0.030,0.735,0.365,0.245)
    draw('b',0.427,0.772,0.125,0.208)
    draw('c',0.574,0.772,0.125,0.208)
    draw('d',0.730,0.772,0.245,0.208)
    draw('e',0.028,0.380,0.370,0.332)
    draw('f',0.420,0.367,0.347,0.394)
    draw('f_upset',0.405,0.655,0.172,0.123)
    draw('g',0.788,0.400,0.195,0.348)
    draw('h',0.027,0.025,0.956,0.312)
    letter('a',0.012,0.983)
    letter('b',0.407,0.983)
    letter('c',0.557,0.983)
    letter('d',0.710,0.983)
    letter('e',0.012,0.733)
    letter('f',0.407,0.770)
    letter('g',0.770,0.763)
    letter('h',0.012,0.353)
  }
  png_path <- panel('Fig2_demo_composite.png')
  pdf_path <- panel('Fig2_demo_composite.pdf')
  grDevices::png(png_path,width=3080,height=2530,res=220,bg='white')
  render()
  grDevices::dev.off()
  grDevices::pdf(pdf_path,width=14,height=11.5,bg='white',onefile=FALSE)
  render()
  grDevices::dev.off()
  invisible(c(png=png_path,pdf=pdf_path))
}
