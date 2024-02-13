plotVafwithGR = function(gr, vafCol = NULL, genes = NULL, top = 10, main_title = NULL,
                         orderByMedian = TRUE, keepGeneOrder = FALSE, flip = FALSE, fn = NULL,
                         gene_fs = 0.8, axis_fs = 0.8, height = 5, width = 5, showN = TRUE, color = NULL){
  
  if(is.null(genes)){
    genes = gr$plotVafIndex %>% unique()
  }
  dat = data.frame(Hugo_Symbol = gr$plotVafIndex, t_vaf = gr$VAF)
  dat = as(dat, "data.table")
  if (is.null(main_title))
  {
    main_title = gr$SYMBOL %>% unique()
  }
  if(!'t_vaf' %in% colnames(dat)){
    if(is.null(vafCol)){
      if(all(c('t_ref_count', 't_alt_count') %in% colnames(dat))){
        message("t_vaf field is missing, but found t_ref_count & t_alt_count columns. Estimating vaf..")
        dat[,t_vaf := as.numeric(as.character(t_alt_count))/(as.numeric(as.character(t_ref_count)) + as.numeric(as.character(t_alt_count)))]
      }else{
        print(colnames(dat))
        stop('t_vaf field is missing. Use vafCol to manually specify vaf column name.')
      }
    }else{
      colnames(dat)[which(colnames(dat) == vafCol)] = 't_vaf'
    }
  }
  
  #dat.genes = data.frame(dat[dat$Hugo_Symbol %in% genes])
  #suppressMessages(datm <- melt(dat.genes[,c('Hugo_Symbol', 't_vaf')]))
  dat.genes = dat[dat$Hugo_Symbol %in% genes]
  datm <- data.table::melt(data = dat.genes[,.(Hugo_Symbol, t_vaf)], id.vars = 'Hugo_Symbol', measure.vars = 't_vaf')
  #remove NA from vcf
  datm = datm[!is.na(value)]
  datm[,value := as.numeric(as.character(value))]
  if(nrow(datm) == 0){
    stop("Nothing to plot.")
  }
  
  #maximum vaf
  if(max(datm$value, na.rm = TRUE) > 1){
    datm$value = datm$value/100
  }
  
  if(keepGeneOrder){
    geneOrder = genes
    datm$Hugo_Symbol = factor(x = datm$Hugo_Symbol, levels = geneOrder)
  }else if(orderByMedian){
    geneOrder = datm[,median(value),Hugo_Symbol][order(V1, decreasing = TRUE)][,Hugo_Symbol]
    datm$Hugo_Symbol = factor(x = datm$Hugo_Symbol, levels = geneOrder)
  }
  
  if(!is.null(color)){
    bcol = color
  }else{
    bcol = c(RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
             RColorBrewer::brewer.pal(n = 8, name = "Accent"))
    if(length(genes) > length(bcol)){
      bcol = sample(x = colors(distinct = TRUE), size = length(genes), replace = FALSE)
    }
  }
  
  if(!is.null(fn)){
    pdf(file = paste0(fn, ".pdf"), width = width, height = height, paper = "special", bg = "white")
  }
  
  if(flip){
    par(mar = c(3, 4, 2, 2))
    b = boxplot(value ~ Hugo_Symbol, data = datm, xaxt="n", boxwex=0.5, outline=FALSE, lty=1, lwd = 1.4, outwex=0,
                staplewex=0, ylim = c(0, 1), axes = FALSE, border = bcol, horizontal = TRUE, ylab = "vaf(%)")
    axis(side = 1, at = seq(0, 1, 0.2), las = 1, font =1, lwd = 1.5, cex.axis = axis_fs)
    axis(side = 2, at = 1:length(b$names), labels = b$names, tick = FALSE, las = 2, font = 3, line = -1, cex.axis = gene_fs)
    if(showN){
      axis(side = 4, at = 1:length(b$names), labels = b$n, font =1, tick = FALSE, line = -1, las = 2, cex.axis = gene_fs)
    }
    abline(v = seq(0, 1, 0.2), h = 1:length(b$names), col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.25), lty = 2)
    
    
    stripchart(value ~ Hugo_Symbol, vertical = FALSE, data = datm,
               method = "jitter", add = TRUE, pch = 16, col = bcol, cex = 0.5)
    
  }else{
    par(mar = c(4, 3, 2, 1), cex.main = 1)
    b = boxplot(value ~ Hugo_Symbol, data = datm, xaxt="n", boxwex=0.5, outline=FALSE, lty=1, lwd = 1.4, outwex=0,
                staplewex=0, ylim = c(0, 1), axes = FALSE, border = bcol,  xlab = "AAChange",
                ylab = "VAF")
    title( main_title, line = 1)
    axis(side = 1, at = 1:length(b$names), labels = b$names, tick = FALSE, las = 2, font = 3, line = -1, cex.axis = gene_fs)
    axis(side = 2, at = seq(0, 1, 0.2), las = 2, cex.axis = axis_fs, lwd = 1.2, font.axis = 2, cex = 1.5, font =1)
    if(showN){
      axis(side = 3, at = 1:length(b$names), labels = b$n, font =1, tick = FALSE,
           line = -1, cex.axis = gene_fs, pos = 0.975)
    }
    abline(h = seq(0, 1, 0.2), v = 1:length(b$names), col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.5), lty = 2)
    
    stripchart(value ~ Hugo_Symbol, vertical = TRUE, data = datm,
               method = "jitter", add = TRUE, pch = 16, col = bcol, cex = 0.5)
    
  }
  
  if(!is.null(fn)){
    dev.off()
  }
}

#' Plotting Distribution of Variant Allele Frequency of Variants
#'@param gr a GRanges object contains annotated variants of patients.
#'@param genes specify genes for which plots has to be generated
#'@param groupByAAchanges specify Amino Acids changes for which plots has to be generated
#'@param bySingleLocus specify the position of the gene products for which plots has to be generated
#'@param title title of the plot
#'@param SNPonly If only plot the SNP.
#'
#'@import data.table
#'@export
#'

plotVAF = function (gr, genes = NULL, groupByAAchanges = NULL, bySingleLocus = NULL, title = NULL, SNPonly = TRUE)
{
  if(SNPonly)
  {
    vartype = getVarType(gr)
    gr = gr[which(vartype == "SNP")]
  }
  
  if ( !is.null(groupByAAchanges) & !is.null(genes) )
  {
    gr = subset(gr, gr$SYMBOL %in% genes) %>% subset(CHANGE %in% groupByAAchanges)
    if (is.null(genes))
    {
      stop('In case specified AAchanges, you must also specify genes.')
    }
    if (length(genes) == 1)
    {
      gr$plotVafIndex = gr$CHANGE
    } else
    {
      gr$plotVafIndex = str_c(gr$SYMBOL ,"_", gr$CHANGE)
    }
  } else if (!is.null(bySingleLocus))
  {
    if (is.null(genes))
    {
      stop('In case specified locus, you must also specify genes.')
    }
    gr = subset(gr, gr$SYMBOL %in% genes) %>% subset(POS %in% bySingleLocus)
    if (length(genes) == 1)
    {
      gr$plotVafIndex = gr$CHANGE
    } else
    {
      gr$plotVafIndex = str_c(gr$SYMBOL ,"_", gr$CHANGE)
    }
  } else if (!is.null(genes))
  {
    gr$plotVafIndex = gr$SYMBOL
    gr = subset(gr, gr$SYMBOL %in% genes)
  } else if (is.null(genes))
  {
    gr$plotVafIndex = gr$SYMBOL
  }
  
  if(length(unique(gr$plotVafIndex)) > 502)
  {
    stop ("Too many unique genes or variants")
  }
  if (is.null(title))
  {
    plotVafwithGR(gr)
  } else
  {
    plotVafwithGR(gr, main_title = title)
  }
}