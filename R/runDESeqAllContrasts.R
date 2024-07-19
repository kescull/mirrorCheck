plot_Volcano <- function(r,title, rowname2symbol = NULL, 
                         p.cutoff = 0.1, fc.cutoff = 2, top.n = 30) {
  padj <- as.data.frame(r) %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(log10padj = -log10(padj),
           diff.expressed = dplyr::if_else(log2FoldChange > fc.cutoff & padj < p.cutoff, "UP",
                                           dplyr::if_else(log2FoldChange < -fc.cutoff & padj < p.cutoff, 
                                            "DOWN","NO"))) %>%
    tibble::rownames_to_column("rowname")
  if (is.null(rowname2symbol)) {
    padj <- padj %>%
      dplyr::mutate(g_symbol = rowname)
  } else {
    padj <- padj  %>%
      dplyr::left_join(rowname2symbol)
  }
  diff.ex <- padj %>% dplyr::filter(diff.expressed != "NO")
  
  sums <- table(diff.ex$diff.expressed)
  maxy <- max(padj$log10padj)
  minx <- min(padj$log2FoldChange)
  maxx <- max(padj$log2FoldChange)
  padj$labels <- ifelse(padj$rowname %in% head(diff.ex[order(diff.ex$padj), 
                                                        "rowname"], top.n), 
                        padj$g_symbol, NA)
  group.colors <- c(DOWN = "#00AFBB",NO = "grey",UP = "#FFDB6D")
  group.labels <- c(DOWN = "Downregulated",NO = "Not significant", UP = "Upregulated")
  
  p <- ggplot2::ggplot(padj, ggplot2::aes(x = log2FoldChange, 
                        y = log10padj, 
                        col = diff.expressed,
                        label = labels)) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = title,
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adj-pvalue")) +
    ggplot2::geom_vline(xintercept = c(-fc.cutoff, fc.cutoff), col = "gray", linetype = 'dashed') +
    ggplot2::geom_hline(yintercept = -log10(p.cutoff), col = "gray", linetype = 'dashed') +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(values = group.colors, 
                       labels = group.labels) +
    ggrepel::geom_text_repel(max.overlaps = Inf, size = 3, color = "black", na.rm = T) +
    ggplot2::annotate("text",x = minx + 2, y = maxy - 1, label = paste("DOWN",sums["DOWN"])) +
    ggplot2::annotate("text",x = maxx - 1, y = maxy - 1, label = paste("UP",sums["UP"]))
}

get_diff_genes_list <- function(res, p.cutoff = 0.1, fc.cutoff = 2) {
  diff.ex <- as.data.frame(res)  %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(log10padj = -log10(padj),
           diff.expressed = dplyr::if_else(log2FoldChange > fc.cutoff & padj < p.cutoff, "UP",
                                           dplyr::if_else(log2FoldChange < -fc.cutoff & padj < p.cutoff, 
                                            "DOWN","NO"))) %>%
    tibble::rownames_to_column("g_symbol") %>% 
    dplyr::filter(diff.expressed != "NO")
  diff.ex$g_symbol
}

plot_heatmap <- function(dg.list,title,dds,annotation.col) {
  if (length(dg.list) < 2) {
    return()
  }
  vsd <- DESeq2::vst(dds)
  df <- as.data.frame(SummarizedExperiment::colData(dds)) %>%
    dplyr::select(dplyr::all_of(annotation.col))
  show_col_names <- T
  if (dim(dds)[2] > 25) {
    show_col_names <- F
  }
  p <- pheatmap::pheatmap(SummarizedExperiment::assay(vsd)[rownames(vsd) %in% dg.list,], 
                          cluster_rows=F, show_rownames=FALSE,
                          show_colnames = show_col_names,
                          cluster_cols=T, annotation_col = df, main = title)
}

get_all_results <- function(i,dds, mode = "lfcShrink", alpha = 0.1) {
  coef = DESeq2::resultsNames(dds)[i+1]
  res <- DESeq2::results(dds, 
                         name = DESeq2::resultsNames(dds)[i+1], 
                         alpha = alpha)
  if (mode == "lfcShrink") {
    res <- DESeq2::lfcShrink(dds, coef=i+1, type="apeglm", res = res)
  }
  resOrdered <- res[order(res$padj),]
}

get_all_MA_plots <- function(res,coef,p.cutoff = 0.1) {
  MA <- DESeq2::plotMA(res, ylim=c(-3,3), main = coef,alpha = p.cutoff)
}

write_diff_genes <- function(res,dg, coef,rowname2symbol = NULL, folder) {
  res.sig <- as.data.frame(res) %>% 
    tibble::rownames_to_column("name") %>% 
    dplyr::filter(name %in% dg)
  if (!is.null(rowname2symbol)) {
    rowname2symbol <- rowname2symbol %>% dplyr::rename(name = rowname)
    res.sig <- res.sig %>% 
      dplyr::left_join(rowname2symbol)
  }
  fn <- file.path(folder,paste0(coef[1],"_DG.csv"))
  write.csv(res.sig,file = fn)
}

# taken from is.integer() documentation
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}

#' Run DESeq for all contrasts
#'
#' Runs [DESeq()] for a `DESeqDataSet` with a simple design. It uses each
#' possible reference level in turn, extracting all contrast results for each.
#'
#' @section Further details: Output includes pdf files with DESeq2 MA plots,
#'   volcano plots and heatmaps of significant genes for each contrast in one
#'   pdf file per reference level, and also a csv of the significant genes for
#'   each contrast (according to the user-specified adj-pvalue and lfc
#'   thresholds). It is possible to run with more complex designs for the first
#'   reference level using the `only.first.iteration` parameter, but this has
#'   not been tested thoroughly and output should be interpreted with care.
#' @param dds Valid `DESeqDataSet` object with a simple design e.g. `design =
#'   ~condition`. If you wish to use a filter such as `edgeR` [filterByExpr()],
#'   please do so prior to running this function. This package is intended to
#'   help assess the usefulness of such filtering.
#' @param folder character string specifying the relative or absolute filepath
#'   of the folder where you want the output files to go. This folder should
#'   exist.
#' @param mode either `"noShrink"` or `"lfcShrink"`. Specifies whether to use
#'   [results()] or [lfcShrink()] (default is `"lfcShrink"`).
#' @param condition character string related to the design formula used when
#'   creating the `DESeqDataSet`. This should be a `colData` variable in the
#'   `DESeqDataSet`; e.g. for `design=~group`, here enter `condition="group"`.
#'   If you use a more complex design, results are not guaranteed. See also
#'   `only.first.iteration`.
#' @param rowname2symbol data.frame or tibble with 2 columns titled `"rownames"`
#'   and `"g_symbol"`. `"rownames"` should match `"rownames(dds)"`. This is used
#'   for labeling the Volcano plot and annotating the tables of differentially
#'   expressed genes
#' @param heatmap.annotation.col character vector specifying which `colData`
#'   variables to use as groups in the heatmap output. If NULL, the heatmap will
#'   use the variable specified in `condition` (default NULL)
#' @param p.cutoff numeric. The significance cut-off for p-adjusted values,
#'   equivalent to parameter alpha in [DESeq2::results()]. (default 0.1)
#' @param fc.cutoff numeric. The log2 fold change significance cut-off; this is
#'   used in the Volcano plots and when printing out significant genes. It is
#'   NOT used to change the null hypothesis in [DESeq2::results()]. (default 2)
#' @param top.n whole number. Number of most highly regulated genes to label
#'   with gene symbols on Volcano plots. (default 30).
#' @param only.first.iteration If True, the function will process all the
#'   contrasts using the first level as the reference level, then exit. This
#'   should produce 1 set of csvs for contrast DGEs and 1 pdf with MA and
#'   Volcano plots. It may be helpful for looking at more complex designs, so
#'   long as `condition` can still be specified as a `colData` variable. But no
#'   promises. (default False)
#' @return None - but creates pdf and csv files in the specified folder.
#' @export
run_DESeq_all_contrasts <- function(dds,folder, 
                                    mode = "lfcShrink", 
                                    condition,
                                    rowname2symbol = NULL,
                                    heatmap.annotation.col = NULL,
                                    p.cutoff = 0.1, 
                                    fc.cutoff = 2, 
                                    top.n = 30,
                                    only.first.iteration = F) {
  if (!mode %in% c("lfcShrink","noShrink")) {
    stop("mode must be lfcShrink or noShrink (default is lfcShrink)")
  }
  
  if (!is.numeric(p.cutoff) | !is.numeric(fc.cutoff) | !is_wholenumber(top.n)) {
    stop(paste(p.cutoff,"p.cutoff and fc.cutoff must be numeric; top.n must be a whole number"))
  }
  if (!condition %in% colnames(SummarizedExperiment::colData(dds))) {
    stop("'condition' should be a variable from colData(dds)")
  }
  groups <- levels(dds[[ {{condition}} ]])
  if(any(!heatmap.annotation.col %in% colnames(SummarizedExperiment::colData(dds)))) {
    stop("heatmap.annotation.col must be a vector of character strings chosen from among those found by calling colData() on the given DESeqDataSet")
  }
  if (!file.exists(file.path(folder))) {
    stop("specified output folder should exist: please make it or check your filepath")
  }
  if (isTRUE(all.equal(DESeq2::design(dds), ~1 ))) {
    stop("DESeq2DataSet needs a design: please enter one (not ~1) and try again")
  }
  if (!is.null(rowname2symbol) & !all(c("rowname","g_symbol") %in% colnames(rowname2symbol))) {
      stop("rowname2symbol must have columns titled `rowname` and `g_symbol`")
  }
  for (g in groups) {
    dds[[ {{condition}} ]] <- relevel(dds[[ {{condition}} ]], ref = g)
    if (length(DESeq2::resultsNames(dds)) < 2) {
      #first time, use first group and run DESeq
      dds <- DESeq2::DESeq(dds)
    } else {
      dds <- DESeq2::nbinomWaldTest(dds)
    }
    
    coef <- DESeq2::resultsNames(dds)[2:length(DESeq2::resultsNames(dds))]
    prefix <- paste0(condition,"_")
    coef <- sub(prefix,"",coef)
    plots.fn <- file.path(folder,paste0("v_",g,"_DGEplots.pdf"))
    
    #get all results and plot
    resList <- lapply(1:length(coef),
                      get_all_results, 
                      dds = dds, mode = mode, alpha = p.cutoff)
    
    pdf(plots.fn, paper = "a4")
    par(mfrow = c(3,2))
    mapply(get_all_MA_plots,resList,coef, MoreArgs = list(p.cutoff = p.cutoff))
    volcList <- mapply(plot_Volcano,resList,coef,
                       MoreArgs = list(rowname2symbol = rowname2symbol,
                                       p.cutoff = p.cutoff, 
                                       fc.cutoff = fc.cutoff, 
                                       top.n = top.n), 
                       SIMPLIFY = F)
    print(volcList)
    dgList <- lapply(resList, get_diff_genes_list, 
                     p.cutoff = p.cutoff, fc.cutoff = fc.cutoff)
    if (!is.null(heatmap.annotation.col)) {
      heatmaps <- mapply(plot_heatmap,dgList,coef,
                         MoreArgs = list(dds = dds,
                                         annotation.col = heatmap.annotation.col), 
                         SIMPLIFY = F)
    } else {
      heatmaps <- mapply(plot_heatmap,dgList,coef,
                         MoreArgs = list(dds = dds,
                                         annotation.col = condition), 
                         SIMPLIFY = F)
    }
    dev.off()
    mapply(write_diff_genes,resList,dgList,coef, 
           MoreArgs = list(rowname2symbol = rowname2symbol, folder = folder))
    if (only.first.iteration) {
      return()
    }
  }  
}
