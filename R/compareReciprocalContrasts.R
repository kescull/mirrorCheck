# newFromList function was sourced from github UpsetR issue 85, solution 
# provided by 'docmanny' in September 2017
newFromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

collapse_reciprocal_results <- function(group1,group2, folder) {
  #get csvs
  fns <- c(paste0(group1,"_vs_",group2,"_DG.csv"),paste0(group2,"_vs_",group1,"_DG.csv"))
  simple.titles <- c(paste0("Against_",group2),paste0("Against_",group1))
  paths <- file.path(folder,fns)
  stopifnot(all(file.exists(paths)))
  
  reciprocal.results <- lapply(paths,read.csv)
  
  names(reciprocal.results) <- simple.titles
  reciprocal.results <- lapply(reciprocal.results, function(x) x %>% 
                                 dplyr::select(-X) %>% 
                                 dplyr::mutate(reg = dplyr::if_else(log2FoldChange > 0, "UP","DOWN")))
  #Venn
  gene.lists <- lapply(reciprocal.results,function(x) x %>% dplyr::pull(name))
  invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
  venn.plot <- VennDiagram::venn.diagram(x = gene.lists, 
                            category.names = names(gene.lists),
                            filename = NULL,
                            imagetype="png",
                            cat.prompts=TRUE,
                            margin=0.2,
                            cat.dist=rep(0.1,length(gene.lists)),
                            ext.dist=rep(-0.1,length(gene.lists)))
  grid::grid.newpage()
  grid::grid.draw(venn.plot)
  
  partitions <- VennDiagram::get.venn.partitions(gene.lists)
  
  #merge results tables
  merged <- reciprocal.results %>% 
    purrr::reduce(dplyr::full_join, by="name", suffix = simple.titles) %>%
    tidyr::unite(combo, tidyselect::starts_with("reg"), remove = F) %>%
    dplyr::mutate(agree = dplyr::if_else(combo == "UP_DOWN" | combo == "DOWN_UP", T,F),
           LFC.diff = rowSums(dplyr::across(tidyselect::starts_with("log2FoldChange"))),
           padj.diff = purrr::reduce(dplyr::across(tidyselect::starts_with("padj")), `-`),
           group = paste(group1,group2,sep=".")) %>%
    dplyr::select(-combo)
  
  overlap <- merged %>% dplyr::filter(!is.na(padj.diff))
  bad <- overlap %>% dplyr::filter(agree == F)
  if (!nrow(bad) == 0) {
    print("Found some that didn't agree on the direction of regulation")
    write.csv(bad,file.path(folder,paste0("direction_error_",group1,"_v_",group2,".csv")))
  }
  if (nrow(overlap) == 0) {
    print(paste("No overlap of regulated genes using different reference level for",
                group1,"and",group2))
  } else {
    write.csv(overlap,file.path(folder,paste0("consensusDG_",group1,"_v_",group2,".csv")))
  }
  
  for.plotting <- merged %>% 
    dplyr::mutate(partition = dplyr::case_when(!is.na(padj.diff) ~ "consensus",
                                 name %in% partitions$..values..[[2]] ~ "group1ref",
                                 name %in% partitions$..values..[[3]] ~ "group2ref"))
}

compare_pairwise_recursively <- function(groups, i, j,folder,res.list=NULL,list.ind = 0) {
  if (is.null(res.list)) {
    res.list <- list()
  }
  if (i > 1) {
    if (j > 0 ) {
      list.ind <- list.ind + 1
      res.list[[list.ind]] <- collapse_reciprocal_results(groups[i],groups[j],folder)
      res.list <- Recall(groups,i,j-1,folder, res.list, list.ind)
    } else {
      res.list <- Recall(groups,i-1,i-2,folder,res.list, list.ind)
    }
  } 
  res.list
}

#' Compare results from reciprocal contrasts
#'
#' This finds all the 'pairs' among the DG csvs created by
#' [run_DESeq_all_contrasts()] and looks for differences in the results found
#' when using each level as the reference level.
#'
#' @section Rationale:
#' 
#'   Let the pair of levels in question be called group1 and group2. Then the
#'   genes found to be upregulated when using group1 as the reference should be
#'   downregulated when group2 is the reference, and vice versa. Also the
#'   magnitude of the log fold change for a gene should be nearly equal for both
#'   sets of results, although with the opposite sign; thus the sum of the two
#'   log fold changes for each gene is expected to be distributed around zero.
#'   Furthermore, there should be a strong consensus set of significantly
#'   up/down-regulated genes i.e. high degree of overlap between the genes seen
#'   using group1 as the reference and those seen using group2. However,
#'   [lfcShrink()] can sometimes produce quite different results for a contrast
#'   depending which level you use as reference, and from our observation is
#'   then more likely to find upregulated genes than down-regulated genes, which
#'   leads to low consensus between the results sets and sum(log fold changes)
#'   density plots which are skewed/shifted to the right. This function produces
#'   1. simple Venn diagrams showing the overlap in results for the pairwise
#'   comparisons, collected in one pdf named Venns.pdf; 
#'   2) diagnostic plots to visualise how well the analysis worked, including 
#'   density plots of the sums of the log fold changes per gene, and stacked bar 
#'   plots summarising the Venn diagrams so that the strength of the consensus 
#'   can be judged readily;
#'   3) csv tables containing only the consensus genes for each pair of levels;
#'   4) Finally, the overlap in the identities of significant genes between
#'   different pairwise contrasts may be of interest. The function uses UpSetR
#'   to plot and print the intersections between the consensus DEG sets from
#'   each pairwise contrast.
#' @param groups character vector specifying the names of the levels in the
#'   factor used for the DESeq analysis by [run_DESeq_all_contrasts()]. i.e. if
#'   [run_DESeq_all_contrasts()] had the parameters `dds = dds, condition = 'set'`,
#'   then here you could use `groups = levels(dds$set)`
#' @param folder character string with the file path to the folder of results
#'   from [run_DESeq_all_contrasts()]. Output figures and tables will also be
#'   deposited here.
#' @return named list of data.frames representing merged tables from each 'pair'
#'   of contrasts. Also creates pdf and csv files in the specified folder.
#' @export
compare_reciprocal_contrasts <- function(groups, folder) {
  # Compare all reciprocal pairs to get consensus DEGs and print Venns
  pdf(file.path(folder,"Venns.pdf"))
  res.list <- compare_pairwise_recursively(levels(as.factor(groups)),
                                           length(groups),
                                           length(groups)-1,
                                           folder)
  dev.off()
  names(res.list) <- lapply(res.list,function (x) x$group[[1]])
  #How many significant genes overall? Make UpSet and intersections (if there is >1 contrast)
  if (length(res.list)>1) {
    listInput <- lapply(res.list, function(x) x %>%
                        dplyr::filter(!is.na(padj.diff)) %>%
                        dplyr::pull(name))
    originalSetNames <- names(listInput)
    setNames <- gsub("\\.","_v_", originalSetNames)
    names(listInput) <- setNames
    intersections <- newFromList(listInput)
    write.csv(intersections, 
              file = file.path(folder,"all_contrasts_DEG_intersections.csv"), 
              na="")
    setNames <- gsub("\\.","\n", originalSetNames)
    names(listInput) <- setNames
  
    upsetPlot <- UpSetR::upset(UpSetR::fromList(listInput), nsets = length(listInput),
                       order.by = "freq",mainbar.y.label = "DEG Intersections", 
                       show.numbers = "no", mb.ratio = c(0.6, 0.4),
                       sets.x.label = "DEGs Per Contrast")
  }  
  # Plot all graphs of summed LFC per gene from each contrast using alternate group as reference level 
  ready.for.lfc.check <- lapply(res.list,function(x) x %>% 
                                  dplyr::filter(!is.na(padj.diff)) %>% 
                                  dplyr::select(LFC.diff,group))
  all <- dplyr::bind_rows(ready.for.lfc.check)
  all$group <- as.factor(all$group) 
  group.levels <- levels(all$group) 
  group.labels <- gsub("\\.","\n", group.levels)
  names(group.labels) <- group.levels
  xlim <- max(abs(min(all$LFC.diff)),abs(max(all$LFC.diff)))
  if (xlim > 0.5) xlim = 0.5
  densityPlot <- ggplot2::ggplot(all, ggplot2::aes(x = LFC.diff, group = group)) +
    ggplot2::geom_density() + 
    ggplot2::facet_wrap(~group, ncol = 3, labeller = ggplot2::labeller(group = group.labels)) +
    ggplot2::xlim(c(-xlim,xlim)) +
    ggplot2::xlab("Sum of the 2 LFCs for the contrast\nwhen using each group as reference level in turn") 
  
  # Plot stacked bar charts of overlap for each contrast
  ready.for.overlap.check <- lapply(res.list,function(x) x %>% 
                                      dplyr::select(partition) %>%
                                      dplyr::mutate(total = dplyr::n()) %>%
                                      dplyr::group_by(partition) %>%
                                      dplyr::summarise(n = dplyr::n(),
                                                total = sum(total)/n,
                                                pc = n/total *100))
  
  all.overlap <- dplyr::bind_rows(ready.for.overlap.check, .id = "group") 
  
  group.labels <- gsub("\\.","\n", all.overlap$group)
  names(group.labels) <- all.overlap$group
  legend.labels <- c("consensus","only with group1 as ref","only with group2 as ref")
  
  barPlot <- ggplot2::ggplot(all.overlap, 
                             ggplot2::aes(x = group, y = pc, fill = partition)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab(ggplot2::element_blank()) +
    ggplot2::ylab("Percentage of differentially expressed genes") +
    ggplot2::scale_x_discrete(labels=group.labels) +
    ggplot2::scale_fill_discrete(labels = legend.labels) + 
    ggplot2::geom_text(ggplot2::aes(label = n),
                       position = ggplot2::position_stack(vjust = 0.5),size = 3) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = ggplot2::element_blank(), 
                                                 reverse=T)) +
    ggplot2::scale_fill_manual(values = c('#EE7733', '#0077BB', '#33BBEE')) +
    ggplot2::coord_flip() 
  
  #Print output graphs
  pdf(file.path(folder,
                "diagnostic plots for contrasts with different ref level.pdf"), 
      paper = "a4")
  print(densityPlot)
  print(barPlot)
  if (length(res.list) > 1) {
    print(upsetPlot)
  }
  dev.off()
  #return data
  res.list
}
