#' @title Fold Change for phyloseq Objects
#' @description Calculate Log10 Fold Change for two-group comparison with phyloseq objects.
#' @param x \code{\link{phyloseq-class}} object or a data matrix 
#'            (features x samples; eg. HITChip taxa vs. samples)
#' @param group Vector with specifying the groups
#' @param sort sort the results
#' @param paired Paired comparison (Default: FALSE)
#' @return Fold change information for two-group comparison.
#' @examples 
#'   data(dietswap)
#'   fc <- check_foldchange(dietswap, "sex")
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
check_foldchange <- function (x, group, sort = FALSE, paired = FALSE) {
  
  # Pick the grouping variable from sample metadata
  g <- group
  if (length(g) == 1) {
    g <- sample_data(x)[[g]]
    
    if (!is.factor(g)) {
      warning(paste("Converting the grouping variable", group, "into a factor."))
      g <- as.factor(g)
    }    
    g <- droplevels(g)
    if (!length(levels(g)) == 2) {
      stop(paste("check_foldchange is currently implemented only for two-group comparisons. The selected variable", group, "has", length(unique(g)), "levels: ", paste(unique(g), collapse = "/")))
    }
  }
  
  if (class(x) == "phyloseq") {    
    x <- abundances(x)
  }
  
  # Calculate fold changes
  if (paired) {
    fc <- apply(x, 1, function (xi) {spl <- split(xi, g);
    log10(mean(spl[[2]] - spl[[1]], na.rm = TRUE))
    })
  } else {
    fc <- apply(x, 1, function (xi) {spl <- split(xi, g);
    log10(mean(spl[[2]], na.rm = T)) - log10(mean(spl[[1]], na.rm = TRUE))
    })
    
  }
  
  # Sort
  if (sort) {
    fc <- sort(fc)
  }
  
  fc
  
}

#' @title Check Wilcoxon
#' @description Calculate Wilcoxon test (with BH correction) for two-group comparison for phyloseq objects.
#' @param x \code{\link{phyloseq-class}} object or a data matrix 
#'            (features x samples; eg. HITChip taxa vs. samples)
#' @param group Vector with specifying the groups
#' @param p.adjust.method p-value correction method for p.adjust function 
#'               (default 'BH'). For other options, see ?p.adjust
#' @param sort sort the results
#' @param paired Paired comparison (Default: FALSE)
#' @return Corrected p-values for two-group comparison.
#' @examples 
#'   #pseq <- download_microbiome("peerj32")$physeq
#'   #pval <- check_wilcoxon(pseq, "gender")
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
check_wilcoxon <- function (x, group, p.adjust.method = "BH", sort = FALSE, paired = FALSE) {
  
  # Also calculate fold changes
  fc <- check_foldchange(x, group, paired = paired)
  
  # Pick the grouping variable from sample metadata
  g <- group
  if (length(g) == 1) {
    g <- sample_data(x)[[g]]
    
    if (!is.factor(g)) {
      warning(paste("Converting the grouping variable", group, "into a factor."))
      g <- as.factor(g)
    }    
    g <- droplevels(g)
    if (!length(levels(g)) == 2) {
      stop(paste("check_wilcoxon is valid only for two-group comparisons. The selected variable", group, "has", length(unique(g)), "levels: ", paste(unique(g), collapse = "/")))
    }
  }
  
  if (class(x) == "phyloseq") {    
    x <- log10(otu_table(x)@.Data)
  }
  
  # Calculate Wilcoxon test with BH p-value correction for gender
  pval <- suppressWarnings(apply(x, 1, function (xi) {wilcox.test(xi ~ g, paired = paired)$p.value}))
  
  # Multiple testing correction
  pval <- p.adjust(pval, method = p.adjust.method)
  
  # Sort p-values
  if (sort) {
    pval <- sort(pval)
  }
  
  data.frame(list(p.value = pval, fold.change.log10 = fc[names(pval)]))
  
  
}



