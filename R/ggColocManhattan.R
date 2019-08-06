#' Creates one Manhattan for two association tests (e.g. GWAS and eQTL)
#' 
#' Creates one Manhattan for two association tests (e.g. GWAS and eQTL) (or any data frame with 
#' chromosome, position, and p-value). This is suitable for visualizing colocalization results.
#' 
#' @param summ.1 A data.frame for the first summary statistics with columns "BP," "CHR," "P," and optionally, "SNP."
#' @param summ.2 A data.frame for the second summary statistics with columns "BP," "CHR," "P," and optionally, "SNP."
#' @param chr A string denoting the column name for the chromosome. Defaults to 
#'   PLINK's "CHR." Said column must be numeric. If you have X, Y, or MT 
#'   chromosomes, be sure to renumber these 23, 24, 25, etc.
#' @param bp A string denoting the column name for the chromosomal position. 
#'   Defaults to PLINK's "BP." Said column must be numeric.
#' @param p A string denoting the column name for the p-value. Defaults to 
#'   PLINK's "P." Said column must be numeric.
#' @param snp A string denoting the column name for the SNP name (rs number). 
#'   Defaults to PLINK's "SNP." Said column should be a character.
#' @param summ.1.name A string for the name of summary statistics 1, "GWAS" by default.
#' @param summ.2.name Similar to \code{summ.1.name}, "eQTL" by default.
#' @param col A character vector indicating which colors to alternate.
#' @param coloc.snp A character of colocalized SNP.
#' @param coloc.gene A character of colocalized gene.
#' @param PP4 A numeric of PP4 value from colocalization analysis to be labeled in the plot.
#' @param logp If TRUE, the -log10 of the p-value is plotted. It isn't very 
#'   useful to plot raw p-values, but plotting the raw value could be useful for
#'   other genome-wide plots, for example, peak heights, bayes factors, test 
#'   statistics, other "scores," etc.
#' @param ... Arguments passed on to other plot/points functions
#'   
#' @return A ggplot object coloc manhattan plot.
#'   
#' @keywords visualization manhattan colocalization
#'   
#' @import utils
#' @import graphics
#' @import stats
#' @import ggplot2
#' @import ggrepel
#' 
#' @examples
#' ggColocManhattan(gwasResults)
#'   
#' @export

ggColocManhattan <- function(summ.1, summ.2, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
                             summ.1.name = "GWAS", summ.2.name = "eQTL",
                             col = c("gray50", "orange"), coloc.snp = "", coloc.gene = "",
                             PP4 = NA, logp = TRUE, ...) {
    # Check for sensible dataset
    ## Make sure you have chr, bp and p columns.
    if (!(chr %in% names(summ.1))) stop(paste("Column", chr, "not found in summ.1!"))
    if (!(bp %in% names(summ.1))) stop(paste("Column", bp, "not found in summ.1!"))
    if (!(p %in% names(summ.1))) stop(paste("Column", p, "not found in summ.1!"))
    if (length(unique(summ.1[[chr]])) > 1) stop("You seem to have more than one chromosome in summ.1!")
    ## warn if you don't have a snp column
    if (!(snp %in% names(summ.1))) warning(paste("No SNP column found in summ.1. OK unless you're trying to highlight."))
    ## make sure chr, bp, and p columns are numeric.
    if (!is.numeric(summ.1[[chr]])) stop(paste(chr, "column in summ.1 should be numeric. Do you have 'X', 'Y', 'MT',etc? If so change to numbers and try again."))
    if (!is.numeric(summ.1[[bp]])) stop(paste(bp, "column in summ.1 should be numeric."))
    if (!is.numeric(summ.1[[p]])) stop(paste(p, "column in summ.1 should be numeric."))
    
    ## Do the same check for summ.2
    if (!(chr %in% names(summ.2))) stop(paste("Column", chr, "not found in summ.2!"))
    if (!(bp %in% names(summ.2))) stop(paste("Column", bp, "not found in summ.2!"))
    if (!(p %in% names(summ.2))) stop(paste("Column", p, "not found in summ.2!"))
    if (length(unique(summ.2[[chr]])) > 1) stop("You seem to have more than one chromosome in summ.2!")
    ## warn if you don't have a snp column
    if (!(snp %in% names(summ.2))) warning(paste("No SNP column found in summ.2. OK unless you're trying to highlight."))
    ## make sure chr, bp, and p columns are numeric.
    if (!is.numeric(summ.2[[chr]])) stop(paste(chr, "column in summ.2 should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(summ.2[[bp]])) stop(paste(bp, "column in summ.2 should be numeric."))
    if (!is.numeric(summ.2[[p]])) stop(paste(p, "column in summ.2 should be numeric."))
    
    # Create a new data frame from the input dataframe
    # If the input data frame has a SNP column, add it to the new data frame you're creating.
    if (!is.null(summ.1[[snp]])) {
        d1 <- data.frame(CHR = summ.1[[chr]], BP = summ.1[[bp]], P = summ.1[[p]], pos  =  NA,
                        index = NA, SNP = summ.1[[snp]], stringsAsFactors = FALSE)
    } else {
        d1 <- data.frame(CHR = summ.1[[chr]], BP = summ.1[[bp]], P = summ.1[[p]], pos  =  NA,
                         index = NA, stringsAsFactors = FALSE)
    }
    
    # Do the same for summ.2
    if (!is.null(summ.2[[snp]])) {
        d2 <- data.frame(CHR = summ.2[[chr]], BP = summ.2[[bp]], P = summ.2[[p]],
                         index = NA, SNP = summ.2[[snp]], stringsAsFactors = FALSE)
    } else {
        d2 <- data.frame(CHR = summ.2[[chr]], BP = summ.2[[bp]], P = summ.2[[p]],
                         index = NA, stringsAsFactors = FALSE)
    }
    
    # Calculate -log(P) for summ.1 and summ.2
    if (logp) {
        d1$logp <- -log10(d1$P)
        d2$logp <- -log10(d2$P)
    } else {
        d1$logp <- d1$P
        d2$logp <- d2$P
    }

    ## Create a data frame combining d1 and d2 to be used in ggplot2
    d1$study.name <- summ.1.name
    d2$study.name <- summ.2.name
    
    labels1 <- as.character(c(0, round(max(d1$logp), 0)))
    labels2 <- as.character(c(0, round(max(d2$logp), 0)))
    
    g <- ggplot() +
        theme_classic(base_size = 12, base_line_size = 1) +
        xlab(paste("Chromosome", unique(d1$CHR))) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5, size = 13),
              axis.line = element_line(color = "grey80"),
              axis.ticks = element_line(colour = "grey80"),
              axis.title.y.left = element_text(color = col[1]),
              axis.title.y.right = element_text(color = col[2])) +
        ggtitle(paste0(coloc.gene, " - ", coloc.snp, " (PP4: ", PP4, ")"))
    
    g <- g +
        geom_point(data = d1, mapping = aes(x = BP, y = rescale(logp)),
                   color = col[1], alpha = 0.4, size = 2, shape = 16) +
        geom_point(data = d2, mapping = aes(x = BP, y = rescale(logp)),
                   color = col[2], alpha = 0.4, size = 2, shape = 16) +
        scale_y_continuous(name = "GWAS -log10(P)",
                           breaks = c(0, 1), labels = labels1, position = "left",
                           sec.axis = sec_axis(trans = ~ ., name = "eQTL -log10(P)",
                                               breaks = c(0, 1),
                                               labels = labels2))
    
    return(g)
}
