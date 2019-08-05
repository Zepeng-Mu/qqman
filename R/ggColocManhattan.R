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
#'   This SNP should be in your dataset.
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
                             col = c("gray80", "gray90"),
                             suggestiveline = -log10(1e-7), genomewideline = NULL, 
                             highlight1 = NULL, highlight2 = NULL,
                             text1 = NULL, text2 = NULL, logp = TRUE, annotatePval = NULL, ...) {
    # Check for sensible dataset
    ## Make sure you have chr, bp and p columns.
    if (!(chr %in% names(summ.1))) stop(paste("Column", chr, "not found in summ.1!"))
    if (!(bp %in% names(summ.1))) stop(paste("Column", bp, "not found in summ.1!"))
    if (!(p %in% names(summ.1))) stop(paste("Column", p, "not found in summ.1!"))
    if (length(unique(summ.1[[chr]] > 1))) stop("You seem to have more than one chromosome in summ.1!")
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
    if (length(unique(summ.2[[chr]] > 1))) stop("You seem to have more than one chromosome in summ.2!")
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
    
    d1$index = rep.int(seq_along(unique(d1$CHR)), times = tapply(d1$SNP, d1$CHR, length))
    d2$index = rep.int(seq_along(unique(d2$CHR)), times = tapply(d2$SNP, d2$CHR, length))
    
    ## Create a data frame combining d1 and d2 to be used in ggplot2
    d1$study.name <- summ.1.name
    d2$study.name <- summ.2.name
    d <- rbind(d1, d2)
    
    g <- ggplot() +
        theme_classic(base_size = 18, base_line_size = 0.5) +
        xlab("Position") +
        scale_y_continuous(name = paste(summ.2.name, expression(-log[10](P))),
                           breaks = , labels = , position = "right") + 
        geom_text(aes(x = min(d$pos), y = max(d$logp), label = "GWAS"),
                  vjust = "inward", hjust = "inward", size = 6, col = "tomato") +
        geom_text(aes(x = max(d$pos), y = max(d$logp), label = "eQTL"),
                  vjust = "inward", hjust = "inward", size = 6, col = "purple") +
        theme(legend.position = "none")
    
    # Add points to the plot
    g <- g +
        geom_point(data = d1, mapping = aes(x = BP, y = rescale(logp)),
                   color = col[1], alpha = 0.3, size = 0.6) +
        scale_y_continuous(name = paste(summ.1.name, expression(-log[10](P))),
                           breaks = , labels = , position = "left") +
        geom_point(data = d2, mapping = aes(x = BP, y = rescale(logp)),
                   color = col[2], alpha = 0.3, size = 0.6) +
        scale_y_continuous(name = paste(summ.2.name, expression(-log[10](P))),
                           breaks = , labels = , position = "right")
    
    # Add gene name labels to the plot
    if (!is.null(text1)) {
        if (any(!(text1$SNP %in% d$SNP))) {
            warning("You're trying to highlight SNPs that don't exist in your results.")
        }
        text1 <- left_join(text1, highlight1, by = "SNP")
        g <- g + geom_text_repel(data = text1, mapping = aes(x = pos, y = logp, label = gene),
                                 color = "grey40", vjust = 0, size = 4,
                                 nudge_y = max(text1$logp) - text1$logp)
    }
    
    return(g)
}
