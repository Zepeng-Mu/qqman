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
#' @param summ.1.name An expression for the name of summary statistics 1, "GWAS -log10(P)" by default.
#' @param summ.2.name Similar to \code{summ.1.name}, "eQTL -log10(P)" by default.
#' @param col A character vector indicating which colors to alternate.
#' @param coloc.snp A character of colocalized SNP.
#' @param coloc.gene A character of colocalized gene.
#' @param PP4 A numeric of PP4 value from colocalization analysis to be labeled in the plot.
#' @param logp If TRUE, the -log10 of the p-value is plotted. It isn't very 
#'   useful to plot raw p-values, but plotting the raw value could be useful for
#'   other genome-wide plots, for example, peak heights, bayes factors, test 
#'   statistics, other "scores," etc.
#' @param size Point size passed to ggplot2
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

ggColocManhattanGrid <- function(
    summ.1, summ.2, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
    summ.1.name = expression(paste("GWAS", -log[10](P))),
    summ.2.name = expression(paste("eQTL", -log[10](P))),
    col = c("gray50", "orange"), coloc.snp = "", coloc.gene = "",
    PP = NA, hypothesis = NA, logp = TRUE, size = 2, lead.snp = NULL, r2 = NULL, ...
) {
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
    
    # Map SNP r2 to discrete colors
    # This is inspired by locuscomparer package
    snp.col <- as.character(cut(r2, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                labels=c('blue4','skyblue','darkgreen','orange','red'),
                                include.lowest=TRUE))
    snp.col[which(d1$BP == lead.snp)] <- "purple"
    
    # Create shapes for lead SNP and other SNPs
    snp.shape <- ifelse(d1$BP == lead.snp, 18, 20)
    
    # Creat size for lead SNP and other SNPs
    snp.size <- ifelse(d1$BP == lead.snp, 4, 2)
    
    g1 <- ggplot(d1, aes(x = BP / 1e6, y = logp)) +
        theme_classic(base_size = 12, base_line_size = 1) +
        geom_point(col = snp.col, shape = snp.shape, size = snp.size) +
        geom_point(aes(x = d1$BP[which(d1$BP == lead.snp)] / 1e6,
                       y = d1$logp[which(d1$BP == lead.snp)]),
                   col = "purple", shape = 18, size = 4) +
        xlab("") + # do not show xlab for d1 data, which is above d2 data
        ylab(summ.1.name)
    
    g2 <- ggplot(d2, aes(x = BP / 1e6, y = logp)) +
        theme_classic(base_size = 12, base_line_size = 1) +
        geom_point(col = snp.col, shape = snp.shape, size = snp.size) +
        geom_point(aes(x = d2$BP[which(d2$BP == lead.snp)] / 1e6,
                       y = d2$logp[which(d2$BP == lead.snp)]),
                   col = "purple", shape = 18, size = 4) +
        xlab(paste("Chromosome", unique(d2$CHR), "(Mb)")) +
        ylab(summ.2.name)
    
    # Create a data frame for scatter plot between d1 and d2
    d3 <- data.frame(pval1 = d1$logp, pval2 = d2$logp,
                     stringsAsFactors = F)
    
    g3 <- ggplot(d3, aes(pval1, pval2)) +
        theme_classic(base_size = 12, base_line_size = 1) +
        geom_point(col = snp.col, shape = snp.shape, size = snp.size) +
        geom_point(aes(x = d3$pval1[which(d1$BP == lead.snp)],
                       y = d3$pval2[which(d1$BP == lead.snp)]),
                   col = "purple", shape = 18, size = 4) +
        xlab(summ.1.name) +
        ylab(summ.2.name)
    
    # Make legend, also from locuscomparer
    legend_box <- data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
    g3 <- ggdraw(g3) +
        geom_rect(data = legend_box,
                  aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
                  color = "black",
                  fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
        draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
        draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) +
        draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) +
        draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
        draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 10)
    
    title <- ggdraw() + 
        draw_label(paste0(coloc.gene, " - ", coloc.snp, " (PP"%&%hypothesis%&%": ", PP, ")"),
                   x = 0,
                   hjust = 0) +
        theme(plot.margin = margin(0, 0, 0, 7))
    
    # Combine the three subplots
    p1 <- plot_grid(g1, g2, align = "v", nrow = 2)
    p2 <- plot_grid(g3, p1)
    p3 <- plot_grid(title, p2, ncol = 1,
                    rel_heights = c(0.1, 1))
    
    return(p3)
}
