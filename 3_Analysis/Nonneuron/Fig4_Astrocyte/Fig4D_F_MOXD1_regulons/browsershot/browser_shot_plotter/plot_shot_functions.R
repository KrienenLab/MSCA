## Assemble multiple tracks of publication-quality browser shots using 
## BRGenomics and ggplot2

suppressPackageStartupMessages({
    require(dplyr)
    require(ggplot2)
    require(patchwork) 
    require(BRGenomics)
})

shot_scalebar <- function(region, plot = TRUE, linewidth = 1) {
    ## A track with only a short, labeled scale bar
    
    # get bar size
    bar_len <- 10^(floor(log10(width(region))))
    if (bar_len > 0.5*width(region))
        bar_len <- bar_len / 10
    if (bar_len > 0.25*width(region))
        bar_len <- bar_len / 2
    
    # bar label
    bar_lab <- ifelse(bar_len < 1000, paste(bar_len, "bp"),
                      paste(bar_len/1e3, "kb"))
    
    ggplot() + 
        geom_segment(aes(x = 0, xend = bar_len, y = 0, yend = 0),
                     color = "black", size = linewidth) +
        scale_x_continuous(limits = c(0, width(region)),
                           breaks = c(0, 0.5*bar_len, bar_len), 
                           labels = c("", bar_lab, ""),
                           expand = expansion(mult = c(0.02, 0))) +
        scale_y_continuous(limits = c(0,0), expand = c(0,0)) +
        labs(x = NULL, y = NULL) +
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(color = "black", vjust = 1,
                                         margin = unit(c(0,0,0,0), "points")),
              axis.ticks = element_blank(),
              plot.margin = unit(c(0,0,0,0), "points")) -> p
    
    if (plot) p else list(p)
}



shot_stranded <- function(grl, region, binsize, ylim = NULL, bin_FUN = sum,
                          expand_ranges = TRUE, smooth = FALSE, plot = TRUE,
                          ncores = getOption("mc.cores", 2L)) {
    ## Vertical line plot both strands of stranded data
    
    # list unlisted data
    if (!is.list(grl) && !is(grl, "GRangesList")) {
        sname <- deparse(substitute(grl))
        grl <- list(grl)
        names(grl) <- sname
    }
    
    if (is.null(names(grl))) {
        if (length(grl) > 1L) {
            names(grl) <- paste(deparse(substitute(grl)), "_", seq_along(grl))
        } else {
            names(grl) <- deparse(substitute(grl))
        }
    }
    
    # if ylim set, pre-filter
    if (!is.null(ylim)) {
        if (min(ylim) >= 0L) grl <- lapply(grl, subset, strand == "+")
        if (max(ylim) <= 0L) grl <- lapply(grl, subset, strand == "-")
    }
    
    smooth_fun <- if (smooth) {
        function(signal) smooth.spline(seq_along(signal), signal)$y
    } else identity
    
    # to get counts on each strand, safest to always do in the "same direction"
    # -> split each GRanges by strand and get counts on unstranded region
    strand(region) <- "*"
    map_fun <- function(gr, nm) {
        lapply(split(gr, strand(gr))[1:2], function(x) {
            getCountsByPositions(x, region, binsize = binsize, FUN = bin_FUN, 
                                 expand_ranges = expand_ranges, ncores = 1)
        }) %>%
            lapply(as.vector) %>%
            lapply(smooth_fun) %>%
            (function(x) data.frame(x = (binsize*seq_along(x[[1]])) - binsize,
                                    sig_p =  x[[1]],
                                    sig_m = -x[[2]],
                                    sample = nm))
    }
    mcMap(map_fun, grl, names(grl), mc.cores = ncores) %>%
        c(make.row.names = FALSE) %>%
        do.call(rbind, .) %>% 
        mutate(sample = factor(sample, levels = unique(sample))) -> 
        df
    
    # get y-axis scales
    if (is.null(ylim))
        c(df$sig_p, df$sig_m) %>% 
        pretty %>% 
        (function(y)
            if (min(y) < 0L && max(y) > 0L)
                c(min(y), 0L, max(y)) else range(y)
        ) ->
        ylim
        
    # saturate off-scaled data (for user supplied ylims)
    #   -> doing it this way allows the axis line to be capped
    df[c("sig_p", "sig_m")] <- lapply(df[c("sig_p", "sig_m")], function(y) {
        y[y > max(ylim)] <- max(ylim)
        y[y < min(ylim)] <- min(ylim)
        y
    })
    
    plotfun_str <- function(df) {
        ggplot(df) + 
            facet_grid(rows = vars(sample)) +
            ggplot2::geom_segment(aes(x = x, xend = x, y = 0, yend = sig_p),
                                  color = "red") + 
            ggplot2::geom_segment(aes(x = x, xend = x, y = 0, yend = sig_m),
                                  color = "blue") +
            scale_y_continuous(limits = range(ylim), breaks = ylim, 
                               labels = ylim, expand = c(0, 0)) + 
            scale_x_continuous(limits = c(0, width(region)),
                               expand = expansion(mult = c(0.02, 0))) +
            labs(x = NULL, y = NULL) +
            theme_classic() +
            theme(panel.grid = element_blank(),
                  panel.border = element_blank(),
                  axis.text.y = element_text(color = "black", size = 6),
                  axis.text.x = element_blank(),
                  axis.ticks.y = element_line(color = "black", size = 0.25),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank(),
                  axis.line.y = element_line(color = "black", size = 0.25),
                  strip.text.y = element_text(angle = 0, hjust = 0),
                  strip.background = element_blank(),
                  plot.margin = unit(c(0,0,3,0), "points"))
    }
    
    plots <- lapply(split(df, df$sample), plotfun_str)
    
    if (!plot)
        return(plots)
    do.call(wrap_plots, c(list(plots), ncol = 1))
}



shot_unstranded <- function(grl, region, binsize, ylim = NULL, bin_FUN = sum,
                            expand_ranges = TRUE, smooth = FALSE, plot = TRUE,
                            ncores = getOption("mc.cores", 2L)) {
    ## Vertical line plot of unstranded data
    
    # list unlisted data
    if (!is.list(grl) && !is(grl, "GRangesList")) {
        sname <- deparse(substitute(grl))
        grl <- list(grl)
        names(grl) <- sname
    }
    
    if (is.null(names(grl))) {
        if (length(grl) > 1L) {
            names(grl) <- paste(deparse(substitute(grl)), "_", seq_along(grl))
        } else {
            names(grl) <- deparse(substitute(grl))
        }
    }
    
    smooth_fun <- if (smooth) {
        function(signal) smooth.spline(seq_along(signal), signal)$y
    } else identity
    
    # get counts
    strand(region) <- "*"
    getCountsByPositions(grl, region, binsize = binsize, FUN = bin_FUN,
                         expand_ranges = expand_ranges, ncores = ncores) %>%
        lapply(as.vector) %>%
        lapply(smooth_fun) %>%
        Map(function(i, nm) 
            data.frame(x = (binsize*seq_along(i)) - binsize, 
                       signal = i, 
                       sample = nm),
            ., names(.)) %>% 
        c(make.row.names = FALSE) %>%
        do.call(rbind, .) %>% 
        mutate(sample = factor(sample, levels = unique(sample))) -> 
        df
    
    # get y-axis scales
    if (is.null(ylim))
        ylim <- c(0, max(pretty(df$signal)))
    
    # saturate off-scaled data (for user supplied ylims)
    #   -> doing it this way allows the axis line to be capped
    df$signal[df$signal > max(ylim)] <- max(ylim)
    df$signal[df$signal < min(ylim)] <- min(ylim)
    
    plotfun_unstr <- function(df) {
        ggplot(df) + 
            facet_grid(rows = vars(sample)) +
            ggplot2::geom_segment(aes(x = x, xend = x, y = 0, yend = signal),
                                  color = "black") + 
            scale_y_continuous(limits = range(ylim), breaks = ylim, 
                               labels = ylim, expand = c(0, 0)) +
            scale_x_continuous(limits = c(0, width(region)),
                               expand = expansion(mult = c(0.02, 0))) +
            labs(x = NULL, y = NULL) +
            theme_classic() +
            theme(panel.grid = element_blank(),
                  panel.border = element_blank(),
                  axis.text.y = element_text(color = "black", size = 6),
                  axis.text.x = element_blank(),
                  axis.ticks.y = element_line(color = "black", size = 0.25),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank(),
                  axis.line.y = element_line(color = "black", size = 0.25),
                  strip.text.y = element_text(angle = 0, hjust = 0),
                  strip.background = element_blank(),
                  plot.margin = unit(c(0,0,3,0), "points"))
    }
    
    plots <- lapply(split(df, df$sample), plotfun_unstr)
    
    if (!plot)
        return(plots)
    do.call(wrap_plots, c(list(plots), ncol = 1))
}


shot_genearrow <- function(region, genelist, gene_names, plot = TRUE) {
    ## Plot simple arrows showing the positions and directions of ranges 
    ## (e.g. genes)
    
    # attach gene_names, and subset genelist
    genelist$gene_names <- gene_names
    gl <- sort(subsetByOverlaps(genelist, region, ignore.strand = TRUE))
    
    # make dataframe for gene models
    df <- data.frame(start = start(gl) - start(region),
                     end = end(gl) - start(region),
                     tx = gl$gene_names)
    
    arrow_ends <- ifelse(strand(gl) == "+", "last", "first")
    
    # due to using scale_x_cont., have to trim out-of-range transcripts;
    # -> but don't want them to have arrows
    # -> will also include unstranded ranges here
    idx.long <- unique(c(
        which(df$start < 0L & arrow_ends == "first"),
        which(df$end > width(region) & arrow_ends == "last"),
        which(strand(gl) == "*")
    ))
    df$start[df$start < 0L] <- 0L
    df$end[df$end > width(region)] <- width(region)
    df.long <- df[idx.long, ]
    if (length(idx.long) > 0L) {
        df <- df[-idx.long, , drop = FALSE]
        arrow_ends <- arrow_ends[-idx.long]
    }
    
    # make and store ggplot
    ggplot() + 
        ggplot2::geom_segment(
            aes(x = start, xend = end, y = tx, yend = tx),
            data = df,
            arrow = if (length(arrow_ends) > 0L) {
                arrow(ends = arrow_ends, type = "closed", 
                      length = unit(4, "points"))
            } else NULL
        ) + 
        ggplot2::geom_segment(
            aes(x = start, xend = end, y = tx, yend = tx),
            data = df.long
        ) + 
        scale_x_continuous(limits = c(0L, width(region)),
                           expand = expansion(mult = c(0.02, 0))) +
        labs(x = NULL, y = NULL) +
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.text.y = element_text(size = 5),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = unit(c(0,0,0,0), "points")) -> p
    
    if (plot) p else list(p)
}


shot_genemodel <- function(region, txdb, plot = TRUE) {
    ## Plot full gene models from a TxDb object
    
    # get transcripts in regions; also use for mapping "TXID" to "TXNAME"
    transcripts(txdb) %>% 
        findOverlapPairs(region, ignore.strand = TRUE) %>%
        pintersect %>% 
        subset(hit) %>%
        (function(gr) data.frame(start = start(gr) - start(region), 
                                 end = end(gr) - start(region),
                                 tx_id = gr$tx_id, tx = gr$tx_name)) -> txs
    # functions for tx features
    # -> Intersect GRangesList ranges over plot region, and return GRanges
    break_grl <- function(grl) {
        findOverlapPairs(grl, region, ignore.strand = TRUE) %>% 
            pintersect %>% 
            endoapply(subset, hit)
    }
    # -> get tx names from GRangesList and add to GRanges mcols
    add_names_grl <- function(grl) {
        txnames <- sapply(names(grl), function(i) txs$tx[txs$tx_id == i])
        mendoapply(function(x, nm) { x$tx_id <- nm; x }, 
                   grl, txnames)
    }
    # -> from GRanges or GRangesList, get dataframe for plotting
    make_df <- function(gr) data.frame(start = start(gr) - start(region),
                                       end = end(gr) - start(region),
                                       tx = as.character(gr$tx_id))
    grl2df <- function(x) make_df(unlist(add_names_grl(x)))
    
    # store UTRs for trimming exons
    fiveUTRsByTranscript(txdb) %>% 
        break_grl %>%
        endoapply(function(x) {
            # oddly there are transcripts with multiple 5' UTRs
            if (length(x) < 1L) 
                return(x)
            if (levels(droplevels(strand(x))) == "+") 
                return(x[which.min(start(x))])
            x[which.max(end(x))]
        }) -> fiveUTR_grl
    threeUTRsByTranscript(txdb) %>% break_grl -> threeUTR_grl
    
    # plot data for UTRs
    fiveUTR <- grl2df(fiveUTR_grl)
    threeUTR <- grl2df(threeUTR_grl)
    
    # plot data for exons - must be trimmed to not overlap UTRs
    trim_overlap <- function(grl1, grl2) {
        nm <- names(grl1)[names(grl1) %in% names(grl2)]
        if (length(nm) > 0L)
            grl1[nm] <- mendoapply(GenomicRanges::setdiff, 
                                   grl1[nm], grl2[nm])
        grl1[lengths(grl1) > 0L]
    }
    exonsBy(txdb, by = "tx") %>% 
        break_grl %>%
        trim_overlap(fiveUTR_grl) %>% 
        trim_overlap(threeUTR_grl) %>%
        grl2df -> exons
    
    # plot data for directional arrows over introns; first arrow spacing
    arrow_spacer <- floor(width(region) / 20) # distance b/w arrows
    
    # get introns for arrow positions and directions
    intronsByTranscript(txdb) %>% break_grl %>%
        endoapply(subset, width > arrow_spacer) %>%
        endoapply(function(x) unlist(tile(x, width = arrow_spacer))) %>%
        resize(2, fix = "center") %>% .[lengths(.) > 0] %>% 
        add_names_grl %>% unlist -> arrows_gr
    
    # arrow positions & directions, and add a slight aesthetic shift
    intron_arrows <- make_df(arrows_gr)
    arrow_ends <- ifelse(strand(arrows_gr) == "-", "first", "last")
    if (length(arrows_gr) > 0L)  
        arrows_gr <- shift(arrows_gr, ifelse(strand(arrows_gr) == "-", 
                                             -0.25 * arrow_spacer, 
                                             0.25 * arrow_spacer))
    
    ggplot() +
        ggplot2::geom_segment( # blank spacers to help set sizes
            mapping = if (nrow(txs)==0L) NULL else
                aes(x = 0, xend = width(region), y = tx, yend = tx),
            data = txs, color = "white", size = 10
        ) +
        ggplot2::geom_segment( # blank spacers
            mapping = if (nrow(exons)==0L) NULL else
                aes(x = start, xend = end, y = tx, yend = tx),
            data = exons, color = "white", size = 5
        ) +
        ggplot2::geom_segment( # plot transcript (intron) lines
            mapping = if (nrow(txs)==0L) NULL else
                aes(x = start, xend = end + 1, y = tx, yend = tx),
            data = txs, size = 0.5
        ) +
        ggplot2::geom_segment( # arrows in introns
            mapping = if (nrow(intron_arrows)==0L) NULL else
                aes(x = start, xend = end, y = tx, yend = tx),
            data = intron_arrows, size = 0.25,
            arrow = if (nrow(intron_arrows)==0L) NULL else
                arrow(ends = arrow_ends, length = unit(4, "points"))
        ) +
        ggplot2::geom_segment( # exon blocks
            mapping = if (nrow(exons)==0L) NULL else
                aes(x = start, xend = end + 1, y = tx, yend = tx),
            data = exons, size = 2.5
        ) +
        ggplot2::geom_segment( # 5' UTR blocks
            mapping = if (nrow(fiveUTR)==0L) NULL else
                aes(x = start, xend = end + 1, y = tx, yend = tx),
            data = fiveUTR, size = 1.25
        ) +
        ggplot2::geom_segment( # 3' UTR blocks
            mapping = if (nrow(threeUTR)==0L) NULL else
                aes(x = start, xend = end + 1, y = tx, yend = tx),
            data = threeUTR, size = 1.25
        ) +
        scale_x_continuous(limits = c(0, width(region)),
                           expand = expansion(mult = c(0.02, 0))) +
        labs(x = NULL, y = NULL) +
        theme_classic() +
        theme(axis.line = element_blank(),
              axis.text.y = element_text(size = 5),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              legend.position = "none",
              plot.margin = unit(c(0,0,0,0), "points")) -> p
    
    if (plot) p else list(p)
}