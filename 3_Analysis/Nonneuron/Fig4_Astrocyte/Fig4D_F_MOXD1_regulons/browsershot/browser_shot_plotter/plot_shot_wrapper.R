

plot_shot <- function(..., region, binsize = NULL, nbins = 250L, 
                      ylim = NULL, expand_ranges = TRUE, bin_FUN = sum, 
                      annotations = NULL, pad_left = 0, pad_right = 0, 
                      model_height = 1, gene_names = NULL, smooth = FALSE,
                      ncores = getOption("mc.cores", 2L)) {
    
    ## for each argument in the functions plot_stranded/plot_unstranded, this
    ## function will make a list with the same name as that argument;
    ##   each element of that list will have the name of one of the dataset
    ##   groups
    
    ## i.e. for the argument ylim, if the dataset groups are PROcap, PROseq, and 
    ## GC, an object named "ylim" will be created, with the format:
    ##   ylim = list(PROcap = c(0, 10), PROseq = c(0, 5), GC = c(0, 100))
    
    # -------------------------------------------------- #
    # Format dataset list, and set all names
    # -------------------------------------------------- #
    
    # each element of grl is a dataset group of one or more GRanges that will 
    # be group-scaled 
    grl <- list(...)
    
    # get any names provided in function call
    call_args <- as.list(match.call())[-1][seq_along(grl)]
    arg_names <- names(call_args)
    
    # for any unnamed arguments, get the names of the objects
    noname <- nchar(arg_names) < 1L
    arg_names[noname] <- as.character(call_args)[noname]
    
    # for any non-list elements of grl, first put them into a list and name the
    # object within
    if (!all(is_list <- vapply(grl, is.list, logical(1))))
        grl[!is_list] <- Map(function(x, nm) setNames(list(x), nm),
                             grl[!is_list], arg_names[!is_list])
    
    # finally set any list names that weren't given
    names(grl)[noname] <- arg_names[noname]
    
    # check for unnamed datasets (GRanges)
    if (any(vapply(grl, function(x) any(is.null(names(x))), logical(1))))
        warning("One or more datasets are not named, most likely due to a list
                with unnamed elements")

    # -------------------------------------------------- #
    # Get region, annotation plot, and scalebar plot
    # -------------------------------------------------- #
    
    # apply region padding
    ranges(region) <- IRanges(start(region) - pad_left, 
                              end(region) + pad_right)
    
    # get annotation plot
    if (is(annotations, "TxDb")) {
        pmod <- shot_genemodel(region = region, txdb = annotations, 
                               plot = FALSE)
    } else {
        pmod <- shot_genearrow(region = region, genelist = annotations, 
                               gene_names = gene_names, plot = FALSE)
    }
    
    # get scalebar
    pbar <- shot_scalebar(region, plot = FALSE)
    
    # -------------------------------------------------- #
    # Determine data plotting functions -> plot_FUN
    # -------------------------------------------------- #
    
    # function to test strandedness for individual GRanges objects
    test_gr_stranded <- function(x) {
        strands <- levels(droplevels(strand(x)))
        if (!"*" %in% strands)
            return(TRUE)
        if (length(strands) > 1)
            stop("Cannot mix stranded and unstranded ranges 
                 within the same GRanges object")
        return(FALSE)
    }
    
    # list of plotting functions
    plot_FUN <- lapply(grl, function(x) {
        stranded <- unique(vapply(x, test_gr_stranded, logical(1)))
        if (length(stranded) > 1)
            stop("Individual lists of GRanges data cannot 
                 mix stranded and unstranded data")
        ifelse(stranded, shot_stranded, shot_unstranded) 
    })
    
    # -------------------------------------------------- #
    # Parse arguments for data plots: variable args
    # -------------------------------------------------- #
    
    # make lists of all args for shot_stranded/shot_unstranded, 
    #   (each list with names = names(grl))
    
    # for args that can vary by dataset group
    arg_merge <- function(arg) {
        argname <- deparse(substitute(arg))
        
        if (!is.list(arg)) {
            # if single value given (or default kept), use for all
            value <- rep(list(arg), length(grl))
            names(value) <- names(grl)
            
        } else {
            # list of values can vary by dataset group
            if (is.null(names(arg))) {
                if (length(arg) != length(grl)) 
                    stop(sprintf("%s not named & not the same length as inputs", 
                                 argname))
                names(arg) <- names(grl)
            }
            # initialize with default values
            default <- formals(plot_shot)[[argname]]
            value <- rep(list(default), length(grl))
            names(value) <- names(grl)
            
            # fill in any values different from the defaults
            value[names(arg)] <- arg
        }
        value
    }
    
    ylim <- arg_merge(ylim)
    bin_FUN <- arg_merge(bin_FUN)
    expand_ranges <- arg_merge(expand_ranges)
    smooth <- arg_merge(smooth)
    
    # -------------------------------------------------- #
    # Parse arguments for data plots: global args
    # -------------------------------------------------- #
    
    # for args that cannot vary across dataset groups
    check_args <- vapply(
        list(region, binsize, nbins, pad_left, pad_right, model_height),
        function(x) is.list(x) || length(x) > 1, 
        logical(1)
    )
    if (any(check_args))
        stop("values for region, binsize, nbins, pad_left, pad_right, 
             and model_height must have a length of 1 (i.e. they cannot 
             vary across plotting groups)")
    
    if (is.null(binsize))
        binsize <- as.integer(width(region) / nbins)
    
    # get lists of args to pass to shot_stranded/shot_unstranded
    region <- rep(list(region), length(grl))
    binsize <- rep(list(binsize), length(grl))
    plot <- rep(list(FALSE), length(grl))
    ncores <- rep(list(ncores), length(grl))
    names(ncores) <- names(plot) <- names(binsize) <- names(region) <- names(grl)
    
    # -------------------------------------------------- #
    # Get data plots, and assemble final plot
    # -------------------------------------------------- #
    
    pdat <- lapply(names(grl), function(x) {
        FUN <- plot_FUN[[x]]
        argnames <- names(formals(FUN))
        args <- lapply(argnames, function(nm) get(nm)[[x]])
        names(args) <- argnames
        do.call(FUN, args)
    })
    
    list(pmod, pbar) %>%
        c(pdat) %>%
        unlist(recursive = FALSE) %>%
        unname -> all_plots
    
    plot_heights <- c(model_height, 0.1, rep(1, length(all_plots) - 2))
    
    wrap_plots(all_plots, ncol = 1, heights = plot_heights)
}