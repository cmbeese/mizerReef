# Overwrite mizer's emptyParams() function to include a slot for vulnerability
# in the rates function

#' Set up empty parameters for a mizerReef model
#'
#' @return An object of type \linkS4class{MizerParams}
#'
#' @export
#' @family functions for setting up models

emptyReefParams <- function(species_params,
                            gear_params = data.frame(),
                            no_w = 100,
                            min_w = 0.001,
                            # w_full = NA,
                            max_w = NA,
                            min_w_pp = 1e-12) {
    assert_that(is.data.frame(species_params),
                is.data.frame(gear_params),
                no_w > 10)

    ## Set defaults ----
    if (is.na(min_w_pp)) min_w_pp <- 1e-12
    species_params <- set_species_param_default(species_params, "w_min", min_w)
    min_w <- min(species_params$w_min)

    species_params <- validSpeciesParams(species_params)
    gear_params <- validGearParams(gear_params, species_params)

    if (is.na(max_w)) {
        max_w <- max(species_params$w_max)
    } else {
        if (max(species_params$w_max) > max_w * (1 + 1e-6)) { # The fudge factor
            # is there to avoid false alerts due to rounding errors.
            too_large <- species_params$species[max_w < species_params$w_max]
            stop("Some of your species have an maximum size larger than max_w: ",
                 toString(too_large))
        }
    }

    # Set up grids ----
    # The following code anticipates that in future we might allow the user to
    # specify a grid with a non-constant log spacing. But we comment this out
    # for now because of the fft.
    # if (missing(w_full)) {
    # set up logarithmic grids
    dx <- log10(max_w / min_w) / (no_w - 1)
    # Community grid
    w <- 10^(seq(from = log10(min_w), by = dx, length.out = no_w))
    # dw[i] = w[i+1] - w[i]. Following formula works also for last entry dw[no_w]
    dw <- (10^dx - 1) * w
    # To avoid issues due to numerical imprecision
    min_w <- w[1]

    # For fft methods we need a constant log bin size throughout.
    # Therefore we use as many steps as are necessary so that the first size
    # class includes min_w_pp.
    if (min_w_pp >= min_w) {
        stop("min_w_pp must be larger than min_w")
    }
    x_pp <- rev(seq(from = log10(min_w),
                    to = log10(min_w_pp),
                    by = -dx)) - dx
    w_full <- c(10^x_pp, w)
    # If min_w_pp happened to lie exactly on a grid point, we now added
    # one grid point too much which we need to remove again
    if (w_full[2] == min_w_pp) {
        w_full <- w_full[2:length(w_full)]
    }
    no_w_full <- length(w_full)
    dw_full <- (10^dx - 1) * w_full
    # } else {
    #     # use supplied w_full
    #     no_w_full <- length(w_full) - 1
    #     dw_full <- diff(w_full)
    #     w_full <- w_full[seq_along(dw_full)]
    #     # Check that sizes are increasing
    #     if (any(dw_full <= 0)) {
    #         stop("w_full must be increasing.")
    #     }
    #     w_min_idx <- match(min_w, w_full)
    #     if (is.na(w_min_idx)) {
    #         stop("w_min must be contained in w_full.")
    #     }
    #     w <- w_full[w_min_idx:no_w_full]
    #     dw <- dw_full[w_min_idx:no_w_full]
    #     no_w <- length(w)
    #     min_w_pp <- w_full[1]
    # }

    # Basic arrays for templates ----
    no_sp <- nrow(species_params)
    species_names <- as.character(species_params$species)
    gear_names <- unique(gear_params$gear)
    mat1 <- array(0, dim = c(no_sp, no_w),
                  dimnames = list(sp = species_names, w = signif(w, 3)))
    ft_pred_kernel <- array(NA, dim = c(no_sp, no_w_full),
                            dimnames = list(sp = species_names, k = 1:no_w_full))
    ft_mask <- t(sapply(species_params$w_max, function(x) w_full < x))

    selectivity <- array(0, dim = c(length(gear_names), no_sp, no_w),
                         dimnames = list(gear = gear_names, sp = species_names,
                                         w = signif(w, 3)))
    catchability <- array(0, dim = c(length(gear_names), no_sp),
                          dimnames = list(gear = gear_names, sp = species_names))
    default_effort <- ifelse(defaults_edition() < 2, 0, 1)
    initial_effort <- rep(default_effort, length(gear_names))
    names(initial_effort) <- gear_names

    interaction <- array(1, dim = c(no_sp, no_sp),
                         dimnames = list(predator = species_names,
                                         prey = species_names))

    vec1 <- as.numeric(rep(NA, no_w_full))
    names(vec1) <- signif(w_full, 3)

    w_min_idx <- get_w_min_idx(species_params, w)

    # Colour and linetype scales ----
    # for use in plots
    # Colour-blind-friendly palettes
    # From http://dr-k-lo.blogspot.co.uk/2013/07/a-color-blind-friendly-palette-for-r.html
    # cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00",
    #                 "#CC79A7", "#F0E442")
    # From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    # cbbPalette <- c("#E69F00", "#56B4E9", "#009E73",
    #                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # Random palette gemerated pm https://medialab.github.io/iwanthue/
    colour_palette <- c("#815f00",
                        "#6237e2",
                        "#8da600",
                        "#de53ff",
                        "#0e4300",
                        "#430079",
                        "#6caa72",
                        "#ee0053",
                        "#007957",
                        "#b42979",
                        "#142300",
                        "#a08dfb",
                        "#644500",
                        "#04004c",
                        "#b79955",
                        "#0060a8",
                        "#dc8852",
                        "#007ca9",
                        "#ab003c",
                        "#9796d9",
                        "#472c00",
                        "#b492b0",
                        "#140000",
                        "#dc8488",
                        "#005c67",
                        "#5c585a")
    # type_palette <- c("solid", "dashed", "dotdash", "longdash",
    #                   "twodash")
    type_palette <- c("solid")

    if ("linecolour" %in% names(species_params)) {
        linecolour <- species_params$linecolour
        # If any NA's first fill them with unused colours
        linecolour[is.na(linecolour)] <-
            setdiff(colour_palette, linecolour)[1:sum(is.na(linecolour))]
        # if there are still NAs, start from beginning of palette again
        linecolour[is.na(linecolour)] <-
            colour_palette[1:sum(is.na(linecolour))]
    } else {
        linecolour <- rep(colour_palette, length.out = no_sp)
    }
    names(linecolour) <- as.character(species_names)
    linecolour <- c(linecolour, "Resource" = "green", "Total" = "black",
                    "Background" = "grey", "Fishing" = "red",
                    "External" = "grey")

    if ("linetype" %in% names(species_params)) {
        linetype <- species_params$linetype
        linetype[is.na(linetype)] <- "solid"
    } else {
        linetype <- rep(type_palette, length.out = no_sp)
    }
    names(linetype) <- as.character(species_names)
    linetype <- c(linetype, "Resource" = "solid", "Total" = "solid",
                  "Background" = "solid", "Fishing" = "solid",
                  "External" = "solid")

    # Make object ----
    # Should Z0, rrPP and ccPP have names (species names etc)?
    params <- new(
        "MizerParams",
        metadata = list(),
        mizer_version = packageVersion("mizer"),
        extensions = vector(mode = "character"),
        time_created = lubridate::now(),
        time_modified = lubridate::now(),
        w = w,
        dw = dw,
        w_full = w_full,
        dw_full = dw_full,
        w_min_idx = w_min_idx,
        maturity = mat1,
        psi = mat1,
        initial_n = mat1,
        intake_max = mat1,
        search_vol = mat1,
        metab = mat1,
        mu_b = mat1,
        ft_pred_kernel_e = ft_pred_kernel,
        ft_pred_kernel_p = ft_pred_kernel,
        pred_kernel = array(),
        gear_params = gear_params,
        selectivity = selectivity,
        catchability = catchability,
        initial_effort = initial_effort,
        rr_pp = vec1,
        cc_pp = vec1,
        sc = w,
        initial_n_pp = vec1,
        species_params = species_params,
        interaction = interaction,
        other_dynamics = list(),
        other_encounter = list(),
        other_mort = list(),
        rates_funcs = list(
            Rates = "mizerRates",
            Vulnerable = "reefVulnerable",
            Encounter = "mizerEncounter",
            FeedingLevel = "mizerFeedingLevel",
            EReproAndGrowth = "mizerEReproAndGrowth",
            PredRate = "mizerPredRate",
            PredMort = "mizerPredMort",
            FMort = "mizerFMort",
            Mort = "mizerMort",
            ERepro = "mizerERepro",
            EGrowth = "mizerEGrowth",
            ResourceMort = "mizerResourceMort",
            RDI = "mizerRDI",
            RDD = "BevertonHoltRDD"
        ),
        resource_dynamics = "resource_semichemostat",
        other_params = list(),
        initial_n_other = list(),
        A = as.numeric(rep(NA, no_sp)),
        linecolour = linecolour,
        linetype = linetype,
        ft_mask = ft_mask
    )

    return(params)
}

# environment(emptyReefParams) <- asNamespace("mizer")
# utils::assignInNamespace("emptyParams", emptyReefParams, ns = "mizer")
