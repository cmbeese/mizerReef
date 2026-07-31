#' Upgrade a mizerReef params object to the current layout
#'
#' Called automatically by [mizer::validParams()] when an object created with
#' an older version of mizerReef is loaded. Migrates data stored in the old
#' custom S4 slots (`refuge_params`, `algae_params`, `detritus_params`) to the
#' `other_params` sub-lists used since version 2.1.0.
#'
#' @param object A `mizerReef` params object (possibly from an older version).
#' @param ... Unused.
#' @return An upgraded `mizerReef` object.
#' @importFrom methods .hasSlot slot
#' @exportS3Method utils::upgrade
upgrade.mizerReef <- function(object, ...) {
    # Migrate old S4 slot data to other_params sub-lists if needed.
    # The old class had three extra slots: refuge_params, algae_params,
    # and detritus_params. In mizerReef >= 2.1.0 these are stored as named
    # lists inside other_params instead.
    #
    # Because the slots no longer exist in the class definition, attempting to
    # access them via @ would raise an error. We use tryCatch to detect the
    # old layout structurally and migrate it safely.
    tryCatch({
        old_refuge <- .hasSlot(object, "refuge_params")
        old_algae  <- .hasSlot(object, "algae_params")
        old_det    <- .hasSlot(object, "detritus_params")
        if (old_refuge && is.null(object@other_params$refuge_params)) {
            object@other_params$refuge_params  <- slot(object, "refuge_params")
        }
        if (old_algae && is.null(object@other_params$algae_params)) {
            object@other_params$algae_params   <- slot(object, "algae_params")
        }
        if (old_det && is.null(object@other_params$detritus_params)) {
            object@other_params$detritus_params <- slot(object, "detritus_params")
        }
    }, error = function(e) NULL)  # silently ignore if slots do not exist
    object
}

#' Upgrade a mizerReef 1.x params object to the current (2.1.0+) layout
#'
#' `mizerReef` 2.x replaced the direct rate-function overrides and flat
#' `other_params` layout used by mizerReef 1.x (the version accompanying the
#' original thesis, on this package's `main` branch) with S3 dispatch on a
#' `mizerReef` marker class and a nested `other_params` layout. Objects
#' created with mizerReef 1.x are never automatically upgraded by
#' [mizer::validParams()], because they were never registered as a
#' `"mizerReef"` extension in the first place -- there is nothing for the
#' automatic [upgrade.mizerReef()] hook to trigger on. Call this function
#' once, by hand, on any params object created with mizerReef 1.x to bring
#' it up to the current structural layout.
#'
#' This is a **structural** migration only. Tuned values that came from your
#' original model -- `algae_growth`, `detritus_external`, capacities, the
#' refuge method and its parameters, species parameters -- are carried over
#' as they were, not reset to package defaults and not silently re-tuned.
#' In particular, the refuge profile (`refuge`/`bin.id`/`refuge_lengths`) is
#' recomputed from your preserved method/method_params/`a_bar`/`b_bar`
#' using the current [getRefuge()] logic, since that is a deterministic
#' function of already-preserved inputs, not a recalibration.
#'
#' Because your original algae/detritus tuning is carried over unchanged, an
#' upgraded object is not guaranteed to be at a genuine joint steady state
#' under the current algae/detritus mechanism (see `vignette("extension_mechanism")`
#' and the algae/detritus sections of [setAlgaeParams()]/[setDetritusParams()]
#' for what changed). Run [reefSteady()] on the result if you want a fresh
#' steady state under the current mechanism.
#'
#' @param params A `MizerParams` object created with mizerReef 1.x (i.e.
#'   from code on this package's `main`/thesis branch, not a `mizerReef`
#'   object already created with this version of the package).
#' @return A `mizerReef` object with the current structural layout.
#' @seealso [upgrade.mizerReef()] for the automatic, incremental micro-version
#'   upgrades that run every time [mizer::validParams()] is called on an
#'   already-registered `mizerReef` object -- this function is for the
#'   one-off 1.x -> 2.x jump only.
#' @export
upgradeReefParams <- function(params) {
    assert_that(is(params, "MizerParams"))

    if (is(params, "mizerReef")) {
        message(
            "This params object is already a mizerReef object created ",
            "with a current version of the package. Nothing to upgrade -- ",
            "future micro-version upgrades happen automatically via ",
            "validParams()."
        )
        return(mizer::validParams(params))
    }

    # Let mizer's own upgrade machinery handle any core mizer-level version
    # drift first (e.g. an old `mizer_version` stamp, a missing `Diffusion`
    # rates_funcs entry). Verified this leaves every reef-specific field
    # untouched -- mizer's upgrade code has no knowledge of mizerReef.
    params <- suppressWarnings(mizer::validParams(params))

    op <- params@other_params

    # Structural signature of a genuine mizerReef 1.x object: the live
    # algae/detritus dynamics components exist (`other_params$algae`/
    # `$detritus`, set up by mizer::setComponent() in newReefParams() since
    # 1.x), but none of the 2.1.0+ nesting/echo structure does.
    is_legacy <- !is.null(op$algae) && is.null(op$algae_params) &&
        !is.null(op$refuge_params) && is.null(op$refuge_params$method_params)

    if (!is_legacy) {
        stop(
            "This does not look like a mizerReef 1.x params object (the ",
            "structure this function expects to migrate). If it was ",
            "created with a recent version of mizerReef, you don't need ",
            "to upgrade it -- validParams() keeps it current automatically."
        )
    }

    # --- species_params: bad_pred -> blocked_pred ---------------------------
    sp <- params@species_params
    if ("bad_pred" %in% names(sp)) {
        if ("blocked_pred" %in% names(sp)) {
            sp$bad_pred <- NULL
        } else {
            names(sp)[names(sp) == "bad_pred"] <- "blocked_pred"
        }
    }
    params@species_params <- sp

    # --- rates_funcs: drop the old direct overrides --------------------------
    # v2 gets reef behaviour via S3 dispatch (project*.mizerReef methods)
    # instead of these direct name overrides -- see vignette("extension_mechanism").
    stock_rates <- list(
        Rates = "mizerRates", Encounter = "mizerEncounter",
        FeedingLevel = "mizerFeedingLevel", PredMort = "mizerPredMort",
        Mort = "mizerMort"
    )
    reset_rates <- character(0)
    for (nm in names(stock_rates)) {
        if (!identical(params@rates_funcs[[nm]], stock_rates[[nm]])) {
            reset_rates <- c(reset_rates, nm)
            params@rates_funcs[[nm]] <- stock_rates[[nm]]
        }
    }

    # --- refuge: rebuild the nested refuge_params structure -------------------
    old_refuge <- op$refuge_params  # 1x6 data.frame in 1.x
    # setRefuge() builds this up as a list (params@other_params$refuge_params$x
    # <- ...); the 1.x object instead has a 1-row data.frame there, which
    # breaks list-shaped assignment (e.g. a multi-row method_params can't be
    # assigned into a single data.frame cell). Clear it to a list first.
    params@other_params$refuge_params <- list()
    params <- setRefuge(
        params,
        method = old_refuge$method,
        method_params = op$method_params,
        a_bar = old_refuge$a_bar, b_bar = old_refuge$b_bar,
        w_settle = old_refuge$w_settle, max_protect = old_refuge$max_protect,
        tau = old_refuge$tau,
        use_dummy_fish_bins = TRUE
    )
    params <- getRefuge(params, use_dummy_fish_bins = TRUE)
    params@other_params$refuge_params$degrade <- isTRUE(op$degrade)
    if (isTRUE(op$degrade)) {
        params@other_params$refuge_params$t_bleach   <- op$t_bleach
        params@other_params$refuge_params$trajectory <- op$trajectory
        params@other_params$refuge_params$deg_scale  <- op$deg_scale
        warning(
            "This object had habitat degradation enabled (`degrade = TRUE`). ",
            "Its degradation trajectory settings (t_bleach/trajectory/",
            "deg_scale) have been relocated into other_params$refuge_params, ",
            "but this path is not exercised by either bundled example model ",
            "and deserves a manual check -- verify with plotDegScale()/",
            "plotDegradationScale() that the migrated trajectory still ",
            "looks right before relying on it.",
            call. = FALSE
        )
    }

    # --- algae: preserve tuned values, populate both the live component and
    #     the (structurally required, but not authoritative) config echo ----
    algae_capacity <- if (is.null(op$algae_capacity)) 1 else op$algae_capacity
    use_UR_cc <- isTRUE(op$carry_capacity)
    params <- setAlgaeParams(
        params,
        algae_growth_initial = op$algae$growth,
        algae_capacity = algae_capacity,
        use_UR_cc = use_UR_cc
    )
    params@other_params$algae <- list(
        rho = op$algae$rho,
        capacity = algae_capacity,
        growth = op$algae$growth
    )
    # setAlgaeParams() doesn't populate the echo's `rho` (only
    # newReefParams()'s own rho-calculation step does that on a fresh
    # build) -- add it here so the echo matches the live component.
    params@other_params$algae_params$rho <- op$algae$rho

    # --- detritus: same treatment ---------------------------------------------
    detritus_capacity <- if (is.null(op$detritus_capacity)) 1 else op$detritus_capacity
    params <- setDetritusParams(
        params,
        detritus_capacity = detritus_capacity,
        sen_decomp = op$detritus$sen_decomp,
        ext_decomp = op$detritus$ext_decomp,
        external = op$detritus$external,
        use_UR_cc = use_UR_cc
    )
    params@other_params$detritus <- list(
        rho = op$detritus$rho,
        capacity = detritus_capacity,
        sen_decomp = op$detritus$sen_decomp,
        ext_decomp = op$detritus$ext_decomp,
        external = op$detritus$external
    )
    params@other_params$detritus_params$rho <- op$detritus$rho

    # --- external mortality / senescence / new_refuge: unchanged names --------
    params@other_params$ext_mort_params  <- op$ext_mort_params
    params@other_params$include_sen_mort <- isTRUE(op$include_sen_mort)
    params@other_params$new_refuge       <- isTRUE(op$new_refuge)

    # --- drop the old flat top-level fields now superseded by the nested/
    #     renamed locations above (would otherwise sit alongside the new
    #     structure as dead, increasingly stale duplicates) ------------------
    stale_top_level <- c(
        "method_params", "refuge", "bin.id", "refuge_lengths",
        "carry_capacity", "algae_capacity", "detritus_capacity",
        "sen_decomp", "ext_decomp", "degrade", "t_bleach", "trajectory",
        "deg_scale",
        # stale bookkeeping copies of the pre-tuning input, not the live
        # tuned value
        "initial_algae_growth", "initial_d_external"
    )
    for (nm in stale_top_level) {
        params@other_params[[nm]] <- NULL
    }

    # --- register the extension chain and promote to the mizerReef class ------
    params@extensions <- mizer::getRegisteredExtensions()
    params <- mizer::coerceToExtensionClass(params)
    params <- mizer::validParams(params)

    if (length(reset_rates) > 0) {
        warning(
            "The old direct rate-function overrides (",
            paste(unlist(stock_rates[reset_rates]), collapse = ", "),
            ") have been reset to mizer's stock rate functions. Reef ",
            "behaviour is now applied through mizerReef's S3 dispatch ",
            "methods instead (see vignette(\"extension_mechanism\")) -- this ",
            "is a real behavioural change (it picks up every fix made to ",
            "those methods since mizerReef 1.x), not just a rename.",
            call. = FALSE
        )
    }
    warning(
        "This object was migrated from mizerReef 1.x. Its algae/detritus ",
        "tuning (algae_growth, detritus external flux, capacities) was ",
        "carried over as-is from the old model and was NOT re-tuned -- it ",
        "is not guaranteed to be at a genuine joint steady state under the ",
        "current algae/detritus mechanism. Run reefSteady() on the result ",
        "if you want a fresh steady state under the current mechanism.",
        call. = FALSE
    )

    params
}
