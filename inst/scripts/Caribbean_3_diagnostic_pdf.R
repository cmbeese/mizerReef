# Simple diagnostic PDF for the bundled caribbean_3_model: biomass match,
# diets, growth curves, spectra. Visual grammar (target-vs-model points,
# one diagnostic per page) borrowed from
# caribbean-10-calibration/generate_diagnostic_pdfs.R, trimmed down --
# caribbean_3 has no soft-target scoreboard or flat-vs-complex comparison.
# Run from the repo root: `Rscript inst/scripts/Caribbean_3_diagnostic_pdf.R`.
suppressMessages({
  library(mizer)
  library(mizerReef)
  library(ggplot2)
})

e <- new.env()
data("caribbean_3_model", envir = e)
params <- e$caribbean_3_model

out_dir <- "artifacts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_file <- file.path(out_dir, paste0("caribbean_3_model-diagnostic-", timestamp, ".pdf"))

target3 <- c(predators = 107, herbivores = 34, inverts = 40)

plot_biomass_match <- function(params, target3) {
  b <- suppressWarnings(getBiomass(params) |> tapply(params@species_params$species, sum))
  df <- data.frame(
    species = names(target3),
    preferred = as.numeric(target3),
    model = as.numeric(b[names(target3)])
  )
  ggplot(df, aes(x = reorder(species, preferred))) +
    geom_point(aes(y = preferred, shape = "FORCE-observed target", color = species),
               fill = "white", size = 4, stroke = 1.2) +
    geom_point(aes(y = model, shape = "Model estimate", color = species), size = 4) +
    scale_shape_manual(name = NULL, values = c("FORCE-observed target" = 23, "Model estimate" = 16),
                        guide = guide_legend(override.aes = list(colour = "grey30", fill = "grey30"))) +
    scale_color_manual(name = "Species", values = params@linecolour,
                        guide = guide_legend(override.aes = list(shape = 16, size = 3))) +
    theme_bw() +
    labs(title = "Biomass match: FORCE-observed targets", x = "Species", y = "Biomass [g/m^2]") +
    theme(legend.position = "right")
}

pdf(out_file, width = 11, height = 8.5, onefile = TRUE)

p_match <- tryCatch(plot_biomass_match(params, target3), error = function(e) NULL)
if (!is.null(p_match)) print(p_match)

p_diet <- tryCatch(plotDiet(params, wlim = c(1, NA)), error = function(e) NULL)
if (!is.null(p_diet)) print(p_diet)

p_growth <- tryCatch({
  gdat <- plotGrowthCurves(params, return_data = TRUE)
  gdat$Size_g <- gdat[["Size [g]"]]
  ggplot(gdat, aes(x = Age, y = Size_g, color = Legend, linetype = Legend)) +
    geom_line(linewidth = 0.4) +
    scale_y_log10() +
    facet_wrap(~Species, scales = "free_y") +
    theme_bw() +
    labs(y = "Size [g]", title = "Growth curves by species (model vs von Bertalanffy)") +
    theme(legend.position = "bottom")
}, error = function(e) NULL)
if (!is.null(p_growth)) print(p_growth)

p_spectra <- tryCatch(plotSpectra(params, total = TRUE, biomass = TRUE), error = function(e) NULL)
if (!is.null(p_spectra)) print(p_spectra)

dev.off()
cat("Wrote diagnostic PDF to", out_file, "\n")
