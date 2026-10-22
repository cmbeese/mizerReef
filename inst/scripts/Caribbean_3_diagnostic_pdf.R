# Diagnostic PDF for a caribbean_3 model: targets, diets, growth, erepro,
# stability, complex vs flat reef and refuge use, one diagnostic per page.
# Page layout and visual grammar (target-vs-model points, literature brackets
# as shaded bands) borrowed from caribbean-10-calibration's
# generate_diagnostic_pdfs.R.
#
# Run from the repo root. With no argument it checks the bundled
# caribbean_3_model; give a saved params file (.rds) to check that instead:
#   Rscript inst/scripts/Caribbean_3_diagnostic_pdf.R
#   Rscript inst/scripts/Caribbean_3_diagnostic_pdf.R path/to/model.rds
# The PDF goes to artifacts/ (git-ignored). Takes under a minute.
suppressMessages({
  library(mizer)
  library(mizerReef)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  model_name <- tools::file_path_sans_ext(basename(args[1]))
  params <- readRDS(args[1])
} else {
  model_name <- "caribbean_3_model"
  params <- mizerReef::caribbean_3_model
}
sp <- params@species_params
colours <- params@linecolour

out_dir <- "artifacts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_file <- file.path(out_dir, paste0(model_name, "-diagnostic-", timestamp, ".pdf"))

# Targets ---------------------------------------------------------------------
# FORCE survey biomasses (Bonaire), above 3.125 g
target_biomass <- c(predators = 107, herbivores = 34, inverts = 40)
# Rogers et al. (2018): recruitment set to give about 10 small fish (< 4 cm)
# per m^2 in each fish guild
target_small_fish <- c(predators = 10, herbivores = 10)
# Rogers et al. (2018) Table S2: predators' fraction of feeding time on each
# prey group, used here as adult diet shares
target_diet <- c(predators = 0.33, herbivores = 0.17, inverts = 0.5)
# Algae: 150-550 g wet weight/m^2 (midpoint 350); detritus: residence time
# 2.8-4.9 h at Palmyra Atoll (Max et al. 2013, midpoint 3.85)
algae_bracket <- c(lower = 150, preferred = 350, upper = 550)
detritus_bracket <- c(lower = 2.8, preferred = 3.85, upper = 4.9)

# Helpers ---------------------------------------------------------------------
numbers_between <- function(p, w_lo, w_hi) {
  sel <- p@w >= w_lo & p@w < w_hi
  rowSums(p@initial_n[, sel, drop = FALSE] *
            rep(p@dw[sel], each = nrow(p@initial_n)))
}
small_fish <- function(p) numbers_between(p, 10^-1.5, 1.6)[names(target_small_fish)]

# Share of predators, herbivores and invertebrates in what predators above
# their maturity size eat (plankton left out, as in the tutorial's targets)
adult_predator_diet <- function(p) {
  diet <- getDiet(p, proportion = FALSE)
  adult <- p@w >= p@species_params["predators", "w_mat"]
  eaten <- colSums(diet["predators", adult, ] *
                     p@initial_n["predators", adult] * p@dw[adult])
  list(prey = eaten[names(target_diet)] / sum(eaten[names(target_diet)]),
       plankton = eaten[["Resource"]] / sum(eaten))
}

# Refuge holes and refuge-using fish in each hole size class. A class is
# contested when there are fewer holes than max_protect times the fish
refuge_table <- function(p) {
  rp <- p@other_params$refuge_params
  md <- rp$method_params
  w_lo <- pmax(rp$a_bar * md[[1]]^rp$b_bar, rp$w_settle)
  w_hi <- rp$a_bar * md[[2]]^rp$b_bar
  users <- p@species_params$species[p@species_params$refuge_user %in% TRUE]
  vulnerable <- getVulnerable(p)
  rows <- lapply(seq_len(nrow(md)), function(k) {
    wi <- which(p@w >= w_lo[k] & p@w < w_hi[k])
    data.frame(bin = paste0(md[[1]][k], "-", md[[2]][k], " cm"),
               holes = md[[3]][k],
               fish = sum(numbers_between(p, w_lo[k], w_hi[k])[users]),
               protected = if (length(wi)) 1 - mean(vulnerable[users[1], wi]) else NA)
  })
  df <- do.call(rbind, rows)
  df$bin <- factor(df$bin, levels = df$bin)
  df$contested <- df$holes < rp$max_protect * df$fish
  df
}

# Open reef: recruits arrive at the complex reef's rate on both reefs (as in
# Rogers et al. 2018, which had constant recruitment); everything else runs
# freely until it settles
hold_recruitment <- function(p, rdd) {
  p@species_params$constant_reproduction <- rdd
  p@rates_funcs$RDD <- "constantRDD"
  p
}
# Settle a reef without changing any parameter (unlike reefSteady(), which
# refits reproduction and detritus input). If it settles onto a limit cycle
# rather than a steady state, follow the cycle for ten periods (at least 50
# years) and average over that.
settle_reef <- function(p) {
  s <- projectUntilSettled(p, t_max = 1000, t_save = 10,
                           progress_bar = FALSE, info_level = 0)
  conv <- attr(s, "convergence")
  t_follow <- if (conv$attractor == "limit_cycle") max(50, 10 * conv$period) else 5
  sim <- project(finalParams(s), t_max = t_follow, t_save = 0.5,
                 progress_bar = FALSE)
  times <- getTimes(sim)[-1]
  fish <- c("predators", "herbivores")
  list(conv = conv, sim = sim,
       biomass = getBiomass(sim)[-1, sp$species, drop = FALSE],
       productivity = getProductivity(sim, time_range = range(times))[, fish, drop = FALSE])
}
run_open_reef <- function(p) {
  rdd <- getRDD(p)
  list(complex = settle_reef(hold_recruitment(p, rdd)),
       flat = settle_reef(hold_recruitment(
         newRefuge(p, new_method = "noncomplex", info_level = 0), rdd)))
}
describe_settling <- function(conv) {
  if (conv$attractor == "limit_cycle") {
    sprintf("limit cycle, period %.1f y, biomass amplitude %.0f%%",
            conv$period, 100 * conv$amplitude)
  } else {
    paste0(gsub("_", " ", conv$attractor), " (", conv$termination, ")")
  }
}
# Model and von Bertalanffy size of the fish at a given age
size_at_age <- function(p, age) {
  g <- plotGrowthCurves(p, return_data = TRUE)
  pick <- function(legend) {
    sapply(c("predators", "herbivores"), function(s) {
      d <- g[g$Species == s & g$Legend == legend, ]
      d[["Size [g]"]][which.min(abs(d$Age - age))]
    })
  }
  list(model = pick("model"), vb = pick("von Bertalanffy"))
}

# Shaded-band target plot, as in the caribbean-10 biomass match
plot_target_match <- function(df, title, y_lab, subtitle = NULL) {
  ggplot(df, aes(x = name)) +
    geom_linerange(data = subset(df, !is.na(lower)),
                   aes(ymin = lower, ymax = upper, colour = colour),
                   linewidth = 9, lineend = "butt", alpha = 0.28,
                   show.legend = FALSE) +
    geom_point(aes(y = target, shape = "Target", colour = colour),
               fill = "white", size = 4, stroke = 1.2) +
    geom_point(aes(y = model, shape = "Model", colour = colour), size = 4) +
    scale_shape_manual(name = NULL, values = c(Target = 23, Model = 16),
                       guide = guide_legend(override.aes = list(
                         colour = "grey30", fill = "grey30"))) +
    scale_colour_identity() +
    facet_wrap(~panel, scales = "free") +
    theme_bw() +
    labs(title = title, subtitle = subtitle, x = NULL, y = y_lab)
}

# Calculations ----------------------------------------------------------------
message("Projecting 100 years with Beverton-Holt reproduction live")
sim <- project(params, t_max = 100, t_save = 1, progress_bar = FALSE)
b_sim <- getBiomass(sim)[, sp$species]
drift <- max(abs(b_sim[nrow(b_sim), ] / b_sim[1, ] - 1))

message("Settling the complex and flat open reef")
reefs <- run_open_reef(params)
prod_total <- lapply(reefs, function(r) rowSums(r$productivity))
prod_ratio <- mean(prod_total$complex) / mean(prod_total$flat)
prod_ratio_range <- mean(prod_total$complex) / rev(range(prod_total$flat))

biomass <- getBiomass(params, use_cutoff = TRUE)
diet <- adult_predator_diet(params)
refuge <- refuge_table(params)
age_input <- setNames(as.numeric(sp$age_mat), sp$species)
age_model <- age_mat(params)

pdf(out_file, width = 11, height = 8.5, onefile = TRUE)

# Page 1: summary -------------------------------------------------------------
fmt <- function(x, digits = 3) paste(signif(x, digits), collapse = " / ")
pct_off <- function(model, target) sprintf("%+.1f%%", 100 * (model / target - 1))
lines <- c(
  paste("Model:", model_name, " | generated", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "Groups: predators / herbivores / inverts",
  "",
  "TARGETS                                   target        model         off by",
  sprintf("  Biomass > 3.125 g [g/m^2]               %-13s %-13s %s",
          fmt(target_biomass), fmt(biomass[names(target_biomass)]),
          paste(pct_off(biomass[names(target_biomass)], target_biomass), collapse = " ")),
  sprintf("  Small fish 1.1-4 cm [per m^2] (pred/herb) %-11s %-13s %s",
          fmt(target_small_fish), fmt(small_fish(params)),
          paste(pct_off(small_fish(params), target_small_fish), collapse = " ")),
  sprintf("  Adult predator diet (pred/herb/inv)     %-13s %-13s plankton %.1f%%",
          fmt(target_diet, 2), fmt(diet$prey, 2), 100 * diet$plankton),
  sprintf("  Predator age at maturity [y]            %-13s %-13s",
          fmt(age_input[["predators"]]), fmt(age_model[["predators"]])),
  sprintf("  Size at age 5 [g] (pred/herb), vs von Bertalanffy %s vs %s",
          fmt(size_at_age(params, 5)$model), fmt(size_at_age(params, 5)$vb)),
  sprintf("  Algae [g WW/m^2] / detritus life [h]    %-13s %-13s",
          paste(algae_bracket[["preferred"]], "/", detritus_bracket[["preferred"]]),
          paste(signif(algae_biomass(params), 3), "/",
                signif(detritus_lifetime(params) * 24 * 365, 3))),
  "",
  "REPRODUCTION AND STABILITY",
  paste("  erepro:", fmt(sp$erepro), if (any(sp$erepro > 1)) " <-- above 1" else " (all below 1)"),
  paste("  Reproduction level:", fmt(getReproductionLevel(params), 2)),
  paste("  Recruits per m^2 per year:", fmt(getRDD(params))),
  paste("  Largest biomass change over 100 y, Beverton-Holt live:", signif(drift, 2)),
  "",
  "COMPLEX VS FLAT REEF (open reef: recruitment held at the complex reef's rate; no parameter refitted)",
  paste("  Complex reef settles to:", describe_settling(reefs$complex$conv)),
  paste("  Flat reef settles to:   ", describe_settling(reefs$flat$conv)),
  paste0("  Fish productivity complex / flat: ", round(prod_ratio, 2),
         " (averaged over the flat reef's cycle; ", round(prod_ratio_range[1], 1),
         "-", round(prod_ratio_range[2], 1), " across the cycle)"),
  paste("  Biomass complex (mean):", fmt(colMeans(reefs$complex$biomass))),
  paste("  Biomass flat (mean):   ", fmt(colMeans(reefs$flat$biomass)),
        " | predators range", fmt(range(reefs$flat$biomass[, "predators"]))),
  "",
  "REFUGE",
  paste0("  Hole size classes with fewer holes than fish: ", sum(refuge$contested),
         " of ", nrow(refuge),
         if (any(refuge$contested))
           paste0(" (", paste(refuge$bin[refuge$contested], collapse = ", "), ")") else ""),
  "",
  "SETTINGS",
  paste0("  Refuge: ", params@other_params$refuge_params$method,
         ", from ", signif(params@other_params$refuge_params$w_settle, 3), " g",
         ", max_protect ", params@other_params$refuge_params$max_protect),
  paste0("  Plankton up to ", signif(params@resource_params$w_pp_cutoff, 3),
         " g, kappa ", signif(params@resource_params$kappa, 4)),
  paste0("  Algae production ", params@other_params$algae$growth,
         " g/m^2/y, capacity ", signif(params@other_params$algae$capacity, 4),
         " g/m^2; detritus capacity ", params@other_params$detritus$capacity,
         ", external input ", signif(params@other_params$detritus$external, 4)),
  paste("  Predator interaction (pred/herb/inv):", fmt(params@interaction["predators", ])),
  paste("  ks:", fmt(sp$ks), " | gamma:", fmt(sp$gamma), " | h:", fmt(sp$h)),
  paste("  Satiation:", paste(sp$satiation, collapse = " / "),
        " (FALSE = no maximum intake, so h has no effect)")
)
old_mar <- par("mar")
par(mar = c(1, 1, 1, 1))
plot.new()
text(0, 1, paste(lines, collapse = "\n"), adj = c(0, 1), cex = 0.8, family = "mono")
par(mar = old_mar)

# Page 2: plot(params): predation mortality, feeding level, spectra ----------
# plot() draws to the device itself; printing it again would add a page
tryCatch(plot(params), error = function(e) message("plot(params) failed: ", conditionMessage(e)))

# Page 3: biomass and small fish vs targets -----------------------------------
df_targets <- rbind(
  data.frame(panel = "Biomass above 3.125 g [g/m^2] (FORCE)",
             name = names(target_biomass), target = target_biomass,
             model = biomass[names(target_biomass)], lower = NA, upper = NA,
             colour = colours[names(target_biomass)]),
  data.frame(panel = "Small fish 1.1-4 cm [per m^2] (Rogers et al. 2018: ~10)",
             name = names(target_small_fish), target = target_small_fish,
             model = small_fish(params), lower = NA, upper = NA,
             colour = colours[names(target_small_fish)]))
print(plot_target_match(df_targets, "Biomass and small fish vs targets", NULL))

# Page 4: algae and detritus vs literature brackets ---------------------------
df_resources <- data.frame(
  panel = c("Algae standing biomass [g WW/m^2]", "Detritus residence time [hours]"),
  name = c("algae", "detritus"),
  target = c(algae_bracket[["preferred"]], detritus_bracket[["preferred"]]),
  model = c(algae_biomass(params), detritus_lifetime(params) * 24 * 365),
  lower = c(algae_bracket[["lower"]], detritus_bracket[["lower"]]),
  upper = c(algae_bracket[["upper"]], detritus_bracket[["upper"]]),
  colour = colours[c("algae", "detritus")])
print(plot_target_match(
  df_resources, "Algae and detritus vs literature brackets", NULL,
  "Shaded band = literature range: algae 150-550 g WW/m^2; detritus 2.8-4.9 h (Max et al. 2013)"))

# Page 5: diets by size --------------------------------------------------------
print(plotDiet(params) + ggtitle("Diet by body size"))

# Page 6: adult predator diet vs Rogers et al. (2018) -------------------------
df_diet <- data.frame(
  prey = factor(rep(names(target_diet), 2), levels = names(target_diet)),
  source = rep(c("Rogers et al. (2018) weights", "Model"), each = 3),
  share = c(target_diet, diet$prey))
print(ggplot(df_diet, aes(x = prey, y = share, fill = prey, alpha = source)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7,
                 colour = "grey30") +
        scale_fill_manual(values = colours[names(target_diet)], guide = "none") +
        scale_alpha_manual(values = c(Model = 1, "Rogers et al. (2018) weights" = 0.35),
                           name = NULL) +
        scale_y_continuous(labels = scales::percent) +
        theme_bw() +
        labs(title = "Adult predator diet vs Rogers et al. (2018)",
             subtitle = sprintf("Predators above maturity size; plankton is %.1f%% of their intake and left out of the shares",
                                100 * diet$plankton),
             x = "Prey", y = "Share of fish and invertebrate prey eaten"))

# Page 7: growth curves --------------------------------------------------------
gdat <- plotGrowthCurves(params, return_data = TRUE)
gdat$Size_g <- gdat[["Size [g]"]]
print(ggplot(gdat, aes(x = Age, y = Size_g, colour = Legend, linetype = Legend)) +
        geom_line(linewidth = 0.5) +
        scale_y_log10() +
        facet_wrap(~Species, scales = "free_y") +
        theme_bw() +
        labs(y = "Size [g]", title = "Growth curves (model vs von Bertalanffy)") +
        theme(legend.position = "bottom"))

# Page 8: age at maturity ------------------------------------------------------
df_age <- data.frame(species = sp$species, input = age_input[sp$species],
                     model = age_model[sp$species])
print(ggplot(df_age, aes(x = species, colour = species)) +
        geom_segment(aes(xend = species, y = input, yend = model),
                     linewidth = 0.7, alpha = 0.6, show.legend = FALSE,
                     na.rm = TRUE) +
        geom_point(aes(y = input, shape = "Species table (input)"), size = 4,
                   na.rm = TRUE) +
        geom_point(aes(y = model, shape = "Model (realised)"), size = 4) +
        scale_shape_manual(name = NULL, values = c("Species table (input)" = 17,
                                                   "Model (realised)" = 16)) +
        scale_colour_manual(values = colours, guide = "none") +
        theme_bw() +
        labs(title = "Age at maturity: species table vs model",
             subtitle = "Only the predators' age at maturity is a calibration target; inverts have none in the table",
             x = NULL, y = "Age at maturity [years]"))

# Page 9: erepro ---------------------------------------------------------------
df_erepro <- data.frame(species = sp$species, erepro = sp$erepro)
print(ggplot(df_erepro, aes(x = reorder(species, erepro), y = erepro, colour = species)) +
        geom_point(size = 4) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick") +
        scale_colour_manual(values = colours, guide = "none") +
        scale_y_log10() +
        coord_flip() +
        theme_bw() +
        labs(title = "Reproductive efficiency (erepro)",
             subtitle = paste0("Dashed line = 1; above it a species would turn more than all its reproductive energy into eggs\n",
                               "Reproduction level: ", fmt(getReproductionLevel(params), 2)),
             x = NULL, y = "erepro (log scale)"))

# Page 10: stability -----------------------------------------------------------
print(plotBiomass(sim) +
        ggtitle("100 years with Beverton-Holt reproduction live",
                subtitle = paste("Largest biomass change:", signif(drift, 2))))

# Page 11: complex vs flat reef ------------------------------------------------
# Means over the averaging window (over the cycle, if there is one)
fish <- c("predators", "herbivores")
df_reefs <- do.call(rbind, lapply(c("Flat", "Complex"), function(reef) {
  r <- reefs[[tolower(reef)]]
  rbind(data.frame(reef = reef, species = fish, measure = "Biomass [g/m^2]",
                   value = colMeans(r$biomass[, fish, drop = FALSE])),
        data.frame(reef = reef, species = fish,
                   measure = "Productivity [g/m^2/year]",
                   value = colMeans(r$productivity)))
}))
df_reefs$reef <- factor(df_reefs$reef, levels = c("Flat", "Complex"))
df_reefs$species <- factor(df_reefs$species, levels = fish)
print(ggplot(df_reefs, aes(x = reef, y = value, fill = species)) +
        geom_col(colour = "black", width = 0.7) +
        facet_wrap(~measure, scales = "free_y") +
        scale_fill_manual(values = colours[fish], name = NULL) +
        theme_classic() +
        labs(title = paste0("Complex vs flat reef: fish productivity ",
                            round(prod_ratio, 1), " times higher on the complex reef"),
             subtitle = paste0("Open reef (recruitment held at the complex reef's rate). ",
                               "Flat reef: ", describe_settling(reefs$flat$conv),
                               "; bars are means over ",
                               round(max(getTimes(reefs$flat$sim))), " years"),
             x = NULL, y = NULL))

# Page 12: the flat reef over the averaging window -------------------------
print(plotBiomass(reefs$flat$sim) +
        ggtitle("Flat reef after settling",
                subtitle = describe_settling(reefs$flat$conv)))

# Page 13: refuge profile ------------------------------------------------------
print(plotRefugeProfile(params) + ggtitle("Refuge profile"))

# Page 14: refuge competition --------------------------------------------------
# Points rather than bars: on a log axis bars have no honest baseline
print(ggplot(refuge, aes(x = bin)) +
        geom_segment(aes(xend = bin, y = holes, yend = fish), colour = "grey70") +
        geom_point(aes(y = holes, colour = "Holes"), size = 4, shape = 15) +
        geom_point(aes(y = fish, colour = "Refuge-using fish"), size = 4) +
        geom_text(aes(y = pmax(holes, fish) * 1.8,
                      label = paste0(round(100 * protected), "%")), size = 3.5) +
        scale_y_log10() +
        scale_colour_manual(values = c(Holes = "grey40", "Refuge-using fish" = "#2C7FB8"),
                            name = NULL) +
        theme_bw() +
        labs(title = "Refuge competition: holes vs fish in each hole size class",
             subtitle = paste0("Label = share of fish protected (at most max_protect = ",
                               params@other_params$refuge_params$max_protect,
                               "); a class is contested when there are fewer holes than fish"),
             x = "Hole size class", y = "Number per m^2 (log scale)") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)))

dev.off()
cat("Wrote diagnostic PDF to", out_file, "\n")
