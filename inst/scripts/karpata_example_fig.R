library(mizer)
library(mizerExperimental)
library(assertthat)
library(mizerReef)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(knitr)
library(patchwork)
library(here)

karpata_model <- caribbean_10_model # karpata_model was renamed to caribbean_10_model


# # Steady State
# ## Species parameter table
# ```{r tabParams}
# #| label: tbl-params
# #| tbl-cap: Species parameters
# sp %>%
#     select(biomass_observed, biomass_cutoff,
#            l_max, #l_mat, 
#            age_mat, a, b, beta) %>%
#     kable(col.names = c("Biomass (g/m^2)", "Cutoff Size (g)",
#                         "Maximum Length (cm)", #"Maturation Length (cm)", 
#                         "Maturation Age (yr)", 
#                         "Length to Weight: a", "Length to Weight: b",
#                         "Predator-Prey Mass Ratio (PPMR)"))
# ```

non_complex <- newRefuge(karpata_model, new_method = "noncomplex")

non_complex <- non_complex |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()


# Fisheries productivity and biomass comparison, L = 20 and L = 40
all_biom11 <- plot2TotalBiomass(non_complex, karpata_model,
                                name1 = "Flat", 
                                name2 = "Complex",
                                stack = TRUE) +
    scale_alpha_manual(values = c(1,1)) + 
    labs(x="", fill ="") +
    guides(alpha = "none") +
    theme_classic()

all_prod11 <- plot2Productivity(non_complex, karpata_model,
                                name1 = "Flat", 
                                name2 = "Complex",
                                stack = TRUE) +
    scale_alpha_manual(values = c(1,1)) + 
    labs(x="", fill ="") +
    guides(alpha = "none") +
    theme_classic()


sum11 <- ggarrange(all_biom11, all_prod11, 
                   ncol = 1, nrow = 2,
                   common.legend = TRUE,
                   legend = "right")

sum12 <- ggarrange(all_biom11, all_prod11, 
                   ncol = 2, nrow = 1,
                   common.legend = TRUE,
                   legend = "top")

# sum12 is the horizontal (side-by-side) layout - not used in the
# manuscript, kept here for reference/comparison only.
ggsave("karpata_example_horizontal.png",
       plot = sum12,
       height = 14,
       width = 20,
       unit = "cm",
       dpi = 600,
       bg = "white")

# sum11 is the vertical (stacked) layout used in the manuscript - saved
# last/under the canonical filename so it's what actually ends up at
# vignettes/figures/karpata_example.png.
ggsave("karpata_example.png",
       plot = sum11,
       height = 20,
       width = 12,
       unit = "cm",
       dpi = 600,
       bg = "white")
