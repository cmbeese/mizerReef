# Refuge Holes -----------------------------------------------------------------
library(tidyverse)
library(here)
library(RColorBrewer)
library(ggpubr)

# Read in refuge hole data
holes <- read.csv(here("data/hole_density.csv"))

# Save size range as a factor
holes$range <- factor(holes$range,
    levels = c(
        "0-5", "5-10",
        "10-15", "15-20",
        "20-25", "25-30",
        "30-35", "35-40",
        "40-45", "45-50", "50+"
    )
)

# Sum over size ranges
hole_sums <- holes %>%
    group_by(reefName, range) %>%
    summarise(
        total_in_range = sum(hole.density),
        density = total_in_range / 120
    )

# Prepare start/end bins
start_L <- seq(0, 45, 5)
end_L <- seq(5, 50, 5)

# Get all reef names
reef_names <- unique(hole_sums$reefName)

# Create lists to store data frames and plots
refuge_dfs <- list()
plots <- list()

for (reef in reef_names) {
    dens <- hole_sums$density[hole_sums$reefName == reef]
    # Remove 50+ bin by merging with previous
    if (length(dens) == 11) {
        dens[10] <- dens[11]
        dens <- dens[-11]
    }
    df <- data.frame(start_L, end_L, refuge_density = dens)
    refuge_dfs[[reef]] <- df
    # Save as CSV and RDA
    write.csv(df, file = paste0("data/", reef, "_refuge.csv"), row.names = FALSE)
    save(df, file = paste0("data/", reef, "_refuge.rda"))
    # Plot
    p <- holes %>%
        mutate(hole.density = hole.density / 10) %>%
        mutate(range = paste(range, "cm")) %>%
        filter(reefName == reef) %>%
        ggplot(aes(x = range, y = hole.density, fill = range)) +
        scale_y_log10() +
        geom_boxplot() +
        scale_fill_brewer(palette = "RdYlBu") +
        labs(
            x = "Refuge Size",
            title = reef,
            y = expression("Log"[10] * "(Refuge Density) (no./m"^2 * ")")
        ) +
        guides(fill = "none") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plots[[reef]] <- p
    ggsave(paste0("figures/", reef, "_holes.png"), plot = p, width = 9, height = 7, units = "cm", bg = "white")
}

# Optionally, arrange all plots together
# library(ggpubr)
# all_holes <- ggarrange(plotlist = plots, ncol = 2)
# print(all_holes)
