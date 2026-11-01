# Plots how the non-NA fold_change cells of merged_fc/merged_fold_change.RData are
# distributed, saved as merged_fc/fold_change_distribution.png: their log2 values,
# the number of cells per key and per pair, and their position in the matrix. NaN
# counts as NA, as is.na() has it. Cells of 0 and Inf have no place on the log2 axis
# but count as cells everywhere else.
# Run from the package root: Rscript --vanilla merged_fc/plot_distribution.R

library(ggplot2)
library(patchwork)

invisible(loadNamespace("SummarizedExperiment"))

ink <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted <- "#898781"
surface <- "#fcfcfb"
grid <- "#e1e0d9"
baseline <- "#c3c2b7"
mark <- "#2a78d6"
value_limit <- 6

env <- new.env()
load("merged_fc/merged_fold_change.RData", envir = env)
fold_change <- SummarizedExperiment::assay(env$merged_fc, "fold_change")
sources <- env$merged_fc$source
set_of_pair <- ifelse(startsWith(sources, "MTBLS"), "mztabm_fc", "mwtab_fc")
has_value <- !is.na(fold_change)
per_key <- rowSums(has_value)
per_pair <- colSums(has_value)

theme_panel <- theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    plot.background = element_rect(fill = surface, color = NA),
    panel.background = element_rect(fill = surface, color = NA),
    plot.title = element_text(color = ink, face = "bold", size = 12),
    plot.subtitle = element_text(color = ink_secondary, size = 9.5,
                                 margin = margin(b = 8)),
    axis.title = element_text(color = ink_secondary, size = 9.5),
    axis.text = element_text(color = ink_muted, size = 9),
    panel.grid.major = element_line(color = grid, linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  )

# Log2 values -------------------------------------------------------------------

positive <- fold_change[is.finite(fold_change) & fold_change > 0]
log_values <- log2(positive)
outside <- sum(abs(log_values) > value_limit)
bins <- hist(log_values[abs(log_values) <= value_limit],
             breaks = seq(-value_limit, value_limit, by = 0.1), plot = FALSE)
n_down <- sum(log_values <= -1)
n_up <- sum(log_values >= 1)
n_within <- sum(abs(log_values) < 1)
region_label <- function(count) {
  sprintf("%s cells\n%.1f%%", format(count, big.mark = ","),
          100 * count / length(log_values))
}
top <- max(bins$counts)

value_panel <- ggplot(data.frame(x = bins$mids, n = bins$counts), aes(x, n)) +
  geom_col(width = 0.085, fill = mark) +
  geom_vline(xintercept = c(-1, 1), color = ink_muted, linewidth = 0.3) +
  annotate("text", x = -1.15, y = top * 0.82, hjust = 1, size = 3.2, color = ink,
           lineheight = 0.95,
           label = paste0("fold change <= 0.5\n", region_label(n_down))) +
  annotate("text", x = 1.15, y = top * 0.82, hjust = 0, size = 3.2, color = ink,
           lineheight = 0.95,
           label = paste0("fold change >= 2\n", region_label(n_up))) +
  scale_x_continuous(breaks = seq(-value_limit, value_limit, by = 2),
                     limits = c(-value_limit, value_limit), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)),
                     labels = function(y) format(y, big.mark = ",", trim = TRUE)) +
  labs(
    title = "Most fold changes sit close to 1",
    subtitle = sprintf(
      "Median %.2f, and %s of %s values (%.1f%%) lie within 2-fold. %d values with |log2| > %d are not drawn.",
      2^stats::median(log_values), format(n_within, big.mark = ","),
      format(length(log_values), big.mark = ","),
      100 * n_within / length(log_values), outside, value_limit
    ),
    x = "log2 fold change", y = "Cells per 0.1 bin"
  ) +
  theme_panel +
  theme(panel.grid.major.x = element_blank())

# Cells per key and per pair ----------------------------------------------------

count_panel <- function(counts, breaks, labels, title, subtitle, x) {
  tallied <- as.integer(table(cut(counts, breaks, labels = labels)))
  stopifnot(sum(tallied) == length(counts))
  ggplot(data.frame(bin = factor(labels, levels = labels), n = tallied),
         aes(bin, n)) +
    geom_col(width = 0.72, fill = mark) +
    geom_text(aes(label = n), vjust = -0.5, size = 3.1, color = ink) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(title = title, subtitle = subtitle, x = x, y = NULL) +
    theme_panel +
    theme(panel.grid.major = element_blank(), axis.text.y = element_blank())
}

key_panel <- count_panel(
  per_key, c(-1, 0, 1, 2, 5, 10, 20, 50, Inf),
  c("0", "1", "2", "3-5", "6-10", "11-20", "21-50", "51+"),
  "Keys are measured in few pairs",
  sprintf("Median %d of %d pairs, at most %d. %d keys have no value at all.",
          as.integer(stats::median(per_key)), ncol(fold_change), max(per_key),
          sum(per_key == 0)),
  "Non-NA cells per key"
)

pair_panel <- count_panel(
  per_pair, c(-1, 0, 10, 25, 50, 100, 200, Inf),
  c("0", "1-10", "11-25", "26-50", "51-100", "101-200", "201+"),
  "Pairs differ widely in coverage",
  sprintf("Median %d of %s keys, at most %d. %d pairs have no value at all.",
          as.integer(stats::median(per_pair)),
          format(nrow(fold_change), big.mark = ","), max(per_pair),
          sum(per_pair == 0)),
  "Non-NA cells per pair"
)

# Position in the matrix --------------------------------------------------------

# Keys run from most to least covered, and pairs likewise within each set.
key_rank <- rank(-per_key, ties.method = "first")
pair_rank <- rank(-per_pair, ties.method = "first")
cells <- which(has_value, arr.ind = TRUE)
pairs_of_set <- table(set_of_pair)
valued_of_set <- table(set_of_pair[per_pair > 0])
map_data <- data.frame(
  key = key_rank[cells[, 1]],
  pair = pair_rank[cells[, 2]],
  set = sprintf("%s\n%d of %d pairs", set_of_pair[cells[, 2]],
                valued_of_set[set_of_pair[cells[, 2]]],
                pairs_of_set[set_of_pair[cells[, 2]]])
)
stopifnot(nrow(map_data) == sum(has_value))
# Pair positions restart in each set, so that each facet starts at its densest.
map_data$pair <- stats::ave(map_data$pair, map_data$set, FUN = function(pair) {
  match(pair, sort(unique(pair)))
})

map_panel <- ggplot(map_data, aes(pair, key)) +
  geom_tile(width = 1, height = 1, fill = mark) +
  facet_grid(~set, scales = "free_x", space = "free_x") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0), limits = c(nrow(fold_change) + 0.5, 0.5)) +
  labs(
    title = sprintf("%s of %s cells hold a value (%.2f%%)",
                    format(sum(has_value), big.mark = ","),
                    format(length(has_value), big.mark = ","),
                    100 * mean(has_value)),
    subtitle = sprintf(
      "One mark per non-NA cell. Sorted by coverage, so the %d empty keys fill the bottom.",
      sum(per_key == 0)
    ),
    x = "Pairs with values, most covered first", y = "Keys, most covered first"
  ) +
  theme_panel +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = baseline, fill = NA, linewidth = 0.3),
    strip.text = element_text(color = ink, size = 9, hjust = 0, lineheight = 0.95),
    strip.clip = "off",
    panel.spacing.x = unit(8, "pt")
  )

figure <- value_panel + key_panel + pair_panel + map_panel +
  plot_layout(design = "AD\nBD\nCD", widths = c(1.15, 1)) &
  theme(plot.background = element_rect(fill = surface, color = NA))

ragg::agg_png("merged_fc/fold_change_distribution.png", width = 13, height = 10,
              units = "in", res = 200, background = surface)
print(figure)
invisible(grDevices::dev.off())

message(sprintf(
  "%d non-NA cells: %d positive finite (%d beyond the axis), %d zero, %d Inf; %d NaN excluded",
  sum(has_value), length(log_values), outside, sum(fold_change == 0, na.rm = TRUE),
  sum(fold_change == Inf, na.rm = TRUE), sum(is.nan(fold_change))
))
