# Create one coefficient figure summarizing the main empirical results.

suppressPackageStartupMessages({
  library(cowplot)
  library(data.table)
  library(ggplot2)
})

project_dir <- normalizePath(
  Sys.getenv("THESIS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
core_dir <- file.path(project_dir, "tables", "chapter5_core")
robustness_dir <- file.path(project_dir, "tables", "chapter5_robustness")
figures_dir <- file.path(project_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
figure_path <- file.path(
  figures_dir,
  "results_chapter_summary.png"
)

read_core <- function(filename) {
  fread(file.path(core_dir, filename))
}

read_robustness <- function(filename) {
  fread(file.path(robustness_dir, filename))
}

standardize_models <- function(x) {
  fifelse(
    x %chin% c("Pooled OLS", "Pooled OLS + wave FE"), "Pooled OLS",
    fifelse(
      x %chin% c("Random effects", "Random effects + wave FE"), "Random effects",
      fifelse(x == "Fixed effects + wave FE", "Fixed effects", x)
    )
  )
}

make_panel <- function(data, panel, category, model, estimate, std_error) {
  data.table(
    panel = panel,
    category = category,
    model = model,
    estimate = estimate,
    std_error = std_error
  )
}

linear <- read_core("table_panel_benchmarks_numeric.csv")
linear_panel <- make_panel(
  linear,
  "A. Linear panel estimates",
  linear$outcome_label,
  standardize_models(linear$model),
  linear$estimate,
  linear$std_error
)

plr <- read_core("table_dml_plr_learner_comparison_numeric.csv")
plr_panel <- make_panel(
  plr,
  "B. PLR-DML estimates",
  plr$outcome_label,
  plr$learner,
  plr$estimate,
  plr$std_error
)

irm <- read_core("table_dml_irm_learner_comparison_numeric.csv")
irm_panel <- make_panel(
  irm,
  "C. IRM-DML estimates",
  irm$outcome_label,
  irm$learner,
  irm$estimate,
  irm$std_error
)

extended <- read_robustness("robustness_extended_waves_numeric.csv")
ordered_recent <- extended[
  specification %chin% c("lo_rich", "lo_common")
]
ordered_recent[, category := fifelse(
  specification == "lo_rich",
  "Waves L-O, rich controls\n(N = 4,316)",
  "Waves L-O, common controls\n(N = 5,001)"
)]

all_wave <- read_robustness(
  "robustness_extended_waves_paired_comparison_numeric.csv"
)
all_wave_plr <- all_wave[panel == "PLR"]
ordered_all_wave <- rbindlist(list(
  data.table(
    model = "Fixed effects",
    estimate = all_wave_plr$fe_estimate[1L],
    std_error = all_wave_plr$fe_std_error[1L]
  ),
  all_wave_plr[, .(
    model = learner,
    estimate = dml_estimate,
    std_error = dml_std_error
  )]
))
ordered_all_wave[, category :=
  "Waves A-O, common controls\n(N = 33,184)"]

weekend <- read_robustness(
  "robustness_weekend_combined_comparison_numeric.csv"
)
weekend_plr <- weekend[panel == "PLR"]
ordered_weekend <- rbindlist(list(
  data.table(
    model = "Fixed effects",
    estimate = weekend_plr$fe_estimate[1L],
    std_error = weekend_plr$fe_std_error[1L]
  ),
  weekend_plr[, .(
    model = learner,
    estimate = dml_estimate,
    std_error = dml_std_error
  )]
))
ordered_weekend[, category :=
  "Weekend social-media use,\nWaves L-O, rich controls\n(N = 4,331)"]

ordered_robustness <- rbindlist(list(
  ordered_recent[, .(category, model, estimate, std_error)],
  ordered_all_wave[, .(category, model, estimate, std_error)],
  ordered_weekend[, .(category, model, estimate, std_error)]
), use.names = TRUE)
ordered_panel <- make_panel(
  ordered_robustness,
  "D. PLR robustness (ordered)",
  ordered_robustness$category,
  ordered_robustness$model,
  ordered_robustness$estimate,
  ordered_robustness$std_error
)

recent_irm_comparison <- read_core(
  "table_fe_irm_comparison_numeric.csv"
)
recent_binary_fe <- data.table(
  category = "Weekday high use, Waves L-O\n(N = 4,316)",
  model = "Fixed effects",
  estimate = recent_irm_comparison$fe_estimate[1L],
  std_error = sqrt(recent_irm_comparison$fe_variance[1L])
)
recent_binary_dml <- copy(irm[outcome == "life_dissatisfaction"])
recent_binary_dml[, category := "Weekday high use, Waves L-O\n(N = 4,316)"]
recent_binary_dml <- recent_binary_dml[, .(
  category,
  model = learner,
  estimate,
  std_error
)]

weekend_irm <- weekend[panel == "IRM"]
weekend_binary <- rbindlist(list(
  data.table(
    category = "Weekend high use, Waves L-O\n(N = 4,331)",
    model = "Fixed effects",
    estimate = weekend_irm$fe_estimate[1L],
    std_error = weekend_irm$fe_std_error[1L]
  ),
  weekend_irm[, .(
    category = "Weekend high use, Waves L-O\n(N = 4,331)",
    model = learner,
    estimate = dml_estimate,
    std_error = dml_std_error
  )]
))

all_wave_irm <- all_wave[panel == "IRM"]
all_wave_binary <- rbindlist(list(
  data.table(
    category = "Weekday high use, Waves A-O\n(N = 33,184)",
    model = "Fixed effects",
    estimate = all_wave_irm$fe_estimate[1L],
    std_error = all_wave_irm$fe_std_error[1L]
  ),
  all_wave_irm[, .(
    category = "Weekday high use, Waves A-O\n(N = 33,184)",
    model = learner,
    estimate = dml_estimate,
    std_error = dml_std_error
  )]
))

binary_robustness <- rbindlist(list(
  recent_binary_fe,
  recent_binary_dml,
  weekend_binary,
  all_wave_binary
), use.names = TRUE)
binary_panel <- make_panel(
  binary_robustness,
  "E. IRM robustness (binary)",
  binary_robustness$category,
  binary_robustness$model,
  binary_robustness$estimate,
  binary_robustness$std_error
)

plot_data <- rbindlist(list(
  linear_panel,
  plr_panel,
  irm_panel,
  ordered_panel,
  binary_panel
), use.names = TRUE)

outcome_order <- c(
  "Life dissatisfaction",
  "Loneliness",
  "School-work dissatisfaction",
  "School dissatisfaction"
)
ordered_order <- c(
  "Waves L-O, rich controls\n(N = 4,316)",
  "Waves L-O, common controls\n(N = 5,001)",
  "Waves A-O, common controls\n(N = 33,184)",
  "Weekend social-media use,\nWaves L-O, rich controls\n(N = 4,331)"
)
binary_order <- c(
  "Weekday high use, Waves L-O\n(N = 4,316)",
  "Weekend high use, Waves L-O\n(N = 4,331)",
  "Weekday high use, Waves A-O\n(N = 33,184)"
)
plot_data[, category := factor(
  category,
  levels = rev(c(outcome_order, ordered_order, binary_order))
)]
plot_data[, panel := factor(
  panel,
  levels = c(
    "A. Linear panel estimates",
    "B. PLR-DML estimates",
    "C. IRM-DML estimates",
    "D. PLR robustness (ordered)",
    "E. IRM robustness (binary)"
  )
)]
plot_data[, model := factor(
  model,
  levels = c(
    "Pooled OLS",
    "Random effects",
    "Fixed effects",
    "Elastic net",
    "Lasso",
    "Random forest",
    "Neural net"
  )
)]
plot_data[, `:=`(
  lower = estimate - 1.96 * std_error,
  upper = estimate + 1.96 * std_error
)]

palette <- c(
  "Pooled OLS" = "#B58CC8",
  "Random effects" = "#91A7D6",
  "Fixed effects" = "#D681A4",
  "Elastic net" = "#7F8FD2",
  "Lasso" = "#72B6B2",
  "Random forest" = "#E5A083",
  "Neural net" = "#CDB96F"
)
shapes <- c(
  "Pooled OLS" = 15,
  "Random effects" = 17,
  "Fixed effects" = 18,
  "Elastic net" = 16,
  "Lasso" = 1,
  "Random forest" = 2,
  "Neural net" = 0
)
dodge <- position_dodge(width = 0.66)

build_panel_plot <- function(panel_name, show_x_title = FALSE) {
  panel_data <- plot_data[panel == panel_name]
  ggplot(
    panel_data,
    aes(
      x = estimate,
      y = category,
      colour = model,
      shape = model
    )
  ) +
    geom_vline(
      xintercept = 0,
      colour = "#8D85A5",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_errorbar(
      aes(xmin = lower, xmax = upper),
      orientation = "y",
      width = 0.18,
      linewidth = 0.58,
      position = dodge
    ) +
    geom_point(
      size = 2.55,
      stroke = 0.95,
      position = dodge
    ) +
    facet_wrap(vars(panel), scales = "free") +
    scale_colour_manual(values = palette, drop = FALSE) +
    scale_shape_manual(values = shapes, drop = FALSE) +
    labs(
      x = if (show_x_title) {
        "Estimated effect (95% confidence interval)"
      } else {
        NULL
      },
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    theme_minimal(base_family = "serif", base_size = 10.5) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        colour = "#EEE9F1",
        linewidth = 0.4
      ),
      strip.text = element_text(
        face = "bold",
        colour = "#5A5165",
        size = 10.2,
        margin = margin(5, 5, 5, 5)
      ),
      strip.background = element_rect(
        fill = "#F2EAF7",
        colour = NA
      ),
      axis.text.y = element_text(
        colour = "#403A47",
        size = 8.7,
        lineheight = 0.94
      ),
      axis.text.x = element_text(
        colour = "#5A5165",
        size = 8.3
      ),
      axis.title.x = element_text(
        face = "bold",
        colour = "#4B4354",
        margin = margin(t = 6)
      ),
      legend.position = "none",
      plot.margin = margin(4, 8, 2, 4)
    )
}

extract_legend <- function(data) {
  legend_plot <- ggplot(
    droplevels(as.data.frame(data)),
    aes(
      x = estimate,
      y = category,
      colour = model,
      shape = model
    )
  ) +
    geom_point(size = 2.55) +
    scale_colour_manual(values = palette, drop = TRUE) +
    scale_shape_manual(values = shapes, drop = TRUE) +
    labs(colour = NULL, shape = NULL) +
    guides(
      colour = guide_legend(nrow = 2, byrow = TRUE),
      shape = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    theme_void(base_family = "serif") +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.margin = margin(0, 0, 0, 0)
    )

  legend_gtable <- ggplotGrob(legend_plot)
  legend_gtable$grobs[[
    which(legend_gtable$layout$name == "guide-box-bottom")
  ]]
}

main_legend <- extract_legend(
  plot_data[panel %chin% c(
    "A. Linear panel estimates",
    "B. PLR-DML estimates"
  )]
)
robustness_legend <- extract_legend(
  plot_data[panel %chin% c(
    "C. IRM-DML estimates",
    "D. PLR robustness (ordered)",
    "E. IRM robustness (binary)"
  )]
)

panel_a <- build_panel_plot("A. Linear panel estimates")
panel_b <- build_panel_plot("B. PLR-DML estimates")
panel_c <- build_panel_plot("C. IRM-DML estimates")
panel_d <- build_panel_plot("D. PLR robustness (ordered)")
panel_e <- build_panel_plot(
  "E. IRM robustness (binary)",
  show_x_title = TRUE
)

summary_plot <- plot_grid(
  panel_a,
  panel_b,
  main_legend,
  ncol = 1,
  rel_heights = c(1, 1, 0.12)
)

robustness_plot <- plot_grid(
  panel_c,
  panel_d,
  panel_e,
  robustness_legend,
  ncol = 1,
  rel_heights = c(1, 1, 0.78, 0.12)
)

ggsave(
  figure_path,
  summary_plot,
  width = 9.6,
  height = 10.2,
  units = "in",
  dpi = 320,
  bg = "white"
)

robustness_figure_path <- file.path(
  figures_dir,
  "results_chapter_summary_robustness.png"
)
ggsave(
  robustness_figure_path,
  robustness_plot,
  width = 9.6,
  height = 11.2,
  units = "in",
  dpi = 320,
  bg = "white"
)

cat("Wrote", figure_path, "\n")
cat("Wrote", robustness_figure_path, "\n")
