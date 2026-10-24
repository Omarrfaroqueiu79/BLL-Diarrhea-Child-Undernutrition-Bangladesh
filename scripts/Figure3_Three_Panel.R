# ==============================================================================
# 13. FIGURE 3: 3-Panel Malnutrition Matrix (Collision-Free Labels)
# ==============================================================================
library(survey)
library(ggplot2)
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

create_prev_panel <- function(outcome_var, panel_title, y_max) {
  form <- as.formula(paste("~", outcome_var))
  raw_res <- svyby(form, ~joint_exp, svy_primary, svymean, na.rm = TRUE)
  
  vals <- as.numeric(coef(raw_res)) * 100
  ses  <- as.numeric(SE(raw_res)) * 100
  
  prev_data <- data.frame(
    Exposure_Group = factor(
      c("Low BLL (<3.5)\nNo Diarrhea\n(Reference)", 
        "Low BLL (<3.5)\n+ Diarrhea", 
        "High BLL (\u22653.5)\nNo Diarrhea", 
        "High BLL (\u22653.5)\n+ Diarrhea"),
      levels = c("Low BLL (<3.5)\nNo Diarrhea\n(Reference)", 
                 "Low BLL (<3.5)\n+ Diarrhea", 
                 "High BLL (\u22653.5)\nNo Diarrhea", 
                 "High BLL (\u22653.5)\n+ Diarrhea")
    ),
    pct    = vals,
    se_pct = ses,
    ci_hi  = vals + 1.96 * ses,
    ci_lo  = pmax(0, vals - 1.96 * ses),
    stringsAsFactors = FALSE
  )
  
  ggplot(prev_data, aes(x = Exposure_Group, y = pct, fill = Exposure_Group)) +
    geom_col(width = 0.58, show.legend = FALSE, color = "black", linewidth = 0.5) +
    # Error bars
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), 
                  width = 0.18, linewidth = 0.75, color = "#1E293B") +
    # Text labels placed CLEARLY ABOVE the top whisker (ZERO collision)
    geom_text(aes(y = 5, label = sprintf("%.1f%%", pct)), 
              vjust = 0, fontface = "bold", size = 3.8, color = "#0F172A") +
    scale_fill_manual(values = c("#A6C8E0", "#6BAED6", "#FDBB84", "#E34A33")) +
    coord_cartesian(ylim = c(0, y_max)) +
    theme_classic(base_size = 11) +
    labs(title = panel_title, x = NULL, y = "Survey-Weighted Prevalence (%) [95% CI]") +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text.x  = element_text(face = "bold", color = "black", size = 9),
      axis.text.y  = element_text(color = "black", size = 10),
      axis.title.y = element_text(face = "bold", size = 10),
      plot.margin  = margin(t = 16, r = 8, b = 6, l = 8)
    )
}

# Panels with expanded y-axis limits to accommodate labels cleanly
panel_A <- create_prev_panel("stunted_num", "A. Stunting (HAZ < -2 SD)", 42)
panel_B <- create_prev_panel("wasted_num",  "B. Wasting (WHZ < -2 SD)", 28)
panel_C <- create_prev_panel("underwt_num", "C. Underweight (WAZ < -2 SD)", 44)

# Combine and save
p_fig3 <- panel_A + panel_B + panel_C + plot_layout(ncol = 3)
ggsave("K:/BGD_2025_MICS7_v01_M/SPSS_Outputs/Figure3_ThreePanel_Malnutrition.png", p_fig3, width = 15.0, height = 5.5, dpi=300)
ggsave("K:/BGD_2025_MICS7_v01_M/SPSS_Outputs/Figure3_ThreePanel_Malnutrition.tif", p_fig3, width = 15.0, height = 5.5, dpi=300, compression="lzw")