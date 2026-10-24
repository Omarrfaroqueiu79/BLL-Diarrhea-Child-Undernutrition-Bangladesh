# ==============================================================================
# SUPPLEMENTARY ANALYSES & REVIEWER REVISIONS SCRIPT
# Addresses comments regarding: Missing data flow, Splines, Hemoglobin, and 5ug/dL
# ==============================================================================

library(rms)
library(survey)
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(officer)
library(flextable)

# ------------------------------------------------------------------------------
# TASK A: LOAD RAW DATA & GENERATE STROBE ATTRITION COUNTS
# ------------------------------------------------------------------------------
cat("\n--- GENERATING STROBE ATTRITION COUNTS ---\n")

ch <- read_sav("K:/BGD_2025_MICS7_v01_M/ch.sav")
hh <- read_sav("K:/BGD_2025_MICS7_v01_M/hh.sav")

cat("Initial Under-5 Children in Survey:", nrow(ch), "\n")

# Ensure merge keys are numeric 
ch <- ch %>% mutate(HH1 = as.numeric(HH1), HH2 = as.numeric(HH2))
hh <- hh %>% mutate(HH1 = as.numeric(HH1), HH2 = as.numeric(HH2))

# Merge datasets
cols_to_keep <- c(setdiff(names(hh), names(ch)), "HH1", "HH2")
hh_subset <- hh[, cols_to_keep]
df <- left_join(ch, hh_subset, by = c("HH1", "HH2"))

# Sequential Exclusions & Counts
df <- df %>% filter(CAGE >= 12 & CAGE <= 59)
cat("Eligible by Age (12-59m):", nrow(df), "\n")

df <- df %>% filter(C_labresult == 1 & !is.na(C_Lead) & C_Lead < 900 & C_Lead > 0)
cat("Valid BLL Measurement:", nrow(df), "\n")

df <- df %>% filter(!is.na(DA1) & DA1 %in% c(1, 2))
cat("Valid Diarrhea Recall:", nrow(df), "\n")

df <- df %>% filter(HAZFLAG == 0 & WHZFLAG == 0 & WAZFLAG == 0)
cat("Final Analytical Sample (Valid Anthropometry):", nrow(df), "\n")
cat("-> Use these exact numbers to draw your Flowchart!\n")

# Variable Derivations for Models
df <- df %>%
  mutate(
    log2_bll = log2(as.numeric(C_Lead)),
    diarrhea = ifelse(DA1 == 1, 1, 0),
    stunted_num = ifelse(HAZ2 < -2.0, 1, 0)
  )

# Convert categorical variables to factors
cat_vars <- c("HL4", "melevel", "helevel", "windex5", "HH6", "HH7", "HC6A", "WS1", "WS11", "HW1")
df[cat_vars] <- lapply(df[cat_vars], haven::as_factor)

# Recreate Survey Design for the Final Sample
covars_m1 <- "CAGE + HL4 + C_Hb + melevel + helevel + windex5 + HH6 + HH7 + HC6A + WS1 + WS11 + HW1"
options(survey.lonely.psu = "adjust")
svy_primary <- svydesign(id = ~HH1, strata = ~HH6, weights = ~chblweight, data = df)

# ------------------------------------------------------------------------------
# TASK B: HEMOGLOBIN SENSITIVITY (Model 1 without C_Hb)
# ------------------------------------------------------------------------------
cat("\n--- SENSITIVITY: EXCLUDING HEMOGLOBIN ---\n")
covars_no_hb <- "CAGE + HL4 + melevel + helevel + windex5 + HH6 + HH7 + HC6A + WS1 + WS11 + HW1"
fit_haz_no_hb <- svyglm(as.formula(paste("HAZ2 ~ log2_bll * diarrhea +", covars_no_hb)), design = svy_primary)

cat("HAZ Beta (Lead) WITHOUT Hemoglobin Adjustment:\n")
print(summary(fit_haz_no_hb)$coefficients["log2_bll", ])

# ------------------------------------------------------------------------------
# TASK C: 5 µg/dL CATEGORICAL SENSITIVITY (Table S5)
# ------------------------------------------------------------------------------
cat("\n--- SENSITIVITY: 5 ug/dL BLL CUTOFF ---\n")
df_sens <- svy_primary$variables
df_sens <- df_sens %>%
  mutate(
    bll_high_5 = ifelse(C_Lead >= 5.0, 1, 0),
    joint_exp_5 = case_when(
      bll_high_5 == 0 & diarrhea == 0 ~ "00_LowBLL_NoDiarrhea",
      bll_high_5 == 0 & diarrhea == 1 ~ "01_LowBLL_Diarrhea",
      bll_high_5 == 1 & diarrhea == 0 ~ "10_HighBLL_NoDiarrhea",
      bll_high_5 == 1 & diarrhea == 1 ~ "11_HighBLL_Diarrhea"
    ),
    joint_exp_5 = factor(joint_exp_5, levels = c("00_LowBLL_NoDiarrhea", "01_LowBLL_Diarrhea", 
                                                 "10_HighBLL_NoDiarrhea", "11_HighBLL_Diarrhea"))
  )

svy_sens <- svydesign(id = ~HH1, strata = ~HH6, weights = ~chblweight, data = df_sens)
fit_pois_5 <- svyglm(as.formula(paste("stunted_num ~ joint_exp_5 +", covars_m1)), 
                     family = quasipoisson(link="log"), design = svy_sens)

cat("Stunting aPR for High BLL (>=5) + Diarrhea vs Low BLL + No Diarrhea:\n")
apr_11 <- exp(summary(fit_pois_5)$coefficients["joint_exp_511_HighBLL_Diarrhea", "Estimate"])
cat("aPR =", round(apr_11, 2), "\n")

# ------------------------------------------------------------------------------
# TASK D: NON-LINEARITY / RESTRICTED CUBIC SPLINES (Updated Style + TIFF)
# ------------------------------------------------------------------------------
cat("\n--- GENERATING SPLINE PLOT (FIGURE S2) ---\n")

# Subset only the 6 variables we need and drop NAs
df_spline <- df_sens %>%
  select(HAZ2, log2_bll, CAGE, HL4, C_Hb, windex5) %>%
  drop_na()

dd <- datadist(df_spline)
options(datadist="dd")

# Fit restricted cubic spline with 3 knots
spline_fit <- ols(HAZ2 ~ rcs(log2_bll, 3) + CAGE + HL4 + C_Hb + windex5, data = df_spline)

cat("\nSpline Non-Linearity p-value (If 'Non-Linear' p > 0.05, it is perfectly linear):\n")
print(anova(spline_fit))

# Extract prediction data manually to avoid the ugly default 'rms' formatting
pred_data <- as.data.frame(Predict(spline_fit, log2_bll))

# Create a highly polished, publication-ready ggplot
p_spline <- ggplot(pred_data, aes(x = log2_bll, y = yhat)) +
  # Confidence interval ribbon (light blue, semi-transparent)
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#A6C8E0", alpha = 0.5) +
  # Main spline prediction line (dark blue, thick)
  geom_line(color = "#1E3A8A", linewidth = 1.2) +
  theme_classic(base_size = 14) +
  labs(
    x = "Blood Lead Level (log2 µg/dL)", 
    y = "Predicted Height-for-Age Z-score (HAZ)", 
    title = "Restricted Cubic Spline: BLL vs HAZ"
  ) +
  theme(
    # Center the title, make it bold
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    # Make axis titles and text bold and crisp
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(color = "black", size = 11, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  )

# 1. Save as PNG
png("K:/BGD_2025_MICS7_v01_M/SPSS_Outputs/FigureS2_Splines.png", width=2200, height=1600, res=300)
print(p_spline)
dev.off()

# 2. Save as TIFF (High resolution with Journal-standard LZW compression)
tiff("K:/BGD_2025_MICS7_v01_M/SPSS_Outputs/FigureS2_Splines.tiff", width=2200, height=1600, res=300, compression="lzw")
print(p_spline)
dev.off()

cat("-> Beautiful FigureS2_Splines saved as BOTH PNG and TIFF to SPSS_Outputs.\n")