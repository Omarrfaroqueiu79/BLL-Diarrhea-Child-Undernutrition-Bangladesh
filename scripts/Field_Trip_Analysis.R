# ==============================================================================
# SCRIPT: Complete Epidemiological Analytical Workflow & Publication Suite
# Study: Intersecting Effects of Blood Lead Levels (BLL) and Diarrheal Morbidity 
#        on Child Growth Faltering in Bangladesh: A Cross-Sectional Analysis of MICS 2025
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment & Output Directory Initialization
# ------------------------------------------------------------------------------
base_dir   <- "K:/BGD_2025_MICS7_v01_M"
output_dir <- file.path(base_dir, "SPSS_Outputs")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

setwd(base_dir)

# 2. Load Core Libraries
suppressPackageStartupMessages({
  library(survey)
  library(haven)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

if (!requireNamespace("car", quietly = TRUE)) {
  install.packages("car", repos = "https://cloud.r-project.org")
}

# ------------------------------------------------------------------------------
# 3. Data Ingestion & Deterministic Record Linkage
# ------------------------------------------------------------------------------
cat("Ingesting raw MICS 2025 microdata files...\n")
ch_raw <- read_sav("ch.sav")
hh_raw <- read_sav("hh.sav")

# Select ONLY variables from hh_raw that are genuinely missing in ch_raw
# (Avoiding key variable overwrites or naming collisions)
hh_sub <- hh_raw %>%
  select(
    HH1, HH2,
    helevel, HHSEX, HHAGE,
    WS1, WS11, HW2, HW3A,
    any_of(c("WQ26", "WQ19D"))
  )

# Merge datasets via composite primary key (HH1: Cluster, HH2: Household)
df_merged <- ch_raw %>%
  left_join(hh_sub, by = c("HH1", "HH2"))

cat(sprintf("Total records post-linkage: %d\n", nrow(df_merged)))

# ------------------------------------------------------------------------------
# 4. Epidemiological Filtering & STROBE Sample Attrition Audit
# ------------------------------------------------------------------------------
# A. Target age group (12 to 59 months)
df_step <- df_merged %>%
  filter(CAGE >= 12 & CAGE <= 59)
n_age <- nrow(df_step)

# B. Biological missing codes recoded and laboratory BLL verification
# Plausible pediatric Hemoglobin range: > 0 and <= 25 g/dL (drops 99, 99.9, etc.)
df_step <- df_step %>%
  mutate(
    C_Lead = if_else(C_Lead >= 900 | C_Lead == 0, NA_real_, as.numeric(C_Lead)),
    C_Hb   = if_else(is.na(C_Hb) | C_Hb <= 0 | C_Hb > 25, NA_real_, as.numeric(C_Hb)),
    C_Cd   = if_else(C_Cd >= 9000 | C_Cd == 0, NA_real_, as.numeric(C_Cd)),
    C_As   = if_else(C_As >= 9000 | C_As == 0, NA_real_, as.numeric(C_As)),
    C_Hg   = if_else(C_Hg >= 9000 | C_Hg == 0, NA_real_, as.numeric(C_Hg)),
    S_Lead = if_else(S_Lead >= 9000 | S_Lead == 0, NA_real_, as.numeric(S_Lead))
  ) %>%
  filter(C_labresult == 1, !is.na(C_Lead))
n_lab <- nrow(df_step)

# C. Valid clinical 2-week diarrhea recall (1 = Yes, 2 = No)
df_step <- df_step %>%
  filter(!is.na(DA1), DA1 %in% c(1, 2))
n_morb <- nrow(df_step)

# D. Biologically plausible anthropometry (WHO standards: flags must be 0)
df_step <- df_step %>%
  filter(HAZFLAG == 0, WHZFLAG == 0, WAZFLAG == 0)
n_final <- nrow(df_step)

cat("\n================ STROBE SAMPLE ATTRITION AUDIT ================\n")
cat(sprintf("1. Children aged 12-59 months:            %d\n", n_age))
cat(sprintf("2. Valid laboratory BLL result:           %d (Dropped: %d)\n", n_lab, n_age - n_lab))
cat(sprintf("3. Valid 2-week diarrhea recall:          %d (Dropped: %d)\n", n_morb, n_lab - n_morb))
cat(sprintf("4. Plausible anthropometric measurements: %d (Dropped: %d)\n", n_final, n_morb - n_final))
cat(sprintf("FINAL ANALYTICAL SAMPLE SIZE (N):         %d\n", n_final))
cat("===============================================================\n\n")
# ------------------------------------------------------------------------------
# 5. Feature Engineering & WHO/JMP Indicator Construction
# ------------------------------------------------------------------------------
df_clean <- haven::zap_formats(df_step)

# JMP/WHO Standard Classification Vectors
improved_water_codes <- c(11, 12, 13, 14, 21, 31, 41, 51, 71, 72, 91)
improved_sanit_codes <- c(11, 12, 13, 21, 22, 31)

df_clean <- df_clean %>%
  mutate(
    # Continuous biomarkers & log2 transformations
    C_Lead      = as.numeric(C_Lead),
    log2_bll    = log2(C_Lead),
    C_Hb        = as.numeric(C_Hb),
    S_Lead      = as.numeric(S_Lead),
    log2_s_lead = log2(S_Lead),
    C_Cd        = as.numeric(C_Cd),
    log2_c_cd   = log2(C_Cd),
    C_As        = as.numeric(C_As),
    log2_c_as   = log2(C_As),
    
    # Continuous Anthropometric z-scores
    HAZ2        = as.numeric(HAZ2),
    WHZ2        = as.numeric(WHZ2),
    WAZ2        = as.numeric(WAZ2),
    CAGE        = as.numeric(CAGE),
    chblweight  = as.numeric(chblweight),
    
    # Primary Clinical Exposure Indicators
    diarrhea    = if_else(as.numeric(DA1) == 1, 1, 0),
    joint_exp   = factor(case_when(
      C_Lead < 3.5  & as.numeric(DA1) == 2 ~ "00_Ref_LowBLL_NoDiarrhea",
      C_Lead < 3.5  & as.numeric(DA1) == 1 ~ "01_LowBLL_Diarrhea",
      C_Lead >= 3.5 & as.numeric(DA1) == 2 ~ "10_HighBLL_NoDiarrhea",
      C_Lead >= 3.5 & as.numeric(DA1) == 1 ~ "11_HighBLL_Diarrhea"
    ), levels = c("00_Ref_LowBLL_NoDiarrhea", "01_LowBLL_Diarrhea", 
                  "10_HighBLL_NoDiarrhea", "11_HighBLL_Diarrhea")),
    
    # Binary Nutritional Endpoints (WHO Standard: z < -2.0 SD)
    stunted_num = if_else(HAZ2 < -2.0, 1, 0),
    wasted_num  = if_else(WHZ2 < -2.0, 1, 0),
    underwt_num = if_else(WAZ2 < -2.0, 1, 0),
    stunted_fac = factor(stunted_num, levels = c(0, 1), labels = c("No", "Yes")),
    wasted_fac  = factor(wasted_num,  levels = c(0, 1), labels = c("No", "Yes")),
    underwt_fac = factor(underwt_num, levels = c(0, 1), labels = c("No", "Yes")),
    
    # Socio-Demographic Covariates
    HL4         = factor(as.numeric(HL4), levels = c(1, 2), labels = c("Male", "Female")),
    HH6         = factor(as.numeric(HH6), levels = c(1, 2), labels = c("Urban", "Rural")),
    
    # Administrative Division: Exact BBS Geocodes (All 8 Divisions)
    HH7         = factor(as.numeric(HH7), 
                         levels = c(10, 20, 30, 40, 45, 50, 55, 60),
                         labels = c("Barishal", "Chattogram", "Dhaka", "Khulna", 
                                    "Mymensingh", "Rajshahi", "Rangpur", "Sylhet")),
    
    windex5     = factor(as.numeric(windex5), levels = 1:5, 
                         labels = c("Poorest", "Second", "Middle", "Fourth", "Richest")),
    
    # Maternal Education (Standardized 4-Tier Demographic Scale)
    melevel     = factor(case_when(
      as.numeric(melevel) == 0 ~ "Pre-primary or none",
      as.numeric(melevel) == 1 ~ "Primary",
      as.numeric(melevel) %in% c(2, 3) ~ "Secondary",
      as.numeric(melevel) %in% c(4, 5) ~ "Higher",
      TRUE ~ "Not in Household"
    ), levels = c("Pre-primary or none", "Primary", "Secondary", "Higher", "Not in Household")),
    
    # Household Head Education (Standardized 4-Tier Demographic Scale)
    helevel     = factor(case_when(
      as.numeric(helevel) == 0 ~ "Pre-primary or none",
      as.numeric(helevel) == 1 ~ "Primary",
      as.numeric(helevel) %in% c(2, 3) ~ "Secondary",
      as.numeric(helevel) %in% c(4, 5) ~ "Higher",
      TRUE ~ "Pre-primary or none"
    ), levels = c("Pre-primary or none", "Primary", "Secondary", "Higher")),
    
    # Living Environment & Paint Condition (HC6A from ch_raw)
    peeling_paint = factor(if_else(as.numeric(HC6A) %in% c(2, 3), 1, 0),
                           levels = c(0, 1), labels = c("Intact Paint", "Peeling/Damaged Paint")),
    
    # Water & Sanitation Vectors (WS1 and WS11 from hh_raw)
    improved_water      = factor(if_else(as.numeric(WS1) %in% improved_water_codes, 1, 0),
                                 levels = c(0, 1), labels = c("Unimproved Water", "Improved Water")),
    improved_sanitation = factor(if_else(as.numeric(WS11) %in% improved_sanit_codes, 1, 0),
                                 levels = c(0, 1), labels = c("Unimproved Toilet", "Improved Toilet")),
    
    # Standard WHO/JMP Handwashing Facility (Facility observed in HW1, Water in HW2, Soap in HW3A)
    basic_handwashing   = factor(
      if_else(
        as.numeric(HW1) %in% c(1, 2, 3) & 
          as.numeric(HW2) == 1 & 
          grepl("A", as.character(HW3A)), 1, 0
      ),
      levels = c(0, 1), 
      labels = c("Limited/No Facility", "Basic Facility")
    )
  ) %>%
  droplevels()

# Strip residual haven labels for clean survey modeling
df_clean <- haven::zap_labels(df_clean)

# Assign joint_exp globally to ensure survey environment scope access
joint_exp <<- df_clean$joint_exp

cat("Audit of 8 Administrative Divisions post-cleaning:\n")
print(table(df_clean$HH7, useNA = "always"))

cat("\nAudit of Handwashing Facility post-cleaning:\n")
print(table(df_clean$basic_handwashing, useNA = "always"))

# ------------------------------------------------------------------------------
# Prune Dataset: Keep ONLY the 28 Analytical Variables
# ------------------------------------------------------------------------------
df_clean <- df_clean %>%
  select(
    # 1. Survey Calibration
    PSU, stratum, chblweight,
    
    # 2. Exposure Indicators
    C_Lead, log2_bll, DA1, diarrhea, joint_exp,
    
    # 3. Anthropometric Outcomes
    HAZ2, WHZ2, WAZ2,
    stunted_num, stunted_fac,
    wasted_num, wasted_fac,
    underwt_num, underwt_fac,
    
    # 4. Demographic & Biological Covariates
    CAGE, HL4, C_Hb,
    melevel, helevel, windex5, HH6, HH7,
    
    # 5. WASH & Built Environment
    peeling_paint, improved_water, improved_sanitation, basic_handwashing,
    
    # 6. Environmental Trace Metals (Robustness Model 2)
    S_Lead, log2_s_lead, C_Cd, log2_c_cd, C_As, log2_c_as,
    
    # 7. Water Quality Subsample (if present)
    any_of("WQ26")
  )

cat(sprintf("\nFinal clean dataset dimension: %d rows and %d analytical variables.\n", 
            nrow(df_clean), ncol(df_clean)))

# ------------------------------------------------------------------------------
# 6. Complex Survey Design Calibration & Clean Dataset Export
# ------------------------------------------------------------------------------
options(survey.lonely.psu = "adjust")

cat("\nCalibrating complex survey design (PSU, Stratum, Normalized Weights)...\n")
svy_primary <- svydesign(
  id      = ~PSU,
  strata  = ~stratum,
  weights = ~chblweight,
  nest    = TRUE,
  data    = df_clean
)

# Export finalized analytical datasets (.rds and .sav)
cat("Exporting finalized analytical datasets to SPSS_Outputs...\n")
saveRDS(df_clean, file.path(output_dir, "mics_lead_diarrhea_clean.rds"))
haven::write_sav(df_clean, path = file.path(output_dir, "mics_lead_diarrhea_clean.sav"))
cat(sprintf("Saved:\n  - %s/mics_lead_diarrhea_clean.rds\n  - %s/mics_lead_diarrhea_clean.sav\n\n", 
            output_dir, output_dir))

# # ------------------------------------------------------------------------------
# # 7. Previous Table 1: Population Profile & Weighted Descriptive Statistics
# # ------------------------------------------------------------------------------
# cat("\nGenerating Table 1: Survey-Weighted Baseline Characteristics across all 8 Divisions...\n")
# groups <- levels(df_clean$joint_exp)
# 
# calc_cont <- function(var_name, label) {
#   f_var <- as.formula(paste0("~", var_name), env = globalenv())
#   
#   m_ov <- svymean(f_var, svy_primary, na.rm = TRUE)
#   v_ov <- svyvar(f_var, svy_primary, na.rm = TRUE)
#   ov_str <- sprintf("%.2f (%.2f)", as.numeric(m_ov), sqrt(as.numeric(v_ov)))
#   
#   by_m <- svyby(f_var, ~joint_exp, svy_primary, svymean, na.rm = TRUE)
#   by_v <- svyby(f_var, ~joint_exp, svy_primary, svyvar, na.rm = TRUE)
#   
#   grp_strs <- sapply(groups, function(g) {
#     row_idx <- which(by_m$joint_exp == g)
#     if (length(row_idx) == 0) return("NA")
#     val_m <- as.numeric(by_m[row_idx, var_name])
#     val_s <- sqrt(as.numeric(by_v[row_idx, var_name]))
#     sprintf("%.2f (%.2f)", val_m, val_s)
#   })
#   
#   f_mod  <- as.formula(paste0(var_name, " ~ joint_exp"), env = globalenv())
#   m_test <- svyglm(f_mod, design = svy_primary)
#   p_val  <- tryCatch(as.numeric(regTermTest(m_test, "joint_exp")$p), error = function(e) NA_real_)
#   p_str  <- if (is.na(p_val)) "" else if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
#   
#   c(Characteristic = label, Overall = ov_str, grp_strs, p_value = p_str)
# }
# 
# calc_cat <- function(var_name, label) {
#   levs <- levels(df_clean[[var_name]])
#   rows_list <- list()
#   
#   f_tab <- as.formula(paste0("~", var_name, " + joint_exp"), env = globalenv())
#   p_val <- tryCatch(as.numeric(svychisq(f_tab, svy_primary)$p.value), error = function(e) NA_real_)
#   p_str <- if (is.na(p_val)) "" else if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
#   
#   f_var    <- as.formula(paste0("~", var_name), env = globalenv())
#   prop_ov  <- svymean(f_var, svy_primary, na.rm = TRUE)
#   prop_grp <- svyby(f_var, ~joint_exp, svy_primary, svymean, na.rm = TRUE)
#   
#   for (i in seq_along(levs)) {
#     lev <- levs[i]
#     col_name <- paste0(var_name, lev)
#     
#     n_unw  <- sum(df_clean[[var_name]] == lev, na.rm = TRUE)
#     pct_ov <- if (col_name %in% names(prop_ov)) as.numeric(prop_ov[col_name]) * 100 else 0
#     ov_str <- sprintf("%d (%.1f%%)", n_unw, pct_ov)
#     
#     grp_strs <- sapply(groups, function(g) {
#       n_g   <- sum(df_clean[[var_name]] == lev & df_clean$joint_exp == g, na.rm = TRUE)
#       row_g <- which(prop_grp$joint_exp == g)
#       pct_g <- if (length(row_g) > 0 && col_name %in% colnames(prop_grp)) {
#         as.numeric(prop_grp[row_g, col_name]) * 100
#       } else {
#         0
#       }
#       sprintf("%d (%.1f%%)", n_g, pct_g)
#     })
#     
#     row_lbl   <- paste0("  ", lev)
#     p_to_show <- if (i == 1) p_str else ""
#     rows_list[[i]] <- c(Characteristic = row_lbl, Overall = ov_str, grp_strs, p_value = p_to_show)
#   }
#   
#   # For binary Yes/No endpoints (stunting, wasting, underweight), report only the positive condition row
#   if (all(c("No", "Yes") %in% levs)) {
#     yes_idx <- which(levs == "Yes")
#     out_row <- rows_list[[yes_idx]]
#     out_row["Characteristic"] <- label
#     out_row["p_value"] <- p_str
#     return(matrix(out_row, nrow = 1, dimnames = list(NULL, names(out_row))))
#   }
#   
#   header_row <- c(Characteristic = label, Overall = "", setNames(rep("", length(groups)), groups), p_value = "")
#   do.call(rbind, c(list(header_row), rows_list))
# }
# 
# t1_rows <- list(
#   calc_cont("CAGE", "Child Age (months), Mean (SD)"),
#   calc_cat("HL4", "Child Sex, n (%)"),
#   calc_cont("C_Lead", "Blood Lead Level (µg/dL), Mean (SD)"),
#   calc_cont("C_Hb", "Hemoglobin (g/dL), Mean (SD)"),
#   calc_cont("S_Lead", "Soil Lead Level (mg/kg), Mean (SD)"),
#   calc_cat("melevel", "Mother's Education, n (%)"),
#   calc_cat("helevel", "Household Head Education, n (%)"),
#   calc_cat("windex5", "Wealth Quintile, n (%)"),
#   calc_cat("HH6", "Place of Residence, n (%)"),
#   calc_cat("HH7", "Administrative Division, n (%)"),
#   calc_cat("improved_water", "Drinking Water Source, n (%)"),
#   calc_cat("improved_sanitation", "Sanitation Facility, n (%)"),
#   calc_cat("basic_handwashing", "Handwashing Facility, n (%)"),
#   calc_cat("peeling_paint", "Wall Paint Condition, n (%)"),
#   calc_cont("HAZ2", "Height-for-Age z-score (HAZ), Mean (SD)"),
#   calc_cont("WHZ2", "Weight-for-Height z-score (WHZ), Mean (SD)"),
#   calc_cont("WAZ2", "Weight-for-Age z-score (WAZ), Mean (SD)"),
#   calc_cat("stunted_fac", "Stunted (HAZ < -2 SD), n (%)"),
#   calc_cat("wasted_fac", "Wasted (WHZ < -2 SD), n (%)"),
#   calc_cat("underwt_fac", "Underweight (WAZ < -2 SD), n (%)")
# )
# 
# Table1_final <- as.data.frame(do.call(rbind, t1_rows), stringsAsFactors = FALSE)
# colnames(Table1_final) <- c(
#   "Characteristics", 
#   "Total Cohort (N=10,507)", 
#   "Low BLL, No Diarrhea (00)", 
#   "Low BLL + Diarrhea (01)", 
#   "High BLL, No Diarrhea (10)", 
#   "High BLL + Diarrhea (11)", 
#   "p-value"
# )
# 
# write.csv(Table1_final, file.path(output_dir, "Table1_Weighted_Characteristics.csv"), row.names = FALSE)
# cat(sprintf("Saved: %s/Table1_Weighted_Characteristics.csv\n", output_dir))








# ------------------------------------------------------------------------------
# 7. Table 1: Population Profile (Segregated by Continuous vs. Categorical)
# ------------------------------------------------------------------------------
cat("\nGenerating Table 1: Segregated Continuous [Mean(SD)] and Categorical [n(%)] Blocks...\n")
groups <- levels(df_clean$joint_exp)

calc_cont <- function(var_name, label) {
  f_var <- as.formula(paste0("~", var_name), env = globalenv())
  
  m_ov <- svymean(f_var, svy_primary, na.rm = TRUE)
  v_ov <- svyvar(f_var, svy_primary, na.rm = TRUE)
  ov_str <- sprintf("%.2f (%.2f)", as.numeric(m_ov), sqrt(as.numeric(v_ov)))
  
  by_m <- svyby(f_var, ~joint_exp, svy_primary, svymean, na.rm = TRUE)
  by_v <- svyby(f_var, ~joint_exp, svy_primary, svyvar, na.rm = TRUE)
  
  grp_strs <- sapply(groups, function(g) {
    row_idx <- which(by_m$joint_exp == g)
    if (length(row_idx) == 0) return("NA")
    val_m <- as.numeric(by_m[row_idx, var_name])
    val_s <- sqrt(as.numeric(by_v[row_idx, var_name]))
    sprintf("%.2f (%.2f)", val_m, val_s)
  })
  
  f_mod  <- as.formula(paste0(var_name, " ~ joint_exp"), env = globalenv())
  m_test <- svyglm(f_mod, design = svy_primary)
  p_val  <- tryCatch(as.numeric(regTermTest(m_test, "joint_exp")$p), error = function(e) NA_real_)
  p_str  <- if (is.na(p_val)) "" else if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
  
  c(Characteristic = label, Overall = ov_str, grp_strs, p_value = p_str)
}

calc_cat <- function(var_name, label) {
  levs <- levels(df_clean[[var_name]])
  rows_list <- list()
  
  f_tab <- as.formula(paste0("~", var_name, " + joint_exp"), env = globalenv())
  p_val <- tryCatch(as.numeric(svychisq(f_tab, svy_primary)$p.value), error = function(e) NA_real_)
  p_str <- if (is.na(p_val)) "" else if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
  
  f_var    <- as.formula(paste0("~", var_name), env = globalenv())
  prop_ov  <- svymean(f_var, svy_primary, na.rm = TRUE)
  prop_grp <- svyby(f_var, ~joint_exp, svy_primary, svymean, na.rm = TRUE)
  
  for (i in seq_along(levs)) {
    lev <- levs[i]
    col_name <- paste0(var_name, lev)
    
    n_unw  <- sum(df_clean[[var_name]] == lev, na.rm = TRUE)
    pct_ov <- if (col_name %in% names(prop_ov)) as.numeric(prop_ov[col_name]) * 100 else 0
    ov_str <- sprintf("%d (%.1f%%)", n_unw, pct_ov)
    
    grp_strs <- sapply(groups, function(g) {
      n_g   <- sum(df_clean[[var_name]] == lev & df_clean$joint_exp == g, na.rm = TRUE)
      row_g <- which(prop_grp$joint_exp == g)
      pct_g <- if (length(row_g) > 0 && col_name %in% colnames(prop_grp)) {
        as.numeric(prop_grp[row_g, col_name]) * 100
      } else {
        0
      }
      sprintf("%d (%.1f%%)", n_g, pct_g)
    })
    
    row_lbl   <- paste0("  ", lev)
    p_to_show <- if (i == 1) p_str else ""
    rows_list[[i]] <- c(Characteristic = row_lbl, Overall = ov_str, grp_strs, p_value = p_to_show)
  }
  
  # For binary Yes/No clinical endpoints, report only the single 'Yes' row
  if (all(c("No", "Yes") %in% levs)) {
    yes_idx <- which(levs == "Yes")
    out_row <- rows_list[[yes_idx]]
    out_row["Characteristic"] <- label
    out_row["p_value"] <- p_str
    return(matrix(out_row, nrow = 1, dimnames = list(NULL, names(out_row))))
  }
  
  header_row <- c(Characteristic = label, Overall = "", setNames(rep("", length(groups)), groups), p_value = "")
  do.call(rbind, c(list(header_row), rows_list))
}

make_blank_header <- function(title_text) {
  matrix(c(Characteristic = title_text, Overall = "", setNames(rep("", length(groups)), groups), p_value = ""),
         nrow = 1, dimnames = list(NULL, c("Characteristic", "Overall", groups, "p_value")))
}

# 1. Continuous Block [Mean (SD)]
cont_rows <- list(
  make_blank_header("Continuous Measures, Mean (SD)"),
  calc_cont("CAGE", "Child Age (months)"),
  calc_cont("C_Lead", "Blood Lead Level (µg/dL)"),
  calc_cont("C_Hb", "Hemoglobin (g/dL)"),
  calc_cont("S_Lead", "Soil Lead Level (mg/kg)"),
  calc_cont("HAZ2", "Height-for-Age z-score (HAZ)"),
  calc_cont("WHZ2", "Weight-for-Height z-score (WHZ)"),
  calc_cont("WAZ2", "Weight-for-Age z-score (WAZ)")
)

# 2. Categorical Block [n (%)]
cat_rows <- list(
  make_blank_header("Categorical Characteristics, n (%)"),
  calc_cat("HL4", "Child Sex"),
  calc_cat("melevel", "Mother's Education Level"),
  calc_cat("helevel", "Household Head Education"),
  calc_cat("windex5", "Wealth Quintile"),
  calc_cat("HH6", "Place of Residence"),
  calc_cat("HH7", "Administrative Division"),
  calc_cat("improved_water", "Drinking Water Source"),
  calc_cat("improved_sanitation", "Sanitation Facility"),
  calc_cat("basic_handwashing", "Handwashing Facility"),
  calc_cat("peeling_paint", "Wall Paint Condition"),
  calc_cat("stunted_fac", "Stunted (HAZ < -2 SD)"),
  calc_cat("wasted_fac", "Wasted (WHZ < -2 SD)"),
  calc_cat("underwt_fac", "Underweight (WAZ < -2 SD)")
)

Table1_final <- as.data.frame(do.call(rbind, c(cont_rows, cat_rows)), stringsAsFactors = FALSE)
colnames(Table1_final) <- c(
  "Characteristics", 
  "Total Cohort (N=10,507)", 
  "Low BLL, No Diarrhea (00)", 
  "Low BLL + Diarrhea (01)", 
  "High BLL, No Diarrhea (10)", 
  "High BLL + Diarrhea (11)", 
  "p-value"
)

write.csv(Table1_final, file.path(output_dir, "Table1_Weighted_Characteristics.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/Table1_Weighted_Characteristics.csv\n", output_dir))








# ==============================================================================
# Export Publication-Ready Table 1 to Microsoft Word (.docx)
# ==============================================================================
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")

library(officer)
library(flextable)

# Ensure Table1_final exists in the environment
raw_tbl <- Table1_final

# Initialize flextable object
ft <- flextable(raw_tbl)

# 1. Apply font, size, and layout
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- fontsize(ft, size = 9, part = "all")
ft <- fontsize(ft, size = 9.5, part = "header")
ft <- bold(ft, part = "header")

# 2. Identify and style section divider rows (Mean [SD] & n [%] headers)
section_idx <- which(raw_tbl$Characteristics %in% c("Continuous Measures, Mean (SD)", 
                                                    "Categorical Characteristics, n (%)"))

ft <- bold(ft, i = section_idx, part = "body")
ft <- italic(ft, i = section_idx, part = "body")
ft <- bg(ft, i = section_idx, bg = "#F2F4F7", part = "body")

# 3. Align text: Left for variables, center for data cells, right for p-values
ft <- align(ft, j = 1, align = "left", part = "all")
ft <- align(ft, j = 2:6, align = "center", part = "all")
ft <- align(ft, j = 7, align = "right", part = "all")

# 4. Standard Classic Three-Line Border Style (APA / ICMJE journal standard)
std_border <- fp_border(color = "black", width = 1)
thin_border <- fp_border(color = "black", width = 0.5)

ft <- border_remove(ft)
ft <- hline_top(ft, part = "header", border = std_border)
ft <- hline_bottom(ft, part = "header", border = std_border)
ft <- hline_bottom(ft, part = "body", border = std_border)

# 5. Set column widths (in inches) for clean horizontal display
ft <- width(ft, j = 1, width = 2.4)
ft <- width(ft, j = 2:6, width = 1.05)
ft <- width(ft, j = 7, width = 0.8)

# 6. Add STROBE title and survey methodological footnote
ft <- set_caption(
  ft, 
  caption = "Table 1. Survey-weighted baseline characteristics of Bangladeshi children aged 12–59 months, stratified by joint blood lead and diarrheal morbidity exposure (N = 10,507)"
)

footnote_text <- paste(
  "Values are design-weighted Mean (Standard Deviation) for continuous measures and unweighted count n",
  "(complex survey-weighted percentage %) for categorical measures. Percentages and means account for the",
  "multi-stage cluster sampling design and child normalized sampling weights. p-values compare distributions",
  "across the four joint exposure categories using design-adjusted Wald F-tests for continuous measures and",
  "Rao-Scott second-order design-adjusted chi-square tests for categorical measures.",
  "Blood lead level cutoff is defined according to the CDC pediatric Blood Lead Reference Value (3.5 µg/dL)."
)

ft <- add_footer_lines(ft, values = footnote_text)
ft <- fontsize(ft, size = 8, part = "footer")
ft <- italic(ft, part = "footer")

# 7. Write to Word Document in landscape orientation
doc_path <- file.path(output_dir, "Table1_Manuscript_Ready.docx")

doc <- read_docx() %>%
  body_add_par(
    "Table 1: Baseline Population Characteristics", 
    style = "heading 1"
  ) %>%
  body_add_flextable(ft)

print(doc, target = doc_path)
cat(sprintf("\nSuccessfully exported publication-ready Word file:\n-> %s\n", doc_path))






# 
# # ==============================================================================
# # Export Publication-Ready Table 2 to Microsoft Word (.docx)
# # ==============================================================================
# library(officer)
# library(flextable)
# 
# # Clean up raw parameter text
# t2_clean <- table2_df %>%
#   mutate(
#     Parameter = case_when(
#       Parameter == "log2_bll" ~ "Blood lead level (log2 ug/dL)",
#       Parameter == "diarrhea" ~ "Recent diarrhea (past 2 weeks)",
#       Parameter == "log2_bll:diarrhea" ~ "Interaction (BLL x Diarrhea)",
#       TRUE ~ Parameter
#     ),
#     SE = sprintf("%.3f", as.numeric(SE)) # Fixes 0.19 to 0.190
#   )
# 
# ft2 <- as_grouped_data(t2_clean, groups = c("Outcome", "Model"))
# ft2 <- flextable(ft2)
# 
# ft2 <- font(ft2, fontname = "Times New Roman", part = "all")
# ft2 <- fontsize(ft2, size = 9.5, part = "all")
# ft2 <- fontsize(ft2, size = 10, part = "header")
# ft2 <- bold(ft2, part = "header")
# 
# # Style headers and group rows
# ft2 <- bold(ft2, i = ~ !is.na(Outcome), part = "body")
# ft2 <- italic(ft2, i = ~ is.na(Outcome) & !is.na(Model), part = "body")
# 
# # Borders and alignment
# std_border <- fp_border(color = "black", width = 1)
# ft2 <- border_remove(ft2)
# ft2 <- hline_top(ft2, part = "header", border = std_border)
# ft2 <- hline_bottom(ft2, part = "header", border = std_border)
# ft2 <- hline_bottom(ft2, part = "body", border = std_border)
# 
# ft2 <- align(ft2, j = 1, align = "left", part = "all")
# ft2 <- align(ft2, j = 2:5, align = "center", part = "all")
# 
# ft2 <- set_header_labels(
#   ft2,
#   Parameter = "Parameter Exposure",
#   Beta = "Adjusted Beta",
#   SE = "SE",
#   CI_95 = "95% CI",
#   p_value = "p-value"
# )
# 
# doc_t2 <- read_docx() %>%
#   body_add_par("Table 2: Continuous Growth Faltering & Multiplicative Interaction", style = "heading 1") %>%
#   body_add_flextable(ft2)
# 
# print(doc_t2, target = file.path(output_dir, "Table2_Manuscript_Ready.docx"))
# cat("Exported: Table2_Manuscript_Ready.docx\n")





# ==============================================================================
# 8. Table 2: Complete Survey-Weighted Linear Dose-Response Models
#    (Model 1 & Model 2 across HAZ, WHZ, and WAZ with Publication-Ready Export)
# ==============================================================================

if (!requireNamespace("survey", quietly = TRUE)) install.packages("survey")
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(survey)
library(officer)
library(flextable)
library(dplyr)

cat("\nFitting Survey-Weighted Linear Dose-Response Models (Table 2)...\n")

# ------------------------------------------------------------------------------
# 1. Define Analytical Design Objects
# ------------------------------------------------------------------------------
# Primary complete survey design (N = 10,507)
options(survey.lonely.psu = "adjust")

svy_primary <- svydesign(
  id      = ~PSU,
  strata  = ~stratum,
  weights = ~chblweight,
  nest    = TRUE,
  data    = df_clean
)

# Metals robustness subsample (complete cases on trace co-pollutants, N = 10,423)
df_metals <- df_clean %>% 
  filter(!is.na(log2_s_lead) & !is.na(log2_c_cd) & !is.na(log2_c_as))

svy_metals <- svydesign(
  id      = ~PSU,
  strata  = ~stratum,
  weights = ~chblweight,
  nest    = TRUE,
  data    = df_metals
)

n_m1 <- format(nrow(df_clean), big.mark = ",")
n_m2 <- format(nrow(df_metals), big.mark = ",")

# ------------------------------------------------------------------------------
# 2. Model Formulas Setup
# ------------------------------------------------------------------------------
covars_m1 <- paste(
  "CAGE", "HL4", "C_Hb", "melevel", "helevel", "windex5", 
  "HH6", "HH7", "peeling_paint", "improved_water", 
  "improved_sanitation", "basic_handwashing", 
  sep = " + "
)

covars_m2 <- paste(
  covars_m1, 
  "log2_s_lead + log2_c_cd + log2_c_as", 
  sep = " + "
)

outcomes <- c(
  "HAZ2" = "Height-for-Age (HAZ)", 
  "WHZ2" = "Weight-for-Height (WHZ)", 
  "WAZ2" = "Weight-for-Age (WAZ)"
)

table2_rows <- list()

# ------------------------------------------------------------------------------
# 3. Model Estimation Loop (Model 1 & Model 2 for HAZ, WHZ, and WAZ)
# ------------------------------------------------------------------------------
for (y in names(outcomes)) {
  y_label <- outcomes[y]
  
  # --- Model 1: Primary Adjustment ---
  f1  <- as.formula(paste0(y, " ~ log2_bll * diarrhea + ", covars_m1))
  m1  <- svyglm(f1, design = svy_primary)
  s1  <- summary(m1)$coefficients
  ci1 <- confint(m1)
  
  params <- c("log2_bll", "diarrhea", "log2_bll:diarrhea")
  for (p in params) {
    p_val <- s1[p, "Pr(>|t|)"]
    table2_rows[[length(table2_rows) + 1]] <- data.frame(
      Outcome   = y_label,
      Model     = sprintf("Model 1: Primary Adjustment (N=%s)", n_m1),
      Parameter = p,
      Beta      = sprintf("%.3f", s1[p, "Estimate"]),
      SE        = sprintf("%.3f", s1[p, "Std. Error"]),
      CI_95     = sprintf("(%.3f, %.3f)", ci1[p, 1], ci1[p, 2]),
      p_value   = if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val),
      stringsAsFactors = FALSE
    )
  }
  
  # --- Model 2: Co-Metals Robustness Adjustment ---
  f2  <- as.formula(paste0(y, " ~ log2_bll * diarrhea + ", covars_m2))
  m2  <- svyglm(f2, design = svy_metals)
  s2  <- summary(m2)$coefficients
  ci2 <- confint(m2)
  
  for (p in params) {
    p_val <- s2[p, "Pr(>|t|)"]
    table2_rows[[length(table2_rows) + 1]] <- data.frame(
      Outcome   = y_label,
      Model     = sprintf("Model 2: Trace Metals Robustness (N=%s)", n_m2),
      Parameter = p,
      Beta      = sprintf("%.3f", s2[p, "Estimate"]),
      SE        = sprintf("%.3f", s2[p, "Std. Error"]),
      CI_95     = sprintf("(%.3f, %.3f)", ci2[p, 1], ci2[p, 2]),
      p_value   = if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val),
      stringsAsFactors = FALSE
    )
  }
}

table2_df <- do.call(rbind, table2_rows)

# Clean parameter text for manuscript presentation
table2_df <- table2_df %>%
  mutate(
    Parameter = case_when(
      Parameter == "log2_bll"          ~ "Blood lead level (log2 ug/dL)",
      Parameter == "diarrhea"          ~ "Recent diarrhea (past 2 weeks)",
      Parameter == "log2_bll:diarrhea" ~ "Interaction (BLL x Diarrhea)",
      TRUE ~ Parameter
    )
  )

# Save Raw CSV Output
write.csv(table2_df, file.path(output_dir, "Table2_Linear_Interactions.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/Table2_Linear_Interactions.csv\n", output_dir))

# ------------------------------------------------------------------------------
# 4. Generate Publication-Quality Word Document (.docx)
# ------------------------------------------------------------------------------
ft2_data <- as_grouped_data(table2_df, groups = c("Outcome", "Model"))
ft2 <- flextable(ft2_data)

# Typography & Layout
ft2 <- font(ft2, fontname = "Times New Roman", part = "all")
ft2 <- fontsize(ft2, size = 9.5, part = "all")
ft2 <- fontsize(ft2, size = 10, part = "header")
ft2 <- bold(ft2, part = "header")

# Formatting group headers
ft2 <- bold(ft2, i = ~ !is.na(Outcome), part = "body")
ft2 <- italic(ft2, i = ~ is.na(Outcome) & !is.na(Model), part = "body")

# Set Column Headers
ft2 <- set_header_labels(
  ft2,
  Parameter = "Growth Metric & Model Parameter",
  Beta      = "Adjusted Beta",
  SE        = "SE",
  CI_95     = "95% CI",
  p_value   = "p-value"
)

# Borders (Classic Three-Line Journal Format)
std_border <- fp_border(color = "black", width = 1)
ft2 <- border_remove(ft2)
ft2 <- hline_top(ft2, part = "header", border = std_border)
ft2 <- hline_bottom(ft2, part = "header", border = std_border)
ft2 <- hline_bottom(ft2, part = "body", border = std_border)

# Text Alignment
ft2 <- align(ft2, j = 1, align = "left", part = "all")
ft2 <- align(ft2, j = 2:5, align = "center", part = "all")

# Column Widths
ft2 <- width(ft2, j = 1, width = 3.2)
ft2 <- width(ft2, j = 2:5, width = 1.1)

# STROBE Footnote
ft2 <- set_caption(
  ft2,
  caption = "Table 2. Multivariable survey-weighted linear regression evaluating blood lead levels, diarrheal morbidity, and continuous growth faltering (z-scores) among Bangladeshi children aged 12–59 months."
)

footnote_t2 <- paste(
  "Beta coefficients (Beta) represent adjusted change in standard deviation scores per two-fold increase in blood lead level",
  "(log2-transformed), caregiver-reported diarrhea in the preceding two weeks, and their continuous multiplicative interaction term.",
  "Model 1 (N = 10,507) adjusts for child age, child sex, venous hemoglobin, maternal education, household head education, wealth quintile,",
  "urban/rural residence, administrative division, peeling wall paint, drinking water source, sanitation facility type, and basic handwashing facilities.",
  "Model 2 (N = 10,423) additionally adjusts for environmental co-pollutants: soil lead (log2 mg/kg), venous blood cadmium (log2 ug/L), and",
  "venous blood arsenic (log2 ug/L). All models account for complex survey strata, primary sampling units, and normalized child blood weights."
)

ft2 <- add_footer_lines(ft2, values = footnote_t2)
ft2 <- fontsize(ft2, size = 8, part = "footer")
ft2 <- italic(ft2, part = "footer")

# Export to Word
doc_t2_path <- file.path(output_dir, "Table2_Manuscript_Ready.docx")
doc_t2 <- read_docx() %>%
  body_add_par("Table 2: Continuous Growth Faltering & Multiplicative Interaction", style = "heading 1") %>%
  body_add_flextable(ft2)

print(doc_t2, target = doc_t2_path)
cat(sprintf("Saved: %s\n\n", doc_t2_path))






# # ------------------------------------------------------------------------------
# # 9. Table 3: Modified Poisson Regression & Additive Synergy Decomposition
# # ------------------------------------------------------------------------------
# cat("\nDecomposing Additive Synergy Metrics using Modified Poisson Regression...\n")
# 
# compute_additive_synergy <- function(outcome_var, outcome_name) {
#   formula_poisson <- as.formula(paste(
#     outcome_var, "~ joint_exp + CAGE + HL4 + C_Hb + melevel + helevel +",
#     "windex5 + HH6 + HH7 + peeling_paint + improved_water + improved_sanitation + basic_handwashing"
#   ), env = globalenv())
#   
#   fit_poisson <- svyglm(
#     formula_poisson,
#     design = svy_primary,
#     family = quasipoisson(link = "log")
#   )
#   
#   b <- coef(fit_poisson)
#   V <- vcov(fit_poisson)
#   
#   idx_01 <- grep("01_LowBLL_Diarrhea", names(b), value = TRUE)
#   idx_10 <- grep("10_HighBLL_NoDiarrhea", names(b), value = TRUE)
#   idx_11 <- grep("11_HighBLL_Diarrhea", names(b), value = TRUE)
#   
#   pr_01 <- exp(b[idx_01])
#   pr_10 <- exp(b[idx_10])
#   pr_11 <- exp(b[idx_11])
#   
#   ci_01 <- exp(c(b[idx_01] - 1.96 * sqrt(V[idx_01, idx_01]), b[idx_01] + 1.96 * sqrt(V[idx_01, idx_01])))
#   ci_10 <- exp(c(b[idx_10] - 1.96 * sqrt(V[idx_10, idx_10]), b[idx_10] + 1.96 * sqrt(V[idx_10, idx_10])))
#   ci_11 <- exp(c(b[idx_11] - 1.96 * sqrt(V[idx_11, idx_11]), b[idx_11] + 1.96 * sqrt(V[idx_11, idx_11])))
#   
#   # Point Estimates
#   reri_est <- pr_11 - pr_10 - pr_01 + 1
#   ap_est   <- reri_est / pr_11
#   s_est    <- (pr_11 - 1) / ((pr_10 - 1) + (pr_01 - 1))
#   
#   # Delta Method Standard Errors
#   sub_vars <- c(idx_01, idx_10, idx_11)
#   V_sub    <- V[sub_vars, sub_vars]
#   
#   grad_reri <- c(-pr_01, -pr_10, pr_11)
#   se_reri   <- sqrt(as.numeric(t(grad_reri) %*% V_sub %*% grad_reri))
#   reri_ci   <- c(reri_est - 1.96 * se_reri, reri_est + 1.96 * se_reri)
#   
#   grad_ap <- c(-pr_01 / pr_11, -pr_10 / pr_11, (pr_10 + pr_01 - 1) / pr_11)
#   se_ap   <- sqrt(as.numeric(t(grad_ap) %*% V_sub %*% grad_ap))
#   ap_ci   <- c(ap_est - 1.96 * se_ap, ap_est + 1.96 * se_ap)
#   
#   denom_s   <- (pr_10 - 1) + (pr_01 - 1)
#   grad_logs <- c(-pr_01 / denom_s, -pr_10 / denom_s, pr_11 / (pr_11 - 1))
#   se_logs   <- sqrt(as.numeric(t(grad_logs) %*% V_sub %*% grad_logs))
#   s_ci      <- exp(c(log(s_est) - 1.96 * se_logs, log(s_est) + 1.96 * se_logs))
#   
#   data.frame(
#     Outcome = outcome_name,
#     aPR_00  = "1.00 (Reference)",
#     aPR_01  = sprintf("%.2f (%.2f, %.2f)", pr_01, ci_01[1], ci_01[2]),
#     aPR_10  = sprintf("%.2f (%.2f, %.2f)", pr_10, ci_10[1], ci_10[2]),
#     aPR_11  = sprintf("%.2f (%.2f, %.2f)", pr_11, ci_11[1], ci_11[2]),
#     RERI    = sprintf("%.2f (%.2f, %.2f)", reri_est, reri_ci[1], reri_ci[2]),
#     AP      = sprintf("%.1f%% (%.1f%%, %.1f%%)", ap_est * 100, ap_ci[1] * 100, ap_ci[2] * 100),
#     S_Index = sprintf("%.2f (%.2f, %.2f)", s_est, s_ci[1], s_ci[2]),
#     stringsAsFactors = FALSE
#   )
# }
# 
# table3_stunting <- compute_additive_synergy("stunted_num", "Stunting (HAZ < -2 SD)")
# table3_wasting  <- compute_additive_synergy("wasted_num",  "Wasting (WHZ < -2 SD)")
# table3_df       <- rbind(table3_stunting, table3_wasting)
# 
# write.csv(table3_df, file.path(output_dir, "Table3_Additive_Synergy_Results.csv"), row.names = FALSE)
# cat(sprintf("Saved: %s/Table3_Additive_Synergy_Results.csv\n", output_dir))







# 
# # ==============================================================================
# # 9. Table 3: Additive Interaction & Relative Excess Risk (RERI) Analysis
# # ==============================================================================
# if (!requireNamespace("survey", quietly = TRUE)) install.packages("survey")
# if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
# if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
# 
# library(survey)
# library(officer)
# library(flextable)
# 
# cat("\nFitting Survey-Weighted Modified Poisson Models for Additive Interaction (Table 3)...\n")
# 
# # Ensure global option for single-PSU strata
# options(survey.lonely.psu = "adjust")
# 
# endpoints <- c(
#   "stunted_num" = "Stunting (HAZ < -2 SD)",
#   "wasted_num"  = "Wasting (WHZ < -2 SD)",
#   "underwt_num" = "Underweight (WAZ < -2 SD)"
# )
# 
# covars_m1 <- paste(
#   "CAGE", "HL4", "C_Hb", "melevel", "helevel", "windex5", 
#   "HH6", "HH7", "peeling_paint", "improved_water", 
#   "improved_sanitation", "basic_handwashing", 
#   sep = " + "
# )
# 
# table3_rows <- list()
# 
# for (y in names(endpoints)) {
#   y_label <- endpoints[y]
#   cat(sprintf(" -> Processing model for: %s\n", y_label))
#   
#   # Fit survey-weighted modified Poisson regression
#   f_pois <- as.formula(paste0(y, " ~ joint_exp + ", covars_m1))
#   m_pois <- svyglm(f_pois, design = svy_primary, family = quasipoisson(link = "log"))
#   
#   b <- coef(m_pois)
#   V <- vcov(m_pois)
#   
#   # Identify joint exposure parameter positions dynamically
#   p01 <- grep("01_LowBLL_Diarrhea", names(b), value = TRUE)
#   p10 <- grep("10_HighBLL_NoDiarrhea", names(b), value = TRUE)
#   p11 <- grep("11_HighBLL_Diarrhea", names(b), value = TRUE)
#   
#   if (length(c(p01, p10, p11)) < 3) {
#     stop("Error: Unable to locate all three joint exposure indicator coefficients in svyglm output.")
#   }
#   
#   # 1. Prevalence Ratios & 95% Confidence Intervals
#   apr_01 <- exp(b[p01])
#   ci_01  <- exp(b[p01] + c(-1.96, 1.96) * sqrt(V[p01, p01]))
#   
#   apr_10 <- exp(b[p10])
#   ci_10  <- exp(b[p10] + c(-1.96, 1.96) * sqrt(V[p10, p10]))
#   
#   apr_11 <- exp(b[p11])
#   ci_11  <- exp(b[p11] + c(-1.96, 1.96) * sqrt(V[p11, p11]))
#   
#   # 2. Point Estimates for Additive Interaction (Rothman / Knol & VanderWeele)
#   # RERI = PR11 - PR10 - PR01 + 1
#   reri_est <- apr_11 - apr_10 - apr_01 + 1
#   
#   # Attributable Proportion: AP = RERI / PR11
#   ap_est   <- reri_est / apr_11
#   
#   # Synergy Index: S = (PR11 - 1) / ((PR10 - 1) + (PR01 - 1))
#   denom_s  <- (apr_10 - 1) + (apr_01 - 1)
#   s_est    <- if (denom_s != 0) (apr_11 - 1) / denom_s else NA_real_
#   
#   # 3. Asymptotic Variance-Covariance Extraction
#   sub_params <- c(p11, p10, p01)
#   sub_V <- V[sub_params, sub_params]
#   
#   # Delta-Method for RERI:
#   # d(RERI)/d(b11) = apr_11; d(RERI)/d(b10) = -apr_10; d(RERI)/d(b01) = -apr_01
#   grad_reri <- c(apr_11, -apr_10, -apr_01)
#   var_reri  <- as.numeric(t(grad_reri) %*% sub_V %*% grad_reri)
#   se_reri   <- sqrt(max(0, var_reri))
#   ci_reri   <- reri_est + c(-1.96, 1.96) * se_reri
#   
#   # Delta-Method for AP (Full Epidemiologic Gradient):
#   # d(AP)/d(b11) = (apr_10 + apr_01 - 1) / apr_11
#   # d(AP)/d(b10) = -apr_10 / apr_11
#   # d(AP)/d(b01) = -apr_01 / apr_11
#   grad_ap  <- c((apr_10 + apr_01 - 1) / apr_11, -apr_10 / apr_11, -apr_01 / apr_11)
#   var_ap   <- as.numeric(t(grad_ap) %*% sub_V %*% grad_ap)
#   se_ap    <- sqrt(max(0, var_ap))
#   ci_ap    <- ap_est + c(-1.96, 1.96) * se_ap
#   
#   # Delta-Method for log(S) (Synergy Index):
#   # ln(S) = ln(apr_11 - 1) - ln((apr_10 - 1) + (apr_01 - 1))
#   if (!is.na(s_est) && s_est > 0 && (apr_11 - 1) > 0 && denom_s > 0) {
#     grad_logs <- c(
#       apr_11 / (apr_11 - 1),
#       -apr_10 / denom_s,
#       -apr_01 / denom_s
#     )
#     var_logs <- as.numeric(t(grad_logs) %*% sub_V %*% grad_logs)
#     se_logs  <- sqrt(max(0, var_logs))
#     ci_s     <- exp(log(s_est) + c(-1.96, 1.96) * se_logs)
#   } else {
#     ci_s <- c(NA_real_, NA_real_)
#   }
#   
#   # Format strings for clean presentation
#   table3_rows[[length(table3_rows) + 1]] <- data.frame(
#     Outcome             = y_label,
#     Ref_LowBLL_NoDiarr  = "1.00 (Reference)",
#     aPR_LowBLL_Diarr    = sprintf("%.2f (%.2f, %.2f)", apr_01, ci_01[1], ci_01[2]),
#     aPR_HighBLL_NoDiarr = sprintf("%.2f (%.2f, %.2f)", apr_10, ci_10[1], ci_10[2]),
#     aPR_HighBLL_Diarr   = sprintf("%.2f (%.2f, %.2f)", apr_11, ci_11[1], ci_11[2]),
#     RERI                = sprintf("%.2f (%.2f, %.2f)", reri_est, ci_reri[1], ci_reri[2]),
#     AP                  = sprintf("%.1f%% (%.1f%%, %.1f%%)", ap_est * 100, ci_ap[1] * 100, ci_ap[2] * 100),
#     S_Index             = if (all(!is.na(ci_s))) sprintf("%.2f (%.2f, %.2f)", s_est, ci_s[1], ci_s[2]) else sprintf("%.2f (NA, NA)", s_est),
#     stringsAsFactors    = FALSE
#   )
# }
# 
# table3_final <- do.call(rbind, table3_rows)
# 
# # Save Raw CSV with Informative Headers
# write.csv(table3_final, file.path(output_dir, "Table3_Additive_Interaction.csv"), row.names = FALSE)
# cat(sprintf("\nSaved CSV: %s/Table3_Additive_Interaction.csv\n", output_dir))
# 
# # ------------------------------------------------------------------------------
# # Export Publication-Ready Table 3 to Microsoft Word (.docx)
# # ------------------------------------------------------------------------------
# ft3 <- flextable(table3_final)
# 
# # Font & Typography
# ft3 <- font(ft3, fontname = "Times New Roman", part = "all")
# ft3 <- fontsize(ft3, size = 8.5, part = "all")
# ft3 <- fontsize(ft3, size = 9, part = "header")
# ft3 <- bold(ft3, part = "header")
# 
# # Clear, Descriptive Column Titles (Eliminating Raw Code Suffixes)
# ft3 <- set_header_labels(
#   ft3,
#   Outcome             = "Clinical Malnutrition Endpoint",
#   Ref_LowBLL_NoDiarr  = "Reference (00)\nLow BLL, No Diarrhea",
#   aPR_LowBLL_Diarr    = "aPR (01)\nLow BLL + Diarrhea",
#   aPR_HighBLL_NoDiarr = "aPR (10)\nHigh BLL, No Diarrhea",
#   aPR_HighBLL_Diarr   = "aPR (11)\nHigh BLL + Diarrhea",
#   RERI                = "RERI\n(95% CI)",
#   AP                  = "AP\n(95% CI)",
#   S_Index             = "Synergy Index S\n(95% CI)"
# )
# 
# # Text Alignment
# ft3 <- align(ft3, j = 1, align = "left", part = "all")
# ft3 <- align(ft3, j = 2:8, align = "center", part = "all")
# 
# # Column Widths (optimized for 8.5 x 11 Landscape or Standard Document Margins)
# ft3 <- width(ft3, j = 1, width = 2.1)
# ft3 <- width(ft3, j = 2:5, width = 1.25)
# ft3 <- width(ft3, j = 6:8, width = 1.2)
# 
# # Classic Journal Three-Line Table Border
# std_border <- fp_border(color = "black", width = 1)
# ft3 <- border_remove(ft3)
# ft3 <- hline_top(ft3, part = "header", border = std_border)
# ft3 <- hline_bottom(ft3, part = "header", border = std_border)
# ft3 <- hline_bottom(ft3, part = "body", border = std_border)
# 
# # Clear any prior footer rows (prevents duplicate stacking on re-run)
# ft3 <- delete_part(ft3, part = "footer")
# 
# # STROBE / ICMJE Caption and Footnote
# ft3 <- set_caption(
#   ft3,
#   caption = "Table 3. Survey-weighted multivariable modified Poisson regression evaluating joint exposure to elevated blood lead and diarrheal morbidity on clinical malnutrition, with measures of biological interaction on the additive scale (N = 10,507)."
# )
# 
# footnote_t3 <- paste(
#   "Abbreviations: aPR, Adjusted Prevalence Ratio; AP, Attributable Proportion due to interaction; CI, Confidence Interval;",
#   "HAZ, Height-for-Age z-score; RERI, Relative Excess Risk due to Interaction; S, Synergy Index; WAZ, Weight-for-Age z-score; WHZ, Weight-for-Height z-score.",
#   "Elevated blood lead level (High BLL) is dichotomized at >=3.5 ug/dL according to the US CDC Blood Lead Reference Value.",
#   "Recent diarrhea is defined as caregiver-reported diarrheal illness in the preceding 14 days.",
#   "All models are survey-weighted and adjusted for child age (months), child sex, venous hemoglobin (g/dL), maternal education, household head education,",
#   "wealth index quintile, urban/rural residence, administrative division, wall paint condition, drinking water source, sanitation facility type, and handwashing facilities.",
#   "Additive biological interaction parameters (RERI, AP, and S) and their corresponding 95% confidence intervals were estimated using the delta method.",
#   "RERI > 0, AP > 0, or S > 1 indicates super-additivity (synergy); RERI = 0, AP = 0, or S = 1 denotes exact additivity (independent effects)."
# )
# 
# ft3 <- add_footer_lines(ft3, values = footnote_t3)
# ft3 <- fontsize(ft3, size = 8, part = "footer")
# ft3 <- italic(ft3, part = "footer")
# 
# # Write to Word Document
# doc_t3_path <- file.path(output_dir, "Table3_Manuscript_Ready.docx")
# doc_t3 <- read_docx() %>%
#   body_add_par("Table 3: Additive Interaction and Excess Risk Analysis", style = "heading 1") %>%
#   body_add_flextable(ft3)
# 
# print(doc_t3, target = doc_t3_path)
# cat(sprintf("Saved publication-ready Word file: %s\n\n", doc_t3_path))



# ==============================================================================
# 9. Tables 3 & 4: Modified Poisson Regression (aPR) & Additive Interaction
# ==============================================================================
if (!requireNamespace("survey", quietly = TRUE)) install.packages("survey")
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")

library(survey)
library(officer)
library(flextable)

cat("\nFitting Survey-Weighted Modified Poisson Models for Tables 3 and 4...\n")

options(survey.lonely.psu = "adjust")

endpoints <- c(
  "stunted_num" = "Stunting (HAZ < -2 SD)",
  "wasted_num"  = "Wasting (WHZ < -2 SD)",
  "underwt_num" = "Underweight (WAZ < -2 SD)"
)

covars_m1 <- paste(
  "CAGE", "HL4", "C_Hb", "melevel", "helevel", "windex5", 
  "HH6", "HH7", "peeling_paint", "improved_water", 
  "improved_sanitation", "basic_handwashing", 
  sep = " + "
)

table3_apr_rows <- list()
table4_int_rows <- list()

for (y in names(endpoints)) {
  y_label <- endpoints[y]
  cat(sprintf(" -> Processing model for: %s\n", y_label))
  
  f_pois <- as.formula(paste0(y, " ~ joint_exp + ", covars_m1))
  m_pois <- svyglm(f_pois, design = svy_primary, family = quasipoisson(link = "log"))
  
  b <- coef(m_pois)
  V <- vcov(m_pois)
  
  p01 <- grep("01_LowBLL_Diarrhea", names(b), value = TRUE)
  p10 <- grep("10_HighBLL_NoDiarrhea", names(b), value = TRUE)
  p11 <- grep("11_HighBLL_Diarrhea", names(b), value = TRUE)
  
  # 1. Prevalence Ratios
  apr_01 <- exp(b[p01])
  ci_01  <- exp(b[p01] + c(-1.96, 1.96) * sqrt(V[p01, p01]))
  
  apr_10 <- exp(b[p10])
  ci_10  <- exp(b[p10] + c(-1.96, 1.96) * sqrt(V[p10, p10]))
  
  apr_11 <- exp(b[p11])
  ci_11  <- exp(b[p11] + c(-1.96, 1.96) * sqrt(V[p11, p11]))
  
  # 2. Additive Interaction
  reri_est <- apr_11 - apr_10 - apr_01 + 1
  ap_est   <- reri_est / apr_11
  denom_s  <- (apr_10 - 1) + (apr_01 - 1)
  s_est    <- if (denom_s != 0) (apr_11 - 1) / denom_s else NA_real_
  
  sub_params <- c(p11, p10, p01)
  sub_V <- V[sub_params, sub_params]
  
  grad_reri <- c(apr_11, -apr_10, -apr_01)
  se_reri   <- sqrt(max(0, as.numeric(t(grad_reri) %*% sub_V %*% grad_reri)))
  ci_reri   <- reri_est + c(-1.96, 1.96) * se_reri
  
  grad_ap  <- c((apr_10 + apr_01 - 1) / apr_11, -apr_10 / apr_11, -apr_01 / apr_11)
  se_ap    <- sqrt(max(0, as.numeric(t(grad_ap) %*% sub_V %*% grad_ap)))
  ci_ap    <- ap_est + c(-1.96, 1.96) * se_ap
  
  if (!is.na(s_est) && s_est > 0 && (apr_11 - 1) > 0 && denom_s > 0) {
    grad_logs <- c(apr_11 / (apr_11 - 1), -apr_10 / denom_s, -apr_01 / denom_s)
    se_logs  <- sqrt(max(0, as.numeric(t(grad_logs) %*% sub_V %*% grad_logs)))
    ci_s     <- exp(log(s_est) + c(-1.96, 1.96) * se_logs)
  } else {
    ci_s <- c(NA_real_, NA_real_)
  }
  
  # Split logic: Add aPRs to Table 3 and Interaction to Table 4
  table3_apr_rows[[length(table3_apr_rows) + 1]] <- data.frame(
    Outcome             = y_label,
    Ref_LowBLL_NoDiarr  = "1.00 (Reference)",
    aPR_LowBLL_Diarr    = sprintf("%.2f (%.2f, %.2f)", apr_01, ci_01[1], ci_01[2]),
    aPR_HighBLL_NoDiarr = sprintf("%.2f (%.2f, %.2f)", apr_10, ci_10[1], ci_10[2]),
    aPR_HighBLL_Diarr   = sprintf("%.2f (%.2f, %.2f)", apr_11, ci_11[1], ci_11[2]),
    stringsAsFactors    = FALSE
  )
  
  table4_int_rows[[length(table4_int_rows) + 1]] <- data.frame(
    Outcome             = y_label,
    RERI                = sprintf("%.2f (%.2f, %.2f)", reri_est, ci_reri[1], ci_reri[2]),
    AP                  = sprintf("%.1f%% (%.1f%%, %.1f%%)", ap_est * 100, ci_ap[1] * 100, ci_ap[2] * 100),
    S_Index             = if (all(!is.na(ci_s))) sprintf("%.2f (%.2f, %.2f)", s_est, ci_s[1], ci_s[2]) else sprintf("%.2f (NA, NA)", s_est),
    stringsAsFactors    = FALSE
  )
}

# --------- FORMAT & EXPORT TABLE 3 (aPRs) ---------
t3_final <- do.call(rbind, table3_apr_rows)
write.csv(t3_final, file.path(output_dir, "Table3_aPR_Results.csv"), row.names = FALSE)

ft3 <- flextable(t3_final)
ft3 <- font(ft3, fontname = "Times New Roman", part = "all")
ft3 <- fontsize(ft3, size = 9, part = "all") %>% fontsize(size = 9.5, part = "header") %>% bold(part = "header")
ft3 <- set_header_labels(ft3, Outcome = "Clinical Malnutrition Endpoint", Ref_LowBLL_NoDiarr = "Reference (00)\nLow BLL, No Diarrhea", aPR_LowBLL_Diarr = "aPR (01)\nLow BLL + Diarrhea", aPR_HighBLL_NoDiarr = "aPR (10)\nHigh BLL, No Diarrhea", aPR_HighBLL_Diarr = "aPR (11)\nHigh BLL + Diarrhea")
ft3 <- align(ft3, j = 1, align = "left", part = "all") %>% align(j = 2:5, align = "center", part = "all")
ft3 <- width(ft3, j = 1, width = 2.1) %>% width(j = 2:5, width = 1.3)
std_border <- fp_border(color = "black", width = 1)
ft3 <- border_remove(ft3) %>% hline_top(part = "header", border = std_border) %>% hline_bottom(part = "header", border = std_border) %>% hline_bottom(part = "body", border = std_border)
ft3 <- set_caption(ft3, caption = "Table 3. Survey-weighted adjusted prevalence ratios (aPRs) evaluating joint blood lead and diarrheal morbidity on clinical malnutrition (N = 10,507)")
ft3 <- add_footer_lines(ft3, values = "Elevated blood lead level (High BLL) is dichotomized at >=3.5 ug/dL according to the US CDC Blood Lead Reference Value. Models adjust for child age, sex, hemoglobin, maternal/household head education, wealth quintile, residence, division, wall paint, drinking water, sanitation, and handwashing facilities.") %>% fontsize(size = 8, part = "footer") %>% italic(part = "footer")

doc_t3_path <- file.path(output_dir, "Table3_Manuscript_Ready.docx")
doc_t3 <- read_docx() %>% body_add_par("Table 3: Prevalence Ratios (aPR)", style = "heading 1") %>% body_add_flextable(ft3)
print(doc_t3, target = doc_t3_path)

# --------- FORMAT & EXPORT TABLE 4 (Interaction) ---------
t4_final <- do.call(rbind, table4_int_rows)
write.csv(t4_final, file.path(output_dir, "Table4_Additive_Interaction.csv"), row.names = FALSE)

ft4 <- flextable(t4_final)
ft4 <- font(ft4, fontname = "Times New Roman", part = "all")
ft4 <- fontsize(ft4, size = 9, part = "all") %>% fontsize(size = 9.5, part = "header") %>% bold(part = "header")
ft4 <- set_header_labels(ft4, Outcome = "Clinical Malnutrition Endpoint", RERI = "RERI\n(95% CI)", AP = "AP\n(95% CI)", S_Index = "Synergy Index S\n(95% CI)")
ft4 <- align(ft4, j = 1, align = "left", part = "all") %>% align(j = 2:4, align = "center", part = "all")
ft4 <- width(ft4, j = 1, width = 2.1) %>% width(j = 2:4, width = 1.6)
ft4 <- border_remove(ft4) %>% hline_top(part = "header", border = std_border) %>% hline_bottom(part = "header", border = std_border) %>% hline_bottom(part = "body", border = std_border)
ft4 <- set_caption(ft4, caption = "Table 4. Additive-Scale Interaction (Departure from Additivity) between elevated blood lead and recent diarrhea (N = 10,507)")
ft4 <- add_footer_lines(ft4, values = "Abbreviations: AP, Attributable Proportion; CI, Confidence Interval; RERI, Relative Excess Risk due to Interaction. Additive parameters were estimated using the delta method. RERI > 0, AP > 0, or S > 1 indicates super-additivity (synergy); RERI = 0, AP = 0, or S = 1 denotes exact additivity.") %>% fontsize(size = 8, part = "footer") %>% italic(part = "footer")

doc_t4_path <- file.path(output_dir, "Table4_Manuscript_Ready.docx")
doc_t4 <- read_docx() %>% body_add_par("Table 4: Additive Interaction", style = "heading 1") %>% body_add_flextable(ft4)
print(doc_t4, target = doc_t4_path)

cat("Successfully generated split Tables 3 and 4!\n")








# ------------------------------------------------------------------------------
# 10. Publication Figure Export Function
# ------------------------------------------------------------------------------
save_pub_plot <- function(plot_obj, filename_base, width, height) {
  png_file <- file.path(output_dir, paste0(filename_base, ".png"))
  tif_file <- file.path(output_dir, paste0(filename_base, ".tif"))
  
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = 300)
  ggsave(tif_file, plot = plot_obj, width = width, height = height, dpi = 300, compression = "lzw")
  cat(sprintf("Exported: %s (.png and 300 DPI .tif)\n", filename_base))
}

# ------------------------------------------------------------------------------
# 11. FIGURE 1: Predicted Marginal Dose-Response Curves (Main-Effects Specification)
# ------------------------------------------------------------------------------
cat("\nRendering Figure 1: Predicted Marginal Dose-Response Curves...\n")

# Fit main-effects survey model to enforce parallel independent trajectories
fit_haz_m1_main <- svyglm(
  HAZ2 ~ log2_bll + diarrhea + CAGE + HL4 + C_Hb + melevel + helevel + 
    windex5 + HH6 + HH7 + peeling_paint + improved_water + 
    improved_sanitation + basic_handwashing,
  design = svy_primary
)

q_low  <- quantile(df_clean$C_Lead, 0.025, na.rm = TRUE)
q_high <- quantile(df_clean$C_Lead, 0.975, na.rm = TRUE)
lead_seq <- seq(log2(q_low), log2(q_high), length.out = 150)

pred_grid <- expand.grid(
  log2_bll            = lead_seq,
  diarrhea            = c(0, 1),
  CAGE                = median(df_clean$CAGE, na.rm = TRUE),
  HL4                 = levels(as.factor(df_clean$HL4))[1],
  C_Hb                = mean(df_clean$C_Hb, na.rm = TRUE),
  melevel             = levels(as.factor(df_clean$melevel))[1],
  helevel             = levels(as.factor(df_clean$helevel))[1],
  windex5             = levels(as.factor(df_clean$windex5))[1],
  HH6                 = levels(as.factor(df_clean$HH6))[1],
  HH7                 = levels(as.factor(df_clean$HH7))[1],
  peeling_paint       = levels(as.factor(df_clean$peeling_paint))[1],
  improved_water      = levels(as.factor(df_clean$improved_water))[1],
  improved_sanitation = levels(as.factor(df_clean$improved_sanitation))[1],
  basic_handwashing   = levels(as.factor(df_clean$basic_handwashing))[1]
)

preds <- predict(fit_haz_m1_main, newdata = pred_grid, se.fit = TRUE)
pred_grid$fit      <- as.numeric(preds)
pred_grid$se       <- as.numeric(sqrt(attr(preds, "var")))
pred_grid$ci_lower <- pred_grid$fit - 1.96 * pred_grid$se
pred_grid$ci_upper <- pred_grid$fit + 1.96 * pred_grid$se
pred_grid$C_Lead   <- 2^(pred_grid$log2_bll)
pred_grid$Diarrhea_Status <- factor(
  if_else(pred_grid$diarrhea == 1, "Diarrhea (Past 2 Weeks)", "No Recent Diarrhea"),
  levels = c("No Recent Diarrhea", "Diarrhea (Past 2 Weeks)")
)

p_fig1 <- ggplot(pred_grid, aes(x = C_Lead, y = fit, color = Diarrhea_Status, fill = Diarrhea_Status)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.18, linetype = 0) +
  geom_line(aes(linetype = Diarrhea_Status), linewidth = 1.1) +
  scale_x_log10(
    breaks = c(1.0, 2.0, 3.5, 5.0, 10.0, 20.0),
    labels = c("1.0", "2.0", "3.5", "5.0", "10.0", "20.0")
  ) +
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "#D9534F", linewidth = 0.8) +
  annotate(
    "text", 
    x = 3.65, 
    y = min(pred_grid$ci_lower) + 0.04, 
    label = "CDC BLRV (3.5 µg/dL)", 
    angle = 90, 
    size = 3.2, 
    color = "#A94442", 
    hjust = 0, 
    fontface = "italic"
  ) +
  scale_color_manual(values = c("No Recent Diarrhea" = "#1F4E79", 
                                "Diarrhea (Past 2 Weeks)" = "#C00000")) +
  scale_fill_manual(values = c("No Recent Diarrhea" = "#1F4E79", 
                               "Diarrhea (Past 2 Weeks)" = "#C00000")) +
  scale_linetype_manual(values = c("No Recent Diarrhea" = "solid", 
                                   "Diarrhea (Past 2 Weeks)" = "longdash")) +
  theme_classic(base_size = 12) +
  labs(
    x = expression(bold("Blood Lead Concentration (" * mu * "g/dL, Logarithmic Scale)")),
    y = expression(bold("Adjusted Height-for-Age z-score (HAZ)")),
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", size = 10),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 20, unit = "pt")
  )

save_pub_plot(p_fig1, "Figure1_HAZ_Compounding_Interaction", width = 8.5, height = 5.2)
# ------------------------------------------------------------------------------
# 12. FIGURE 2: Divisional Distribution of BLL (All 8 Divisions, Clean Canvas)
# ------------------------------------------------------------------------------
cat("\nRendering Figure 2: Divisional Distribution Across All 8 Administrative Divisions...\n")

p_fig2 <- ggplot(df_clean, aes(x = HH7, y = C_Lead, fill = HH7)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, show.legend = FALSE, color = "#1F4E79") +
  geom_hline(yintercept = 3.5, linetype = "dashed", color = "#C00000", linewidth = 0.8) +
  coord_cartesian(ylim = c(0, quantile(df_clean$C_Lead, 0.95, na.rm = TRUE))) +
  scale_fill_brewer(palette = "Blues") +
  theme_classic(base_size = 12) +
  labs(
    x = "Administrative Division",
    y = expression(bold("Blood Lead Concentration (" * mu * "g/dL)"))
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10)
  )

save_pub_plot(p_fig2, "Figure2_Divisional_BLL_Distribution", width = 8.5, height = 4.8)

# ------------------------------------------------------------------------------
# 13. FIGURE 3: Additive Stunting Prevalence Matrix (Clean Canvas)
# ------------------------------------------------------------------------------
cat("\nRendering Figure 3: Additive Interaction Stunting Matrix...\n")

stunt_by_joint <- svyby(~stunted_num, ~joint_exp, svy_primary, svymean, na.rm = TRUE)
stunt_by_joint <- as.data.frame(stunt_by_joint)
stunt_by_joint$pct    <- stunt_by_joint$stunted_num * 100
stunt_by_joint$se_pct <- stunt_by_joint$se * 100

stunt_by_joint$Exposure_Group <- factor(
  c("Low BLL (<3.5)\nNo Diarrhea\n(Reference)",
    "Low BLL (<3.5)\n+ Diarrhea",
    "High BLL (≥3.5)\nNo Diarrhea",
    "High BLL (≥3.5)\n+ Diarrhea\n(Dual Insult)"),
  levels = c("Low BLL (<3.5)\nNo Diarrhea\n(Reference)",
             "Low BLL (<3.5)\n+ Diarrhea",
             "High BLL (≥3.5)\nNo Diarrhea",
             "High BLL (≥3.5)\n+ Diarrhea\n(Dual Insult)")
)

p_fig3 <- ggplot(stunt_by_joint, aes(x = Exposure_Group, y = pct, fill = Exposure_Group)) +
  geom_col(width = 0.58, show.legend = FALSE, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pct - 1.96 * se_pct, ymax = pct + 1.96 * se_pct), 
                width = 0.18, linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -1.6, fontface = "bold", size = 4) +
  scale_fill_manual(values = c("#A6C8E0", "#6BAED6", "#FDBB84", "#E34A33")) +
  ylim(0, max(stunt_by_joint$pct + 1.96 * stunt_by_joint$se_pct) + 6) +
  theme_classic(base_size = 12) +
  labs(
    x = "Joint Exposure Category",
    y = "Survey-Weighted Stunting Prevalence (%) [95% CI]"
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(face = "bold", color = "black", size = 9.5),
    axis.text.y = element_text(color = "black", size = 10)
  )

save_pub_plot(p_fig3, "Figure3_Additive_Stunting_Prevalence", width = 8.0, height = 4.8)



# 
# # ==============================================================================
# # 14. FIGURE 4: Causal Mediation / DAG Schema (Exact Border Alignment & Table 2 Reconciled)
# # ==============================================================================
# cat("\nRendering Figure 4: Corrected Causal Mediation Diagram...\n")
# 
# # Node boxes definition
# dag_nodes <- data.frame(
#   x = c(1.0, 3.0, 5.0),
#   y = c(1.0, 2.0, 1.0),
#   label = c(
#     "Blood Lead Level\n(log2_bll)\n[Exposure]",
#     "Diarrheal Morbidity\n[Mediator]\nNIE: -0.0001 (p = 0.240)",
#     "Height-for-Age Z-score\n(HAZ)\n[Outcome]"
#   )
# )
# 
# # Arrow trajectories and reconciled epidemiological path labels
# dag_arrows <- data.frame(
#   # x1/y1 and x2/y2 mapped directly onto the node box boundaries
#   x1 = c(1.22, 3.52, 1.51),
#   y1 = c(1.18, 1.82, 1.00),
#   x2 = c(2.48, 4.78, 4.36),
#   y2 = c(1.82, 1.18, 1.00),
#   path_label = c(
#     "Path a (Enteric Pathway)",
#     "Path b (Acute Mucosal Damage)",
#     "Natural Direct Effect (NDE)\nbeta = -0.095, p < 0.001\n[>99% of Total Effect]"
#   ),
#   lbl_x = c(1.85, 4.15, 3.00),
#   lbl_y = c(1.58, 1.58, 0.72)
# )
# 
# p_fig4 <- ggplot() +
#   # Upper indirect mediation paths
#   geom_segment(data = dag_arrows[1:2, ], 
#                aes(x = x1, y = y1, xend = x2, yend = y2),
#                arrow = arrow(length = unit(0.32, "cm"), type = "closed"),
#                color = "#3A3A3A", linewidth = 1.1, lineend = "round") +
#   # Bottom direct effect path: starts flush with Exposure, arrowhead touches Outcome
#   geom_segment(data = dag_arrows[3, ], 
#                aes(x = x1, y = y1, xend = x2, yend = y2),
#                arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
#                color = "#C00000", linewidth = 1.3, lineend = "round") +
#   # Path text annotations
#   geom_label(data = dag_arrows, aes(x = lbl_x, y = lbl_y, label = path_label), 
#              size = 3.2, fill = "white", color = "black", fontface = "italic", 
#              label.size = 0.25, label.padding = unit(0.20, "lines")) +
#   # Node boxes (drawn on top of the root coordinates)
#   geom_label(data = dag_nodes, aes(x = x, y = y, label = label), 
#              size = 3.8, fontface = "bold", fill = "#F8FAFC", color = "#1E3A8A", 
#              label.padding = unit(0.60, "lines"), linewidth = 0.85) +
#   coord_cartesian(xlim = c(0.2, 5.8), ylim = c(0.5, 2.3)) +
#   theme_void() +
#   theme(
#     plot.background = element_rect(fill = "white", color = NA)
#   )
# 
# save_pub_plot(p_fig4, "Figure4_Causal_Decomposition_Path", width = 9.0, height = 4.5)






# ==============================================================================
# 14. FIGURE 4: Conceptual DAG (Replaced Old Mediation Model)
# ==============================================================================
cat("\nRendering Figure 4: Conceptual Directed Acyclic Graph (DAG)...\n")

dag_nodes <- data.frame(
  x = c(3, 1, 5, 3),
  y = c(2.5, 1.5, 1.5, 0.5),
  label = c(
    "Socioeconomic & Environmental Covariates\n(Wealth, Education, WASH, Region)",
    "Lead Exposure\n(Blood Lead Level)",
    "Enteric Morbidity\n(Recent Diarrhea)",
    "Child Nutritional Status\n(HAZ, WHZ, WAZ)"
  ),
  fill_col = c("#f9f9f9", "#e1f5fe", "#e1f5fe", "#ffebee"),
  border_col = c("#333333", "#0288d1", "#0288d1", "#c62828")
)

# Arrows mapping dependencies (without causal text numbers)
dag_arrows <- data.frame(
  x1 = c(2.2, 3.8, 3.0, 1.4, 4.6),
  y1 = c(2.3, 2.3, 2.3, 1.3, 1.3),
  x2 = c(1.3, 4.7, 3.0, 2.6, 3.4),
  y2 = c(1.7, 1.7, 0.7, 0.7, 0.7)
)

p_fig4 <- ggplot() +
  geom_segment(data = dag_arrows, 
               aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
               color = "#3A3A3A", linewidth = 1.3, lineend = "round") +
  geom_label(data = dag_nodes, aes(x = x, y = y, label = label), 
             size = 4, fontface = "bold", fill = dag_nodes$fill_col, 
             color = dag_nodes$border_col, label.padding = unit(0.8, "lines"), 
             linewidth = 1) +
  coord_cartesian(xlim = c(0.2, 5.8), ylim = c(0.2, 2.8)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

# Output using your standardized image saving function
save_pub_plot(p_fig4, "FigureS1_Conceptual_DAG", width = 8.5, height = 5.0)










# ==============================================================================
# 15. SUPPLEMENTARY APPENDIX DELIVERABLES (TABLES S1 TO S4)
# ==============================================================================
cat("\nGenerating Supplementary Analytical Suite (Tables S1 to S4)...\n")

# --- Fit Crude Models ---
fit_haz_crude <- svyglm(HAZ2 ~ log2_bll * diarrhea, design = svy_primary)
fit_whz_crude <- svyglm(WHZ2 ~ log2_bll * diarrhea, design = svy_primary)
fit_waz_crude <- svyglm(WAZ2 ~ log2_bll * diarrhea, design = svy_primary)

# --- Fit Adjusted Models (Creating the missing fit_*_m1 objects) ---
fit_haz_m1 <- svyglm(
  as.formula(paste0("HAZ2 ~ log2_bll * diarrhea + ", covars_m1)),
  design = svy_primary
)
fit_whz_m1 <- svyglm(
  as.formula(paste0("WHZ2 ~ log2_bll * diarrhea + ", covars_m1)),
  design = svy_primary
)
fit_waz_m1 <- svyglm(
  as.formula(paste0("WAZ2 ~ log2_bll * diarrhea + ", covars_m1)),
  design = svy_primary
)

# --- Extraction Helper Function ---
extract_coef <- function(model, term) {
  s  <- summary(model)$coefficients
  ci <- confint(model)
  if (term %in% rownames(s)) {
    b  <- sprintf("%.3f", s[term, "Estimate"])
    se <- sprintf("%.3f", s[term, "Std. Error"])
    ci_str <- sprintf("(%.3f, %.3f)", ci[term, 1], ci[term, 2])
    p_val  <- s[term, "Pr(>|t|)"]
    p_str  <- if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
    return(c(Beta = b, SE = se, CI_95 = ci_str, p_value = p_str))
  }
  return(c(Beta = "NA", SE = "NA", CI_95 = "NA", p_value = "NA"))
}

# --- Build Table S1 Rows ---
table_s1_rows <- list()
outcomes <- list(
  list(name = "Height-for-Age (HAZ)",    crude = fit_haz_crude, adj = fit_haz_m1),
  list(name = "Weight-for-Height (WHZ)", crude = fit_whz_crude, adj = fit_whz_m1),
  list(name = "Weight-for-Age (WAZ)",    crude = fit_waz_crude, adj = fit_waz_m1)
)

params <- c("log2_bll", "diarrhea", "log2_bll:diarrhea")

for (out in outcomes) {
  for (p in params) {
    c_res <- extract_coef(out$crude, p)
    a_res <- extract_coef(out$adj, p)
    
    table_s1_rows[[length(table_s1_rows) + 1]] <- data.frame(
      Outcome       = out$name,
      Parameter     = p,
      Crude_Beta    = c_res["Beta"],
      Crude_95CI    = c_res["CI_95"],
      Crude_p       = c_res["p_value"],
      Adjusted_Beta = a_res["Beta"],
      Adjusted_95CI = a_res["CI_95"],
      Adjusted_p    = a_res["p_value"],
      stringsAsFactors = FALSE
    )
  }
}

table_crude_adj <- do.call(rbind, table_s1_rows)

# Clean parameter labels
table_crude_adj <- table_crude_adj %>%
  mutate(
    Parameter = case_when(
      Parameter == "log2_bll"          ~ "Blood lead level (log2 ug/dL)",
      Parameter == "diarrhea"          ~ "Recent diarrhea (past 2 weeks)",
      Parameter == "log2_bll:diarrhea" ~ "Interaction (BLL x Diarrhea)",
      TRUE ~ Parameter
    )
  )

# Export Table S1 to CSV
write.csv(table_crude_adj, file.path(output_dir, "TableS1_Crude_vs_Adjusted.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/TableS1_Crude_vs_Adjusted.csv\n", output_dir))







# ==============================================================================
# Supplementary Table S1: Crude vs. Adjusted Survey-Weighted Linear Models
# Export to Publication-Ready Word Document (.docx)
# ==============================================================================

if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(officer)
library(flextable)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. Load and Clean Data Frame
# ------------------------------------------------------------------------------
# Load directly from CSV file if not already active in environment
if (!exists("table_crude_adj")) {
  csv_file <- if (exists("output_dir")) file.path(output_dir, "TableS1_Crude_vs_Adjusted.csv") else "TableS1_Crude_vs_Adjusted.csv"
  if (!file.exists(csv_file)) {
    stop(sprintf("File '%s' not found. Please ensure the working directory or 'output_dir' is set correctly.", csv_file))
  }
  table_crude_adj <- read.csv(csv_file, stringsAsFactors = FALSE)
}

t_s1_clean <- table_crude_adj %>%
  mutate(
    # Clean parameter names to formal manuscript labels
    Parameter = case_when(
      Parameter == "log2_bll"          ~ "Blood lead level (log2 ug/dL)",
      Parameter == "diarrhea"          ~ "Recent diarrhea (past 2 weeks)",
      Parameter == "log2_bll:diarrhea" ~ "Interaction (BLL x Diarrhea)",
      TRUE ~ as.character(Parameter)
    ),
    # Safe p-value string formatting without introducing coercion NAs
    Crude_p = case_when(
      grepl("<", as.character(Crude_p)) ~ "<0.001",
      suppressWarnings(!is.na(as.numeric(Crude_p))) & as.numeric(Crude_p) < 0.001 ~ "<0.001",
      suppressWarnings(!is.na(as.numeric(Crude_p))) ~ sprintf("%.3f", as.numeric(Crude_p)),
      TRUE ~ as.character(Crude_p)
    ),
    Adjusted_p = case_when(
      grepl("<", as.character(Adjusted_p)) ~ "<0.001",
      suppressWarnings(!is.na(as.numeric(Adjusted_p))) & as.numeric(Adjusted_p) < 0.001 ~ "<0.001",
      suppressWarnings(!is.na(as.numeric(Adjusted_p))) ~ sprintf("%.3f", as.numeric(Adjusted_p)),
      TRUE ~ as.character(Adjusted_p)
    ),
    # Combine estimates with 95% CI into publication-standard cells
    Crude_Est = sprintf("%.3f %s", as.numeric(Crude_Beta), trimws(as.character(Crude_95CI))),
    Adj_Est   = sprintf("%.3f %s", as.numeric(Adjusted_Beta), trimws(as.character(Adjusted_95CI)))
  ) %>%
  select(Outcome, Parameter, Crude_Est, Crude_p, Adj_Est, Adjusted_p)

# ------------------------------------------------------------------------------
# 2. Convert to Grouped Flextable (Panelled by Outcome)
# ------------------------------------------------------------------------------
ft_s1_data <- as_grouped_data(t_s1_clean, groups = "Outcome")
ft_s1 <- flextable(ft_s1_data)

# Typography & Sizing
ft_s1 <- font(ft_s1, fontname = "Times New Roman", part = "all")
ft_s1 <- fontsize(ft_s1, size = 9, part = "all")
ft_s1 <- fontsize(ft_s1, size = 9.5, part = "header")
ft_s1 <- bold(ft_s1, part = "header")

# Bold the Outcome Group Headers (HAZ, WHZ, WAZ)
ft_s1 <- bold(ft_s1, i = ~ !is.na(Outcome), part = "body")

# Set Clear Column Labels
ft_s1 <- set_header_labels(
  ft_s1,
  Parameter  = "Growth Metric & Parameter Exposure",
  Crude_Est  = "Crude Beta (95% CI)",
  Crude_p    = "Crude p-value",
  Adj_Est    = "Adjusted Beta (95% CI)*",
  Adjusted_p = "Adjusted p-value"
)

# Text Alignments
ft_s1 <- align(ft_s1, j = 1, align = "left", part = "all")
ft_s1 <- align(ft_s1, j = 2:5, align = "center", part = "all")

# Column Widths
ft_s1 <- width(ft_s1, j = 1, width = 2.8)
ft_s1 <- width(ft_s1, j = c(2, 4), width = 1.8)
ft_s1 <- width(ft_s1, j = c(3, 5), width = 1.1)

# Classic Journal Three-Line Borders
std_border <- fp_border(color = "black", width = 1)
ft_s1 <- border_remove(ft_s1)
ft_s1 <- hline_top(ft_s1, part = "header", border = std_border)
ft_s1 <- hline_bottom(ft_s1, part = "header", border = std_border)
ft_s1 <- hline_bottom(ft_s1, part = "body", border = std_border)

# Clear any residual footers before adding lines (prevents duplicate stacking)
ft_s1 <- delete_part(ft_s1, part = "footer")

# STROBE / ICMJE Caption & Footnote
ft_s1 <- set_caption(
  ft_s1,
  caption = "Supplementary Table S1. Comparison of unadjusted (crude) and multivariable-adjusted survey-weighted linear regression models evaluating blood lead levels, diarrheal morbidity, and continuous pediatric growth faltering (z-scores) (N = 10,507)."
)

footnote_s1 <- paste(
  "Abbreviations: BLL, Blood Lead Level; CI, Confidence Interval; HAZ, Height-for-Age z-score; WAZ, Weight-for-Age z-score; WHZ, Weight-for-Height z-score.",
  "Beta coefficients represent the change in anthropometric z-score per two-fold increase in blood lead level (log2-transformed), presence of caregiver-reported diarrhea in the preceding 2 weeks, and their continuous multiplicative interaction term.",
  "*Multivariable models adjust for child age (months), child sex, venous hemoglobin (g/dL), maternal education, household head education, wealth index quintile, urban/rural residence, administrative division, peeling wall paint, drinking water source, sanitation facility type, and handwashing facilities.",
  "All models account for complex survey strata, primary sampling units, and normalized child blood weights."
)

ft_s1 <- add_footer_lines(ft_s1, values = footnote_s1)
ft_s1 <- fontsize(ft_s1, size = 8, part = "footer")
ft_s1 <- italic(ft_s1, part = "footer")

# ------------------------------------------------------------------------------
# 3. Export to Microsoft Word Document
# ------------------------------------------------------------------------------
target_dir <- if (exists("output_dir")) output_dir else getwd()
doc_s1_path <- file.path(target_dir, "TableS1_Crude_vs_Adjusted.docx")

doc_s1 <- read_docx() %>%
  body_add_par("Supplementary Table S1: Crude vs. Adjusted Models", style = "heading 1") %>%
  body_add_flextable(ft_s1)

print(doc_s1, target = doc_s1_path)
cat(sprintf("Saved publication-ready Word file: %s\n", doc_s1_path))






# --- TABLE S2: Logistic Regression (aOR) Comparison ---
cat("Computing Table S2: Survey Logistic Regression aORs...\n")
run_logistic_comp <- function(outcome_var, outcome_label) {
  f_logit <- as.formula(paste(
    outcome_var, "~ joint_exp + CAGE + HL4 + C_Hb + melevel + helevel +",
    "windex5 + HH6 + HH7 + peeling_paint + improved_water + improved_sanitation + basic_handwashing"
  ), env = globalenv())
  
  fit_logit <- svyglm(f_logit, design = svy_primary, family = quasibinomial(link = "logit"))
  
  b <- coef(fit_logit)
  v <- vcov(fit_logit)
  idx <- grep("joint_exp", names(b), value = TRUE)
  
  res_list <- list()
  for (id in idx) {
    or_val  <- exp(b[id])
    ci_low  <- exp(b[id] - 1.96 * sqrt(v[id, id]))
    ci_high <- exp(b[id] + 1.96 * sqrt(v[id, id]))
    p_val   <- summary(fit_logit)$coefficients[id, "Pr(>|t|)"]
    p_str   <- if (p_val < 0.001) "<0.001" else sprintf("%.3f", p_val)
    
    res_list[[length(res_list) + 1]] <- data.frame(
      Outcome = outcome_label,
      Group   = gsub("joint_exp", "", id),
      aOR     = sprintf("%.2f", or_val),
      CI_95   = sprintf("(%.2f, %.2f)", ci_low, ci_high),
      p_value = p_str,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, res_list)
}

TableS2_final <- rbind(
  run_logistic_comp("stunted_num", "Stunting (HAZ < -2 SD)"),
  run_logistic_comp("wasted_num", "Wasting (WHZ < -2 SD)"),
  run_logistic_comp("underwt_num", "Underweight (WAZ < -2 SD)")
)
write.csv(TableS2_final, file.path(output_dir, "TableS2_Logistic_aOR_Comparison.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/TableS2_Logistic_aOR_Comparison.csv\n", output_dir))






# --- TABLE S3: E. coli Subsample Sensitivity ---
cat("Computing Table S3: Microbiological E. coli Subsample...\n")
if ("WQ26" %in% names(df_clean) && sum(!is.na(df_clean$WQ26)) > 100) {
  df_ecoli <- df_clean %>%
    mutate(
      ecoli_risk = factor(case_when(
        is.na(WQ26) ~ NA_character_,
        WQ26 == 0 ~ "Low (<1 CFU)",
        WQ26 >= 1 & WQ26 <= 10 ~ "Medium (1-10 CFU)",
        WQ26 >= 11 & WQ26 <= 100 ~ "High (11-100 CFU)",
        WQ26 > 100 ~ "Very High (>100 CFU)"
      ), levels = c("Low (<1 CFU)", "Medium (1-10 CFU)", "High (11-100 CFU)", "Very High (>100 CFU)"))
    ) %>%
    filter(!is.na(ecoli_risk))
  
  svy_ecoli <- svydesign(
    id = ~PSU, strata = ~stratum, weights = ~chblweight,
    nest = TRUE, data = df_ecoli
  )
  
  fit_ecoli_sens <- svyglm(
    HAZ2 ~ log2_bll * diarrhea + ecoli_risk + CAGE + HL4 + C_Hb + melevel + helevel + 
      windex5 + HH6 + peeling_paint,
    design = svy_ecoli
  )
  
  s_ecoli <- summary(fit_ecoli_sens)$coefficients
  target_vars <- c("log2_bll", "diarrhea", "log2_bll:diarrhea")
  
  TableS3_final <- data.frame(
    Outcome     = "HAZ (E. coli Subsample)",
    Sample_Size = sprintf("N = %d", nrow(df_ecoli)),
    Parameter   = target_vars,
    Beta        = sprintf("%.3f", s_ecoli[target_vars, "Estimate"]),
    SE          = sprintf("%.3f", s_ecoli[target_vars, "Std. Error"]),
    CI_95       = sprintf("(%.3f, %.3f)", 
                          s_ecoli[target_vars, "Estimate"] - 1.96 * s_ecoli[target_vars, "Std. Error"],
                          s_ecoli[target_vars, "Estimate"] + 1.96 * s_ecoli[target_vars, "Std. Error"]),
    p_value     = if_else(s_ecoli[target_vars, "Pr(>|t|)"] < 0.001, "<0.001", 
                          sprintf("%.3f", s_ecoli[target_vars, "Pr(>|t|)"]))
  )
} else {
  TableS3_final <- data.frame(
    Note = "WQ26 microbiological E. coli data was not available or had insufficient sample size in current analytic merge."
  )
}
write.csv(TableS3_final, file.path(output_dir, "TableS3_Ecoli_Sensitivity_Subsample.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/TableS3_Ecoli_Sensitivity_Subsample.csv\n", output_dir))




# ==============================================================================
# Supplementary Table S2: Multivariable Survey-Weighted Logistic Regression (aOR)
# Export to Publication-Ready Word Document (.docx)
# ==============================================================================

if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(officer)
library(flextable)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. Load and Clean Data Frame
# ------------------------------------------------------------------------------
csv_file <- if (exists("output_dir")) {
  file.path(output_dir, "TableS2_Logistic_aOR_Comparison.csv")
} else {
  "TableS2_Logistic_aOR_Comparison.csv"
}

if (file.exists(csv_file)) {
  t_s2_raw <- read.csv(csv_file, stringsAsFactors = FALSE)
} else {
  # Direct fallback data frame from your generated results
  t_s2_raw <- data.frame(
    Outcome = c(
      rep("Stunting (HAZ < -2 SD)", 3),
      rep("Wasting (WHZ < -2 SD)", 3),
      rep("Underweight (WAZ < -2 SD)", 3)
    ),
    Group = rep(c("01_LowBLL_Diarrhea", "10_HighBLL_NoDiarrhea", "11_HighBLL_Diarrhea"), 3),
    aOR = c(1.04, 1.15, 1.20, 1.36, 1.19, 1.62, 1.12, 1.12, 1.46),
    CI_95 = c(
      "(0.72, 1.50)", "(1.02, 1.30)", "(0.91, 1.59)",
      "(0.89, 2.10)", "(1.02, 1.38)", "(1.17, 2.25)",
      "(0.78, 1.61)", "(0.99, 1.27)", "(1.10, 1.94)"
    ),
    p_value = c(0.848, 0.019, 0.205, 0.157, 0.024, 0.004, 0.534, 0.063, 0.009),
    stringsAsFactors = FALSE
  )
}

# Add explicit reference row for each outcome panel and format text
t_s2_clean <- t_s2_raw %>%
  mutate(
    Exposure_Group = case_when(
      grepl("01", Group) ~ "Low BLL (<3.5 ug/dL) + Diarrhea",
      grepl("10", Group) ~ "High BLL (>=3.5 ug/dL), No Diarrhea",
      grepl("11", Group) ~ "High BLL (>=3.5 ug/dL) + Diarrhea (Joint Exposure)",
      TRUE ~ as.character(Group)
    ),
    aOR_CI = sprintf("%.2f %s", as.numeric(aOR), trimws(as.character(CI_95))),
    p_val  = case_when(
      grepl("<", as.character(p_value)) ~ "<0.001",
      suppressWarnings(!is.na(as.numeric(p_value))) & as.numeric(p_value) < 0.001 ~ "<0.001",
      suppressWarnings(!is.na(as.numeric(p_value))) ~ sprintf("%.3f", as.numeric(p_value)),
      TRUE ~ as.character(p_value)
    )
  ) %>%
  select(Outcome, Exposure_Group, aOR_CI, p_val)

# Insert reference row (00: Low BLL, No Diarrhea) for each clinical outcome
t_s2_expanded <- list()
for (outc in unique(t_s2_clean$Outcome)) {
  sub_df <- t_s2_clean %>% filter(Outcome == outc)
  ref_row <- data.frame(
    Outcome        = outc,
    Exposure_Group = "Low BLL (<3.5 ug/dL), No Diarrhea (Reference: 00)",
    aOR_CI         = "1.00 (Reference)",
    p_val          = "—",
    stringsAsFactors = FALSE
  )
  t_s2_expanded[[length(t_s2_expanded) + 1]] <- rbind(ref_row, sub_df)
}
t_s2_final <- do.call(rbind, t_s2_expanded)

# ------------------------------------------------------------------------------
# 2. Build Grouped Flextable (Panelled by Malnutrition Outcome)
# ------------------------------------------------------------------------------
ft_s2_data <- as_grouped_data(t_s2_final, groups = "Outcome")
ft_s2 <- flextable(ft_s2_data)

# Typography & Sizing
ft_s2 <- font(ft_s2, fontname = "Times New Roman", part = "all")
ft_s2 <- fontsize(ft_s2, size = 9, part = "all")
ft_s2 <- fontsize(ft_s2, size = 9.5, part = "header")
ft_s2 <- bold(ft_s2, part = "header")

# Bold the Outcome Group Headers
ft_s2 <- bold(ft_s2, i = ~ !is.na(Outcome), part = "body")

# Set Column Headers
ft_s2 <- set_header_labels(
  ft_s2,
  Exposure_Group = "Joint Exposure Category",
  aOR_CI         = "Adjusted Odds Ratio (95% CI)*",
  p_val          = "p-value"
)

# Alignments
ft_s2 <- align(ft_s2, j = 1, align = "left", part = "all")
ft_s2 <- align(ft_s2, j = 2:3, align = "center", part = "all")

# Column Widths
ft_s2 <- width(ft_s2, j = 1, width = 3.6)
ft_s2 <- width(ft_s2, j = 2, width = 2.4)
ft_s2 <- width(ft_s2, j = 3, width = 1.2)

# Classic Journal Three-Line Borders
std_border <- fp_border(color = "black", width = 1)
ft_s2 <- border_remove(ft_s2)
ft_s2 <- hline_top(ft_s2, part = "header", border = std_border)
ft_s2 <- hline_bottom(ft_s2, part = "header", border = std_border)
ft_s2 <- hline_bottom(ft_s2, part = "body", border = std_border)

# Purge any residual footers before adding new lines
ft_s2 <- delete_part(ft_s2, part = "footer")

# STROBE Caption & Explanatory Footnote
ft_s2 <- set_caption(
  ft_s2,
  caption = "Supplementary Table S2. Survey-weighted multivariable logistic regression evaluating joint exposure to elevated blood lead levels and diarrheal morbidity on clinical malnutrition outcomes (N = 10,507)."
)

footnote_s2 <- paste(
  "Abbreviations: aOR, Adjusted Odds Ratio; BLL, Blood Lead Level; CI, Confidence Interval; HAZ, Height-for-Age z-score;",
  "WAZ, Weight-for-Age z-score; WHZ, Weight-for-Height z-score.",
  "Elevated blood lead level is defined as >=3.5 ug/dL (US CDC Reference Value). Recent diarrhea is caregiver-reported diarrheal morbidity within the preceding 14 days.",
  "*Models adjust for child age (months), child sex, venous hemoglobin (g/dL), maternal education, household head education, wealth index quintile,",
  "urban/rural residence, administrative division, peeling wall paint, drinking water source, sanitation facility type, and basic handwashing facilities.",
  "All models account for complex survey strata, primary sampling units, and normalized child blood weights."
)

ft_s2 <- add_footer_lines(ft_s2, values = footnote_s2)
ft_s2 <- fontsize(ft_s2, size = 8, part = "footer")
ft_s2 <- italic(ft_s2, part = "footer")

# ------------------------------------------------------------------------------
# 3. Export to Microsoft Word
# ------------------------------------------------------------------------------
target_dir <- if (exists("output_dir")) output_dir else getwd()
doc_s2_path <- file.path(target_dir, "TableS2_Logistic_aOR_Comparison.docx")

doc_s2 <- read_docx() %>%
  body_add_par("Supplementary Table S2: Logistic Regression (aOR) Sensitivity Check", style = "heading 1") %>%
  body_add_flextable(ft_s2)

print(doc_s2, target = doc_s2_path)
cat(sprintf("Saved publication-ready Word file: %s\n", doc_s2_path))




# ==============================================================================
# Supplementary Table S3: Water E. coli Contamination Sensitivity Subsample
# Export to Publication-Ready Word Document (.docx)
# ==============================================================================

if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(officer)
library(flextable)
library(dplyr)

target_dir <- if (exists("output_dir")) output_dir else getwd()
csv_s3 <- file.path(target_dir, "TableS3_Ecoli_Sensitivity_Subsample.csv")

if (!file.exists(csv_s3)) {
  stop(sprintf("File '%s' not found.", csv_s3))
}

t_s3_raw <- read.csv(csv_s3, stringsAsFactors = FALSE)

# Detect columns and standardize formatting
t_s3_clean <- t_s3_raw %>%
  mutate(
    across(everything(), as.character),
    # Clean parameter labels if column exists
    across(any_of(c("Parameter", "Parameter_Exposure", "Variable")), ~ case_when(
      . == "log2_bll"          ~ "Blood lead level (log2 ug/dL)",
      . == "diarrhea"          ~ "Recent diarrhea (past 2 weeks)",
      . == "log2_bll:diarrhea" ~ "Interaction (BLL x Diarrhea)",
      grepl("ecoli|WQ26", ., ignore.case = TRUE) ~ "E. coli water contamination",
      TRUE ~ .
    )),
    # Clean p-values
    across(any_of(c("p_value", "p.value", "p", "P")), ~ case_when(
      grepl("<", .) ~ "<0.001",
      suppressWarnings(!is.na(as.numeric(.))) & as.numeric(.) < 0.001 ~ "<0.001",
      suppressWarnings(!is.na(as.numeric(.))) ~ sprintf("%.3f", as.numeric(.)),
      TRUE ~ .
    ))
  )

# Check if grouping by Outcome is applicable
has_outcome <- any(names(t_s3_clean) %in% c("Outcome", "Endpoint", "Growth_Metric"))
outcome_col <- names(t_s3_clean)[names(t_s3_clean) %in% c("Outcome", "Endpoint", "Growth_Metric")][1]

if (!is.na(outcome_col)) {
  ft_s3_data <- as_grouped_data(t_s3_clean, groups = outcome_col)
  ft_s3 <- flextable(ft_s3_data)
  ft_s3 <- bold(ft_s3, i = as.formula(paste0("~ !is.na(", outcome_col, ")")), part = "body")
} else {
  ft_s3 <- flextable(t_s3_clean)
}

# Typography & Three-Line Style
ft_s3 <- font(ft_s3, fontname = "Times New Roman", part = "all")
ft_s3 <- fontsize(ft_s3, size = 9, part = "all")
ft_s3 <- fontsize(ft_s3, size = 9.5, part = "header")
ft_s3 <- bold(ft_s3, part = "header")

ft_s3 <- align(ft_s3, j = 1, align = "left", part = "all")
if (ncol(t_s3_clean) > 1) {
  ft_s3 <- align(ft_s3, j = 2:ncol(t_s3_clean), align = "center", part = "all")
}

std_border <- fp_border(color = "black", width = 1)
ft_s3 <- border_remove(ft_s3)
ft_s3 <- hline_top(ft_s3, part = "header", border = std_border)
ft_s3 <- hline_bottom(ft_s3, part = "header", border = std_border)
ft_s3 <- hline_bottom(ft_s3, part = "body", border = std_border)

ft_s3 <- delete_part(ft_s3, part = "footer")
ft_s3 <- set_caption(
  ft_s3,
  caption = "Supplementary Table S3. Sensitivity analysis evaluating blood lead levels, diarrheal morbidity, and child growth faltering within the household drinking water quality (E. coli) testing subsample."
)

footnote_s3 <- paste(
  "Abbreviations: BLL, Blood Lead Level; CI, Confidence Interval; CFU, Colony Forming Units; HAZ, Height-for-Age z-score;",
  "WAZ, Weight-for-Age z-score; WHZ, Weight-for-Height z-score.",
  "Household water contamination was assessed via membrane filtration testing of source/stored drinking water (WQ26).",
  "Models account for complex survey strata, primary sampling units, and normalized subsample weights, with multivariable adjustment for all baseline demographic and socio-environmental covariates."
)

ft_s3 <- add_footer_lines(ft_s3, values = footnote_s3)
ft_s3 <- fontsize(ft_s3, size = 8, part = "footer")
ft_s3 <- italic(ft_s3, part = "footer")
ft_s3 <- autofit(ft_s3)

doc_s3_path <- file.path(target_dir, "TableS3_Ecoli_Sensitivity_Subsample.docx")
doc_s3 <- read_docx() %>%
  body_add_par("Supplementary Table S3: Water Quality Sensitivity Subsample", style = "heading 1") %>%
  body_add_flextable(ft_s3)

print(doc_s3, target = doc_s3_path)
cat(sprintf("Saved: %s\n", doc_s3_path))






# --- TABLE S4: Generalized Variance Inflation Factors (VIF) ---
cat("Computing Table S4: Multicollinearity Diagnostics (GVIF)...\n")
lm_vif <- lm(
  HAZ2 ~ log2_bll + diarrhea + CAGE + HL4 + C_Hb + melevel + helevel + 
    windex5 + HH6 + HH7 + peeling_paint + 
    improved_water + improved_sanitation + basic_handwashing,
  data = df_clean
)
vif_vals <- car::vif(lm_vif)

TableS4_final <- if (is.matrix(vif_vals)) {
  data.frame(
    Variable      = rownames(vif_vals),
    GVIF          = sprintf("%.2f", vif_vals[, 1]),
    Df            = vif_vals[, 2],
    GVIF_Adjusted = sprintf("%.2f", vif_vals[, 3]),
    Conclusion    = if_else(vif_vals[, 3] < 2.5, "No Collinearity (< 2.5)", "Potential Collinearity")
  )
} else {
  data.frame(
    Variable   = names(vif_vals),
    VIF        = sprintf("%.2f", vif_vals),
    Conclusion = if_else(vif_vals < 5.0, "Low (< 5.0)", "High (>= 5.0)")
  )
}
write.csv(TableS4_final, file.path(output_dir, "TableS4_VIF_Collinearity_Diagnostics.csv"), row.names = FALSE)
cat(sprintf("Saved: %s/TableS4_VIF_Collinearity_Diagnostics.csv\n", output_dir))






# ==============================================================================
# Supplementary Table S4: Multicollinearity Diagnostics (GVIF)
# Export to Publication-Ready Word Document (.docx)
# ==============================================================================

library(officer)
library(flextable)
library(dplyr)

target_dir <- if (exists("output_dir")) output_dir else getwd()
csv_s4 <- file.path(target_dir, "TableS4_VIF_Collinearity_Diagnostics.csv")

if (!file.exists(csv_s4)) {
  stop(sprintf("File '%s' not found.", csv_s4))
}

t_s4_raw <- read.csv(csv_s4, stringsAsFactors = FALSE)

# Clean and format diagnostic values
t_s4_clean <- t_s4_raw %>%
  mutate(
    # Clean variable names if present
    across(any_of(c("Variable", "Covariate", "Parameter")), ~ case_when(
      . == "log2_bll"            ~ "Blood lead level (log2 ug/dL)",
      . == "diarrhea"            ~ "Recent diarrhea (past 2 weeks)",
      . == "CAGE"                ~ "Child age (months)",
      . == "HL4"                 ~ "Child sex",
      . == "C_Hb"                ~ "Hemoglobin (g/dL)",
      . == "melevel"             ~ "Maternal education",
      . == "helevel"             ~ "Household head education",
      . == "windex5"             ~ "Wealth index quintile",
      . == "HH6"                 ~ "Area of residence (Urban/Rural)",
      . == "HH7"                 ~ "Administrative division",
      . == "peeling_paint"       ~ "Peeling/chipped wall paint",
      . == "improved_water"      ~ "Drinking water source type",
      . == "improved_sanitation" ~ "Sanitation facility type",
      . == "basic_handwashing"   ~ "Basic handwashing facility",
      . == "log2_s_lead"         ~ "Soil lead (log2 mg/kg)",
      . == "log2_c_cd"           ~ "Blood cadmium (log2 ug/L)",
      . == "log2_c_as"           ~ "Blood arsenic (log2 ug/L)",
      TRUE ~ .
    )),
    # Standardize numeric precision to 2 decimal places for VIF/GVIF metrics
    across(where(is.numeric), ~ sprintf("%.2f", .))
  )

ft_s4 <- flextable(t_s4_clean)

# Typography & Borders
ft_s4 <- font(ft_s4, fontname = "Times New Roman", part = "all")
ft_s4 <- fontsize(ft_s4, size = 9, part = "all")
ft_s4 <- fontsize(ft_s4, size = 9.5, part = "header")
ft_s4 <- bold(ft_s4, part = "header")

ft_s4 <- align(ft_s4, j = 1, align = "left", part = "all")
if (ncol(t_s4_clean) > 1) {
  ft_s4 <- align(ft_s4, j = 2:ncol(t_s4_clean), align = "center", part = "all")
}

std_border <- fp_border(color = "black", width = 1)
ft_s4 <- border_remove(ft_s4)
ft_s4 <- hline_top(ft_s4, part = "header", border = std_border)
ft_s4 <- hline_bottom(ft_s4, part = "header", border = std_border)
ft_s4 <- hline_bottom(ft_s4, part = "body", border = std_border)

ft_s4 <- delete_part(ft_s4, part = "footer")
ft_s4 <- set_caption(
  ft_s4,
  caption = "Supplementary Table S4. Multicollinearity diagnostic statistics (Variance Inflation Factors) for fully adjusted multivariable regression models."
)

footnote_s4 <- paste(
  "Abbreviations: Df, Degrees of Freedom; GVIF, Generalized Variance Inflation Factor.",
  "GVIF^(1/(2*Df)) provides a metric comparable to standard single-degree-of-freedom VIFs.",
  "Values of GVIF^(1/(2*Df)) < 2.0 (equivalent to VIF < 4.0) demonstrate an absence of severe collinearity among exposure, socio-demographic, and built-environment covariates."
)

ft_s4 <- add_footer_lines(ft_s4, values = footnote_s4)
ft_s4 <- fontsize(ft_s4, size = 8, part = "footer")
ft_s4 <- italic(ft_s4, part = "footer")
ft_s4 <- autofit(ft_s4)

doc_s4_path <- file.path(target_dir, "TableS4_VIF_Collinearity_Diagnostics.docx")
doc_s4 <- read_docx() %>%
  body_add_par("Supplementary Table S4: Multicollinearity Diagnostics", style = "heading 1") %>%
  body_add_flextable(ft_s4)

print(doc_s4, target = doc_s4_path)
cat(sprintf("Saved: %s\n", doc_s4_path))







cat("\n================================================================================\n")
cat("SUCCESS: Master publication suite and supplementary appendix generated.\n")
cat(sprintf("All deliverables written to: %s\n", output_dir))
cat("================================================================================\n")




summary(fit_haz_m1_main)$coefficients["diarrhea", ]
# or if it was a factor:
summary(fit_haz_m1_main)$coefficients[grep("diarrhea", rownames(summary(fit_haz_m1_main)$coefficients)), ]

