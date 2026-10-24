# Blood Lead Exposure, Recent Diarrhea, and Child Undernutrition in Bangladesh

### A National Cross-Sectional Analysis of Joint Associations and Interaction Among Children Aged 12–59 Months

**Authors:**  
- **Omar Faroque** ([GitHub](https://github.com/Omarrfaroqueiu79))  
- **Mst. Homayra Jannaty**  
- **Suman Biswas*** (Corresponding: suman@stat.iu.ac.bd)  

*Department of Statistics and Data Science, Islamic University, Kushtia-7003, Bangladesh*

---

## 📌 Overview

This repository contains the replication code, complex survey design specifications, multivariable regression modeling routines, and publication-ready figures for the research study evaluating individual and joint associations of venous blood lead levels (BLL) and recent diarrheal morbidity with child growth faltering in Bangladesh.

The analysis is based on nationally representative microdata from the **Bangladesh Multiple Indicator Cluster Survey (MICS) 2025** ( = 10,507$).

---

## 📊 Key Findings

1. **High Pediatric Lead Burden:**
   - Survey-weighted mean BLL was **5.22 µg/dL** (SD 3.78).
   - **64.3%** of Bangladeshi children aged 12–59 months exceeded the US CDC Blood Lead Reference Value of **3.5 µg/dL**.
2. **Linear & Ponderal Growth Faltering:**
   - Each doubling of continuous BLL ($\log_2$-transformed) was significantly associated with lower height-for-age z-score (**HAZ**: $\beta = -0.095$, \%\text{ CI}: -0.132\text{ to } -0.059$,  < 0.001$) and lower weight-for-age z-score (**WAZ**: $\beta = -0.072$, \%\text{ CI}: -0.104\text{ to } -0.040$,  < 0.001$).
   - Associations remained robust after extensive adjustment for socioeconomic status, maternal education, WASH infrastructure, child hemoglobin, and environmental co-pollutants (soil lead, blood cadmium, blood arsenic).
3. **Joint Exposure & Additive Interaction:**
   - Children with joint exposure to high BLL ($\ge 3.5\ \mu\text{g/dL}$) and recent diarrhea exhibited the highest prevalence of acute wasting (**aPR = 1.49**, \%\text{ CI}: 1.14–1.95$) and underweight (**aPR = 1.29**, \%\text{ CI}: 1.07–1.56$).
   - Multiplicative and additive interaction metrics (RERI, AP, Synergy Index $) showed no statistically detectable departure from additivity, indicating that concurrent exposure represents a compounding double burden of risk rather than biological synergy.

---

## 📁 Repository Structure

`
├── figures/
│   ├── Figure1_Divisional_BLL_Distribution.png
│   ├── Figure2_HAZ_Compounding_Interaction.png
│   ├── Figure3_ThreePanel_Malnutrition.png
│   ├── FigureS1_STROBE_Flowchart.png
│   ├── FigureS2_Conceptual_DAG.png
│   └── FigureS3_Splines.png
├── tables/
│   ├── Table1_Weighted_Characteristics.csv / .docx
│   ├── Table2_Linear_Interactions.csv / .docx
│   ├── Table3_aPR_Results.csv / .docx
│   ├── Table4_Additive_Interaction.csv / .docx
│   ├── TableS1_Crude_vs_Adjusted.csv / .docx
│   ├── TableS2_Logistic_aOR_Comparison.csv / .docx
│   ├── TableS3_Ecoli_Sensitivity_Subsample.csv / .docx
│   └── TableS4_VIF_Collinearity_Diagnostics.csv / .docx
├── scripts/
│   ├── Field_Trip_Analysis.R       # Primary master survey analysis script
│   ├── Figure3_Three_Panel.R       # Standardized 3-panel malnutrition bar chart
│   ├── Reviewer_Revisions.R        # Sensitivity models and additive metrics
│   └── Fix_Fig3_Collision.R        # Clean label positioning routines
├── manuscript/
│   ├── BLL_manuscript_Draft_Updated_24_10_26.docx  # Full updated manuscript
│   └── Manuscript_Revisions_Master_Action_Plan.docx # Redline revision guide
├── .gitignore
└── README.md
`

---

## 🛠 Software Requirements

The analyses were conducted using **R (version 4.4.2)** with the following core packages:
- survey (v4.4-2) — Complex survey design, weights, clustering, and Taylor linearization
- dplyr (v1.1.4) & 	idyr — Data wrangling and transformation
- ggplot2 & patchwork — High-resolution publication visualizations
- haven — SPSS .sav file importation
- lextable & officer — Publication-ready table generation

---

## 📜 Data Access

The raw microdata analyzed in this study belong to the **Bangladesh Multiple Indicator Cluster Survey (MICS) 2025**, conducted by the Bangladesh Bureau of Statistics (BBS) and UNICEF Bangladesh. The datasets are freely accessible to registered researchers upon submission of a brief research proposal at the official UNICEF MICS portal:  
🔗 [https://mics.unicef.org/surveys](https://mics.unicef.org/surveys)

---

## 📄 License & Citation

This repository is maintained for open and reproducible scientific research. If you use the analytical scripts or figures, please cite the corresponding publication:

> **Faroque O, Jannaty MH, Biswas S.** *Blood Lead Exposure, Recent Diarrhea, and Child Undernutrition in Bangladesh: A National Cross-Sectional Analysis of Joint Associations and Interaction Among Children Aged 12–59 Months.* (Under Review, 2026).
