# Endolaserless Analysis

## Repository Summary

This repository contains analysis code for the **Laserless Study (NCT02976012)**, a prospective randomized trial of aflibercept monotherapy after endolaserless vitrectomy for proliferative diabetic retinopathy (PDR)-related vitreous hemorrhage.

We processed and analyzed ultra-widefield fluorescein angiography (UWFA) images graded by the Wisconsin Reading Center to quantify the non-perfusion index (NPI) over 3 years in eyes randomized to q8-week versus q16-week intravitreal aflibercept injection schedules.

## Key Findings

The analysis code supports the following key findings:

- **NPI progression**: Retinal nonperfusion increased significantly over 3 years despite anti-VEGF therapy
- **Treatment frequency effect**: Eyes treated q16 weeks showed significant NPI progression, while q8-week eyes did not
- **Rescue injections**: The q16 group required more unscheduled (PRN) rescue aflibercept injections, yet still progressed
- **Association with injection burden**: In the q16 group, higher cumulative injections were positively associated with greater NPI burden; no such association was seen in the q8 group

The analysis highlights that more frequent, earlier anti-VEGF dosing may slow progression of retinal ischemia compared to extended dosing, with implications for future long-acting anti-VEGF therapies.

## Repository Structure

```
endolaserless_analysis/
├── scripts/
│   ├── count_prn_injections.R          # PRN injection analysis and counting
│   ├── endolaserless_analysis-2.R      # Main NPI analysis script
│   └── helper_scripts.R                # Utility functions and helper code
├── npi_project/
│   ├── prism_output/                   # Prism-generated figures and visualizations
│   │   ├── Final/                      # Final publication figures
│   │   └── [other analysis plots]
│   └── output/                         # Analysis output files (not committed)
├── README.md                           # This file
└── .gitignore                          # Git ignore rules
```

## Prism Output Figures

All visualizations were created using **GraphPad Prism**. The following key figures are included in this repository:

### Final Publication Figures

**Figure 1. Combined q8 and q16 treatment schedules data shows that NPI increase over time despite intravitreal aflibercept monotherapy.** Bar graphs showing nonperfusion index (NPI) changes over time. Blue squares represent q8 (N=9) group patients. Orange circles represent q16 (N=14) group patients. Statistics using mixed effects model with Dunnett posttest and comparisons to week 4. ***p<0.001, ****p<0.0001. Error bars, mean ± SEM.

**Figure 2. In q16 treatment schedule but not q8, NPI increases significantly from baseline.** (A) Bar graphs showing baseline non-perfusion index (NPI) at week 4 for q8 (N=9) and q16 (N=14) treatment schedules. Unpaired Student's T-test. p=ns. (B) Timeline graph showing raw NPI and (C) bar graph showing the change in NPI from week 4 baseline at each follow-up appointment stratified by treatment group. Statistics using mixed effects model with Tukey posttest and all multiple comparisons. Error bars, mean ± SEM. ****p<0.0001.

**Figure 3. Decreased treatment load and later PRN injections are insufficient to prevent RNP progression.** Box-and-whisker plots displaying the (A) cumulative total and (B) cumulative PRN (non-scheduled) intravitreal aflibercept injections received at each timepoint stratified by q8 (N=9) and q16 (N=14) treatment groups. Box-and-whisker plot shows the median (line inside the box), interquartile range (box), and whiskers representing 1.5 times the interquartile range. Dotted horizontal line in (A) indicates 10 cumulative injections. Statistics using two-way repeated measures ANOVA with (A) Holm-Sidak posttest and comparisons between treatment schedules at each timepoint and (B) Tukey posttest and all multiple comparisons. *p<0.05, **p<0.01, ***p<0.001, ****p<0.0001. (C) Scatter plot showing the total number of injections versus the area under the curve (AUC) of nonperfusion index (NPI) from baseline week 4 to week 152 for each subject stratified by treatment group. Statistics using simple linear regression with p-value representing the result of F-test for non-zero slope. *p<0.05.

### Supplementary Figures

**Supplementary Figure 1. Combined q8 and q16 treatment schedules data shows that NPA increase over time despite intravitreal aflibercept monotherapy.** Bar graphs showing nonperfusion area (NPA) changes over time. Blue squares represent q8 (N=9) group patients. Orange circles represent q16 (N=14) group patients. Statistics using mixed effects model with Dunnett posttest and comparisons to week 4. *p<0.05, **p<0.01. Error bars, mean ± SEM.

## Usage

1. **Setup**: Ensure R is installed with required packages (see script headers for dependencies)
2. **Data Preparation**: Place required data files in the appropriate directories
3. **Run Analysis**: Execute scripts in order:
   - `endolaserless_analysis-2.R` - Main analysis pipeline
   - `count_prn_injections.R` - PRN injection analysis
4. **Output**: Results are saved to `npi_project/output/` directory

## Dependencies

The analysis scripts require several R packages including:
- `tidyverse` for data manipulation
- `lme4` and `lmerTest` for mixed effects modeling
- `ggplot2` and `ggprism` for plotting
- `openxlsx` for Excel file handling
- Additional statistical packages as specified in script headers

## Citation

If you use this analysis code, please cite the original study:
**Laserless Study (NCT02976012)** - A prospective randomized trial of aflibercept monotherapy after endolaserless vitrectomy for proliferative diabetic retinopathy (PDR)-related vitreous hemorrhage.

## License

This repository contains research analysis code. Please contact the authors for usage permissions.
