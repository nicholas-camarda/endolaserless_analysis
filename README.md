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
│   ├── count_prn_injections.R          # PRN injection analysis
│   ├── endolaserless_analysis-2.R      # Main NPI analysis pipeline
│   └── helper_scripts.R                # Utility functions
├── data/
│   └── old/                            # Input data files
├── processed_data/
│   └── npi_project/                    # Processed analysis data
├── output/
│   └── npi_project/
│       └── prism_output/               # Prism-generated figures
└── README.md
```

## Analysis Scripts

### `endolaserless_analysis-2.R`
Main analysis pipeline that:
- Processes Wisconsin Reading Center NPI data
- Performs mixed-effects modeling
- Generates Prism-compatible output files
- Analyzes complete cases and full datasets
- Calculates area-under-curve (AUC) metrics

### `count_prn_injections.R`
Analyzes PRN (pro re nata) injection patterns:
- Tracks rescue injection timing and frequency
- Correlates injection burden with NPI progression
- Generates injection summary statistics

### `helper_scripts.R`
Utility functions for:
- AUC calculations
- Data processing and formatting
- Statistical analysis helpers
- Prism data formatting

## Data Sources

- **Wisconsin Reading Center**: UWFA grading data for NPI quantification
- **Clinical Trial Data**: Injection schedules, PRN injections, and clinical outcomes
- **RedCap Database**: Additional clinical parameters and follow-up data

## Visualization

All statistical visualizations were generated using **GraphPad Prism** for publication-quality figures. The analysis scripts output data in Prism-compatible formats for final figure generation.

## Dependencies

Required R packages:
- tidyverse (dplyr, ggplot2, tidyr, purrr)
- readxl
- openxlsx
- lmerTest
- ggprism
- GetoptLong
- broom.mixed
- rstatix

## Usage

1. Install required R packages
2. Ensure input data files are in the correct locations
3. Run `endolaserless_analysis-2.R` for main analysis
4. Run `count_prn_injections.R` for PRN injection analysis
5. Import generated Excel files into Prism for final visualizations

## Citation

If using this analysis code, please cite the original Laserless Study publication and acknowledge the Wisconsin Reading Center for UWFA grading.

## License

This repository contains analysis code for research purposes. Please contact the study investigators for data access and collaboration opportunities.
