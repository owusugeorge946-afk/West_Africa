# Tano River Basin hydro-climatic reproducibility archive

Reproducibility code for **Remote Sensing and Machine Learning-Based Assessment of Hydro-Climatic Variability and Drought-Flood Hazard in the Transboundary Tano River Basin, West Africa**.

## Scope
Study period: 2000–2024. Core variables: CHIRPS annual rainfall, Rx1day, Rx5day, ETCCDI CDD; Landsat NDVI and LST; DCI and ERPI; validation; Spearman correlation; Random Forest descriptive reconstruction.

## Important revision controls
- CDD = maximum annual run of consecutive days with precipitation < 1 mm, not the annual number of dry days.
- Trend inference uses a Hamed–Rao variance-corrected Mann–Kendall test plus Sen's slope. The implementation is explicit in `02_scripts/02_modified_MK_SenSlope.py`.
- The pooled NDVI trend must be interpreted together with sensor-era sensitivity (2000–2012 vs 2013–2024). It must not be described as robust sensor-independent greening unless harmonized sensor analysis supports that conclusion.
- DCI/ERPI Random Forest models are descriptive reconstructions of constructed indices, not causal driver models.
- Monthly CHIRPS–IMERG validation does not directly validate daily Rx1day, Rx5day, or CDD.
- Spatial maps must be regenerated from analytical rasters. Do not redraw or cosmetically replace scientific maps.

## Required local input
Place the authoritative annual table at:
`01_data/Tano_Final_Annual_Hydroclimatic_Results.csv`

Required columns:
`year,rainfall,rx1day,rx5day,cdd,ndvi,lst`
Optional existing columns may include `dci_entropy,erpi_entropy,landsat_scenes,valid_landsat_pixels`.

For validation scripts, place paired monthly CHIRPS/IMERG values at:
`01_data/validation/chirps_imerg_monthly.csv` with columns `date,chirps,imerg`.

## Run order
1. `python 02_scripts/01_descriptive_statistics.py`
2. `python 02_scripts/02_modified_MK_SenSlope.py`
3. `python 02_scripts/03_NDVI_sensor_sensitivity.py`
4. `python 02_scripts/04_correlation_analysis.py`
5. `python 02_scripts/05_DCI_ERPI.py`
6. `python 02_scripts/06_DCI_NDVI_sensitivity.py`
7. `python 02_scripts/07_random_forest_LOOCV.py`
8. `python 02_scripts/08_validation_metrics.py`
9. Run `02_scripts/10_GEE_spatial_analysis.js` in Google Earth Engine and export the analytical rasters/validation table.
10. Generate manuscript figures only after numerical and spatial validation checks pass.

## Earth Engine
Project: `ee-george946`
Basin asset: `projects/ee-george946/assets/Tano_River`
Selected map years: 2001, 2010, 2020, 2023.

## Reproducibility rule
The repository intentionally does not fabricate missing source observations or manuscript outputs. Scripts fail clearly when required inputs are absent. Commit final numerical outputs only after reproducing them from the authoritative data.

## License/citation
Add the final manuscript citation and Zenodo DOI after the revised archive is deposited.
