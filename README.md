# AO-FACED Manuscript Analysis

MATLAB and Python analysis code used to process and compare adaptive-optics (AO) and system-AO/NoAO FACED microscopy data. The repository contains manuscript figure scripts for 3D bead characterization, in vivo 3D imaging, vessel structure and blood flow, glutamate responses, and deformable-mirror characterization.

## Citation

If you use this code, please cite:

Zhu, J., Natan, R. G., Zhong, J., Kang, I., & Ji, N. (2026). *In vivo aberration measurement and correction for ultrafast FACED two-photon fluorescence microscopy of the brain*. bioRxiv. https://doi.org/10.64898/2026.02.06.704504

## Main workflows

### 1. 3D bead-stack comparison

- `Image_compare_beads.m`: registers system-AO and NoAO fluorescent-bead volumes, segments individual beads, extracts axial intensity profiles, estimates bead FWHM, and maps image-quality changes across the FACED field of view.

### 2. In vivo 3D comparison

1. Run `caiman_register_tiff_stack.py` to motion-correct TIFF stacks with CaImAn.
2. Run `StackCompare_Motion_corr_M507_data1.m` to align system-AO and NoAO volumes, crop matched regions, and prepare structural projections for comparison.

### 3. Vessel images and blood flow

1. Run `VesselFig_Image_compare_VS_v8.m` to generate vessel kymographs, select ROIs, flatten vessel bands, and average them into line profiles.
2. Run `VesselFig_Line_profile_compare.m` to compare system-AO and NoAO plasma intensity and peak-valley distinguishability across ROIs and fields of view.
3. Run `VesselFig_flow_quant_main_v4.m` to quantify flow from accepted system-AO and NoAO kymographs, select an AO shift parameter, save results, and calculate quality-control metrics.
4. Run `VesselFig_plot_flow_compare.m` to compare fit quality and summarize measured velocities.

### 4. Glutamate imaging

- `GlutamateFig_iGlu4FlashFig_v4.m`: performs trial alignment and averaging, extracts ROI traces, and compares system-AO and NoAO glutamate responses.
- `GlutamateFig_iGlu4FlashFig_v5_sta_AllFOV.m`: calculates population-level glutamate-response statistics across fields of view.

## Other analysis

- `LP_contrast.m`: line-profile, peak-valley, and background-corrected contrast measurements.
- `readDM.m`: deformable-mirror maps, sample/system wavefront subtraction, RMS error, tilt montages, and FACED modulation analysis.

## Core helper functions

- `flow_quant_func.m`: LSPIV-based flow quantification and velocity-trace generation.
- `imgflat.m`: interpolation and flattening of kymograph ROIs.
- `SliceXZ.m`: extraction of kymographs along selected vessel paths.
- `VesselLineCst_v2.m`: intensity and line-profile measurements.
- `localPeakValleyDistinguishability_v2.m`: adaptive peak-valley distinguishability calculation.
- `plot_compare_v3.m`: paired population plots and statistics.
- `ReadImageJROI.m`: import of ImageJ ROI files.
- `saveFigs_mod.m`: figure export helper.

## Requirements

- MATLAB with Image Processing Toolbox and Statistics and Machine Learning Toolbox.
- Parallel Computing Toolbox for the parallel shift search in `VesselFig_flow_quant_main_v4.m`.
- Python with CaImAn for `caiman_register_tiff_stack.py`.
- ImageJ/Fiji for creating the ROI ZIP files used by the vessel workflow.
- Additional laboratory helper functions may need to be added to the MATLAB path, including plotting, scale-bar, motion-correction, and ROI utilities referenced by the scripts.

## Data setup

The manuscript scripts contain experiment-specific absolute paths and ROI selections. Before running a script:

1. Replace its input and output paths with local paths.
2. Confirm that the expected MAT, TIFF, FIG, CSV, and ImageJ ROI files are present.
3. Add this repository and any external helper-function directories to the MATLAB path.
4. Review dataset selectors, ROI indices, crop ranges, and display limits near the top of the script.

Large raw and processed imaging datasets are not included in this repository.

## Notes

Many scripts are interactive and may open figures, request ROI selection, or display MATLAB menus. Run scripts section by section when adapting them to a new dataset. Generated MAT, TIFF, PDF, FIG, and CSV outputs are written to the configured experiment folders.
