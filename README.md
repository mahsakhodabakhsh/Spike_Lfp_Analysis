# Spike & LFP Neural Data Analysis

MATLAB pipeline for analyzing spike and local field potential (LFP) recordings from a multi-electrode neurophysiology dataset, covering spike-train statistics, information-theoretic discriminability, population decoding, and spectral/cross-frequency analysis.

## Data

Recordings from 16 electrodes across multiple task conditions (angles), with spike times and LFP for each. Data was provided as part of the **IPM Neural Data Analysis Summer School**.
Link_Data: https://drive.google.com/file/d/1HwVHaaTSrYtI5wJVzQbPxAtan2sUy64Y/view?usp=sharing

## Helper functions

Analysis relies on a set of `ndass_*` helper functions (smoothing, mutual information, ROC, SVM decoding, wavelet transform, PAC, spike-field locking, etc.), written by **Ehsan Rezayat** for the IPM Neural Data Analysis Summer School. algorithms implemented in them are cited in-code to their source papers (e.g. Tort et al., 2010 for phase-amplitude coupling).
Link_Functions: https://drive.google.com/file/d/1gvdD2CFJYnI04hfbkXeW544MvWbZHi9O/view?usp=drive_link

## What the code does

**Spike analysis**
- Loads spike trains for 27 units across 16 stimulus conditions
- Raster plots: per-trial spike rasters for 4 example neurons, per condition
- PSTH: smoothed firing rate over time per condition, per neuron
- Mutual information & ROC: quantifies how well single-neuron activity discriminates the target condition from others, computed in two time windows (1000–1500 ms and 2700–3200 ms) and compared against a within-condition control
- Population SVM decoding: trains a classifier on population activity in sliding time windows to decode stimulus condition, plots decoding performance over time

**LFP analysis**
- Loads LFP from all 16 electrodes
- Preprocessing: 50 Hz notch filter, artifact rejection, z-score normalization
- Multitaper spectral power (Chronux): baseline-normalized power spectra per condition, for the same two time windows as above
- Phase-amplitude coupling (Tort et al., 2010): cross-frequency coupling between low-frequency phase (4–8 Hz) and high-frequency power (30–60 Hz), within and across two electrodes, with shuffle-corrected significance
- Spike-field coherence: phase-locking between spike times and LFP phase across frequency bands, for both time windows

## Run

Requires MATLAB, Chronux toolbox, and the `ndass_*` functions from the summer school materials on the path.
