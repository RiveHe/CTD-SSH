# CTD + Copernicus SSH Bottom-Pressure Reconstruction (HOT Project)

This repository contains scripts to **reconstruct ocean-only bottom pressure** by combining:
- **HOT Project CTD** (deep casts + shallow casts; temperature/salinity profiles → steric height / density structure)
- **Copernicus satellite Sea Surface Height (SSH)**

We compare the reconstructed bottom pressure against observed bottom pressure (when available) to quantify **how deep the CTD needs to be** (minimum cast depth) to best reproduce bottom-pressure variability. The resulting “ocean water pressure” estimate can be used to:
- interpret bottom pressure in both **geophysics and oceanography**
- provide an **oceanographic correction** to **denoise bottom pressure sensors** (remove steric + barotropic components captured by CTD+SSH)

---

## Scientific motivation

Bottom pressure measured on the seafloor includes signals from:
- oceanographic variability (water mass/density changes + sea-level changes)
- geophysical processes (e.g., deformation, loading, etc., depending on application)

This project focuses on reconstructing the **ocean-only pressure contribution** using CTD and SSH, and then evaluating **CTD depth requirements** for accurate reconstruction. This is especially useful when:
- CTD profiles are not always full-depth
- you want a principled way to decide whether a cast is “deep enough” for bottom-pressure studies
- you want to estimate and remove oceanographic variability from BPR observations

---

## Data sources

### HOT Project CTD
- Deep casts and shallow casts
- Variables typically used: Temperature, Salinity, Pressure (and derived density / specific volume anomaly)

### Copernicus SSH
- Gridded satellite altimetry sea surface height (SSH) products
- Used to capture sea-level / barotropic variability that complements CTD steric effects

> **Note:** Raw CTD and Copernicus datasets are not included in this repo. Configure local paths or use the included download scripts (if present).

---

## Method overview

### 1) Download and organize data
- Fetch HOT CTD deep and shallow casts for the target period(s)
- Fetch Copernicus SSH over the HOT region and time range

### 2) CTD processing
- QC and standardize CTD profiles
- Interpolate profiles to a common vertical grid (optional)
- Compute density-related quantities needed for steric integration
- Create “truncated” versions of each profile at multiple candidate maximum depths to test sensitivity

### 3) Compute ocean-only pressure contribution at the bottom
We combine:
- **Steric component** from CTD (density structure)
- **Surface height component** from SSH (sea-level variability)

to estimate the time-varying **ocean water pressure at the seafloor** (oceanographic contribution).

### 4) Depth requirement experiment
For each candidate CTD maximum depth (e.g., 500 m, 1000 m, 2000 m, …):
- reconstruct bottom pressure using the truncated CTD + SSH
- compare to a reference:
  - either **full-depth CTD reconstruction** (if available), or
  - **observed bottom pressure** (if you have BPR at/near HOT)
- quantify performance using metrics such as:
  - correlation (R)
  - RMS error / variance explained
  - spectral coherence (optional)

### 5) Output
- A recommended **minimum CTD depth** for reliable reconstruction at the site/time window
- Time series and figures showing reconstruction skill vs CTD depth
- Reconstructed ocean-only pressure that can be used for **BPR denoising**

---


