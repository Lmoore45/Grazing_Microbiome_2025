# Microbes, Minerals, and Managment: What Controls Mineral-Associated Organic Carbon in Semi-Arid Rangelands?
*Companion code repository for manuscript submission*

## Manuscript Abstract  
*Insert abstract text here once finalized.*  
This repository supports the analyses and figures described in the manuscript, providing reproducible workflows for data processing, sPLS/PLS modeling, and supplementary figure generation using microbial and environmental data.

Due to rancher privacy agreements, precise geographic coordinates used in the manuscript analyses are not included in this repository. The version of the code shared here uses county-level centroids as placeholders.
As a result, certain figures (e.g., climate figures, PLS, RDA) may differ slightly from those in the published manuscript. However, all primary findings and interpreations remain consistent. 

---

## Repository Structure and How to Use 

The code is organized by numbered sections. To fully reproduce the analyses, follow the instructions below for obtaining and using the appropriate input files and scripts.

### 1.0 – Data Preparation  
This section loads metadata and OTU tables.

You will need to download the following input files from the repository or dataset archive:
- `1.1_Metadata.csv` – Core metadata table  
- `1.2_16S_OTU_table.csv.zip` – Rarefied 16S OTU table (zipped)
- `1.3_ITS_OTU_table.txt` – Raw ITS OTU table  

---

### 2.0 – Climate Data  
This section uses an API script to pull climate variables for each site.

Due to rancher privacy agreements, precise site coordinates are not publicly provided. For reproducibility and transparency, the metadata included contains county centroid coordinates for each site, which differ from the actual coordinates used in the original analyses.
The script `2.0_climate_API_pull.R` demonstrates the API call.

---

### 6.1, 6.2 – Sparse Partial Least Squares (sPLS)  
These scripts perform sPLS regression on both 16S and ITS data to identify microbial taxa most predictive of MAOC.

To run these scripts, be sure to also download and include:
- `6.0_VIP.R` – Custom function script used to calculate VIP scores from sPLS models.

---

### 7.0 – Partial Least Squares (PLS)  
Similar to sPLS but using all features (no sparsity constraint). Also requires:
- `6.0_VIP.R` – Shared VIP function used in both sPLS and PLS scripts.
- Note: This analysis originally used precise goegraphic coordinates, but the version shared here substitues county-level centroids to protect rancher privacy. As a result, minor differences may appear in spatially influenced model outputs, though overall patterns remain consisten. 

---

### S1.0, S2.1, S2.2 – Supplementary Figures  
These scripts generate additional exploratory figures and supporting plots used in the supplemental materials of the manuscript.

---

## Software Requirements  
These analyses were developed and tested in R. Major packages used include:
- `tidyverse`
- `mixOmics`
- `pls` 
- `compositions`
- `ggplot2`
- `corrplot`
- `Hmisc`


---

## License  
This repository is released under the [MIT License](LICENSE). You are free to use, modify, and share this code with attribution.

---

## Questions or Contributions?  
For inquiries related to the manuscript or code, please reach out to Laura.Moore@colostate.edu, or open an issue here in the repository.

