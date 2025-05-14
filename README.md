# Microbes, Minerals, and Managment: What Controls Mineral-Associated Organic Carbon in Semi-Arid Rangelands?
*Companion code repository for manuscript submission*

## Manuscript Abstract  
*Insert abstract text here once finalized.*  
This repository supports the analyses and figures described in the manuscript, providing reproducible workflows for data processing, sPLS/PLS modeling, and supplementary figure generation using microbial and environmental data.

---

## Repository Structure and How to Use 

The code is organized by numbered sections. To fully reproduce the analyses, follow the instructions below for obtaining and using the appropriate input files and scripts.

### 1.0 – Data Preparation  
This section loads metadata and OTU tables.

You will need to download the following input files from the repository or dataset archive:
- `1.1_Metadata.csv` – Core metadata table  
- `1.2_16S_OTU_table.csv` – Rarefied 16S OTU table  
- `1.3_ITS_OTU_table.txt` – Raw ITS OTU table  

---

### 2.0 – Climate Data  
This section uses an API script to pull climate variables for each site.

Due to rancher privacy agreements, precise geographic coordinates are not publicly provided.  
The script `2.0_climate_API_pull.R` demonstrates the API call but **will not execute fully without coordinate data**. You may substitute your own coordinates if replicating this analysis elsewhere.

---

### 6.1, 6.2 – Sparse Partial Least Squares (sPLS)  
These scripts perform sPLS regression on both 16S and ITS data to identify microbial taxa most predictive of MAOC.

To run these scripts, be sure to also download and include:
- `6.0_VIP.R` – Custom function script used to calculate VIP scores from sPLS models.

---

### 7.0 – Partial Least Squares (PLS)  
Similar to sPLS but using all features (no sparsity constraint). Also requires:
- `6.0_VIP.R` – Shared VIP function used in both sPLS and PLS scripts.

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

