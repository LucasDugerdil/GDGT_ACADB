# Cofounding factors mitigating the brGDGT-based paleothermometer in drylands

## Overview
This GitHub project is associated to the article "*BrGDGT-based paleothermometer in drylands: the necessity to constrain aridity and salinity as confounding factors to ensure the robustness of calibrations*" published in *Biogeosciences* in 2026 (Dugerdil et al., 2026a).

**Author**: **Lucas Dugerdil**<sup>1,2</sup>

**Affiliations**:
1. Univ. Lyon, ENS de Lyon, Université Lyon 1, CNRS, UMR 5276 LGL-TPE, F-69364, Lyon, France1
2. Université de Montpellier, CNRS, IRD, EPHE, UMR 5554 ISEM, Montpellier, France

**ORCID**: [0000-0003-0266-564X](https://orcid.org/0000-0003-0266-564X)

**Funding**: ANR, Grant [ANR‐22‐CE27‐0018](https://anr.fr/Project-ANR-22-CE27-0018) (STEPABILITY), Sébastien Joannin and ANR, Grant [ANR-20-FRAL-0006](https://anr.fr/Project-ANR-20-FRAL-0006) (KUR(A)GAN), Giulio Palumbi.

**Open Access**:

<table width="100%">
  <tr>
    <td width="33.33%" align="left" valign="middle"><strong>Research article</strong></td>
    <td width="33.33%" align="left" valign="middle"><strong>Published release</strong></td>
    <td width="33.33%" align="left" valign="middle"><strong>Data repository</strong></td>
  </tr>
  <tr>
    <td width="33.33%" align="left" valign="middle">

[![Static Badge](https://img.shields.io/badge/DOI-10.5194%2Fbg--23--1013--2026-yellow)](https://doi.org/10.5194/bg-23-1013-2026)

</td>
    <td width="33.33%" align="left" valign="middle">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18429660.svg)](https://doi.org/10.5281/zenodo.18429660)

</td>
    <td width="33.33%" align="left" valign="middle">

[![Static Badge](https://img.shields.io/badge/DOI-10.1594%2FPANGAEA.983391-green)](https://doi.org/10.1594/PANGAEA.983391)

</td>
  </tr>
</table>

## Description
This R script is the full script develloped in the publication.
It is usefull to verify the replicability of this study.
The script could be modified for calibrations in other study areas. 

## How to install/run the ACADB brGDGT calibrations?
1. Install [R](https://larmarange.github.io/analyse-R/installation-de-R-et-RStudio.html)
2. It is easier to use [Rstudio](https://posit.co/downloads/)
3. Download this GitHub repository from ZIP file (by clicing on the green button `<> Code` beyond.

## To run the full code
	- Many `R` package are necessary to run this script. Check if all required package are already installed in case of error appearing during the script process.
	- Open the `GDGT_ACADB.Rproj` file in Rstudio
	- In the section `#### Select figures / tables to plot ####` change the boolean value to `TRUE` or `FALSE` to plot the selected figure or table
