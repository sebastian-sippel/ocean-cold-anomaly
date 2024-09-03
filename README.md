# ocean-cold-anomaly

## Name
Early 20th century cold bias in observational ocean surface temperature

[![DOI](https://zenodo.org/badge/817366226.svg)](https://zenodo.org/doi/10.5281/zenodo.13646027)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


## Description
This repository contains code to reproduce the findings of the paper 'Early 20th century cold bias in observational ocean surface temperature' (by Sippel, Kent, Meinshausen, Chan, Kadow, Neukom, Fischer, Humphrey, Rohde, de Vries, and Knutti), and this README describes the code and data for reproducibility. 

Code availability: We provide all code to reproduce the main results and all figures of the paper. Either browse/download files individually or clone the repository to your local machine (git clone https://github.com/sebastian-sippel/ocean-cold-anomaly.git).
Data availability: We provide the newly derived reconstructions as .txt files under data/04_final in this repository (described below). We provide all preprocessed intermediate data and post-processed data needed to reproduce the study's results and figures, all of which are available under https://data.iac.ethz.ch/Sippel_et_al_2024_ocean-cold-anomaly/. Due to storage and copyright constraints, original observational datasets and CMIP6 data have to be downloaded from their original sources (given in Data Availability Statement).

This README file contains:

1. Directory Structure 
2. Reproduction of figures
3. Statistical Learning code for reconstruction of global mean surface temperature (GMST).
4. Processing data from original sources

### 1. Directory Structure 
We describe the most important functions and directories.

	/figures
		Contains sub-directories with scripts to reproduce respective figures. We describe how to reproduce in detail below (2. Reproduction of figures).

	/code

		/_convenience
			Contains convenience functions and methods projection and plotting of spatial data.
		_functions_CMIP6.R / 
			Contains all functions to read CMIP6 data. 
		_attribution_hildreth-lu.R
			Attribution method based on Hildreth-Lu regression
		_ridge_fun_v2.R 
			Functions to generate statistical reconstructions.

	/data
		Contains the newly derived reconstructions, and all intermediate data to reproduce the study's results. Due to storage constraints, the data directories /00_ up to /_03 are available under: https://data.iac.ethz.ch/Sippel_et_al_2024_ocean-cold-anomaly/.

		/_final
			This directory contains the new, final reconstructions as .txt files. The reconstructions presented in the main text are labeled as GMST_reconstructions.txt, and the reconstructions without adding bias and uncertainty estimates for training are labeled as GMST_reconstructions_no-training-bias-unc.txt. Both files contain 13 columns with GMST reconstructions based on CRUTEM5, HadSST4, HadSST4-unadj, ClassNMAT, CoastalHybridSST, ERSSTv5, COBE-SST2, and BEST-Land, as described in the paper. For CRUTEM5 and HadSST4, the 2.5th and 97.5 percentiles (shown in Fig. 1) are also given.

		Further directories available from https://data.iac.ethz.ch/Sippel_et_al_2024_ocean-cold-anomaly/:

		/00_DATASET
			Should be set up to contain all original datasets used in the present study. These include the datasets described in the paper and are futher specified below:	
			/cmip6
				All cmip6 data files as described in Supplementary Table 4 for the variables tas and tos.
		/01_processed4train_CMIP6
			This folder contains .RData files, stratified by month for tos and tas, which are used in the statistical model training step.
		/01_processed4train_CRU
			This folder contains .RData files of the HadSST4 and CRUTEM5 observations, which are used as input to the trained statistical models. 
		/02_trained_models
			This folder contains .RData files that contain the trained statistical regression models (as well as evaluation metrics), stratified by predictor setup (tas, tos, tos_MR, etc.) and the reconstruction target metric (GMST, GMLSAT, etc.).
		/03_processed4analysis
			This folder contains .RData files that are ready for the analyses presented in the paper. 
		/03_processedCMIP6_reconstr
			Folder contains CMIP6 reconstructions (based on observational coverage and errors) shown as comparison to observations in Fig. 2 and Fig. 5.
		/03_processedOBS_reconstr
			Folder contains .RData files of the observations-based reconstructions. Final reconstructions as .txt files in /data/_final.
		/03_processed_CMIP6_4evaluation
			Folder contains .RData files of the CMIP6 files for evaluation, shown in Extended Data Figures 2-3 and Supplementary Figure 1.

	/scripts

		Read-in of observational datasets and masks:
			00_readCRU.R
				Processes HadSST4, CRUTEM5 and CLASSNMat observations to be used for reconstructions.

		CMIP6 read-in scripts:
			_00a_process_CMIP6_mon.R
				Processes CMIP6 data with R and CDO into /data/00_DATASET/
				Processes and reads CMIP6 data into /data/01_processed4train_CMIP6/
			_00b_process_CMIP6_GMLSAT-GMSST_v2.R
				Calculate blended data based on cmip6 and extracts target regional means. Amends files in /data/01_processed4train_CMIP6/
			
		Training GMST-reconstruction model scripts:
			01_train_hybrid36.R
			01_train_tas_land.R
			01_train_tas_sea_CLASSNMAT.R
			01_train_tos.R
			01_train_tos_MR.R

		Data processing scripts:
			02b_process4evaluation.R
			02b_process4evaluation_CLASSNMAT.R
			02b_process4evaluation_MR.R
			02b_process4evaluation_hybrid36.R

		Data loading scripts:
			03a_load_global_observations.R
			03b_load_global_obs_reconstruction.R
			03b_load_global_obs_reconstruction_ANCDATA.R
			03c_load_global_proxy_data.R
			03d_read_ocean2k.R

		Master scripts for plotting of figures:
			04a_master_load_reconstructions.R
			04b_master_read_paleo_reconstructions.R
			04c_compute_trends_4paleo-comparison.R

### 2. Reproduction of Figures:
To reproduce figures, first go to https://data.iac.ethz.ch/Sippel_et_al_2024_ocean-cold-anomaly/ and obtain all required data and save it to the /data directory. The above described /data directory structured needs to be preserved. Then, go to /figures. Each folder (e.g., /01_reconstruction) contains an .R file that is executable and produces the respective figure, using code and data from the /code, /scripts, and /data directories. Additionally all raw figures are already contained in each figure directory as .pdf files.

### 3. Statistical Learning Code for reconstruction of global mean surface temperature (GMST).
The statistical learning code for training the GMST reconstructions (and those for other target metrics) in the respective folder for training GMST-reconstruction model scripts: /scripts/01_ Those scripts make use of the pre-processed observational coverage maps and CMIP6 data. The functions for fitting the ridge regression model is available in /code/_ridge_fun_v2.R

### 4. Processing data from original sources
To process all data from their original source requires to download all data from original sources given in the Data Availability Statement. Further processing is possible by running all scripts contained in /scripts in their successive order. 


## Authors and acknowledgment
The Github repository is maintained by the corresponding author (Sebastian Sippel, sebastian.sippel@uni-leipzig.de). All acknowledgements and references are available in the published paper.


## License
This project is licensed under the MIT License. For more details, see the [LICENSE](LICENSE) file.


***
