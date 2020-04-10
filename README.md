# NEFI_microbe
The goal of this project is to forecast microbial relative abundances at NEON sites in space and time. 

All paths for data products are located in the 'paths.r' file. If you are replicating these analyses on your machine you will need to change where the master data directory points to here. This is also setup to work across multiple machines based on the hostname of a particular machine.

*data_construction* 
This folder contains scripts to process pre-release OTU tables from NEON, queries associated meta-data from the NEON API, and extracts spatial products of interest. There is also code to format ITS data and associated meta-data from Tedersoo et al. 2014 Science, which is being used as a spatial prior.

*analysis*
This folder contains scripts to fit dirichlet models to multivariate fungal functional group data.

*uncertainty_products*
We love quantifying and propogating uncertainty. Many spatial products do not come with an associated uncertainty product. When possible we create these by fitting predicted vs. observed models. This allows us to estimate both a value and its uncertainty when extracting data from spatial products.

*functions*
All custom functions used in this project are stored here.

*figures*
Just a place to save figures to. This directory is ignored in git tracking.
