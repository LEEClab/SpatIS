# SpatIS - Spatial Individual Specialization Indices

This repository shares the code and details on the Spatial Individual Specialization Index (*SpatIS*) and the Spatis Individual Complementary Specialization Index (*SpatICS*). Those two indices were developed and first presented in Kerches-Rogeri et al (2020, see reference below).

For more information on how to run it for your data please look at the SpatIS website: https://leeclab.github.io/SpatIS/.  
You will also find further information in the main article and supplementary information of the reference below (Kerches-Rogeri et al. 2020). 

We present here the source code for running *SpatIS* and *SpatICS* using movement data of individuals. We also show VHF movement data for the yellow–shouldered bats *Sturnira lilium* in Southeastern Brazil, an R script for analyzing this data using *SpatIS* and *SpatICS*, and other scripts and data used in Kerches-Rogeri et al (*2020*).

## Code

The scripts in this repository are the following (look at the folder [code](https://github.com/LEEClab/SpatIS/tree/master/code)):
- [spatis_source_code_v1_0.R](https://github.com/LEEClab/SpatIS/blob/master/code/spatis_source_code_v1_0.R): source code with SpatIS and SpatIS_randomize R functions (details within the script); this is the main file for running *SpatIS* and *SpatICS*;
- [individual_specialization_space_use.R](https://github.com/LEEClab/SpatIS/blob/master/code/individual_specialization_space_use.R): script for spatial and individual specialization analyses from Kerches-Rogeri et al. (*2020);
- [test_for_body_size_and_sex.R](https://github.com/LEEClab/SpatIS/blob/master/code/test_for_body_size_and_sex.R): script for testing for body size and sex effects on individual use of space through binomial GLMs;
- [use_foraging_areas_fruit_density.R](https://github.com/LEEClab/SpatIS/blob/master/code/use_foraging_areas_fruit_density.R): script for testing the effects of fruit density on the use of foraging areas y each individual;
- spatis_tutorial (in [pdf](https://github.com/LEEClab/SpatIS/blob/master/spatis_tutorial/spatis_tutorial.pdf) and [Rmarkdown](https://github.com/LEEClab/SpatIS/blob/master/code/spatis_tutorial.Rmd) formats): Tutorial for applying SpatIS for other data sets. It may also be found [here](https://github.com/LEEClab/SpatIS/tree/master/spatis_tutorial).
- scripts and pdf for the other supplementary material from Kerches-Rogeri et al (*2020*). 

## Data

We present here all the data used for analysis by Kerches-Rogeri et al (*2020*):
- [spatis_paper_bat_positions.csv](https://github.com/LEEClab/SpatIS/blob/master/data/spatis_paper_bat_positions.csv): telemetry fixes for 11 individual bats monitored in Sao Carlos, SP, Brazil.
- [spatis_paper_individual_bat_info.csv](https://github.com/LEEClab/SpatIS/blob/master/data/spatis_paper_individual_bat_info.csv): individual information for the tracked bats.
- [proportion_use_landcover_polygons.csv](https://github.com/LEEClab/SpatIS/blob/master/data/proportion_use_landcover_polygons.csv): proportion of use of land cover polygons in the landscape by each individual.
- [proportion_use_foraging_sites.csv](https://github.com/LEEClab/SpatIS/blob/master/data/proportion_use_foraging_sites.csv): proportion of use of foraging sites by each individual.
- [use_foraging_areas_density_fruits.csv](https://github.com/LEEClab/SpatIS/blob/master/data/use_foraging_areas_density_fruits.csv): density of the main fruit genera consumed by bats in each foraging area.

## Citation

If you use *SpatIS*, please refer to

Kerches-Rogeri, P.; Niebuhr, B. B.; Muylaert, R. L., Mello, M. A. R. 2020. [Individual specialization in the use of space of frugivorous bats](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2656.13339). Journal of Animal Ecology, *accepted*. DOI: 10.1111/1365-2656.13339.

You can also cite this repository, which contains both the code and data used in the publication:  
Niebuhr, B. B.; Kerches-Rogeri P.; Muylaert, R. L., Mello, M. A. R. (2020). SpatIS: Spatial Individual Specialization Indices (Version v1.0). Zenodo.
[![DOI](https://zenodo.org/badge/95931986.svg)](https://zenodo.org/badge/latestdoi/95931986)

The movement data of bats presented here may also be found and downloaded from [MoveBank](www.movebank.org), in the study ["Individual specialization in the use of space by frugivorous bats"](https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study577905925).

## Contact

If you have questions or suggestions, do not hesitate to contact us (or open an issue [here](https://github.com/LEEClab/SpatIS/issues)):
+ Patricia Kerches-Rogeri <<parogeri@gmail.com>>  
+ Bernardo Brandão Niebuhr <<bernardo_brandaum@yahoo.com.br>>  
+ Renara L. Muylaert <<renatamuy@gmail.com>>  
+ Marco A. R. Mello <<marmello@usp.br>>
