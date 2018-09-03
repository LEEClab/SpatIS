# SpatIS - Spatial Individual Specialization Index

This repository shares the code and details on the Spatial Individual Specialization Index (SpatIS). More information on how to run it for your data may be found [in this tutorial](https://github.com/LEEClab/SpatIS/tree/master/spatis_tutorial) and within the main article and supplementary information of the reference below (Kerches-Rogeri et al. *in prep*). 

We present here the source code for running SpatIS using movement data of individuals. We also show VHF movement data for the yellow–shouldered bats *Sturnira lilium* in Southeastern Brazil, an R script for analyzing this data using *SpatIS*, and other scripts and data used in Kerches-Rogeri et al (*in prep.*).

## Code

The scripts in this repository are the following (look at the folder [code](https://github.com/LEEClab/SpatIS/tree/master/code):
- [spatis_source_code_v1_0.R](https://github.com/LEEClab/SpatIS/blob/master/code/spatis_source_code.R): source code with SpatIS and SpatIS_randomize R functions (details within the script)
- [individual_specialization_space_use.R](https://github.com/LEEClab/SpatIS/blob/master/code/individual_specialization_space_use.R): script for spatial and individual specialization analyses from Kerches-Rogeri et al (*in prep.*).
- [test_for_body_size_and_sex.R](https://github.com/LEEClab/SpatIS/blob/master/code/test_for_body_size_and_sex.R): script for testing for body size and sex effects on individual use of space
- spatis_tutorial (in [pdf](https://github.com/LEEClab/SpatIS/blob/master/code/spatis_tutorial.pdf) and [markdown](https://github.com/LEEClab/SpatIS/blob/master/code/spatis_tutorial.Rmd) formats): Tutorial for applying SpatIS for other data sets. It may also be found [here](https://github.com/LEEClab/SpatIS/tree/master/spatis_tutoria).

## Data

We present here all the data used for analysis by Kerches-Rogeri et al (*in prep.*):
- [bat_locations_ufscar_VHF.csv](https://github.com/LEEClab/SpatIS/blob/master/data/bat_locations_ufscar_VHF.csv): telemetry fixes for 11 individual bats monitored in Sao Carlos, SP, Brazil.
- [test_for_bodysize_sex_polygons.csv](https://github.com/LEEClab/SpatIS/blob/master/data/test_for_bodysize_sex_polygons.csv): proportion of use of habitat polygons in the landscape by each individual.
- [test_for_bodysize_sex_foraging_sites.csv](https://github.com/LEEClab/SpatIS/blob/master/data/test_for_bodysize_sex_foraging_sites.csv): proportion of use of foraging sites by each individual.

## Citation

If you use *SpatIS*, please refer to

Kerches-Rogeri, P.; Niebuhr, B. B.; Muylaert, R. L., Mello, M. A. R. Individual specialization in the space use of frugivorous bats. *In prep.*

## Contact

If you have questions or suggestions, do not hesitate to contact us (or open an issue [here](https://github.com/LEEClab/SpatIS/issues)):
+ Patricia Kerches-Rogeri <<parogeri@gmail.com>>  
+ Bernardo Brandão Niebuhr <<bernardo_brandaum@yahoo.com.br>>  
+ Renara L. Muylaert <<renatamuy@gmail.com>>  
+ Marco A. R. Mello <<marmello@gmail.com>>
