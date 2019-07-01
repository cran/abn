# abn devel project (PRIVATE)


## structure
* abn/ contains source code
* Use testthat for tests

## Modifications
* bug correction for simulation and prediction
* plot function with mb, fitabn output
* implementation of forumla in: Buildscore/fitabn (modification of rd file) 
* plot function + rd file
* prediction function + rd file
* simulation function + rd file

## To do:
### R Package (release)
* rd file: plot, simulation, prediction
* implementation of test for plot/prediction/simulation functions
* diagnostic plots
* bootstrapping
* Implement parallel computing in fitabn(), buildscorecache() 

### R package (extension):
* Implement new prior MP (multinomial + data separation) 
* Data separation 
* Random effect (already one layer) 
* Multinomial
* Negative binomial, ... 
* Time/Memory resources optimization
* Handling missing data 