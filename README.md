Exploring Causal Estimands
==========================

<img src="https://img.shields.io/badge/Study%20Status-Started-blue.svg" alt="Study Status: Started">

- Analytics use case(s): **Population-Level Estimation**
- Study type: **Methods Research**
- Tags: **-**
- Study lead: **Martijn Schuemie**
- Study lead forums tag: **[schuemie](https://forums.ohdsi.org/u/schuemie)**
- Study start date: **April 1, 2025**
- Study end date: **-**
- Protocol: **-**
- Publications: **-**
- Results explorer: **-**

An exploration and evaluation of various causal models and their estimands using simulation and real data. 

# How to run

1. Follow [these instructions](https://ohdsi.github.io/Hades/rSetup.html) for seting up your R environment, including RTools and Java. 

2. Open your study package in RStudio. Use the following code to install all the dependencies:

	```r
	renv::restore()
	```
	
3. See the Simulations and RealworldExample folders for the R scripts used in this study. Running the RealworldExample experiment from scratch requires a database in [Common Data Model version 5](https://github.com/OHDSI/CommonDataModel).


# License

The contents of this repository are licensed under the Apache License 2.0.