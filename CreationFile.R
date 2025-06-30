# meta -------------------------------------------------------------------------
#
# Creating the alcanR package
#
# Author: Jared Parrish, PhD
# Date: 11/20/2024


# Step 1: Load Libraries and any needed datasets -------------------------------
install.packages("devtools")
 # make sure you have roxygen2 for documentation

# Step 2: Using RStudio create a new package -----------------------------------
# File > New Project > New Directory > R Package.
#  enter the package name that you want and select the directory
#
# Components
#  Description: Metadata about the package
#  NAMESPACE: Controls which functions are exported from package
#  R/: Directory where R script functions are stored.
#  man/: Directory for the documentation files.

# Step 3: Write functions ------------------------------------------------------
# describe the paramaters going into the function and other details

# Step 4: Document functions ---------------------------------------------------
# Run the code below to generate the documentation. Documentation is created
# from the function files tagging and commenting.
devtools::document()

# Step 5: Develop tests --------------------------------------------------------
# install "testthat" package if you don't have it.
# create sub directory tests and inside testthat.
# develop test scripts and place in the testthat sub.
#
# base code for testing
library(survival)
library(survey)
load(".\\data\\AlcanLink22.rda")
t1 <- svydesign(id=~1,strata=~SUD_NEST,fpc=~fpc,weights=~WTANAL,data=AlcanLink22)

source(".\\R\\SvyIncProp.R")

frpt_IP <- svyIncProp(SurveyDesign = t1,
                      event = "FST_RPT_Indicator",
                      eventime = "FST_RPT_CensorAge",
                      nreps = 20)
# run test
devtools::load_all()
devtools::test()

# Step 6: Build and check package ----------------------------------------------
# Requires: rtools
devtools::check()

devtools::build()

# Step 7: Load package ---------------------------------------------------------
install.packages("C:\\Users\\jwparrish\\Documents\\alcanR_0.1.0.tar.gz")

