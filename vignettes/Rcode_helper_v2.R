## NOW FOR MARGINAL DIST CHECK
rm(list=ls());
#library(abn);
if(.Platform$OS.type=="windows"){
setwd("C:/Users/fraser/fraser_local/BayesianNetworks_subproject/BN_R_libraries/active/abn/src");
} else {setwd("~/fraser_local/BayesianNetworks_subproject/BN_R_libraries/active/abn/src");}
source("../R/fitabn.r"); source("../R/bn-internal.r"); source("../R/getmarginal.r");source("../R/unique_nets.r")
source("../R/searchabn.r"); source("../R/bn-internal.r");source("../R/hillsearchabn.r");   
source("../R/tographviz.r");source("../R/fitbn.r");source("../R/searchbn.r");source("../R/hillsearchbn.r");  
source("../R/arcfreq.r");
source("../R/allnodesbn.r");
source("../R/allnodesabn.r");
source("../R/findmostprobablebn.r");
source("../R/getposteriorfeaturesbn.r");
if(.Platform$OS.type=="windows"){dyn.load("abn.dll");} else {dyn.load("abn.so");}
#setwd("/home/fraser/fraser_local/BayesianNetworks_subproject/BN_R_libraries/active/jags/BPHS_examples");

if(.Platform$OS.type=="windows"){
setwd("C:/Users/fraser/fraser_local/BayesianNetworks_subproject/BN_R_libraries/active/vignettes");
} else {setwd("~/fraser_local/BayesianNetworks_subproject/BN_R_libraries/active/abn/vignettes");}

load("../data/var33.RData");
load("../data/pigs1par.RData");
#setwd("/home/fraser/fraser_local/BayesianNetworks_subproject/BN_R_libraries/active/jags/bootstrap_example");

