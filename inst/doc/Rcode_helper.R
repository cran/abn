rm(list=ls());
#library(abn);
if(.Platform$OS.type=="windows"){
setwd("C:/Users/fraser/fraser_local/BayesianNetworks_subproject/BN_R_libraries/abn_v0.3-2/abn/src");
} else {setwd("~/fraser_local/BayesianNetworks_subproject/BN_R_libraries/abn_v0.3-2/abn/src");}
source("../R/fitabn.r"); source("../R/bn-internal.r"); source("../R/getmarginal.r");
source("../R/searchabn.r"); source("../R/bn-internal.r");source("../R/hillsearchabn.r");source("../R/unique_nets.r");       
source("../R/prune_nets.r");source("../R/tographviz.r");source("../R/fitbn.r");source("../R/searchbn.r");source("../R/hillsearchbn.r");     
if(.Platform$OS.type=="windows"){dyn.load("abn.dll");} else {dyn.load("abn.so");}
# now specify the model DAG
if(.Platform$OS.type=="windows"){
setwd("C:/Users/fraser/fraser_local/BayesianNetworks_subproject/BN_R_libraries/abn_v0.3-2/abn/inst/doc");
} else {setwd("~/fraser_local/BayesianNetworks_subproject/BN_R_libraries/abn_v0.3-2/abn/inst/doc");}

load("../../data/var33.RData");



rm(list=ls());
#library(abn);
if(.Platform$OS.type=="windows"){
setwd("C:/Users/fraser/fraser_local/BayesianNetworks_subproject/BN_R_libraries/abn_v0.3-2/abn/src");
} else {setwd("~/fraser_local/BayesianNetworks_subproject/BN_R_libraries/abn_v0.3-2/abn/src");}
source("../R/fitabn.r"); source("../R/bn-internal.r"); source("../R/getmarginal.r");
source("../R/searchabn.r"); source("../R/bn-internal.r");source("../R/hillsearchabn.r");   
source("../R/tographviz.r");source("../R/fitbn.r");source("../R/searchbn.r");source("../R/hillsearchbn.r");     
if(.Platform$OS.type=="windows"){dyn.unload("abn.dll");} else {dyn.unload("abn.so");}
