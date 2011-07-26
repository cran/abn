## fitbn.R --- 
## Author          : Fraser Lewis
## Created On      : Sun May 13:43 2010
## Last Modified By: Fraser Lewis
## Last Modified On: Sun May 13:43 2010
## Update Count    : 0
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Fraser Lewis
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

prunenets<-function(dags.list,threshold.freq){
            ## step 1. get the number of times each individual arc appears
            numvars<-dim(dags.list[[1]])[1];
            tot<-matrix(rep(0,numvars*numvars),ncol=numvars);## this will contain the total freq of each arc
            for(i in 1:length(dags.list)){tot<-tot+dags.list[[i]];}
            ## step 2. prune each network             
            trimmed<-dags.list;## a copy just to avoid addition mallocing - hopefully. Each entry is overwritten
            for(i in 1:length(dags.list)){
            best<-dags.list[[i]];
            best2<-best*tot;## this removes from tot any arcs which were not in best (NOTE: best is contained in tot)
            best.trimmed<-apply(best2,c(1,2),FUN=trim,threshold=threshold.freq); 
            trimmed[[i]]<-best.trimmed;
            }
            return(list(arcs.sum=tot,dags=trimmed)); 
}

