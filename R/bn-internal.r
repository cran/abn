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

## Coerce observed data.frame into factors then put the levels into integers
makeintofactors<-function(data.df){
     ##coerce data.frame into factors
        if(!is.data.frame(data.df)){stop("data must be in a data.frame");}
        for(i in 1:dim(data.df)[2]){##for each variable
          if(!is.factor(data.df[,i])){## if not already a factor then coerce
                cat("warning: one or more variables coerced into factors\n");
                if(length(which(is.na(data.df[,i])==T))>0){stop("missing values must be removed");}
                data.df[,i]<-as.factor(data.df[,i]);}
                ## now coerce factor levels into integers
                data.df[,i]<-as.integer(data.df[,i]);
                }
          return(data.df);
          }

## check matrix passed as DAG is valid and if so coerces to integer. Use a square format, rows are children and cols parents 
##    0 - 00001000110000, node 0 has parents comprising of node 4,8,9
##    1 - 00010010011001
##    2 - 01010100001000
##    etc
makeintobinarymatrix<-function(dag.df,ban){
        if(!is.matrix(dag.df)){stop("DAG must be a matrix");}
        if(is.null(colnames(dag.df)) || is.null(rownames(dag.df))){stop("must have both row and column names set in DAG");}
        for(i in 1:dim(dag.df)[1]){##for each variable
          if(length(which(dag.df[i,]==1)) + length(which(dag.df[i,]==0))!=dim(dag.df)[2]){## if not already a factor then coerce
                stop("DAG definition invalid\n");}}
        if(dim(dag.df)[1]!=dim(dag.df)[2]){stop("matrix must be square");}
        for(i in 1:dim(dag.df)[1]){##for each variable
                                      if(dag.df[i,i]!=0 && !ban){stop("child cannot be its own parent - diagonal entries must be zero");}}
        ##end of checks so now coerce 
        local.df<-matrix(integer(dim(dag.df)[1]*dim(dag.df)[2]),ncol=dim(dag.df)[2]);
         for(i in 1:dim(dag.df)[1]){##for each variable   
                ## now coerce into ints
                local.df[i,]<-as.integer(dag.df[i,]);}
        maxparents<-max(apply(local.df,1,sum));
          return(list(dag=local.df,maxparents=maxparents));
          }

## trivial function for use in apply within prune.nets
trim<-function(val,threshold){if(val<threshold){return(0);} else {return(1);}}


