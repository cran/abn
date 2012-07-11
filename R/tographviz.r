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

tographviz<-function(dag,data.df,outfile){
            if(!is.matrix(dag)){stop("must be a matrix");}
            if(is.null(rownames(dag)) || is.null(colnames(dag))){stop("names must be set");}
            ## create header part
            cat("digraph dag {","\n\n",file=outfile,append=FALSE);
            for(i in 1:length(colnames(dag))){
                                 if(is.factor(data.df[,i])){cat(paste("\"",colnames(dag)[i],"\"[shape=square];\n",sep=""),file=outfile,append=TRUE);
                       } else {cat(paste("\"",colnames(dag)[i],"\"[shape=oval];\n",sep=""),file=outfile,append=TRUE);}
            }
            cat("\n\n\n",file=outfile,append=TRUE)
    
    for(i in colnames(dag)){##for each variable
             children<-which(dag[,i]==1);##get row with children
             if(length(children)>=1){##if have at least one child
             child.nom<-rownames(dag)[children];
               for(j in child.nom){cat("\"",i,"\"","->","\"",j,"\";","\n",sep="",file=outfile,append=TRUE);
                 }
                }
            }

            ## footer part
            cat("\n}\n",file=outfile,append=TRUE);
            
}


