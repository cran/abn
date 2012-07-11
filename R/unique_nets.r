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

uniquenets<-function(dags.list){  ## pass a list of matrices and return a vector of the indices of unique matrices
 compactnets<-rep("",length(dags.list));
 for(i in 1:length(dags.list)){
 tmp<-paste(as.vector(dags.list[[i]]),collapse=",");
 compactnets[i]<-paste(as.vector(tmp),collapse=",");}
 compactnets.unq<-unique(compactnets);
 unq.indexes<-rep(-1,length(compactnets.unq));
 for(i in 1:length(compactnets.unq)){## for each unique network
 thisone<-which(compactnets==compactnets.unq[i])[1];##take first hit only
 unq.indexes[i]<-thisone;
 }
 return(unq.indexes);
}

