/***************************************************************************
 *   Copyright (C) 2006 by F. I. Lewis   
 *   fraser.lewis@ed.ac.uk   
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/** *********************************************************************** 
 * definitions of structures
 *
 ***************************************************************************/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct diskdatamatrix_struct diskdatamatrix;

/** designed to hold integer data and column names **/
struct datamatrix_struct {
      int **data;/** the observed (multinomial) categories for each variable, each row a data point */
      gsl_matrix *datamatrix;
      int numDataPts;/** total number of datapoints**/
      int numVars;/** total number of variables/nodes **/
     /* char **namesVars;*//** array of strings denoting node/variable names, in order of data columns */
      int *numVarlevels;/** number of unique categories per variable */
     /*double *weights;*//** hold a double value which is the weight of each observed data point */
     /* double *xmin;
      double *xmax;*/
      int numparams;
      gsl_vector *priormean;
      gsl_vector *priorsd;
      gsl_vector *gslvec1;
      gsl_vector *gslvec2;
      /*double relerr;*/
      gsl_vector *Y;
  
};

typedef struct datamatrix_struct datamatrix;


/** designed to hold a network definition **/
struct network_struct {
      int **defn;/** each row a variable and each col the parents of the variable, indexes from 0*/
      int numNodes;/** total number of variables/nodes **/
      /*char **namesNodes;*//** array of strings denoting node/variable names, in order of data columns */
      double **nodesparameters;/** node->parentcombinationindex->dirichlet param **/
      int **nodesparameters_lookup;/** node->parentcombinationindex->parentcombination **/
      int *numNodeLevels; /** will hold the total number of parent combinations each node has **/
      int maxparents;/** max number of parents allowed per node */
      int maxParentCombinations;/** maximum number of parent combination which could occur */
      double networkScore;
      int **banlist;
};

typedef struct network_struct network;


/** designed to hold a network definition **/
struct cycle_struct {
     unsigned int *isactive;
     unsigned int *incomingedges;
     unsigned int **graph;
     
};

typedef struct cycle_struct cycle;

/** designed to hold a network definition **/
struct storage_struct {
     int *parentindexes;
     double *n_ij;
     int *multipliers;/** note this is from crossmultiply() **/
     unsigned int **order;/** from generate_random_dag() */
     unsigned int *indexes;/** from generate_random_dag() */

};

typedef struct storage_struct storage;

 struct rparams 
     {double a; 
      double b;};

struct fnparams
       {
         gsl_vector *Y;
	 gsl_vector *vectmp1;
	 gsl_vector *vectmp2;
	 gsl_vector *vectmp1long;
	 gsl_vector *vectmp2long;
	  gsl_vector *vectmp3long;
	 gsl_vector *term1;
	 gsl_vector *term2;
	 gsl_vector *term3;
       gsl_matrix *X;
  gsl_matrix *mattmp1;
  gsl_matrix *mattmp2;
	 gsl_vector *priormean;
	 gsl_vector *priorsd;
	 gsl_vector *betafull;
	 double betafixed;
	 int betaindex;
	 gsl_vector *dgvalues;
	 gsl_matrix *hessgvalues;
	 
	};
	
