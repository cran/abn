#!/bin/bash
# init_order_par.R

Rout="Rout"
Rfiles="R"

echo

# Standard syntax.
for a in 2 3 4 5 6
do
  rm -f $Rfiles/"torun.r"
  echo "max.par <- $a" >> $Rfiles/torun.r
  cat $Rfiles/init_order_par.R >> $Rfiles/torun.r
  R CMD BATCH $Rfiles/torun.r $Rout/outfile$a.txt
  echo -n "$a" 
done  

rm -f $Rfiles/"torun.r"

echo; echo
