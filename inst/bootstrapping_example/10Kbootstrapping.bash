#!/bin/bash
# 10KBootstrapping.R

Rout="Rout"
Rfiles="Jags"

echo

# Standard syntax.
for a in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
do
  rm -f $Rfiles/"torun.r"
  echo "index <- $a" >> $Rfiles/torun.r
  cat $Rfiles/10KBootstrapping.R >> $Rfiles/torun.r
  R CMD BATCH $Rfiles/torun.r $Rout/out10Kboot$a.txt
  echo -n "$a" 
done  

rm -f $Rfiles/"torun.r"

echo; echo
