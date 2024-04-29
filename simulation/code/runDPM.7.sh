#!/bin/bash

i=3061
imax=3571
resultsPath='resultsDir'
while [ $i != $imax ]
do
	./revisit_optimize_multidrug_responses_trimrates_strategies_index1 2 two_drug_argumentlist1.csv $i 10000 $resultsPath/param_ALLDRUG_ $resultsPath/stopt_ALLDRUG_ $resultsPath/dosage_ALLDRUG_ $resultsPath/pop_ALLDRUG_
	i=`expr $i + 1`
done
