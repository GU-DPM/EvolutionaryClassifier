#!/bin/bash

i=0;
while [ $i != 4038 ]; do
 qsub -N "p$i" -q yeangque -v id=$i revisit_optimize_multidrug_responses_trimrates_strategies_index1_twodrugs.sh 
 i=`expr $i + 1`
done
