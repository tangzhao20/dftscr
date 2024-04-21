#!/bin/bash
#
# This scripts prepare the job directorys and files for the AFM simulation.
# The afm.py script should be run first, and then the PARSEC input file should be modified as needed.
#
# Input: afm.in

if [ ! -f "afm.in" ]; then
    echo "Error: afm.in doesn't exist"
    exit
fi

parallel=$(awk -v p="parallel" '$1==p { print $2 }' afm.in)

i=1
k=1
manualname="manual_1_1.dat"
while [ -f $manualname ]; do
    for (( j=1 ; j<=$parallel ; j++ )) do
        dirname="seq_${i}_${j}"
        mkdir $dirname
        cd $dirname
        cp ../manual_${i}_${j}.dat manual.dat
        cat ../parsec.in.head > parsec.in
        STEP=`awk -v j=$j 'NR==j {print}' ../steps.dat`
        echo "movement_num  $((STEP-1))" >> parsec.in
        cat ../parsec_st_${i}_${j}.dat >> parsec.in
        cat ../job.sh | sed "s/%%jobname%%/a_${i}_${j}/g" > job.sh
        if [ "$1" == "sbatch" ]; then
            printf "$k $i $j " >> ../sbatch.log
            sbatch job.sh | tail -n1 >> ../sbatch.log
            tail -n1 ../sbatch.log
        fi
        cd ..
        k=$((k+1))
    done
    i=$((i+1))
    manualname="manual_${i}_1.dat"
done
