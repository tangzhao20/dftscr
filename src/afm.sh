#!/bin/bash
#
# This scripts prepare the job directorys and files for the AFM simulation.
# The afm.py script should be run first, and then the PARSEC input file should be modified as needed.
#
# Input: afm.in
#
if [ "$1" == "remove" ]; then
    rm manual_*.dat parsec_st_*.dat steps.dat toten.dat sbatch.log
    exit
fi

if [ "$1" == "seq" ]; then

    if [ ! -f "afm.in" ]; then
        echo "Error: afm.in doesn't exist"
        exit
    fi
    parallel=$(awk -v p="parallel" '$1==p { print $2 }' afm.in)
    lfdet=0
    lfdet=$(awk 'NF==1 && $1=="fdet" { print 1; exit }' afm.in)
    lfdet=$(awk '$1=="fdet" && tolower($2) ~ /^(true|\.true\.)$/ { print 1; exit }' afm.in)

    # if FDET, create the spot dir
    if [[ "$lfdet" == 1 ]]; then
        manualname="manual_${i}_1.dat"
        for (( j=1 ; j<=$parallel ; j++ )) do
            dirname="spot"
            mkdir $dirname
            cd $dirname
            cat ../parsec.in.head > parsec.in
            cat ../parsec_st_spot.dat >> parsec.in
            cat ../job.sh | sed "s/%%jobname%%/a_${i}_${j}/g" > job.sh
            cd ..
        done
    fi

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
            cd ..
            k=$((k+1))
        done
        i=$((i+1))
        manualname="manual_${i}_1.dat"
    done
fi

if [ "$1" == "sbatch" ]; then
    i=1
    k=1
    dirname1="seq_1_1"
    while [ -d $dirname1 ]; do
        for (( j=1 ; j<=$parallel ; j++ )) do
            dirname="seq_${i}_${j}"
            cd $dirname
            printf "$k $i $j " >> ../sbatch.log
            sbatch job.sh | tail -n1 >> ../sbatch.log
            tail -n1 ../sbatch.log
            cd ..
            k=$((k+1))
        done
        i=$((i+1))
        dirname1="seq_${i}_1"
    done
fi
