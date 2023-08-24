#!/bin/bash

# This scripts prepare the job directorys and files for the AFM simulation.
#
# dftscr_dir/src should be added to path

if [ ! -f "afm.in" ]; then
    echo "Error: afm.in doesn't exist"
    stop
fi

afm.py

parallel=$(awk -v p="parallel" '$1==p { print $2 }' afm.in)
for (( i=1 ; i<=3 ; i++ )) do
    for (( j=1 ; j<=$parallel ; j++ )) do
        dirname="seq_${i}_${j}"
        mkdir $dirname
        cd $dirname
        mv ../manual_${i}_${j}.dat manual.dat
        cat ../parsec.in.head > parsec.in
        STEP=`awk -v j=$j 'NR==j {print}' ../steps.dat`
        echo "movement_num  $((STEP-1))" >> parsec.in
        cat ../parsec_st.dat >> parsec.in
        cat ../job.sh | sed "s/%%jobname%%/a_${i}_${j}/g" > job.sh
        cd ..
    done
done
