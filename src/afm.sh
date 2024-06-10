#!/bin/bash
#
# This scripts prepare the job directorys and files for the AFM simulation.
# The afm.py script should be run first, and then the PARSEC input file should be modified as needed.
#
# Input: afm.in
#
if [[ "$1" == "remove" ]]; then
    rm manual_*.dat parsec_st_*.dat steps.dat toten.dat sbatch.log
    exit
fi

if [[ "$1" == "seq" || "$1" == "sbatch" || "$1" == "check" ]]; then
    # read afm.in
    if [[ ! -f "afm.in" ]]; then
        echo "Error: afm.in doesn't exist"
        exit
    fi
    sed 's/[#!].*//' afm.in > afm.in.tmp

    parallel=$(awk -v flag="parallel" '$1==flag { print $2 }' afm.in.tmp)
    lfdet=$(awk '$1=="fdet" && ( NF==1 || tolower($2) ~ /^(true|\.true\.)$/ ) { print 1; found=1; exit } END { if (!found) print 0 }' afm.in.tmp)

    z_min=$(awk -v flag="z_range" '$1==flag { print $2 }' afm.in.tmp)
    z_max=$(awk -v flag="z_range" '$1==flag { print $3 }' afm.in.tmp)
    z_spacing="0.3"
    z_spacing=$(awk -v flag="z_spacing" '$1==flag { print $2; found=1; exit } END { if (!found) print 0.3 } ' afm.in.tmp)
    Nz=$(echo "($z_max - $z_min) / $z_spacing + 0.000001" | bc | cut -d '.' -f1)

    rm afm.in.tmp
fi

if [[ "$1" == "seq" ]]; then
    # if FDET, create the spot dir
    if [[ "$lfdet" == 1 ]]; then
        dirname="spot"
        mkdir $dirname
        cd $dirname
        cat ../parsec.in.head > parsec.in
        echo "" >> parsec.in
        cat ../parsec_st_spot.dat >> parsec.in
        cat ../job.sh | sed "s/%%jobname%%/a_spot/g" > job.sh
        cd ..
    fi

    k=1
    for (( i=1 ; i<=$Nz ; i++ )) do
        for (( j=1 ; j<=$parallel ; j++ )) do
            dirname="seq_${i}_${j}"
            mkdir $dirname
            cd $dirname
            cp ../manual_${i}_${j}.dat manual.dat
            cat ../parsec.in.head > parsec.in
            echo "" >> parsec.in
            echo "#---------output from afm.sh----------" >> parsec.in
            echo "Potential_Field: .TRUE." >> parsec.in
            echo "Potential_Field_Name: s_pot.dat" >> parsec.in
            echo "Kinetic_Energy_Functional: pb" >> parsec.in
            echo "" >> parsec.in
            echo "Minimization manual" >> parsec.in
            STEP=`awk -v j=$j 'NR==j {print}' ../steps.dat`
            echo "movement_num  $((STEP-1))" >> parsec.in
            echo "" >> parsec.in
            cat ../parsec_st_${i}_${j}.dat >> parsec.in
            cat ../job.sh | sed "s/%%jobname%%/a_${i}_${j}/g" > job.sh
            ln -s ../spot/pot.dat s_pot.dat
            cd ..
            k=$((k+1))
        done
    done
fi

if [[ "$1" == "sbatch" ]]; then
    k=1
    for (( i=1 ; i<=$Nz ; i++ )) do
        for (( j=1 ; j<=$parallel ; j++ )) do
            dirname="seq_${i}_${j}"
            cd $dirname
            printf "$k $i $j " >> ../sbatch.log
            sbatch job.sh | tail -n1 >> ../sbatch.log
            tail -n1 ../sbatch.log
            cd ..
            k=$((k+1))
        done
    done
fi

if [[ "$1" == "check" ]]; then
    k=1
    for (( i=1 ; i<=$Nz ; i++ )) do
        for (( j=1 ; j<=$parallel ; j++ )) do
            dirname="seq_${i}_${j}"
            cd $dirname
            if [[ ! -f "parsec.out" ]]; then
                echo "$dirname : parsec.out not found"
                cd ..
                continue
            fi
            line4=$(tail -n 4 parsec.out | head -n 1)
            if [[ $line4 != " Current date/time:"* ]]; then
                echo "$dirname : the calculation did not finish"
            fi
            cd ..
            k=$((k+1))
        done
    done
fi
