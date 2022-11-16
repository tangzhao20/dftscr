# submit the jobs of multiple structures
# pospath.sh $1
# $1 is the number of structures between _i and _f, same as pospath.py

N=$1

for ((i=1;i<=$((3*N+1));i++)) 
do
    if [ -d "$i" ]; then
        continue
    fi
    mkdir $i
    cd $i
    POSCARname="POSCAR.$i"
    cp ../$POSCARname POSCAR
    cp ../job.sh job.sh
    sed -i "s/%%jobname%%/$i/g" job.sh
    printf "$i " >> ../sbatch.log
    sbatch job.sh >> ../sbatch.log
    tail -n1 ../sbatch.log
    cd ..
done
