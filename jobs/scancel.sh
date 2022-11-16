# scancel.sh $1 $2 (from $1 to $2)
# or scancel.sh $1 (only $1)

# read and modify sbatch.log


if [ -f "sbatch.log.old" ]; then
    rm sbatch.log.old
fi
cp sbatch.log sbatch.log.old

if [ -z "$1" ]; then
    echo "Need more arguments."
    exit
elif [ -z "$2" ]; then
    i1=$1
    i2=$1
else
    i1=$1
    i2=$2
fi

while IFS= read -r var
do
    HEADER=`echo "$var" | cut -d ' ' -f 2`
    JOB=`echo "$var" | cut -d ' ' -f 1`
    if [ "$HEADER" == "Submitted" ] && [ "$JOB" -ge "$i1" ] && [ "$JOB" -le "$i2" ]; then
        JOBID=`echo "$var" | cut -d ' ' -f 5`
        CLUSTER=`echo "$var" | cut -d ' ' -f 8`
        echo "Canceling job $JOBID in cluster $CLUSTER"
        scancel -M $CLUSTER $JOBID
        sed -i "/$var/d" sbatch.log
        rm -r pos_$((JOB))
    fi
done < sbatch.log
