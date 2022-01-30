#!/bin/bash

while getopts ":r:f:t:m:s:b:" opt
do
    case $opt in
        r)
        dirL=$OPTARG;
        ;;
        f)
        from=$OPTARG;
        ;;
        t)
        to=$OPTARG;
        ;;
        m)
        model=$OPTARG;
        ;;
        s)
        sex=$OPTARG;
        ;;
        b)
        inBed=$OPTARG;
        ;;
        ?)
        echo "未知参数$OPTARG"
        exit 1;;
    esac
done

if [ -z $dirL ]; then
    echo "-r is necessary in this program! It is the root dir: (-r /data/gestation)";
    exit;
fi;

if [ -z $from ]; then
    echo "-f is necessary in this program! It is the number of the firstly dir: (-f 14)";
    exit;
fi;

if [ -z $to ]; then
    echo "-t is necessary in this program! It is the number of the lastly dir: (-t 20)";
    exit;
fi;

if [ -z $model ]; then
    echo "-m is necessary in this program! It is the working model: (-m c/m/a)";
    echo "c: childer"
    echo "m: maternal"
    echo "a: all"
    exit;
fi;

if [[ $model != "c" && $model != "m" && $model != "a" ]]; then
    echo "-m is necessary in this program! It is the working model: (-m c/m/a)";
    echo "c: childer";
    echo "m: maternal";
    echo "a: all";
    exit;
fi;

if [ -z $sex ]; then
    echo "-s is necessary in this program: (-s B/G)";
    exit;
fi;

if [[ ${sex^^} != "G" && ${sex^^} != "B" ]]; then
    echo "-s is necessary in this program: (-s B/G)";
    exit;
fi;

if [ -z $inBed ]; then
    echo "-b is necessary in this program! It is the target range for screening: (-b ...bed)";
    exit;
fi;

if [ ! -f $inBed ]; then
    echo "The file: "$inBed" is not exists.";
    exit;
fi;


tmp_fifofile="$dirL/$.fifo"
mkfifo $tmp_fifofile;
exec 6<>$tmp_fifofile;
rm $tmp_fifofile;

for i in {1..3}
do
    echo
done >&6

#call snp for chromosomes
for i in `seq $from $to`;
do
    read -u6
    {
        echo "Processing bam files in $dirL/${sex^^}$i........!";
        bash $dirL/process.sh -d ${sex^^}$i -r $dirL -t $model -b $inBed;
        echo >&6;   
    } &
done;
wait
exec 6>&-
echo "All processes have finished!";



