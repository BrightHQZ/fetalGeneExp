#!/bin/bash

while getopts ":d:r:t:b:" opt
do
    case $opt in
        d)
        dirL=$OPTARG;
        ;;
        r)
        rootDir=$OPTARG;
        ;;
        t)
        inType=$OPTARG;
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
  echo "-d is necessary in this program! It is the work dir: (-d B14)";
  exit;
fi;

if [ -z $rootDir ]; then
  echo "-r is necessary in this program! It is the root dir: (-r /data/gestation)";
  exit;
fi;

if [ -z $inType ]; then
  echo "-t is necessary in this program! It is the work for childer, maternal or all: (-t c/m/a)";
  exit;
fi;

if [[ $inType != "c" && $inType != "m" && $inType != "a" ]]; then
  echo "The input value of -t should in the three char 'c','m' or 'a'";
  exit;
fi;

if [ ! -d "$rootDir/$dirL/childer/" ]; then
  mkdir $rootDir/$dirL/childer/
fi;

if [ ! -d "$rootDir/$dirL/maternal/" ]; then
  mkdir $rootDir/$dirL/maternal/
fi;

if [ -d $rootDir/$dirL ]; then
  echo "The bam files of in $rootDir/$dirL will processed............";
  cd $rootDir/$dirL
else
  echo "The dir is not exist: $rootDir/$dirL";
  exit;
fi; 

for file in `ls *.bam | sed s/.bam//`
do
  echo $rootDir/$dirL/$file".bam is processing";
  start=$(date +%s)
  if [ ! -f $rootDir/$dirL/childer/$dirL"_finalBamFile_childer.bam" ] && [ ! -f $rootDir/$dirL/childer/$file".bam" ] && [[ $inType == "c"  || $inType == "a" ]]; then 
    time /data/bioTools/bin/sambamba markdup -r $rootDir/$dirL/$file".bam" $rootDir/$dirL/$file"_filter.bam" -t 8 --tmpdir=$rootDir/$dirL/tmp
    time samtools view -@ 8 -q 20 -f 0 -f 16 -h $rootDir/$dirL/$file"_filter.bam" | awk '{if ($0~/^@/) print($0); else if (length($10) >= 120 && length($10) < 156 ) print($0); }' > $rootDir/$dirL/childer/$file".sam"
    time samtools sort -@ 8 $rootDir/$dirL/childer/$file".sam" > $rootDir/$dirL/childer/$file".bam"
    rm $rootDir/$dirL/childer/$file".sam" 
    rm $rootDir/$dirL/$file"_filter"*.*
  fi;
  if [ ! -f $rootDir/$dirL/maternal/$dirL"_finalBamFile_maternal.bam" ] && [ ! -f $rootDir/$dirL/maternal/$file".bam" ] && [[ "$inType" == "m" || "$inType" == "a" ]]; then
    time /data/bioTools/bin/sambamba markdup -r $rootDir/$dirL/$file".bam" $rootDir/$dirL/$file"_filter.bam" -t 8 --tmpdir=$rootDir/$dirL/tmp
    time samtools view -@ 8 -q 20 -f 0 -f 16 -h $rootDir/$dirL/$file"_filter.bam" | awk '{if ($0~/^@/) print($0); else if (length($10) >= 156 && length($10) < 200 && $3 != "chrY") print($0); }' > $rootDir/$dirL/maternal/$file".sam"
    time samtools sort -@ 8 $rootDir/$dirL/maternal/$file".sam" > $rootDir/$dirL/maternal/$file".bam"
    rm $rootDir/$dirL/maternal/$file".sam" 
    rm $rootDir/$dirL/$file"_filter"*.*
  fi;
  end=$(date +%s)
  take=$(( end - start ))
  echo $rootDir/$dirL/$file" is finished, time is ${take} seconds";
done 

rm -rf $rootDir/$dirL/tmp

cd $rootDir;

if [[ $inType == "c" || $inType == "a" ]]; then
  if [ ! -f $rootDir/$dirL/childer/$dirL"_finalBamFile_childer.bam" ]; then 
    time samtools merge $rootDir/$dirL/childer/$dirL"_finalBamFile_childer.bam" $rootDir/$dirL/childer/*.bam
    time ls $rootDir/$dirL/childer/C*.bam | xargs rm
    #time samtools sort -@ 8  $rootDir/$dirL/childer/$dirL"_finalBamFile_childer.bam" >  $rootDir/$dirL/childer/$dirL"_finalBamFile_childer_sort.bam"
    time samtools index -@ 8  $rootDir/$dirL/childer/$dirL"_finalBamFile_childer.bam"
  fi;
  time samtools bedcov $inBed $rootDir/$dirL/childer/$dirL"_finalBamFile_childer.bam" | awk '{print $0."\t"sprintf("%.3f",$5/($3-$2))}' > $rootDir/$dirL/childer/$dirL"_childer_TSS_P_S.bed"

  if [ ! -d "$rootDir/$dirL/childer/bed/" ]; then
    mkdir $rootDir/$dirL/childer/bed/
  fi;

  time perl $rootDir/splitDeepBEDToGenes.pl --bed $rootDir/$dirL/childer/$dirL"_childer_TSS_P_S.bed" --out $rootDir/$dirL/childer/bed
  rm -rf $rootDir/$dirL/childer/tmp/
fi;

if [[ $inType == "m" || $inType == "a" ]]; then
  if [ ! -f $rootDir/$dirL/maternal/$dirL"_finalBamFile_maternal.bam" ]; then
    time samtools merge $rootDir/$dirL/maternal/$dirL"_finalBamFile_maternal.bam" $rootDir/$dirL/maternal/*.bam
    time ls $rootDir/$dirL/maternal/C*.bam | xargs rm
    #time samtools sort -@ 8  $rootDir/$dirL/maternal/$dirL"_finalBamFile_maternal.bam" >  $rootDir/$dirL/maternal/$dirL"_finalBamFile_maternal_sort.bam"
    time samtools index -@ 8  $rootDir/$dirL/maternal/$dirL"_finalBamFile_maternal.bam"
  fi;
  time samtools bedcov $inBed $rootDir/$dirL/maternal/$dirL"_finalBamFile_maternal.bam" | awk '{print $0."\t"sprintf("%.3f",$5/($3-$2))}' > $rootDir/$dirL/maternal/$dirL"_maternal_TSS_P_S.bed"

  if [ ! -d "$rootDir/$dirL/maternal/bed/" ]; then
    mkdir $rootDir/$dirL/maternal/bed/
  fi;

  time perl $rootDir/splitDeepBEDToGenes.pl --bed $rootDir/$dirL/maternal/$dirL"_maternal_TSS_P_S.bed" --out $rootDir/$dirL/maternal/bed
  rm -rf $rootDir/$dirL/maternal/tmp/
fi;


#/data/bioTools/bin/bedtools/bedtools genomecov -bga -split -ibam /data/gestation/$dirL/childer/$dirL"_finalBamFile_childer.bam" > /data/gestation/$dirL/childer/$dirL"_childer.bed"
#/data/bioTools/bin/bedtools/bedtools genomecov -bga -split -ibam /data/gestation/$dirL/maternal/$dirL"_finalBamFile_maternal.bam" > /data/gestation/$dirL/maternal/$dirL"_maternal.bed"
#/data/bioTools/bin/bedtools/bedtools intersect -a /data/gestation/TSS_protein.bed -b /data/gestation/$dirL/childer/$dirL"_childer.bed" -wa -wb | awk '{ print($6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$5) }' > /data/gestation/$dirL/childer/$dirL"_childer_TSS_P_S.bed"
#/data/bioTools/bin/bedtools/bedtools intersect -a /data/gestation/TSS_protein.bed -b /data/gestation/$dirL/maternal/$dirL"_maternal.bed" -wa -wb | awk '{ print($6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$5) }' > /data/gestation/$dirL/maternal/$dirL"_maternal_TSS_P_S.bed"


#/data/bioTools/bin/bedtools/bedtools merge -d -5 -i /data/gestation/$dirL/childer/$dirL"_childer_filter.bed" -c 4 -o max,count | awk '{if($5 > 1) print($0"\t"$3-$2)}' | awk '$6 > 50' > /data/gestation/$dirL/childer/$dirL"_childer_filter_merged.bed"
#/data/bioTools/bin/bedtools/bedtools genomecov -bga -split -ibam /data/gestation/$dirL/maternal/$dirL"_finalBamFile_childer.bam" | awk '$4 > 9' > /data/gestation/$dirL/maternal/$dirL"_maternal_filter.bed"
#/data/bioTools/bin/bedtools/bedtools merge -d -5 -i /data/gestation/$dirL/maternal/$dirL"_maternal_filter.bed" -c 4 -o max,count | awk '{if($5 > 1) print($0"\t"$3-$2)}' | awk '$6 > 50' > /data/gestation/$dirL/maternal/$dirL"_maternal_filter_merged.bed"

#samtools view -@ 8 -q 20 -f 0 -f 16 -h test.bam | awk '{if ($0~/^@/) print($0); else if (length($10) >= 120 && length($10) < 156 ) print($0); }' > test_childer.sam
#samtools view -@ 8 -q 20 -f 0 -f 16 -h test.bam | awk '{if ($0~/^@/) print($0); else if (length($10) >= 156 && length($10) < 200 && $3 != "chrY") print($0); }' > test_maternal.sam

#samtools sort -@ 8 test_childer.sam > test_childer.bam
#samtools sort -@ 8 test_maternal.sam > test_maternal.bam

#samtools index test_childer.bam
#samtools index test_maternal.bam

#samtools merge finalBamFile.bam *.bam


