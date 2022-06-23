#!/usr/bin/env bash
# Merged forward and reverse reads

flash=~/programs/FLASH/FLASH-1.2.11-Linux-x86_64/flash

R1=i39_S27_L001_R1_001.fastq
R2=i39_S27_L001_R2_001.fastq
out="i39"

# ids="L422T01 L422T02 L422T03 L422T04 L422T05 L422T06 L422T07 L422T08 L422T09 L422T10 L422T11 L422T12 L422T13 L422T14 undetermined_fastq"
ids="undetermined_S0"
for reads in ${ids[@]}; do
    R1=${reads}/${reads}.R1.fastq.gz
    R2=${reads}/${reads}.R2.fastq.gz
   $flash -t 8 $R1 $R2 -o merged/$reads -M 250
done

