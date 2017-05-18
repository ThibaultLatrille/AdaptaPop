#!/usr/bin/env bash

cd /mnt/sda1/AdaptaPop/data/pb_cleaned_mutsel/
for file in ./*.param
do
find $file -type f -print0 | xargs -0 sed -i 's/pandata\/tlatrill/mnt\/sda1/g'
/home/thibault/Tools/pbmpi2/data/readpb_mpi -x 100 -om "${file%.*}"
done

cd /mnt/sda1/AdaptaPop/data/pb_cleaned_mutselfreeomega/
for file in ./*.param
do
find $file -type f -print0 | xargs -0 sed -i 's/pandata\/tlatrill/mnt\/sda1/g'
/home/thibault/Tools/pbmpi2/data/readpb_mpi -x 100 -om "${file%.*}"
done
