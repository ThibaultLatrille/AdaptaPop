#!/usr/bin/env bash
# cd /pandata/tlatrill/AdaptaPop/data/pb_mutsel/
cd /home/thibault/AdaptaPop/data/pb_mutsel/
find /home/thibault/AdaptaPop/data/pb_mutsel/ -type f -print0 | xargs -0 sed -i 's/pandata\/tlatrill/home\/thibault/g'
for file in ./*.param
do
/home/thibault/pbmpi2/data/readpb_mpi -x 100 -ss "${file%.*}"
/home/thibault/pbmpi2/data/readpb_mpi -x 100 -om "${file%.*}"
done
cd /home/thibault/AdaptaPop/data/pb_mutselfreeomega
find /home/thibault/AdaptaPop/data/pb_mutselfreeomega/ -type f -print0 | xargs -0 sed -i 's/pandata\/tlatrill/home\/thibault/g'
for file in ./*.param
do
/home/thibault/pbmpi2/data/readpb_mpi -x 100 -ss "${file%.*}"
/home/thibault/pbmpi2/data/readpb_mpi -x 100 -om "${file%.*}"
done
