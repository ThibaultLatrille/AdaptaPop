#!/usr/bin/env bash
for EXP in ./*;
do
  for FOLDER in "${EXP}"/*;
  do
    for FASTA in "${FOLDER}"/*.fasta;
    do
      gzip "${FASTA}"
    done
    for ERRORS in "${FOLDER}"/*.errors;
    do
      gzip "${ERRORS}"
    done
  done
done