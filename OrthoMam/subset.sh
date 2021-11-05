#!/usr/bin/env bash
for FILE in ./Experiments/*/sitemutsel_1.run
do 
cp --parents $FILE ..
done
echo "Created folder"
for FOLDER in ../Experiments/*
do 
rm -rf "$FOLDER"/sitemutsel_1.run
done
echo "Removed subfolder"
for FILE in ./Experiments/*/*.run.ci0.05.tsv
do 
cp $FILE .$FILE
done
echo "Copied files"
for FILE in ./Experiments/*/*.run.ci0.0025.tsv
do 
cp $FILE .$FILE
done
echo "Copied files"
