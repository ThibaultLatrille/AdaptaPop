#!/usr/bin/env bash
for file in Experiments/*/*.chain.gz
do
gunzip "${file}"
done
