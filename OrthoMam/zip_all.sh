#!/usr/bin/env bash
for file in Experiments/*/*.chain
do
gzip -f "${file}"
done
