#!/usr/bin/env bash
for param in ./Experiments/*/*.param; do
  sed -i "s/.run	1	1000	/.run	1	2000	/g" ${param}
done
