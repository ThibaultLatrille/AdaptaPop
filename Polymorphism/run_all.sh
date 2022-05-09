#### Local
for FOLDER in ./*; do
  if [ -d "$FOLDER" ]; then
    cd "$FOLDER"
    snakemake -j 8 -k --printshellcmds
    cd ..
  fi
done
