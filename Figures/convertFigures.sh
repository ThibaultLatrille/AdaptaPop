for file in ./*
do
filename=$(basename "$file")
extension="${filename##*.}"
body="${filename%.*}"
echo $body
echo $extension
if  [ "$extension" = "svg" ]
then
  inkscape -D -z --file=$filename --export-pdf=$body.pdf
fi
done

