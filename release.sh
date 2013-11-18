#!/bin/bash
LANG=C
COMMON_DIR="../Common/"
COMMON_FILES=(FindLinearFit.m GetConstants.m GetMaterial.m ImportCSV.m SaveAsCSV.m)
SRC_FILES=(example.m TranslateData.m)
OTHER_FILES=(README LICENSE data.csv)

files=()
for fname in "${COMMON_FILES[@]}" 
do
	files+=("$COMMON_DIR$fname")
done
rm -rf tmp
mkdir -p tmp/hilo
mkdir tmp/hilo/processed
cp ${files[*]} ${SRC_FILES[*]} ${OTHER_FILES[*]} tmp/hilo

version=`git describe`
sed 's/\$VERSION_STRING\$/'"$version"'/' <hilo.m >tmp/hilo/hilo.m

cd tmp
zip -rq "hilo_$version.zip" hilo/* -x "*.DS_Store"
mv *.zip ../
cd ..
rm -rf tmp
