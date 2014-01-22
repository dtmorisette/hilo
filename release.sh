#!/bin/bash
LANG=C
COMMON_DIR="../MatlabUtilities/"
COMMON_FILES=(FindLinearFit.m GetConstants.m GetMaterial.m ImportCSV.m SaveAsCSV.m dydx.m)
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

cd $COMMON_DIR
com_version=`git describe`
cd ../hilo
hilo_version=`git describe`
sed 's/\$VERSION_STRING\$/'"hilo: $hilo_version; common: $com_version"'/' <hilo.m >tmp/hilo/hilo.m

cd tmp
zip -rq "hilo_$hilo_version.zip" hilo/* -x "*.DS_Store"
mv *.zip ../
cd ..
rm -rf tmp
