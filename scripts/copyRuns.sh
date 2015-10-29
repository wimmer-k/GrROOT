#!/bin/bash

if [ $# -lt 2 ]
then
  echo "Not enough input parameters"
  echo "Usage: copyRuns.sh RUNDIR1 [...] OUTPUTDIR"
  exit
fi

SRC=${@:1:$#-1}
DEST=${@: -1}

for FOLDER in $SRC
do
  echo "Copying $FOLDER to $DEST/$FOLDER"
  mkdir $DEST/$FOLDER
  cp $FOLDER/Log $FOLDER/Setpoints $FOLDER/Statistics $DEST/$FOLDER
  gzip -c3 $FOLDER/Global.dat > $DEST/$FOLDER/Global.dat.gz
  gzip -c3 $FOLDER/GlobalRaw.dat > $DEST/$FOLDER/GlobalRaw.dat.gz
  echo "Done copying $FOLDER"
done