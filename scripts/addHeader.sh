#!/bin/bash

if [ $# -lt 2 ]
then
  echo "Not enough input parameters"
  echo "Usage: addHeader.sh HEADER SOURCE [...]"
  exit
fi

HEADER=${@:1:1}
SOURCEFILES=${@:2:$#}

for FILE in $SOURCEFILES
do
  echo "Adding $HEADER at the top of $FILE"
  mv $FILE temp.bak
  cp $HEADER $FILE
  cat temp.bak >> $FILE
  rm temp.bak
  echo "Done modifying $FILE"
done