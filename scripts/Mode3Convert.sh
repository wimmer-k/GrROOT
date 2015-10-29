#!/bin/bash

input=$1
settings=$2
output=$3

tempfile=$(dirname $output)/tmp.$RANDOM.root

GrROOT -i $input -o $tempfile -s $settings
MakeMode2 -i $tempfile -o $output -s $settings
rm -f $tempfile