#!/bin/bash

export PATH=/usr/bin:$PATH

rm extract_data_fullerenes.x

gfortran-mp-7 -o extract_data_fullerenes.x extract_data_fullerenes.f90

echo " "
echo " The compilation was properly performed. "
echo " "

./extract_data_fullerenes.x