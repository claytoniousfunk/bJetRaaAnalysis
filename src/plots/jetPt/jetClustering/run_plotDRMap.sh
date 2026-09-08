#!/bin/bash

for cent in {1..16}
do

    for jet in {1..6}
    do

	root -l -b -q "plotDRMap.C($cent,$jet)"
	       
    done

done 
