#!/bin/bash


for i in {161..170}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


