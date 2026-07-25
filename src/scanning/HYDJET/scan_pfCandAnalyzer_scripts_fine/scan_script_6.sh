#!/bin/bash


for i in {51..60}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


