#!/bin/bash


for i in {11..20}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


