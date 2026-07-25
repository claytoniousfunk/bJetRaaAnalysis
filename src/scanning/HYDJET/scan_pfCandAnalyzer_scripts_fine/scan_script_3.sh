#!/bin/bash


for i in {21..30}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


