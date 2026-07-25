#!/bin/bash


for i in {151..160}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


