#!/bin/bash


for i in {101..110}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


