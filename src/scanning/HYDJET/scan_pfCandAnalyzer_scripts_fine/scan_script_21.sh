#!/bin/bash


for i in {201..210}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


