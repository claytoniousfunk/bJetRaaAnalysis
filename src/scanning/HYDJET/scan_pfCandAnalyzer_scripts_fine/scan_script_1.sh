#!/bin/bash


for i in {1..10}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


