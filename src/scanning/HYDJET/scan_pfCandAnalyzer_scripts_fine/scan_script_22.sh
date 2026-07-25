#!/bin/bash


for i in {211..220}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


