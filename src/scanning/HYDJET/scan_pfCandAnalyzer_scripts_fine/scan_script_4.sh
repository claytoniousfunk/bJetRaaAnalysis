#!/bin/bash


for i in {31..40}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


