#!/bin/bash


for i in {111..120}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


