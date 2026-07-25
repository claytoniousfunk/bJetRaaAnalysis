#!/bin/bash


for i in {251..260}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


