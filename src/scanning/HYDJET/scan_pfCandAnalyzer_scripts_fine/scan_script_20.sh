#!/bin/bash


for i in {191..200}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


