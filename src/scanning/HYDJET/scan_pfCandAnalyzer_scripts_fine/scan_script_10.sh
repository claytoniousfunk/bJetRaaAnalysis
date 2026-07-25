#!/bin/bash


for i in {91..100}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


