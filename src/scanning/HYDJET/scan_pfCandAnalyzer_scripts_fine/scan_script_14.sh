#!/bin/bash


for i in {131..140}

do

    root -l  -b -q "run_pfCandAnalyzer.C($i)"

done

      


