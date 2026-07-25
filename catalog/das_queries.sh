#!/usr/bin/env bash
# Example dasgoclient queries for datasets used in this analysis.
#
# Requires a valid CMS grid proxy (voms-proxy-init --voms cms) and
# dasgoclient on PATH -- available on lxplus / cmsconnect / any CVMFS-mounted
# CMSSW area. These are reference examples to copy/paste, not a script meant
# to be run top to bottom.

# 1. Find the exact registered dataset name -- version/era suffixes drift,
#    so wildcard-search before trusting a hardcoded path.
dasgoclient --query="dataset dataset=/HIHardProbes/HIRun2018A*/AOD"

# 2. Summary stats: event count, file count, total size on disk.
dasgoclient --query="dataset=/HIHardProbes/HIRun2018A-04Apr2019-v1/AOD | grep dataset.nevents, dataset.nfiles, dataset.size"

# 3. Full file list -- what CRAB actually reads as input.
dasgoclient --query="file dataset=/HIHardProbes/HIRun2018A-04Apr2019-v1/AOD"

# 4. Per-file event counts -- cross-check against local forest output counts.
dasgoclient --query="file dataset=/HIHardProbes/HIRun2018A-04Apr2019-v1/AOD | grep file.name, file.nevents"

# 5. Runs covered by the dataset -- intersect against certificationFiles/*.txt
#    to see what fraction is actually certified good-for-physics.
dasgoclient --query="run dataset=/HIHardProbes/HIRun2018A-04Apr2019-v1/AOD"

# 6. Same pattern for MC: just a different primary dataset name and tier.
dasgoclient --query="dataset dataset=/BJet_pThat-15_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8*/*/AODSIM"
