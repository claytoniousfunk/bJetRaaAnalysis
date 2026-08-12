#pragma once
// Thin wrapper kept so existing #include lines keep working.
//
// The mapping itself lives in getCentBinImpl.h and uses half-open [lo,hi)
// intervals. Class count follows the centrality header in scope, as it always did here (4 with centrality_4CentBins.h).
#include "getCentBinImpl.h"

inline int getCentBin(int hiBin){ return getCentBinUpTo(hiBin, NCentralityIndices - 1); }
