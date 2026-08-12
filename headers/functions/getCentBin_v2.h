#pragma once
// Thin wrapper kept so existing #include lines keep working.
//
// The mapping itself lives in getCentBinImpl.h and uses half-open [lo,hi)
// intervals. Keeps its historical 2-class behaviour: only the first two classes of the edge array are used, everything above returns -1.
#include "getCentBinImpl.h"

inline int getCentBin(int hiBin){ return getCentBinUpTo(hiBin, 2); }
