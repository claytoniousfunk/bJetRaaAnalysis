#pragma once
// Thin wrapper kept so existing #include lines keep working.
//
// The mapping itself lives in getCentBinImpl.h and uses half-open [lo,hi)
// intervals. Keeps its historical 16-class behaviour (ultra-fine).
#include "getCentBinImpl.h"

inline int getCentBin(int hiBin){ return getCentBinUpTo(hiBin, 16); }
