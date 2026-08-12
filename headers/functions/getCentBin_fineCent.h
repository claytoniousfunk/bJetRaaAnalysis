#pragma once
// Standalone 11-class scheme with its own hardcoded edges (last class 100-180,
// wider than the others and extending past the usual 160 cut). Unlike the other
// getCentBin_*.h files this one does not read centEdges, so it keeps a local
// edge array rather than including getCentBinImpl.h.
//
// Converted to the same half-open [lo,hi) convention as getCentBinImpl.h.

inline int getCentBin(int hiBin)
{
  static const int fineCentEdges[12] =
    {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 180};
  if(hiBin < fineCentEdges[0]) return -1;
  for(int i = 1; i < 12; i++)
    if(hiBin < fineCentEdges[i]) return i;
  return -1;
}
