#pragma once
//
// Single source of truth for the hiBin -> centrality-class mapping.
//
// Convention: half-open intervals, [ centEdges[i-1], centEdges[i] ).
//
//   * every class is exactly (hi - lo) hiBin units wide
//   * no hiBin value belongs to two classes
//   * no hiBin value inside the range is silently dropped
//
// Returns 1 .. nClasses, or -1 for hiBin below centEdges[0] or at/above
// centEdges[nClasses].
//
// The edges come from whichever centrality header is in scope (centEdges,
// NCentralityIndices). nClasses is passed in by the thin wrapper headers so
// that each of them keeps the class count it has always had -- the wrappers
// used to hand-code both the count AND the interval convention, and it was the
// convention that had drifted apart between them, not the count.
//
// PREVIOUS CONVENTION. The old hand-written versions used two different rules:
//
//   A  hiBin >= lo && hiBin <  hi      (getCentBin_v2, _v3, _fineCent)
//   B  hiBin >  lo && hiBin <= hi      (getCentBin, _V2p4), with the FIRST
//                                       bin inclusive on both ends
//
// They disagree at every boundary, and B's first bin was one unit wider than
// its own label. This header adopts A. Anything produced with B put each
// boundary value in the class below, so ROOT files written before this change
// are offset by one hiBin unit at every boundary relative to files written
// after it. The convention is not recorded inside the output -- histogram
// titles are built from centEdges either way -- so a mixed set cannot be told
// apart after the fact. Regenerate rather than mix.

inline int getCentBinUpTo(int hiBin, int nClasses)
{
  if(hiBin < centEdges[0]) return -1;
  const int nMax = (NCentralityIndices - 1);
  if(nClasses > nMax) nClasses = nMax;
  for(int i = 1; i <= nClasses; i++)
    if(hiBin < centEdges[i]) return i;
  return -1;   // at or above the top edge of the last class
}
