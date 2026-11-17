#pragma once
// ATLAS b-jet cross section in pp at 5.02 TeV (2017, 260 pb^-1), anti-kt R = 0.4,
// |y| < 2.1, d^2 sigma / dpT dy in pb/GeV (per unit y, averaged over |y| < 2.1).
// arXiv:2204.13530v2, Eur. Phys. J. C 83 (2023) 438, Figure 7 left (red points).
// PDF: ~/Documents/nuclear/Literature/HIN-23-008/2204.13530v2.pdf
//
// EXTRACTED FROM THE VECTOR FIGURE, 2026-11-16 (page 19, `mutool draw -F trace`),
// as ATLAS_bJetRAA_2204_13530.h: no table in the paper, no HEPData on disk.
// x axis calibrated on the 20 GeV ticks, y (log) on the decade ticks; value = the
// marker centre, sysHi / sysLo = the red box (relative), stat = the bar where it
// leaves the marker (last three bins; elsewhere hidden, below the marker radius,
// 7%: statBelow). Expect ~0.5% precision on the values.
//
// DEFINITION, differs from ours: b jet = b hadron (pT > 5 GeV) within DeltaR < 0.3;
// measured with muon-tagged jets, unfolded to truth b jets whose pT INCLUDES the
// muon and neutrino momenta, then divided by B(b -> mu X) = 20.6% and corrected
// for muons below 3 GeV. Our gen jets carry no neutrinos, so at fixed pT the
// ATLAS cross section is higher. Sized with Pythia 8 (tools/compareBJetDefinitions_pythia.py,
// 1.6M events, 2026-11-16): ATLAS definition (|y| < 2.1, per unit y) / ours (visible
// jets, |eta| < 1.6, per unit eta) = 1.25 / 1.25 / 1.23 / 1.20 (+- 0.05) at 80-100 /
// 100-120 / 120-150 / 150-200 GeV. Pieces: acceptance 0.90-0.84, neutrinos in every b
// jet x1.23-1.27, muon-jet sample / B instead of all b jets on top of that.

namespace atlasBXsecPP {

const int    nPt = 10;
const double ptEdges[nPt + 1] = {80.7, 90.4, 101.3, 113.4, 127.1, 142.3, 159.4, 178.5, 200.0, 240.0, 290.0};
const double xsec[nPt]  = {230.4, 133.2, 73.41, 42.54, 21.49, 11.33, 6.265, 3.017, 1.385, 0.368};   // pb/GeV
const double sysHi[nPt] = {0.52, 0.27, 0.15, 0.15, 0.20, 0.22, 0.23, 0.22, 0.16, 0.15};             // relative
const double sysLo[nPt] = {0.38, 0.24, 0.14, 0.14, 0.18, 0.19, 0.21, 0.20, 0.15, 0.15};
const double stat[nPt]  = {-1., -1., -1., -1., -1., -1., -1., 0.11, 0.15, 0.27};                    // relative; -1 = hidden
const double statBelow  = 0.07;

}
