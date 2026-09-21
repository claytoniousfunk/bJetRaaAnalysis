// Can a calo jet be matched to a PF jet? Geometry and energy, jet by jet.
//
// pp PYTHIA forest, 42k events, with ak4CaloJetAnalyzer/t and ak4PFJetAnalyzer/t
// written from the SAME events, so the two collections can be compared entry by
// entry with no event matching needed.
//
// Direction of the match: calo -> PF. For every calo jet passing the fiducial
// cuts, find the PF jet with the smallest dR. That is the natural direction for
// this analysis, where the calo spectrum is the measurement and the PF response
// matrix would be the correction; it answers "given a calo jet, is there a PF
// jet to associate it with", not the reverse.
//
// The closest PF jet is the closest one above pfPtFloor, with no pT-balance
// requirement. So a bad match shows up as a large dR or a large dpT rather than
// as a failure, which is the point -- the tails are the measurement here.
//
// Figures, all in figures/jetMatching/:
//   dR_vs_ptCalo.pdf            dR of the closest PF jet vs calo jet pT. Each
//                               pT column is normalized to unit area, so the
//                               map reads as the dR distribution at that pT
//                               and is not swamped by the falling spectrum.
//                               dR axes stop at 0.4-0.5 throughout: the match
//                               is tight enough that a wider axis is blank.
//   dR_dists_ptSlices.pdf       the same information as 1D slices, unit area.
//                               Log y: the quantity of interest is the
//                               mismatch rate in the tail, which is orders of
//                               magnitude below the peak.
//   matchEff_vs_ptCalo.pdf      fraction of calo jets whose closest PF jet is
//                               within dR < 0.1, 0.2, 0.4. Binomial errors.
//   dR_vs_dpt.pdf               dpT = pT(PF) - pT(calo), absolute, vs dR.
//   dR_vs_dptRel.pdf            the same divided by pT(calo), which is the
//                               comparable quantity across the pT range.
//   dptRel_dists_ptSlices.pdf   dpT/pT for geometrically matched jets only
//                               (dR < dRMatch), unit area, per pT slice.
// Occupancy maps use log z; the two-dimensional densities span several decades.
//
// PT CONVENTION, chosen by ptMode:
//   0  forest jtpt for both. In this forest the calo jets receive the same
//      ~1.09 factor as PF instead of an AK4Calo correction, so calo closes on
//      gen at 0.783 (PF 0.997) and dpT/pT peaks at +0.5 -> +0.1, not at zero.
//   1  rawpt for both.
//   2  calo: Spring18_ppRef5TeV_V6 L2Relative AK4Calo applied to rawpt, the
//      same payload pp_scan.C's manual-JEC path uses; PF: forest jtpt, which
//      already closes. L2L3Residual is left off because it is a data-only
//      residual; on this MC it overshoots to 1.05. With this, calo closes at
//      1.014 and dpT/pT centers at -0.01 to -0.03 -- the Gaussian-at-zero that
//      two properly calibrated collections should give. The small negative
//      offset is the 1.014 vs 0.997 closure difference.
//      CAVEAT: the forest only stores calo jets with jtpt >= 30, i.e. raw
//      >~ 27.5 GeV, which the AK4Calo correction moves to >~ 40 GeV. In this
//      mode the 30-50 slice is threshold-sculpted (641 jets, against 13k in
//      mode 0) and should not be read; 50 GeV and up is complete.
// The printout reports the gen closure of the calo pT actually used.
//
// Usage: root -l -b -q 'plotCaloPFMatching_PYTHIA.C'       // forest jtpt
//        root -l -b -q 'plotCaloPFMatching_PYTHIA.C(1)'    // raw pT
//        root -l -b -q 'plotCaloPFMatching_PYTHIA.C(2)'    // calo AK4Calo JEC
// Run from: src/plots/jetPt/jetCollection/

#include "../../../../headers/plotting/plotStyle.h"
#include "../../../../JetEnergyCorrections/JetCorrector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include "TF1.h"
#include "TLine.h"
#include <cstdio>

const char *forestFile = "/home/clayton/Downloads/HiForestAOD_PYTHIA.root";
const char *outDir     = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetMatching/";
const char *caloL2Rel  = "/home/clayton/Analysis/code/bJetRaaAnalysis/JetEnergyCorrections/"
                         "Spring18_ppRef5TeV_V6_DATA_L2Relative_AK4Calo.txt";

// fiducial: the calo jet is the probe, the PF side is left loose on purpose
const double caloPtFloor = 30.;   // below the 50 GeV analysis floor, to see the turn-on of the matching itself
const double caloEtaMax  = 2.0;   // wider than the analysis etaMax = 1.6, so edge effects are visible
const double pfPtFloor   = 10.;   // "closest PF jet" means closest above this
const double dRMatch     = 0.2;   // geometric match for the dpT distributions

const int    nSlice = 5;
const double sliceLo[nSlice] = { 30.,  50.,  80., 120., 200.};
const double sliceHi[nSlice] = { 50.,  80., 120., 200., 400.};

// Okabe-Ito indices, skipping 4 (yellow), which plotStyle.h flags as unreadable
// on white for lines and markers.
const int sliceColor[nSlice] = {0, 1, 2, 3, 6};
const int sliceMark[nSlice]  = {markFilledCircle, markFilledSquare, markFilledDiamond,
                                markOpenCircle, markOpenSquare};

static const int MAXJET = 500;

static double dPhi(double a, double b)
{
  double d = a - b;
  while(d >  TMath::Pi()) d -= TMath::TwoPi();
  while(d < -TMath::Pi()) d += TMath::TwoPi();
  return d;
}

// normalize every x column of a 2D map to unit area over y
static void columnNormalize(TH2D *h)
{
  for(int ix = 1; ix <= h->GetNbinsX(); ix++){
    double s = 0.;
    for(int iy = 1; iy <= h->GetNbinsY(); iy++) s += h->GetBinContent(ix, iy);
    if(s <= 0.) continue;
    for(int iy = 1; iy <= h->GetNbinsY(); iy++){
      h->SetBinContent(ix, iy, h->GetBinContent(ix, iy) / s);
      h->SetBinError  (ix, iy, h->GetBinError  (ix, iy) / s);
    }
  }
}

static void drawMap(TH2D *h, const char *xT, const char *yT, const char *zT,
                    bool logZ, const char *note, const char *outPath)
{
  gStyle->SetPalette(kViridis);
  TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", 800, 700);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.19);
  c->SetTopMargin(0.13);  c->SetBottomMargin(0.13);
  if(logZ) c->SetLogz();

  h->SetTitle(""); h->SetStats(0);
  h->GetXaxis()->SetTitle(xT); h->GetXaxis()->SetTitleSize(0.045); h->GetXaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetTitle(yT); h->GetYaxis()->SetTitleSize(0.045); h->GetYaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetTitleOffset(1.25);
  h->GetZaxis()->SetTitle(zT); h->GetZaxis()->SetTitleSize(0.040); h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleOffset(1.40);
  h->Draw("COLZ");

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
  la.DrawLatex(0.13, 0.935, "PYTHIA pp 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4");
  la.SetTextSize(0.031);
  la.DrawLatex(0.13, 0.895, note);

  c->SaveAs(outPath);
  delete c;
}

void plotCaloPFMatching_PYTHIA(int ptMode = 0)
{
  if(ptMode < 0 || ptMode > 2){ printf("ERROR: ptMode must be 0, 1 or 2\n"); return; }
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(forestFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", forestFile); return; }

  // Get() picks the highest cycle, which is the complete 42000-entry write; the
  // PF directory also holds a truncated t;1 from an earlier pass.
  TTree *tC = (TTree*) f->Get("ak4CaloJetAnalyzer/t");
  TTree *tP = (TTree*) f->Get("ak4PFJetAnalyzer/t");
  if(!tC || !tP){ printf("ERROR: a jet tree is missing\n"); return; }
  if(tC->GetEntries() != tP->GetEntries()){
    printf("ERROR: trees are not the same length (%lld vs %lld) -- not the same events\n",
           tC->GetEntries(), tP->GetEntries());
    return;
  }

  const char *ptLabelV[3] = {"forest-corrected", "raw", "AK4Calo-corrected calo"};
  const char *ptTagV[3]   = {"", "_rawPt", "_caloJEC"};
  const char *ptLabel = ptLabelV[ptMode];
  const char *ptTag   = ptTagV[ptMode];

  std::vector<std::string> jecFiles = {caloL2Rel};
  JetCorrector caloJEC(jecFiles);

  Int_t nC = 0, nP = 0;
  Float_t etaC[MAXJET], phiC[MAXJET], refEtaC[MAXJET], refPhiC[MAXJET], refPtC[MAXJET];
  Float_t etaP[MAXJET], phiP[MAXJET], refEtaP[MAXJET], refPhiP[MAXJET];
  Float_t corrC[MAXJET], rawC[MAXJET], corrP[MAXJET], rawP[MAXJET];
  Float_t manC[MAXJET];   // AK4Calo L2Relative on rawpt, filled per event

  // Each branch is addressed exactly once: a second SetBranchAddress on the
  // same branch silently replaces the first, leaving the earlier array unfilled.
  // ptC/ptP are aliases into whichever convention is selected.
  tC->SetBranchAddress("nref", &nC);
  tC->SetBranchAddress("jteta", etaC);
  tC->SetBranchAddress("jtphi", phiC);
  tC->SetBranchAddress("refeta", refEtaC);
  tC->SetBranchAddress("refphi", refPhiC);
  tC->SetBranchAddress("refpt", refPtC);
  tC->SetBranchAddress("jtpt", corrC);
  tC->SetBranchAddress("rawpt", rawC);

  tP->SetBranchAddress("nref", &nP);
  tP->SetBranchAddress("jteta", etaP);
  tP->SetBranchAddress("jtphi", phiP);
  tP->SetBranchAddress("refeta", refEtaP);
  tP->SetBranchAddress("refphi", refPhiP);
  tP->SetBranchAddress("jtpt", corrP);
  tP->SetBranchAddress("rawpt", rawP);

  Float_t * const ptC = ptMode == 1 ? rawC : (ptMode == 2 ? manC : corrC);
  Float_t * const ptP = ptMode == 1 ? rawP : corrP;   // PF already closes; mode 2 leaves it alone

  double sumClos = 0.; long nClos = 0;   // gen closure of the calo pT in use

  // --- histograms -----------------------------------------------------------
  const double ptLo = caloPtFloor, ptHi = 400.;
  const int    nPt  = (int)((ptHi - ptLo) / 5.);

  // dR ranges stop at 0.5: the matching is so tight that a wider axis is empty
  // canvas. Overflow still counts in the denominators and the printout.
  TH2D *h2_dR_pt   = new TH2D("h2_dR_pt", "", nPt, ptLo, ptHi, 80, 0., 0.4);
  TH2D *h2_dR_dpt  = new TH2D("h2_dR_dpt", "", 50, 0., 0.5, 120, -120., 120.);
  TH2D *h2_dR_dptR = new TH2D("h2_dR_dptR", "", 50, 0., 0.5, 125, -1., 1.5);
  h2_dR_pt->Sumw2(); h2_dR_dpt->Sumw2(); h2_dR_dptR->Sumw2();

  TH1D *hDR[nSlice], *hDPtR[nSlice];
  for(int s = 0; s < nSlice; s++){
    hDR[s]   = new TH1D(Form("hDR_s%d", s),   "", 100, 0., 0.5);
    // symmetric enough about zero for mode 2, wide enough for mode 0's +0.5 peak
    hDPtR[s] = new TH1D(Form("hDPtR_s%d", s), "", 100, -1.0, 1.5);
    hDR[s]->Sumw2(); hDPtR[s]->Sumw2();
  }

  const double dRCut[3] = {0.1, 0.2, 0.4};
  TH1D *hNum[3], *hDen;
  hDen = new TH1D("hDen", "", nPt, ptLo, ptHi);
  for(int i = 0; i < 3; i++) hNum[i] = new TH1D(Form("hNum_%d", i), "", nPt, ptLo, ptHi);

  // per-slice bookkeeping
  long nJet[nSlice] = {0}, nNoPF[nSlice] = {0}, nSameGen[nSlice] = {0},
       nGenKnown[nSlice] = {0}, nSplit[nSlice] = {0}, nMatched[nSlice] = {0};
  double sumDPtR[nSlice] = {0}, sumDPtR2[nSlice] = {0},
         sumRatioCorr[nSlice] = {0}, sumRatioRaw[nSlice] = {0};

  const Long64_t nEv = tC->GetEntries();
  for(Long64_t ie = 0; ie < nEv; ie++){
    tC->GetEntry(ie);
    tP->GetEntry(ie);

    if(ptMode == 2){
      for(int i = 0; i < nC && i < MAXJET; i++){
        caloJEC.SetJetPT(rawC[i]); caloJEC.SetJetEta(etaC[i]); caloJEC.SetJetPhi(phiC[i]);
        manC[i] = caloJEC.GetCorrectedPT();
      }
    }

    for(int i = 0; i < nC; i++){
      if(fabs(etaC[i]) > caloEtaMax) continue;
      if(refPtC[i] > 50.){ sumClos += ptC[i] / refPtC[i]; nClos++; }
      if(ptC[i] < caloPtFloor) continue;

      // closest PF jet above the floor, and how many sit within 0.4
      int    best = -1;
      double bestDR = 1e9;
      int    nNear = 0;
      for(int j = 0; j < nP; j++){
        if(ptP[j] < pfPtFloor) continue;
        double de = etaC[i] - etaP[j], dp = dPhi(phiC[i], phiP[j]);
        double dr = sqrt(de*de + dp*dp);
        if(dr < 0.4) nNear++;
        if(dr < bestDR){ bestDR = dr; best = j; }
      }

      int s = -1;
      for(int k = 0; k < nSlice; k++)
        if(ptC[i] >= sliceLo[k] && ptC[i] < sliceHi[k]){ s = k; break; }
      if(s < 0) continue;
      nJet[s]++;

      hDen->Fill(ptC[i]);

      if(best < 0){ nNoPF[s]++; continue; }   // no PF jet at all: counts against the efficiency

      double dpt  = ptP[best] - ptC[i];
      double dptR = dpt / ptC[i];

      h2_dR_pt->Fill(ptC[i], bestDR);
      h2_dR_dpt ->Fill(bestDR, dpt);
      h2_dR_dptR->Fill(bestDR, dptR);
      hDR[s]->Fill(bestDR);
      for(int k = 0; k < 3; k++) if(bestDR < dRCut[k]) hNum[k]->Fill(ptC[i]);
      if(nNear >= 2) nSplit[s]++;

      if(bestDR < dRMatch){
        nMatched[s]++;
        hDPtR[s]->Fill(dptR);
        sumDPtR[s]  += dptR;
        sumDPtR2[s] += dptR*dptR;
        if(corrC[i] > 0.) sumRatioCorr[s] += corrP[best] / corrC[i];
        if(rawC[i]  > 0.) sumRatioRaw[s]  += rawP[best]  / rawC[i];

        // do the two reco jets point at the same gen jet? refeta/refphi are the
        // forest's own gen match, so equality means a common gen jet.
        if(refEtaC[i] > -900. && refEtaP[best] > -900.){
          nGenKnown[s]++;
          if(fabs(refEtaC[i] - refEtaP[best]) < 1e-4 &&
             fabs(dPhi(refPhiC[i], refPhiP[best])) < 1e-4) nSameGen[s]++;
        }
      }
    }
  }

  // --- printout -------------------------------------------------------------
  printf("\n  calo -> PF matching, PYTHIA pp, %lld events, %s pT\n", nEv, ptLabel);
  printf("  calo probe: pT > %.0f, |eta| < %.1f   PF: closest with pT > %.0f\n\n",
         caloPtFloor, caloEtaMax, pfPtFloor);
  printf("  %-11s %8s %7s | %8s %8s %8s | %8s %9s | %8s %8s %8s\n",
         "pT(calo)", "N", "no PF", "dR<0.1", "dR<0.2", "dR<0.4",
         "medianDR", "split", "<dpT/pT>", "RMS", "sameGen");
  for(int s = 0; s < nSlice; s++){
    if(!nJet[s]) continue;
    double q[1], p[1] = {0.5};
    hDR[s]->GetQuantiles(1, q, p);
    double f1 = 0., f2 = 0., f4 = 0.;
    long n1 = 0, n2 = 0, n4 = 0;
    for(int b = 1; b <= hDR[s]->GetNbinsX(); b++){
      double c = hDR[s]->GetBinContent(b), up = hDR[s]->GetXaxis()->GetBinUpEdge(b);
      if(up <= 0.1 + 1e-9) n1 += (long)c;
      if(up <= 0.2 + 1e-9) n2 += (long)c;
      if(up <= 0.4 + 1e-9) n4 += (long)c;
    }
    f1 = (double)n1/nJet[s]; f2 = (double)n2/nJet[s]; f4 = (double)n4/nJet[s];
    double mean = nMatched[s] ? sumDPtR[s]/nMatched[s] : 0.;
    double rms  = nMatched[s] ? sqrt(TMath::Max(0., sumDPtR2[s]/nMatched[s] - mean*mean)) : 0.;
    printf("  %4.0f-%-6.0f %8ld %6.2f%% | %7.3f%% %7.3f%% %7.3f%% | %8.3f %9.3f | %+8.3f %8.3f %7.2f%%\n",
           sliceLo[s], sliceHi[s], nJet[s], 100.*nNoPF[s]/nJet[s],
           100.*f1, 100.*f2, 100.*f4, q[0], (double)nSplit[s]/nJet[s],
           mean, rms, nGenKnown[s] ? 100.*nSameGen[s]/nGenKnown[s] : 0.);
  }
  printf("\n  'split' is the fraction of calo jets with two or more PF jets within dR < 0.4,\n");
  printf("  i.e. how often the association is ambiguous rather than one-to-one.\n");
  printf("  'sameGen' is, among dR < %.1f matches, the fraction sharing one gen jet.\n", dRMatch);

  printf("\n  mean pT(PF)/pT(calo) for dR < %.1f matches, both conventions\n", dRMatch);
  printf("  %-11s %12s %12s\n", "pT(calo)", "forest jtpt", "raw");
  for(int s = 0; s < nSlice; s++){
    if(!nMatched[s]) continue;
    printf("  %4.0f-%-6.0f %12.4f %12.4f\n", sliceLo[s], sliceHi[s],
           sumRatioCorr[s]/nMatched[s], sumRatioRaw[s]/nMatched[s]);
  }

  // Why the two agree: the correction in this forest is the same for both
  // collections, so it cancels in the ratio. Closing each collection against
  // its own gen jet is what exposes that, and the calo side does not close.
  printf("\n  JEC closure against the forest's own gen match, refpt > 50 GeV, |eta| < %.1f\n", caloEtaMax);
  printf("  %-8s %14s %14s %14s\n", "coll", "<jtpt/rawpt>", "<jtpt/refpt>", "<rawpt/refpt>");
  for(const char *dn : {"ak4CaloJetAnalyzer", "ak4PFJetAnalyzer"}){
    TTree *t = (TTree*) f->Get(Form("%s/t", dn));
    if(!t) continue;
    const char *cut = Form("refpt > 50. && abs(jteta) < %f", caloEtaMax);
    double v[3];
    const char *ex[3] = {"jtpt/rawpt", "jtpt/refpt", "rawpt/refpt"};
    for(int i = 0; i < 3; i++){
      t->Draw(Form("%s>>hTmp(400,0.,4.)", ex[i]), cut, "goff");
      v[i] = ((TH1*) gDirectory->Get("hTmp"))->GetMean();
    }
    printf("  %-8s %14.4f %14.4f %14.4f\n",
           strstr(dn, "Calo") ? "calo" : "PF", v[0], v[1], v[2]);
  }
  printf("  calo pT as used in this run (%s): <pT/refpt> = %.4f\n",
         ptLabel, nClos ? sumClos/nClos : -1.);

  // If both collections are calibrated, dpT/pT is centered on zero with a width
  // set by the two resolutions. Core Gaussian fit, +-1.2 RMS around the peak.
  printf("\n  (pT_PF - pT_calo)/pT_calo for dR < %.1f matches, %s pT\n", dRMatch, ptLabel);
  printf("  %-11s %8s %8s %8s %8s\n", "pT(calo)", "median", "gaus mu", "gaus sig", "RMS");
  for(int s = 0; s < nSlice; s++){
    if(hDPtR[s]->GetEntries() < 50) continue;
    double q, p = 0.5; hDPtR[s]->GetQuantiles(1, &q, &p);
    double pk = hDPtR[s]->GetBinCenter(hDPtR[s]->GetMaximumBin()), r = hDPtR[s]->GetRMS();
    TF1 fg(Form("fg_%d", s), "gaus", pk - 1.2*r, pk + 1.2*r);
    hDPtR[s]->Fit(&fg, "QNR0");
    printf("  %4.0f-%-6.0f %+8.3f %+8.3f %8.3f %8.3f\n", sliceLo[s], sliceHi[s],
           q, fg.GetParameter(1), fg.GetParameter(2), r);
  }

  // --- figures --------------------------------------------------------------
  // log z even after column normalization: the core bin holds a few tens of
  // percent and the tail bins 1e-4, so a linear scale shows the core only.
  columnNormalize(h2_dR_pt);
  drawMap(h2_dR_pt, "calo jet #it{p}_{T} [GeV]", "#Delta#it{R}(calo, closest PF)",
          "fraction of calo jets per column", true,
          Form("closest PF jet, %s #it{p}_{T}, each #it{p}_{T} column normalized", ptLabel),
          Form("%sdR_vs_ptCalo%s.pdf", outDir, ptTag));

  drawMap(h2_dR_dpt, "#Delta#it{R}(calo, closest PF)",
          "#Delta#it{p}_{T} = #it{p}_{T}^{PF} #minus #it{p}_{T}^{calo} [GeV]",
          "calo jets", true,
          Form("%s #it{p}_{T}, calo #it{p}_{T} > %.0f GeV", ptLabel, caloPtFloor),
          Form("%sdR_vs_dpt%s.pdf", outDir, ptTag));

  drawMap(h2_dR_dptR, "#Delta#it{R}(calo, closest PF)",
          "(#it{p}_{T}^{PF} #minus #it{p}_{T}^{calo}) / #it{p}_{T}^{calo}",
          "calo jets", true,
          Form("%s #it{p}_{T}, calo #it{p}_{T} > %.0f GeV", ptLabel, caloPtFloor),
          Form("%sdR_vs_dptRel%s.pdf", outDir, ptTag));

  // dR distributions, unit area, log y for the tail
  {
    TCanvas *c = new TCanvas("c_dR", "", 700, 800);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.12);
    c->SetTopMargin(0.13);  c->SetRightMargin(0.05);
    c->SetLogy();

    TLegend *leg = makeLegend(0.50, 0.55, 0.93, 0.80, 0.033);
    double yMax = 0.;
    for(int s = 0; s < nSlice; s++){
      if(hDR[s]->Integral() > 0.) hDR[s]->Scale(1./hDR[s]->Integral());
      yMax = TMath::Max(yMax, hDR[s]->GetMaximum());
    }
    for(int s = 0; s < nSlice; s++){
      styleH(hDR[s], okabeHex[sliceColor[s]], sliceMark[s], 0.9);
      if(s == 0){
        hDR[s]->SetTitle(""); hDR[s]->SetStats(0);
        hDR[s]->SetMinimum(2e-5); hDR[s]->SetMaximum(yMax * 8.);
        hDR[s]->GetXaxis()->SetTitle("#Delta#it{R}(calo, closest PF)");
        hDR[s]->GetXaxis()->SetTitleSize(0.045); hDR[s]->GetXaxis()->SetLabelSize(0.040);
        hDR[s]->GetYaxis()->SetTitle("fraction of calo jets");
        hDR[s]->GetYaxis()->SetTitleSize(0.045); hDR[s]->GetYaxis()->SetLabelSize(0.040);
        hDR[s]->GetYaxis()->SetTitleOffset(1.45);
        hDR[s]->Draw("E1");
      } else hDR[s]->Draw("E1 same");
      leg->AddEntry(hDR[s], Form("%.0f < #it{p}_{T}^{calo} < %.0f GeV", sliceLo[s], sliceHi[s]), "lp");
    }
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.15, 0.935, "PYTHIA pp 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4");
    la.SetTextSize(0.031);
    la.DrawLatex(0.15, 0.895, Form("calo #rightarrow closest PF jet, %s #it{p}_{T}", ptLabel));

    c->SaveAs(Form("%sdR_dists_ptSlices%s.pdf", outDir, ptTag));
    delete c;
  }

  // matching efficiency vs pT
  {
    TCanvas *c = new TCanvas("c_eff", "", 700, 800);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.12);
    c->SetTopMargin(0.13);  c->SetRightMargin(0.05);

    TLegend *leg = makeLegend(0.50, 0.28, 0.94, 0.48, 0.034);
    const int mk[3] = {markFilledCircle, markFilledSquare, markFilledDiamond};
    for(int k = 0; k < 3; k++){
      TH1D *e = (TH1D*) hNum[k]->Clone(Form("eff_%d", k));
      e->Divide(hNum[k], hDen, 1., 1., "B");
      styleH(e, okabeHex[k], mk[k], 1.0);
      if(k == 0){
        e->SetTitle(""); e->SetStats(0);
        // zoomed: all three curves sit above 0.9, and the only structure is the
        // low pT turn-on of the dR < 0.1 curve
        e->SetMinimum(0.86); e->SetMaximum(1.03);
        e->GetXaxis()->SetTitle("calo jet #it{p}_{T} [GeV]");
        e->GetXaxis()->SetTitleSize(0.045); e->GetXaxis()->SetLabelSize(0.040);
        e->GetYaxis()->SetTitle("fraction of calo jets with a PF jet within #Delta#it{R}");
        e->GetYaxis()->SetTitleSize(0.040); e->GetYaxis()->SetLabelSize(0.040);
        e->GetYaxis()->SetTitleOffset(1.65);
        e->Draw("E1");
      } else e->Draw("E1 same");
      leg->AddEntry(e, Form("#Delta#it{R} < %.1f", dRCut[k]), "lp");
    }
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.15, 0.935, "PYTHIA pp 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4");
    la.SetTextSize(0.031);
    la.DrawLatex(0.15, 0.895, Form("denominator: all calo jets, |#eta| < %.1f, %s #it{p}_{T}", caloEtaMax, ptLabel));

    c->SaveAs(Form("%smatchEff_vs_ptCalo%s.pdf", outDir, ptTag));
    delete c;
  }

  // dpT/pT distributions for geometrically matched jets
  {
    TCanvas *c = new TCanvas("c_dpt", "", 700, 800);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.12);
    c->SetTopMargin(0.13);  c->SetRightMargin(0.05);

    TLegend *leg = makeLegend(0.54, 0.56, 0.94, 0.81, 0.033);
    double yMax = 0.;
    for(int s = 0; s < nSlice; s++){
      if(hDPtR[s]->Integral() > 0.) hDPtR[s]->Scale(1./hDPtR[s]->Integral());
      yMax = TMath::Max(yMax, hDPtR[s]->GetMaximum());
    }
    for(int s = 0; s < nSlice; s++){
      styleH(hDPtR[s], okabeHex[sliceColor[s]], sliceMark[s], 0.9);
      if(s == 0){
        hDPtR[s]->SetTitle(""); hDPtR[s]->SetStats(0);
        hDPtR[s]->SetMinimum(0.); hDPtR[s]->SetMaximum(yMax * 1.55);
        hDPtR[s]->GetXaxis()->SetTitle("(#it{p}_{T}^{PF} #minus #it{p}_{T}^{calo}) / #it{p}_{T}^{calo}");
        hDPtR[s]->GetXaxis()->SetTitleSize(0.045); hDPtR[s]->GetXaxis()->SetLabelSize(0.040);
        hDPtR[s]->GetYaxis()->SetTitle("fraction of matched calo jets");
        hDPtR[s]->GetYaxis()->SetTitleSize(0.045); hDPtR[s]->GetYaxis()->SetLabelSize(0.040);
        hDPtR[s]->GetYaxis()->SetTitleOffset(1.45);
        hDPtR[s]->Draw("E1");
      } else hDPtR[s]->Draw("E1 same");
      leg->AddEntry(hDPtR[s], Form("%.0f < #it{p}_{T}^{calo} < %.0f GeV", sliceLo[s], sliceHi[s]), "lp");
    }
    leg->Draw();

    // where two calibrated collections should peak
    TLine zero; zero.SetLineStyle(7); zero.SetLineColor(kGray + 2);
    zero.DrawLine(0., 0., 0., yMax * 1.55);

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.15, 0.935, "PYTHIA pp 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4");
    la.SetTextSize(0.031);
    la.DrawLatex(0.15, 0.895, Form("#Delta#it{R} < %.1f matches, %s #it{p}_{T}", dRMatch, ptLabel));

    c->SaveAs(Form("%sdptRel_dists_ptSlices%s.pdf", outDir, ptTag));
    delete c;
  }

  printf("\n  figures in %s\n", outDir);
  f->Close();
}
