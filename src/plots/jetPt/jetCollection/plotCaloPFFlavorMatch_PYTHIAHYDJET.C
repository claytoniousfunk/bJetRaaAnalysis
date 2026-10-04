// Flavor assigned to a calo jet vs flavor assigned to its matched PF jet,
// PYTHIA+HYDJET, per centrality class.
//
// PbPb counterpart of plotCaloPFFlavorMatch_PYTHIA.C (pp). That study showed the
// calo jet's own label (refparton_flavor) keeps only ~35-45% of the jets PF
// calls b, which is why the scans can relabel a calo jet with the flavor of the
// nearest PF jet (caloFlavorFromPFMatch). The pp geometry was clean: plain ak4
// on both sides, median dR ~ 0.02. Here the two collections are subtracted
// differently -- akPu4Calo (pileup subtraction) vs akCs4PF (constituent
// subtraction) -- and sit in a HYDJET background, so whether that relabeling is
// sound has to be checked separately, per centrality.
//
// WHAT IS COMPARED. Both labels are exactly what PYTHIAHYDJET_scan.C uses:
//   calo  refparton_flavorForB, -999 (no matched parton) -> 0 -> x
//         (recoJetPartonFlavor's fallback). No bHadronNumber, so never bGS.
//   PF    matchedPartonFlavor, promoted to 17 (bGS) when |flavor| == 5 and
//         bHadronNumber == 2 (PYTHIAHYDJET_scan.C:1529).
// NOTE the PbPb forests carry no jtPartonFlavor on akCs4PF. The PF-match path in
// PYTHIAHYDJET_scan.C and PYTHIAHYDJET_scan_response.C reads jtPartonFlavor, so
// on these forests it prints a warning and falls back to refparton_flavorForB:
// switching caloFlavorFromPFMatch on currently changes nothing. To make it work,
// the branch would have to become matchedPartonFlavor, the PbPb PF definition.
//
// MATCHING copies the scan's caloFlavorByPFMatch: nearest akCs4PF jet in dR, no
// PF pT floor (the forest stores jtpt > 20 anyway), accepted within
// caloPFMatchDR = 0.3. Matched pairs are also checked for sharing a gen jet
// (same refpt and refeta, both > 0), which separates "same physical jet,
// labels disagree" from "the match picked a different jet".
//
// Unweighted, like the pp study: a per-pair agreement in a narrow pT slice
// hardly depends on the pThat weight, and the counts stay readable as counts.
//
// Figures, all in figures/jetMatching/PYTHIAHYDJET/:
//   flavorMap_<class>.pdf          8x8 map per centrality class, calo pT > 30,
//                                  each calo column normalized to 100%
//   flavorMap_<class>_pt80to200.pdf  the same, for 80-200 GeV calo jets, where
//                                  the analysis lives and fakes no longer dominate
//   flavorAgreement_vs_pt_cent.pdf b efficiency (PF b -> calo b), b purity
//                                  (calo b -> PF b) and overall agreement vs
//                                  calo pT, one panel each, four classes overlaid
//   matchQuality_vs_pt_cent.pdf    fraction of calo jets with a PF jet inside
//                                  dR < 0.3, and fraction of those pairs sharing
//                                  one gen jet
//   dR_dists_cent.pdf              dR to the nearest PF jet, calo pT > 50, unit area
//
// Usage: root -l -b -q 'plotCaloPFFlavorMatch_PYTHIAHYDJET.C'
// Run from: src/plots/jetPt/jetCollection/

#include <cmath>
#include "../../../../headers/plotting/plotStyle.h"
#include "../../../../headers/plotting/coarseCent.h"
#include "../../../scanning/JetFlavorIdx.h"
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
#include <cstdio>

// two independent 9k-event PYTHIA+HYDJET DiJet forests, same tree layout
const char *forestFiles[] = {
  "/home/clayton/Downloads/HiForestAOD_PH.root",
  "/home/clayton/Downloads/HiForestAOD_100.root"
};
const int nFiles = sizeof(forestFiles)/sizeof(forestFiles[0]);
const char *outDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetMatching/PYTHIAHYDJET/";

const double caloPtFloor = 30.;
const double caloEtaMax  = 2.0;
const double dRMatch     = 0.3;   // caloPFMatchDR in PYTHIAHYDJET_scan*.C
const double dRPtFloor   = 50.;   // for the dR distributions

const int    nSlice = 6;   // 0 inclusive, 1..5 pT slices, as in the pp study
const double sliceLo[nSlice] = { 30.,  30.,  50.,  80., 120., 200.};
const double sliceHi[nSlice] = {400.,  50.,  80., 120., 200., 400.};

// classes from coarseCent.h, by hiBin (0.5% units)
const int hiBinLo[NCoarse] = {  0,  20,  60, 100};
const int hiBinHi[NCoarse] = { 20,  60, 100, 160};
const char *classHex[NCoarse]  = {"#0072B2", "#D55E00", "#009E73", "#CC79A7"};
const int   classMark[NCoarse] = {markFilledCircle, markFilledSquare, markFilledDiamond, markCross};

const int   nFlav = 8;
const char *flavLabel[nFlav] = {"d", "u", "s", "c", "b", "b_{GS}", "g", "x"};

static const int MAXJET = 500;

static double dPhi(double a, double b)
{
  double d = a - b;
  while(d >  TMath::Pi()) d -= TMath::TwoPi();
  while(d < -TMath::Pi()) d += TMath::TwoPi();
  return d;
}
static int caloJetPartonFlavor(int ref){ return ref < -900 ? 0 : ref; }
static int mergeBGS(int idx){ return idx == kBGSJets ? kBJets : idx; }
static int centClass(int hiBin)
{
  for(int k = 0; k < NCoarse; k++) if(hiBin >= hiBinLo[k] && hiBin < hiBinHi[k]) return k;
  return -1;
}

static void drawFlavorMap(TH2D *h, const char *note1, const char *note2, const char *outPath)
{
  gStyle->SetPalette(kViridis);
  gStyle->SetPaintTextFormat(".1f");

  TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", 800, 780);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.17);
  c->SetTopMargin(0.18);  c->SetBottomMargin(0.12);

  h->SetTitle(""); h->SetStats(0);
  h->GetXaxis()->SetTitle("calo jet flavor (refparton_flavorForB)");
  h->GetYaxis()->SetTitle("PF jet flavor (matchedPartonFlavor)");
  h->GetXaxis()->SetTitleSize(0.045); h->GetXaxis()->SetLabelSize(0.058);
  h->GetYaxis()->SetTitleSize(0.045); h->GetYaxis()->SetLabelSize(0.058);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetZaxis()->SetTitle("percent of calo column");
  h->GetZaxis()->SetTitleSize(0.040); h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleOffset(1.25);
  h->SetMarkerSize(1.25);
  h->SetMinimum(0.); h->SetMaximum(100.);
  h->Draw("COLZ TEXT");

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.037);
  la.DrawLatex(0.13, 0.955, "PYTHIA+HYDJET 5.02 TeV, akPu4Calo #rightarrow akCs4PF");
  la.SetTextSize(0.030);
  la.DrawLatex(0.13, 0.917, note1);
  la.DrawLatex(0.13, 0.883, note2);

  c->SaveAs(outPath);
  delete c;
}

static void normalizeColumns(TH2D *h)
{
  for(int ix = 1; ix <= h->GetNbinsX(); ix++){
    double sum = 0.;
    for(int iy = 1; iy <= h->GetNbinsY(); iy++) sum += h->GetBinContent(ix, iy);
    if(sum <= 0.) continue;
    for(int iy = 1; iy <= h->GetNbinsY(); iy++)
      h->SetBinContent(ix, iy, 100. * h->GetBinContent(ix, iy) / sum);
  }
}

static TH2D* newMap(const char *name)
{
  TH2D *h = new TH2D(name, "", nFlav, 0.5, nFlav + 0.5, nFlav, 0.5, nFlav + 0.5);
  for(int i = 0; i < nFlav; i++){
    h->GetXaxis()->SetBinLabel(i + 1, flavLabel[i]);
    h->GetYaxis()->SetBinLabel(i + 1, flavLabel[i]);
  }
  return h;
}

void plotCaloPFFlavorMatch_PYTHIAHYDJET()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  // [class][slice]
  TH2D *hMap[NCoarse], *hMapMid[NCoarse];
  TH1D *hDR[NCoarse];
  long nCalo[NCoarse][nSlice] = {{0}}, nPair[NCoarse][nSlice] = {{0}};
  long nSameGen[NCoarse][nSlice] = {{0}}, nSame[NCoarse][nSlice] = {{0}};
  long nCaloB[NCoarse][nSlice] = {{0}}, nPFB[NCoarse][nSlice] = {{0}}, nBothB[NCoarse][nSlice] = {{0}};
  long nCaloFake[NCoarse][nSlice] = {{0}};
  long nNoFlav = 0, nEvTot = 0, nEvClass[NCoarse] = {0};
  for(int k = 0; k < NCoarse; k++){
    hMap[k]    = newMap(Form("hMap_%s", coarseTag[k]));
    hMapMid[k] = newMap(Form("hMapMid_%s", coarseTag[k]));
    hDR[k]     = new TH1D(Form("hDR_%s", coarseTag[k]), "", 50, 0., 0.5);
  }

  for(int fi = 0; fi < nFiles; fi++){
    TFile *f = TFile::Open(forestFiles[fi]);
    if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", forestFiles[fi]); return; }
    TTree *tC = (TTree*) f->Get("akPu4CaloJetAnalyzer/t");
    TTree *tP = (TTree*) f->Get("akCs4PFJetAnalyzer/t");
    TTree *tE = (TTree*) f->Get("hiEvtAnalyzer/HiTree");
    if(!tC || !tP || !tE){ printf("ERROR: a tree is missing in %s\n", forestFiles[fi]); return; }
    if(tC->GetEntries() != tP->GetEntries() || tC->GetEntries() != tE->GetEntries()){
      printf("ERROR: trees are not the same length in %s\n", forestFiles[fi]); return;
    }
    if(!tC->GetBranch("refparton_flavorForB") || !tP->GetBranch("matchedPartonFlavor") ||
       !tP->GetBranch("bHadronNumber")){
      printf("ERROR: a flavor branch is missing from %s\n", forestFiles[fi]); return;
    }

    Int_t nC = 0, nP = 0, hiBin = -1;
    Float_t ptC[MAXJET], etaC[MAXJET], phiC[MAXJET], refptC[MAXJET], refetaC[MAXJET];
    Float_t ptP[MAXJET], etaP[MAXJET], phiP[MAXJET], refptP[MAXJET], refetaP[MAXJET];
    Int_t   refFlavC[MAXJET], flavP[MAXJET], bHadP[MAXJET];

    tE->SetBranchAddress("hiBin", &hiBin);
    tC->SetBranchAddress("nref", &nC);
    tC->SetBranchAddress("jtpt", ptC);
    tC->SetBranchAddress("jteta", etaC);
    tC->SetBranchAddress("jtphi", phiC);
    tC->SetBranchAddress("refpt", refptC);
    tC->SetBranchAddress("refeta", refetaC);
    tC->SetBranchAddress("refparton_flavorForB", refFlavC);
    tP->SetBranchAddress("nref", &nP);
    tP->SetBranchAddress("jtpt", ptP);
    tP->SetBranchAddress("jteta", etaP);
    tP->SetBranchAddress("jtphi", phiP);
    tP->SetBranchAddress("refpt", refptP);
    tP->SetBranchAddress("refeta", refetaP);
    tP->SetBranchAddress("matchedPartonFlavor", flavP);
    tP->SetBranchAddress("bHadronNumber", bHadP);

    const Long64_t nEv = tC->GetEntries();
    for(Long64_t ie = 0; ie < nEv; ie++){
      tE->GetEntry(ie); tC->GetEntry(ie); tP->GetEntry(ie);
      nEvTot++;
      int k = centClass(hiBin);
      if(k < 0) continue;
      nEvClass[k]++;

      for(int i = 0; i < nC; i++){
        if(ptC[i] < caloPtFloor || ptC[i] >= sliceHi[0]) continue;
        if(fabs(etaC[i]) > caloEtaMax) continue;

        int best = -1; double bestDR = 1e9;
        for(int j = 0; j < nP; j++){
          double de = etaC[i] - etaP[j], dp = dPhi(phiC[i], phiP[j]);
          double dr = sqrt(de*de + dp*dp);
          if(dr < bestDR){ bestDR = dr; best = j; }
        }
        if(ptC[i] >= dRPtFloor && best >= 0) hDR[k]->Fill(TMath::Min(bestDR, 0.4999));

        const bool matched = (best >= 0 && bestDR < dRMatch);
        const bool caloFake = (refptC[i] <= 0.);
        const bool sameGen  = matched && !caloFake && refptP[best] > 0. &&
                              fabs(refptC[i] - refptP[best]) < 1e-3 &&
                              fabs(refetaC[i] - refetaP[best]) < 1e-3;

        int idxC = -1, idxP = -1;
        if(matched){
          idxC = getFlavorIdx(caloJetPartonFlavor(refFlavC[i]));
          int fp = flavP[best];
          if(abs(fp) == 5 && bHadP[best] == 2) fp = 17;
          idxP = getFlavorIdx(fp);
          if(idxC < 0 || idxP < 0){ nNoFlav++; }
        }

        for(int s = 0; s < nSlice; s++){
          if(ptC[i] < sliceLo[s] || ptC[i] >= sliceHi[s]) continue;
          nCalo[k][s]++;
          if(caloFake) nCaloFake[k][s]++;
          if(!matched || idxC < 0 || idxP < 0) continue;
          nPair[k][s]++;
          if(sameGen) nSameGen[k][s]++;
          if(mergeBGS(idxC) == mergeBGS(idxP)) nSame[k][s]++;
          bool bC = mergeBGS(idxC) == kBJets, bP = mergeBGS(idxP) == kBJets;
          if(bC) nCaloB[k][s]++;
          if(bP) nPFB[k][s]++;
          if(bC && bP) nBothB[k][s]++;
          if(s == 0) hMap[k]->Fill(idxC, idxP);
        }
        if(matched && idxC >= 0 && idxP >= 0 && ptC[i] >= 80. && ptC[i] < 200.)
          hMapMid[k]->Fill(idxC, idxP);
      }
    }
    f->Close();
  }

  // --- printout -------------------------------------------------------------
  printf("\n  calo vs PF flavor, PYTHIA+HYDJET, %ld events (%ld / %ld / %ld / %ld by class)\n",
         nEvTot, nEvClass[0], nEvClass[1], nEvClass[2], nEvClass[3]);
  printf("  calo akPu4Calo pT > %.0f, |eta| < %.1f -> nearest akCs4PF within dR < %.1f\n",
         caloPtFloor, caloEtaMax, dRMatch);
  printf("  %ld matched pairs had an unmappable flavor code\n", nNoFlav);
  for(int k = 0; k < NCoarse; k++){
    printf("\n  %s\n", coarseLabel[k]);
    printf("  %-10s %7s %7s | %7s %8s | %7s | %6s %6s %6s | %7s %8s\n",
           "pT(calo)", "calo", "fake", "match", "sameGen", "agree", "caloB", "pfB", "bothB",
           "b eff", "b pur");
    for(int s = 0; s < nSlice; s++){
      if(!nCalo[k][s]) continue;
      long np = nPair[k][s];
      printf("  %-10s %7ld %6.1f%% | %6.1f%% %7.1f%% | %6.1f%% | %6ld %6ld %6ld | %6.1f%% %7.1f%%\n",
             s == 0 ? "incl" : Form("%.0f-%.0f", sliceLo[s], sliceHi[s]),
             nCalo[k][s], 100.*nCaloFake[k][s]/nCalo[k][s],
             100.*np/nCalo[k][s], np ? 100.*nSameGen[k][s]/np : 0.,
             np ? 100.*nSame[k][s]/np : 0.,
             nCaloB[k][s], nPFB[k][s], nBothB[k][s],
             nPFB[k][s]   ? 100.*nBothB[k][s]/nPFB[k][s]   : 0.,
             nCaloB[k][s] ? 100.*nBothB[k][s]/nCaloB[k][s] : 0.);
    }
  }
  printf("\n  fake    = calo jet with no gen match (refpt <= 0)\n");
  printf("  match   = calo jets with a PF jet inside dR; sameGen = of those, same gen jet\n");
  printf("  b eff   = of PF b (incl. bGS), fraction calo also labels b\n");
  printf("  b pur   = of calo b, fraction PF also labels b (or bGS)\n");

  // calo label of PF b jets, 80-200 GeV, per class
  printf("\n  calo label of the PF b jets (b and bGS rows), 80 < calo pT < 200\n");
  for(int k = 0; k < NCoarse; k++){
    double bTot = 0.;
    for(int ix = 1; ix <= nFlav; ix++)
      bTot += hMapMid[k]->GetBinContent(ix, kBJets) + hMapMid[k]->GetBinContent(ix, kBGSJets);
    printf("  %-7s", coarseLabel[k]);
    for(int ix = 1; ix <= nFlav; ix++){
      double v = hMapMid[k]->GetBinContent(ix, kBJets) + hMapMid[k]->GetBinContent(ix, kBGSJets);
      if(v > 0.) printf(" %s %.1f%% ", flavLabel[ix-1], bTot > 0. ? 100.*v/bTot : 0.);
    }
    printf(" (%.0f jets)\n", bTot);
  }

  // --- flavor maps ------------------------------------------------------------
  for(int k = 0; k < NCoarse; k++){
    normalizeColumns(hMap[k]);
    drawFlavorMap(hMap[k],
                  Form("%s, calo #it{p}_{T} > %.0f GeV, #Delta#it{R} < %.1f, each calo column = 100%%",
                       coarseLabel[k], caloPtFloor, dRMatch),
                  "calo has no bHadronNumber: the b_{GS} column is empty by construction",
                  Form("%sflavorMap_%s.pdf", outDir, coarseTag[k]));
    normalizeColumns(hMapMid[k]);
    drawFlavorMap(hMapMid[k],
                  Form("%s, 80 < calo #it{p}_{T} < 200 GeV, #Delta#it{R} < %.1f, each calo column = 100%%",
                       coarseLabel[k], dRMatch),
                  "calo has no bHadronNumber: the b_{GS} column is empty by construction",
                  Form("%sflavorMap_%s_pt80to200.pdf", outDir, coarseTag[k]));
  }

  // --- vs pT curves -----------------------------------------------------------
  const int nB = nSlice - 1;
  double edges[nB + 1] = {30., 50., 80., 120., 200., 400.};
  auto binomHist = [&](const char *name, int k, long (*num)[nSlice], long (*den)[nSlice]){
    TH1D *h = new TH1D(Form("%s_%s", name, coarseTag[k]), "", nB, edges);
    for(int s = 1; s < nSlice; s++){
      if(!den[k][s]) continue;
      double p = (double)num[k][s]/den[k][s];
      h->SetBinContent(s, p);
      h->SetBinError(s, sqrt(TMath::Max(1e-12, p*(1.-p)/den[k][s])));
    }
    styleH(h, classHex[k], classMark[k], classMark[k] == markFilledDiamond ? 1.4 : 1.1);
    return h;
  };

  auto drawPanels = [&](int nPanel, const char **ytitle, long (*num[])[nSlice], long (*den[])[nSlice],
                        const char *outName, const char *sub){
    TCanvas *c = new TCanvas(Form("c_%s", outName), "", 520*nPanel, 560);
    c->Divide(nPanel, 1, 0.002, 0.);
    for(int p = 0; p < nPanel; p++){
      c->cd(p + 1);
      gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.04);
      gPad->SetTopMargin(0.17);  gPad->SetBottomMargin(0.13);
      TLegend *leg = makeLegend(0.22, 0.20, 0.62, 0.42, 0.040);
      for(int k = 0; k < NCoarse; k++){
        TH1D *h = binomHist(Form("%s_p%d", outName, p), k, num[p], den[p]);
        if(k == 0){
          h->SetMinimum(0.); h->SetMaximum(1.05);
          h->GetXaxis()->SetTitle("calo jet #it{p}_{T} [GeV]");
          h->GetYaxis()->SetTitle(ytitle[p]);
          h->GetXaxis()->SetTitleSize(0.050); h->GetXaxis()->SetLabelSize(0.045);
          h->GetYaxis()->SetTitleSize(0.050); h->GetYaxis()->SetLabelSize(0.045);
          h->GetYaxis()->SetTitleOffset(1.5);
          h->Draw("E1");
        } else h->Draw("E1 same");
        leg->AddEntry(h, coarseLabel[k], "lp");
      }
      if(p == 0) leg->Draw();
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.040);
      la.DrawLatex(0.17, 0.93, "PYTHIA+HYDJET, akPu4Calo #rightarrow akCs4PF");
      la.SetTextSize(0.037);
      la.DrawLatex(0.17, 0.875, sub);
    }
    c->SaveAs(Form("%s%s.pdf", outDir, outName));
    delete c;
  };

  {
    const char *yt[3] = {"PF b #rightarrow calo b (b efficiency)",
                         "calo b #rightarrow PF b (b purity)",
                         "same flavor label (b_{GS} = b)"};
    long (*num[3])[nSlice] = {nBothB, nBothB, nSame};
    long (*den[3])[nSlice] = {nPFB,   nCaloB, nPair};
    drawPanels(3, yt, num, den, "flavorAgreement_vs_pt_cent",
               Form("matched pairs, #Delta#it{R} < %.1f", dRMatch));
  }
  {
    const char *yt[2] = {Form("calo jets with a PF jet in #Delta#it{R} < %.1f", dRMatch),
                         "matched pairs sharing one gen jet"};
    long (*num[2])[nSlice] = {nPair, nSameGen};
    long (*den[2])[nSlice] = {nCalo, nPair};
    drawPanels(2, yt, num, den, "matchQuality_vs_pt_cent", "nearest PF jet, no PF #it{p}_{T} floor");
  }

  // --- dR distributions -------------------------------------------------------
  {
    TCanvas *c = new TCanvas("c_dR", "", 700, 640);
    c->SetLeftMargin(0.15); c->SetTopMargin(0.13); c->SetBottomMargin(0.12);
    c->SetLogy();
    TLegend *leg = makeLegend(0.60, 0.62, 0.92, 0.84, 0.036);
    for(int k = 0; k < NCoarse; k++){
      if(hDR[k]->Integral() > 0) hDR[k]->Scale(1./hDR[k]->Integral());
      styleH(hDR[k], classHex[k], classMark[k], classMark[k] == markFilledDiamond ? 1.3 : 1.0);
      if(k == 0){
        hDR[k]->SetMinimum(1e-5); hDR[k]->SetMaximum(2.);
        hDR[k]->GetXaxis()->SetTitle("#Delta#it{R}(calo, nearest PF)  [last bin = overflow]");
        hDR[k]->GetYaxis()->SetTitle("fraction of calo jets");
        hDR[k]->GetXaxis()->SetTitleSize(0.045); hDR[k]->GetYaxis()->SetTitleSize(0.045);
        hDR[k]->Draw("E1");
      } else hDR[k]->Draw("E1 same");
      leg->AddEntry(hDR[k], coarseLabel[k], "lp");
    }
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.15, 0.945, "PYTHIA+HYDJET 5.02 TeV, akPu4Calo #rightarrow akCs4PF");
    la.SetTextSize(0.031);
    la.DrawLatex(0.15, 0.905, Form("calo #it{p}_{T} > %.0f GeV, |#eta| < %.1f, unit area", dRPtFloor, caloEtaMax));
    c->SaveAs(Form("%sdR_dists_cent.pdf", outDir));
    delete c;
  }

  printf("\n  figures in %s\n", outDir);
}
