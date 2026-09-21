// Flavor assigned to a calo jet vs flavor assigned to its matched PF jet.
//
// Companion to plotCaloPFMatching_PYTHIA.C, which established that the calo ->
// PF geometric match is essentially perfect above 50 GeV (100% within dR < 0.4,
// median dR ~ 0.02, ~100% of pairs sharing one gen jet). So a disagreement in
// the map below is a disagreement between the two flavor DEFINITIONS on one and
// the same physical jet, not a matching failure.
//
// FLAVOR CONVENTIONS, copied from PYTHIA_scan.C:1139-1146 so the map describes
// what the analysis actually does.
//   PF    jetFlavorInt = (int) jtPartonFlavor, then promoted to 17 (bGS) when
//         |flavor| == 5 && bHadronNumber == 2.
//   calo  jetFlavorInt = caloJetPartonFlavor(refparton_flavor), i.e. the same
//         value except that -999 (no matched parton) is mapped to 0 -> xJets.
//         ak4CaloJetAnalyzer/t has no jtPartonFlavor, jtHadronFlavor or
//         bHadronNumber, so bHadronNumber is held at 0 and NO calo jet can ever
//         be labeled bGS. The bGS column is therefore empty by construction,
//         not by measurement, and every PF bGS jet is forced off the diagonal.
// Both sides go through getFlavorIdx() from src/scanning/JetFlavorIdx.h, the
// same mapping the templates use.
//
// Because of that structural asymmetry the agreement numbers are quoted twice:
// strict (8 categories) and merged (bGS folded into b, which is the comparison
// that means something).
//
// Figures, all in figures/jetMatching/:
//   flavorMap_<ptbin>.pdf        8x8 map, calo label on x, PF label on y, one
//                                per pT slice plus inclusive. Each calo COLUMN
//                                is normalized to 100%, so a column reads as
//                                "given a calo jet labeled F, what does PF call
//                                it". Cell text is that percentage.
//   flavorMapCounts_incl.pdf     the inclusive map as raw counts, so the
//                                column normalization cannot hide a column
//                                that is carrying almost no jets.
//   flavorAgreement_vs_pt.pdf    fraction of matched pairs with the same label
//                                vs pT (strict and merged), plus the b-jet
//                                efficiency and purity of the calo label.
//
// pT binning and slicing are on the CALO jet pT, since the calo spectrum is the
// thing being measured. Note the calo energy scale in this forest is 22% low
// (see plotCaloPFMatching_PYTHIA.C), so a calo pT slice does not hold the same
// jets as the PF slice of the same name; that is a scale problem, not a flavor
// problem, and it does not affect the per-pair comparison made here.
//
// Usage: root -l -b -q 'plotCaloPFFlavorMatch_PYTHIA.C'
// Run from: src/plots/jetPt/jetCollection/

#include <cmath>
#include "../../../../headers/plotting/plotStyle.h"
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

const char *forestFile = "/home/clayton/Downloads/HiForestAOD_PYTHIA.root";
const char *outDir     = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/jetMatching/";

const double caloPtFloor = 30.;
const double caloEtaMax  = 2.0;
const double pfPtFloor   = 10.;
const double dRMatch     = 0.2;   // justified by plotCaloPFMatching_PYTHIA.C

// slice 0 is inclusive; 1..5 are the pT slices
const int    nSlice = 6;
const double sliceLo[nSlice] = { 30.,  30.,  50.,  80., 120., 200.};
const double sliceHi[nSlice] = {400.,  50.,  80., 120., 200., 400.};
const char  *sliceTag[nSlice] = {"incl", "pt30to50", "pt50to80", "pt80to120",
                                 "pt120to200", "pt200to400"};

// the eight specific categories getFlavorIdx can return, in enum order
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

// PYTHIA_scan.C:259 -- calo jets with no matched parton land in xJets
static int caloJetPartonFlavor(int refpartonFlavor)
{
  return refpartonFlavor < -900 ? 0 : refpartonFlavor;
}

// fold bGS into b, for the comparison that is not rigged against calo
static int mergeBGS(int idx){ return idx == kBGSJets ? kBJets : idx; }

static void drawFlavorMap(TH2D *h, const char *note, bool asCounts, const char *outPath)
{
  gStyle->SetPalette(kViridis);
  gStyle->SetPaintTextFormat(asCounts ? ".0f" : ".1f");

  TCanvas *c = new TCanvas(Form("c_%s", h->GetName()), "", 800, 780);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.17);
  c->SetTopMargin(0.18);  c->SetBottomMargin(0.12);   // room for three caption lines
  if(asCounts) c->SetLogz();

  h->SetTitle(""); h->SetStats(0);
  h->GetXaxis()->SetTitle("calo jet flavor (refparton_flavor)");
  h->GetYaxis()->SetTitle("PF jet flavor (jtPartonFlavor)");
  h->GetXaxis()->SetTitleSize(0.045); h->GetXaxis()->SetLabelSize(0.058);
  h->GetYaxis()->SetTitleSize(0.045); h->GetYaxis()->SetLabelSize(0.058);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetZaxis()->SetTitle(asCounts ? "matched jet pairs" : "percent of calo column");
  h->GetZaxis()->SetTitleSize(0.040); h->GetZaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetTitleOffset(1.25);
  h->SetMarkerSize(1.25);
  h->Draw("COLZ TEXT");

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.037);
  la.DrawLatex(0.13, 0.955, "PYTHIA pp 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4");
  la.SetTextSize(0.030);
  la.DrawLatex(0.13, 0.917, note);
  la.DrawLatex(0.13, 0.883, "calo has no bHadronNumber: the b_{GS} column is empty by construction");

  c->SaveAs(outPath);
  delete c;
}

void plotCaloPFFlavorMatch_PYTHIA()
{
  initPlotStyle();
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(forestFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", forestFile); return; }

  // Get() takes the highest cycle, the complete 42000-entry write; the PF
  // directory also holds a truncated t;1 from an earlier pass.
  TTree *tC = (TTree*) f->Get("ak4CaloJetAnalyzer/t");
  TTree *tP = (TTree*) f->Get("ak4PFJetAnalyzer/t");
  if(!tC || !tP){ printf("ERROR: a jet tree is missing\n"); return; }
  if(tC->GetEntries() != tP->GetEntries()){
    printf("ERROR: trees are not the same length -- not the same events\n"); return;
  }
  if(!tC->GetBranch("refparton_flavor") || !tP->GetBranch("jtPartonFlavor") ||
     !tP->GetBranch("bHadronNumber")){
    printf("ERROR: a flavor branch is missing from this forest\n"); return;
  }
  if(tC->GetBranch("jtPartonFlavor"))
    printf("NOTE: this calo tree DOES have jtPartonFlavor, unlike the withGS forest.\n"
           "      The map still uses refparton_flavor, which is what the scan uses.\n");

  Int_t   nC = 0, nP = 0;
  Float_t ptC[MAXJET], etaC[MAXJET], phiC[MAXJET];
  Float_t ptP[MAXJET], etaP[MAXJET], phiP[MAXJET];
  Int_t   refFlavC[MAXJET];
  Float_t partFlavP[MAXJET];            // jtPartonFlavor is stored as a float
  Int_t   bHadNP[MAXJET];

  tC->SetBranchAddress("nref", &nC);
  tC->SetBranchAddress("jtpt", ptC);
  tC->SetBranchAddress("jteta", etaC);
  tC->SetBranchAddress("jtphi", phiC);
  tC->SetBranchAddress("refparton_flavor", refFlavC);

  tP->SetBranchAddress("nref", &nP);
  tP->SetBranchAddress("jtpt", ptP);
  tP->SetBranchAddress("jteta", etaP);
  tP->SetBranchAddress("jtphi", phiP);
  tP->SetBranchAddress("jtPartonFlavor", partFlavP);
  tP->SetBranchAddress("bHadronNumber", bHadNP);

  // --- histograms -----------------------------------------------------------
  TH2D *hMap[nSlice];
  for(int s = 0; s < nSlice; s++){
    hMap[s] = new TH2D(Form("hMap_%s", sliceTag[s]), "",
                       nFlav, 0.5, nFlav + 0.5, nFlav, 0.5, nFlav + 0.5);
    for(int i = 0; i < nFlav; i++){
      hMap[s]->GetXaxis()->SetBinLabel(i + 1, flavLabel[i]);
      hMap[s]->GetYaxis()->SetBinLabel(i + 1, flavLabel[i]);
    }
  }

  long nPair[nSlice] = {0}, nSameStrict[nSlice] = {0}, nSameMerged[nSlice] = {0};
  long nCaloB[nSlice] = {0}, nPFB[nSlice] = {0}, nBothB[nSlice] = {0};
  long nUnmatched = 0, nNoFlav = 0;

  const Long64_t nEv = tC->GetEntries();
  for(Long64_t ie = 0; ie < nEv; ie++){
    tC->GetEntry(ie);
    tP->GetEntry(ie);

    for(int i = 0; i < nC; i++){
      if(ptC[i] < caloPtFloor || ptC[i] >= sliceHi[0]) continue;
      if(fabs(etaC[i]) > caloEtaMax) continue;

      int best = -1; double bestDR = 1e9;
      for(int j = 0; j < nP; j++){
        if(ptP[j] < pfPtFloor) continue;
        double de = etaC[i] - etaP[j], dp = dPhi(phiC[i], phiP[j]);
        double dr = sqrt(de*de + dp*dp);
        if(dr < bestDR){ bestDR = dr; best = j; }
      }
      if(best < 0 || bestDR > dRMatch){ nUnmatched++; continue; }

      int idxC = getFlavorIdx(caloJetPartonFlavor(refFlavC[i]));

      int flavP = (int) partFlavP[best];
      if(fabs(flavP) == 5 && bHadNP[best] == 2) flavP = 17;   // bGS, as in the scan
      int idxP = getFlavorIdx(flavP);

      if(idxC < 0 || idxP < 0){ nNoFlav++; continue; }

      for(int s = 0; s < nSlice; s++){
        if(ptC[i] < sliceLo[s] || ptC[i] >= sliceHi[s]) continue;
        hMap[s]->Fill(idxC, idxP);
        nPair[s]++;
        if(idxC == idxP) nSameStrict[s]++;
        if(mergeBGS(idxC) == mergeBGS(idxP)) nSameMerged[s]++;
        bool bC = (mergeBGS(idxC) == kBJets), bP = (mergeBGS(idxP) == kBJets);
        if(bC) nCaloB[s]++;
        if(bP) nPFB[s]++;
        if(bC && bP) nBothB[s]++;
      }
    }
  }

  // --- printout -------------------------------------------------------------
  printf("\n  calo vs PF flavor, PYTHIA pp, %lld events\n", nEv);
  printf("  calo probe: pT > %.0f, |eta| < %.1f, matched to closest PF within dR < %.1f\n",
         caloPtFloor, caloEtaMax, dRMatch);
  printf("  %ld calo jets had no PF jet inside dR, %ld had an unmappable flavor code\n\n",
         nUnmatched, nNoFlav);

  printf("  %-12s %8s | %9s %9s | %7s %7s %7s | %8s %8s\n",
         "pT(calo)", "pairs", "same(8)", "same(bGS=b)", "caloB", "pfB", "bothB",
         "b eff", "b purity");
  for(int s = 0; s < nSlice; s++){
    if(!nPair[s]) continue;
    printf("  %-12s %8ld | %8.2f%% %10.2f%% | %7ld %7ld %7ld | %7.1f%% %8.1f%%\n",
           s == 0 ? "inclusive" : Form("%.0f-%.0f", sliceLo[s], sliceHi[s]),
           nPair[s], 100.*nSameStrict[s]/nPair[s], 100.*nSameMerged[s]/nPair[s],
           nCaloB[s], nPFB[s], nBothB[s],
           nPFB[s]   ? 100.*nBothB[s]/nPFB[s]   : 0.,
           nCaloB[s] ? 100.*nBothB[s]/nCaloB[s] : 0.);
  }
  printf("\n  'b eff'    = of PF jets labeled b (or bGS), the fraction calo also labels b\n");
  printf("  'b purity' = of calo jets labeled b, the fraction PF also labels b (or bGS)\n");

  // where do the calo x jets go? that is the refparton_flavor failure, in full
  printf("\n  PF label of the calo jets that refparton_flavor left unassigned (x), inclusive\n  ");
  double xTot = 0.;
  for(int iy = 1; iy <= nFlav; iy++) xTot += hMap[0]->GetBinContent(kXJets, iy);
  for(int iy = 1; iy <= nFlav; iy++){
    double v = hMap[0]->GetBinContent(kXJets, iy);
    if(v > 0.) printf(" %s %.1f%% ", flavLabel[iy-1], xTot > 0. ? 100.*v/xTot : 0.);
  }
  printf("\n  (%.0f jets, %.1f%% of all matched pairs)\n",
         xTot, nPair[0] ? 100.*xTot/nPair[0] : 0.);

  // the reverse projection for b: where the b jets the calo label loses end up.
  // This is the number that matters for a calo b-jet measurement.
  printf("\n  calo label of the PF b jets (b and bGS rows), inclusive\n  ");
  double bTot = 0.;
  for(int ix = 1; ix <= nFlav; ix++)
    bTot += hMap[0]->GetBinContent(ix, kBJets) + hMap[0]->GetBinContent(ix, kBGSJets);
  for(int ix = 1; ix <= nFlav; ix++){
    double v = hMap[0]->GetBinContent(ix, kBJets) + hMap[0]->GetBinContent(ix, kBGSJets);
    if(v > 0.) printf(" %s %.1f%% ", flavLabel[ix-1], bTot > 0. ? 100.*v/bTot : 0.);
  }
  printf("\n  (%.0f PF b jets matched to a calo jet)\n", bTot);

  // --- figures --------------------------------------------------------------
  // counts version first, before the column normalization overwrites the map
  TH2D *hCounts = (TH2D*) hMap[0]->Clone("hMapCounts_incl");
  drawFlavorMap(hCounts, Form("matched pairs, calo #it{p}_{T} > %.0f GeV, raw counts", caloPtFloor),
                true, Form("%sflavorMapCounts_incl.pdf", outDir));

  for(int s = 0; s < nSlice; s++){
    if(!nPair[s]) continue;
    // normalize each calo column to 100%
    for(int ix = 1; ix <= nFlav; ix++){
      double sum = 0.;
      for(int iy = 1; iy <= nFlav; iy++) sum += hMap[s]->GetBinContent(ix, iy);
      if(sum <= 0.) continue;
      for(int iy = 1; iy <= nFlav; iy++)
        hMap[s]->SetBinContent(ix, iy, 100. * hMap[s]->GetBinContent(ix, iy) / sum);
    }
    hMap[s]->SetMinimum(0.); hMap[s]->SetMaximum(100.);
    TString note = (s == 0)
      ? Form("calo #it{p}_{T} > %.0f GeV, each calo column normalized to 100%%", caloPtFloor)
      : Form("%.0f < calo #it{p}_{T} < %.0f GeV, each calo column normalized to 100%%",
             sliceLo[s], sliceHi[s]);
    drawFlavorMap(hMap[s], note, false, Form("%sflavorMap_%s.pdf", outDir, sliceTag[s]));
  }

  // agreement and b-tag correspondence vs pT
  {
    const int nB = nSlice - 1;
    double edges[nB + 1] = {30., 50., 80., 120., 200., 400.};
    TH1D *hStrict = new TH1D("hStrict", "", nB, edges);
    TH1D *hMerged = new TH1D("hMerged", "", nB, edges);
    TH1D *hBEff   = new TH1D("hBEff",   "", nB, edges);
    TH1D *hBPur   = new TH1D("hBPur",   "", nB, edges);
    for(int s = 1; s < nSlice; s++){
      int b = s;   // slice s covers bin b = s
      auto setBinom = [&](TH1D *h, long num, long den){
        if(!den) return;
        double p = (double)num/den;
        h->SetBinContent(b, p);
        h->SetBinError(b, sqrt(TMath::Max(1e-12, p*(1.-p)/den)));
      };
      setBinom(hStrict, nSameStrict[s], nPair[s]);
      setBinom(hMerged, nSameMerged[s], nPair[s]);
      setBinom(hBEff,   nBothB[s],      nPFB[s]);
      setBinom(hBPur,   nBothB[s],      nCaloB[s]);
    }

    TCanvas *c = new TCanvas("c_agree", "", 700, 800);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.12);
    c->SetTopMargin(0.13);  c->SetRightMargin(0.05);

    // hStrict is not drawn: bGS is rare enough that it lies exactly under
    // hMerged (57.57% vs 57.58% inclusive). Both are in the printout.
    styleH(hMerged, okabeHex[0], markFilledCircle);
    styleH(hBEff,   okabeHex[2], markFilledDiamond, 1.3);
    styleH(hBPur,   okabeHex[6], markOpenCircle);

    hMerged->SetTitle(""); hMerged->SetStats(0);
    // headroom for the legend, which would otherwise sit on the b efficiency
    hMerged->SetMinimum(0.); hMerged->SetMaximum(1.45);
    hMerged->GetXaxis()->SetTitle("calo jet #it{p}_{T} [GeV]");
    hMerged->GetXaxis()->SetTitleSize(0.045); hMerged->GetXaxis()->SetLabelSize(0.040);
    hMerged->GetYaxis()->SetTitle("fraction");
    hMerged->GetYaxis()->SetTitleSize(0.045); hMerged->GetYaxis()->SetLabelSize(0.040);
    hMerged->GetYaxis()->SetTitleOffset(1.45);
    hMerged->Draw("E1");
    hBEff  ->Draw("E1 same");
    hBPur  ->Draw("E1 same");

    // bGS is merged into b for the agreement curve, and the b curves are the
    // purity and the efficiency of the calo label; kept out of the legend by
    // request, documented in the header and printed in the table.
    TLegend *leg = makeLegend(0.52, 0.66, 0.94, 0.855, 0.031);
    leg->AddEntry(hMerged, "same flavor label", "lp");
    leg->AddEntry(hBPur,   "calo b #rightarrow PF b", "lp");
    leg->AddEntry(hBEff,   "PF b #rightarrow calo b", "lp");
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.15, 0.945, "PYTHIA pp 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4");
    la.SetTextSize(0.031);
    la.DrawLatex(0.15, 0.905, Form("matched pairs, #Delta#it{R} < %.1f", dRMatch));

    c->SaveAs(Form("%sflavorAgreement_vs_pt.pdf", outDir));
    delete c;
  }

  printf("\n  figures in %s\n", outDir);
  f->Close();
}
