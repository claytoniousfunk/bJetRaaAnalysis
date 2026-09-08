// X-projection of h_muonDR_inclusiveClosestJet_C* (dR vs. closest-jet pT, all
// muons -- not restricted to the muon's own tagged jet) in three jet pT
// windows, merged into the standard 4 coarse centrality classes.
//
// One canvas per coarse class, three overlaid dR shapes (one per pT window),
// each unit-normalized so the curves compare SHAPE, not rate -- absolute
// yields already differ by two orders of magnitude across the pT windows and
// would dominate an un-normalized overlay.
//
// STATISTICS ARE THIN outside 0-10%, and the legend reports raw entries per
// curve so that is visible rather than hidden: merged-class entry counts run
// from 364 (0-10%, 50-60 GeV) down to 4 (50-80%, 60-80 GeV). Treat the
// peripheral panels as indicative shapes, not measurements.
//
// dR rebinned into 5 uniform 0.1-wide bins -- wider than
// plotFastJetMuonDR_inJetVsInclusive.cc's convention, since the merged
// coarse-class statistics here (as low as 4 entries in a pT window) do not
// support the finer low-dR binning used there.
//
// Usage: root -l -b -q 'plotMuonDR_inclusiveClosestJet_coarseCent.C'
// Run from: src/plots/muJetDr/

#include "../../../headers/functions/divideByBinwidth.h"

const char *inFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-9-1_ultraFineCentBins.root";

const char *histBase = "h_muonDR_inclusiveClosestJet";

const char *outDir = "../../../figures/muJetDr/";

// coarse class -> [first, last] ultra-fine slice index, standard convention
const int NClass = 4;
const int sliceLo[NClass] = { 1,  3,  7, 11};
const int sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const int    NPt = 3;
const double ptLo[NPt] = {50, 60,  80};
const double ptHi[NPt] = {60, 80, 120};
const char  *ptHex[NPt]   = {"#0072B2", "#D55E00", "#009E73"};   // Okabe-Ito
const int    ptMarker[NPt] = {20, 21, 33};                        // circle, square, diamond -- no triangles

const int    NEdge_dR = 6;
double       edge_dR[NEdge_dR] = {0, 0.1, 0.2, 0.3, 0.4, 0.5};

const double plotDRmax = 0.5;

void plotMuonDR_inclusiveClosestJet_coarseCent(){
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(inFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", inFile); return; }

  for(int ci = 0; ci < NClass; ci++){

    // sum the fine slices making up this coarse class
    TH2D *hSum = nullptr;
    for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
      TH2D *h = nullptr;
      f->GetObject(Form("%s_C%d", histBase, si), h);
      if(!h){ printf("WARNING: %s_C%d missing\n", histBase, si); continue; }
      if(!hSum){ hSum = (TH2D*) h->Clone(Form("hSum_%d", ci)); hSum->SetDirectory(nullptr); }
      else hSum->Add(h);
    }
    if(!hSum){ printf("WARNING: no histograms found for class %s, skipping\n", classLabel[ci]); continue; }

    TH1D *proj[NPt];
    int   nEntries[NPt];

    for(int p = 0; p < NPt; p++){
      int b1 = hSum->GetYaxis()->FindBin(ptLo[p] + 1e-6);
      int b2 = hSum->GetYaxis()->FindBin(ptHi[p] - 1e-6);
      TH1D *raw = hSum->ProjectionX(Form("proj_%d_%d", ci, p), b1, b2);
      nEntries[p] = (int) raw->GetEntries();

      TH1D *h = (TH1D*) raw->Rebin(NEdge_dR-1, Form("h_dR_%d_%d", ci, p), edge_dR);
      h->SetDirectory(nullptr);
      if(h->Integral() > 0.) h->Scale(1./h->Integral());
      divideByBinwidth(h);
      proj[p] = h;
    }

    printf("%-8s", classLabel[ci]);
    for(int p = 0; p < NPt; p++) printf("  %.0f-%.0f GeV: %d entries", ptLo[p], ptHi[p], nEntries[p]);
    printf("\n");

    TCanvas *c = new TCanvas(Form("c_%d", ci), "", 700, 700);
    TPad *pad = new TPad("pad", "", 0, 0, 1, 1);
    pad->SetLeftMargin(0.14);
    pad->SetBottomMargin(0.12);
    pad->Draw(); pad->cd();

    double ymax = 0.;
    for(int p = 0; p < NPt; p++) if(proj[p]->GetMaximum() > ymax) ymax = proj[p]->GetMaximum();

    for(int p = 0; p < NPt; p++){
      int col = TColor::GetColor(ptHex[p]);
      proj[p]->SetLineColor(col);
      proj[p]->SetMarkerColor(col);
      proj[p]->SetMarkerStyle(ptMarker[p]);
      proj[p]->SetMarkerSize(1.0);
      proj[p]->SetLineWidth(2);
      proj[p]->SetTitle("");
      proj[p]->GetXaxis()->SetTitle("#it{#Delta}#it{R}(#it{#mu},jet)");
      proj[p]->GetYaxis()->SetTitle("1/#it{N} d#it{N}/d(#it{#Delta}#it{R})");
      proj[p]->GetXaxis()->SetTitleSize(0.045);
      proj[p]->GetYaxis()->SetTitleSize(0.045);
      proj[p]->GetXaxis()->SetLabelSize(0.04);
      proj[p]->GetYaxis()->SetLabelSize(0.04);
      proj[p]->GetXaxis()->SetRangeUser(0., plotDRmax);
      proj[p]->SetMaximum(ymax * 1.6);
      proj[p]->SetMinimum(0.);
      proj[p]->Draw(p == 0 ? "ep" : "ep same");
    }

    TLegend *leg = new TLegend(0.40, 0.60, 0.88, 0.80);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
    for(int p = 0; p < NPt; p++)
      leg->AddEntry(proj[p], Form("%.0f < p_{T}^{jet} < %.0f GeV  (n=%d)", ptLo[p], ptHi[p], nEntries[p]), "lp");
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.032);
    la.DrawLatex(0.17, 0.90, Form("PbPb SingleMuon (5.02 TeV), %s", classLabel[ci]));
    la.DrawLatex(0.17, 0.855, "p_{T}^{#mu} > 15 GeV, |#eta^{#mu}| < 2, |#eta^{jet}| < 1.6");

    TString out = TString(outDir) + Form("muonDR_inclusiveClosestJet_coarseCent_%s.pdf",
                    TString(classLabel[ci]).ReplaceAll("%","pct").ReplaceAll("-","to").Data());
    c->SaveAs(out);
    printf("Saved %s\n\n", out.Data());
  }
}
