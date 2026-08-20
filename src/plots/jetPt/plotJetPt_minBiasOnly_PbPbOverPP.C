// PbPb/pp jet pT shape comparison -- data (MinBias, no stitching, no
// unfolding, no per-event or per-Z normalization) vs the equivalent MC
// (PYTHIA+HYDJET/PYTHIA), one panel per coarse PbPb centrality class.
//
// Each curve (data PbPb per class, data pp, MC PbPb+HYDJET per class, MC pp)
// is independently scaled so its own integral above ptNormMin matches --
// every curve is forced to 1 up there by construction. This sidesteps
// needing a common, trustworthy absolute normalization across four different
// samples (two data, two MC) and asks only "do the shapes agree," which a
// high-pT-normalized ratio answers directly: flat at 1 above ptNormMin,
// deviations below it are shape differences, not a residual normalization
// mismatch.
//
// MC inclusive reco jet pT: matched-only (projected off the 2D
// matchedRecoJetPt-vs-genJetPt response histogram onto its reco-pT axis).
//
// An earlier version of this plot also added each sample's "unmatched" reco
// jets (no gen-level partner) to better match what a data spectrum -- which
// has no truth information to filter on -- actually is. That version used
// PYTHIAHYDJET_response_..._2026-7-28_{even,odd}Events.root for the PbPb+
// HYDJET side, which turned out to have doPThatWeight=false baked in at scan
// time (headers/config/config_PYTHIAHYDJET.h) -- a genuinely unweighted
// sample masquerading as a physical spectrum, which produced a nonsensical
// ratio (confirmed by tracing the w_pthat assignment in
// PYTHIAHYDJET_scan_response.C). The replacement file below is confirmed
// pThat-weighted, but is older and was not built with an "unmatched" reco-
// jet histogram at all -- so this version uses matched-only jets for BOTH
// MC curves. Using unmatched-inclusive for pp but matched-only for PbPb+
// HYDJET would reintroduce exactly the kind of construction asymmetry that
// caused the previous problem, so pp was switched to matched-only too rather
// than left inconsistent with PbPb+HYDJET.
//
// MC sample provenance (as of 2026-08-18):
//   pp   : PYTHIA_DiJet_response_..._2026-8-3.root (bJetRaaAnalysis) --
//          combined file; checked against its own evenEvents/oddEvents split
//          and confirmed to be their exact sum, so used directly.
//   PbPb : bJetMuonTaggingAnalysis/.../PYTHIAHYDJET_response_DiJet_pThat-15_
//          ..._doPThatCorrelationFilter_2026-3-23.root -- read-only source
//          repo (bJetMuonTaggingAnalysis is never modified from here); single
//          file, no even/odd split, no unmatched-jet histogram. Axis is
//          20-500 GeV / 96 bins (still 5 GeV/bin, just starting at 20 instead
//          of 0 like the other MC file) -- confirmed compatible with ptEdge
//          below, which already starts at 20.
// Both files carry genuine per-event cross-section weights (bin content !=
// raw entry count in both); no additional reweighting applied here.
//
// Centrality indexing for the MC files: C0 = combined (unused here), C1-C4 =
// 0-10/10-30/30-50/50-80%, same order as the data side's sliceLo/sliceHi/
// classLabel convention.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_PbPbOverPP.C'
// Run from: src/plots/jetPt/

const char *pbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-8-13_ultraFineCentBins.root";
const char *ppFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root";
const char *ppHistName = "h_inclRecoJetPt";

const char *mcPPFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIA/PYTHIA_DiJet_response_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-8-3.root";
const char *mcPbPbFile =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/response/PYTHIAHYDJET_response_DiJet_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_doPThatCorrelationFilter_2026-3-23.root";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_PbPbOverPP.pdf";

// PbPb coarse class -> [first, last] ultra-fine slice index (data side) /
// _C{1-4} suffix (MC side), same convention as makeFakeJetFile.C:
//   C1 = 0-10%  : slices  1- 2      C3 = 30-50% : slices  7-10
//   C2 = 10-30% : slices  3- 6      C4 = 50-80% : slices 11-16
const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

// Both curves in a pair are scaled so their integral above this threshold matches.
const double ptNormMin = 200.;

const double plotPtMin = 30., plotPtMax = 500.;

// Variable-width binning: 5 GeV steps from 20 GeV (the analysis jetPtCut, where
// both spectra turn on) up to 100 GeV, then the same coarse edges calculateRAA.C
// uses above 100 GeV, where raw statistics run thin. Reused for MC too so both
// sides bin identically.
const int    NEdge = 29;
double       ptEdge[NEdge] = {
  20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,
  110,120,130,140,150,160,180,200,240,280,350,500
};

const char *dataHex = "#0072B2", *mcHex = "#D55E00";
const int   dataMark = 20, mcMark = 21;

static TH1D* buildShape(TH1D *raw, const char *name)
{
  TH1D *h = (TH1D*) raw->Rebin(NEdge-1, name, ptEdge);
  h->SetDirectory(nullptr);
  int b0 = h->FindBin(ptNormMin);
  double norm = h->Integral(b0, h->GetNbinsX()+1);
  if(norm <= 0.){ printf("ERROR: %s has zero integral above %.0f GeV\n", name, ptNormMin); return nullptr; }
  h->Scale(1./norm);
  h->Scale(1., "width");  // genuine dN/dpT; cancels in any ratio of same-binned curves
  return h;
}

// Matched-only reco-level jet pT: projected off the 2D recoPt-vs-genPt
// response onto its reco-pT axis. See the file header for why this doesn't
// also add each sample's "unmatched" component.
static TH1D* buildMatchedRecoPt(TFile *f, const char *suffix, const char *name)
{
  TH2D *matched2D = nullptr;
  f->GetObject(Form("h_matchedRecoJetPt_genJetPt_allJets%s", suffix), matched2D);
  if(!matched2D){ printf("ERROR: missing matched2D hist for suffix '%s'\n", suffix); return nullptr; }

  TH1D *matched1D = (TH1D*) matched2D->ProjectionX(name);
  matched1D->SetDirectory(nullptr);
  return matched1D;
}

void plotJetPt_minBiasOnly_PbPbOverPP()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fP = TFile::Open(pbpbFile);
  TFile *fp = TFile::Open(ppFile);
  if(!fP || fP->IsZombie()){ printf("ERROR: cannot open %s\n", pbpbFile); return; }
  if(!fp || fp->IsZombie()){ printf("ERROR: cannot open %s\n", ppFile);   return; }

  TFile *fMCpp   = TFile::Open(mcPPFile);
  TFile *fMCpbpb = TFile::Open(mcPbPbFile);
  if(!fMCpp || fMCpp->IsZombie()){ printf("ERROR: cannot open %s\n", mcPPFile); return; }
  if(!fMCpbpb || fMCpbpb->IsZombie()){ printf("ERROR: cannot open %s\n", mcPbPbFile); return; }

  // --- data pp reference ---
  TH1D *hpp0 = nullptr;
  fp->GetObject(ppHistName, hpp0);
  if(!hpp0){ printf("ERROR: %s not found in %s\n", ppHistName, ppFile); return; }
  TH1D *hpp = buildShape(hpp0, "hpp");
  if(!hpp) return;

  // --- MC pp reference ---
  TH1D *mcPPraw = buildMatchedRecoPt(fMCpp, "", "mcPPraw");
  if(!mcPPraw) return;
  TH1D *hMCpp = buildShape(mcPPraw, "hMCpp");
  if(!hMCpp) return;

  TH1D *hPbPb[NClass], *hRatio[NClass];
  TH1D *hMCPbPb[NClass], *hMCRatio[NClass];
  for(int ci = 0; ci < NClass; ci++){
    // --- data PbPb, this class ---
    TH1D *sum = nullptr;
    bool ok = true;
    for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
      TH1D *h = nullptr;
      fP->GetObject(Form("h_inclRecoJetPt_C%d", si), h);
      if(!h){ printf("ERROR: h_inclRecoJetPt_C%d not found\n", si); ok = false; break; }
      if(!sum){ sum = (TH1D*) h->Clone(Form("rawSum_%d", ci)); sum->SetDirectory(nullptr); }
      else sum->Add(h);
    }
    if(!ok) return;
    hPbPb[ci] = buildShape(sum, Form("hPbPb_C%d", ci));
    if(!hPbPb[ci]) return;
    hRatio[ci] = (TH1D*) hPbPb[ci]->Clone(Form("hRatio_C%d", ci));
    hRatio[ci]->Divide(hpp);

    // --- MC PbPb+HYDJET, this class ---
    TString suffix = Form("_C%d", ci+1);
    TH1D *mcRaw = buildMatchedRecoPt(fMCpbpb, suffix.Data(), Form("mcPbPb_C%d", ci));
    if(!mcRaw) return;
    hMCPbPb[ci] = buildShape(mcRaw, Form("hMCPbPb_C%d", ci));
    if(!hMCPbPb[ci]) return;
    hMCRatio[ci] = (TH1D*) hMCPbPb[ci]->Clone(Form("hMCRatio_C%d", ci));
    hMCRatio[ci]->Divide(hMCpp);
  }

  printf("\n%-12s", "pT [GeV]");
  for(int ci = 0; ci < NClass; ci++) printf(" %14s %14s", Form("%s data", classLabel[ci]), Form("%s MC", classLabel[ci]));
  printf("\n");
  for(int b = 1; b <= hpp->GetNbinsX(); b++){
    double lo = hpp->GetXaxis()->GetBinLowEdge(b);
    if(lo < plotPtMin || lo >= plotPtMax) continue;
    printf("%4.0f-%-6.0f", lo, hpp->GetXaxis()->GetBinUpEdge(b));
    for(int ci = 0; ci < NClass; ci++)
      printf(" %8.3f+-%.3f %8.3f+-%.3f",
             hRatio[ci]->GetBinContent(b), hRatio[ci]->GetBinError(b),
             hMCRatio[ci]->GetBinContent(b), hMCRatio[ci]->GetBinError(b));
    printf("\n");
  }

  int colData = TColor::GetColor(dataHex), colMC = TColor::GetColor(mcHex);

  double ratioMin[NClass], ratioMax[NClass];
  for(int ci = 0; ci < NClass; ci++){
    ratioMin[ci] = 1e300; ratioMax[ci] = -1e300;
    for(int b = 1; b <= hRatio[ci]->GetNbinsX(); b++){
      double lo = hRatio[ci]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      for(double v : {hRatio[ci]->GetBinContent(b), hMCRatio[ci]->GetBinContent(b)}){
        if(v == 0.) continue;
        if(v < ratioMin[ci]) ratioMin[ci] = v;
        if(v > ratioMax[ci]) ratioMax[ci] = v;
      }
    }
    double pad = (ratioMax[ci] - ratioMin[ci]) * 0.15;
    ratioMin[ci] = TMath::Max(0., ratioMin[ci] - pad);
    ratioMax[ci] += pad;
  }

  TCanvas *c = new TCanvas("cMinBiasRatio_dataVsMC", "", 1200, 1000);
  c->Divide(2, 2, 0.001, 0.001);

  for(int ci = 0; ci < NClass; ci++){
    c->cd(ci + 1);
    gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.09);  gPad->SetBottomMargin(0.14);

    hRatio[ci]->SetLineColor(colData); hRatio[ci]->SetMarkerColor(colData);
    hRatio[ci]->SetMarkerStyle(dataMark); hRatio[ci]->SetMarkerSize(1.0);
    hRatio[ci]->SetLineWidth(2);
    hRatio[ci]->SetTitle("");
    hRatio[ci]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
    hRatio[ci]->GetYaxis()->SetTitle("PbPb / pp");
    hRatio[ci]->GetXaxis()->SetTitleSize(0.048); hRatio[ci]->GetXaxis()->SetLabelSize(0.042);
    hRatio[ci]->GetYaxis()->SetTitleSize(0.048); hRatio[ci]->GetYaxis()->SetLabelSize(0.042);
    hRatio[ci]->GetYaxis()->SetTitleOffset(1.45);
    hRatio[ci]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hRatio[ci]->SetMinimum(ratioMin[ci]); hRatio[ci]->SetMaximum(ratioMax[ci]);
    hRatio[ci]->Draw("ep");

    hMCRatio[ci]->SetLineColor(colMC); hMCRatio[ci]->SetMarkerColor(colMC);
    hMCRatio[ci]->SetMarkerStyle(mcMark); hMCRatio[ci]->SetMarkerSize(1.0);
    hMCRatio[ci]->SetLineWidth(2);
    hMCRatio[ci]->Draw("ep same");

    TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

    TLatex clat; clat.SetNDC(); clat.SetTextSize(0.052); clat.SetTextFont(62);
    clat.SetTextAlign(13);
    clat.DrawLatex(0.16 + 0.03, 1. - 0.09 - 0.03, Form("PbPb %s", classLabel[ci]));

    if(ci == 0){
      TLegend *leg = new TLegend(0.38, 0.68, 0.94, 0.85);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
      leg->AddEntry(hRatio[ci],   "data (PbPb/pp)",            "lp");
      leg->AddEntry(hMCRatio[ci], "MC (PYTHIA+HYDJET/PYTHIA)", "lp");
      leg->Draw();
    }
  }

  c->cd(0);
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.020);
  lat.DrawLatex(0.06, 0.985, "MinBias only, shape only -- each curve normalized above 200 GeV");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
