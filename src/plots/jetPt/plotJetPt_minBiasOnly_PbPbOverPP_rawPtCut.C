// Same as plotJetPt_minBiasOnly_PbPbOverPP.C, but the PbPb side uses
// h_inclRecoJetPt_rawPtCut_C* -- jets whose UNCORRECTED pT clears 40 GeV
// (see PbPb_scan.C's rawJetPtCut), x axis still the JEC-corrected pT. Intended
// to suppress combinatorial jets, which are built from the UE pedestal and so
// sit at low raw pT even where JEC pushes them into the signal region.
//
// The pp side is unchanged (plain h_inclRecoJetPt): pp_scan.C has no raw-pT-cut
// variant, and pp MinBias has no meaningful combinatorial background for such a
// cut to suppress in the first place -- there is no PbPb-style underlying event
// to build a fake jet's pT from. So this compares filtered PbPb against
// unfiltered pp, deliberately: the question is whether the cut brings PbPb's
// shape into line with the (already clean) pp shape, not a like-for-like
// selection on both sides.
//
// Same normalization/binning scheme and the same caveats apply -- see the
// unfiltered version for the full explanation.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_PbPbOverPP_rawPtCut.C'
// Run from: src/plots/jetPt/

const char *pbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-8-13_ultraFineCentBins.root";
const char *ppFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root";
const char *ppHistName   = "h_inclRecoJetPt";
const char *pbpbHistBase = "h_inclRecoJetPt_rawPtCut";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_PbPbOverPP_rawPtCut.pdf";

// PbPb coarse class -> [first, last] ultra-fine slice index (10 hiBin units
// each), same convention as makeFakeJetFile.C:
//   C1 = 0-10%  : slices  1- 2      C3 = 30-50% : slices  7-10
//   C2 = 10-30% : slices  3- 6      C4 = 50-80% : slices 11-16
const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

// Both curves are scaled so their integral above this threshold matches.
const double ptNormMin = 200.;

// Plotted x-range. Kept at 30 GeV to match the unfiltered version for a direct
// visual comparison, even though the raw-pT cut is expected to tame the
// near-threshold peak that motivated excluding 20-30 GeV there.
const double plotPtMin = 30.;

// Fixed, not autoscaled -- see the printed NOTE check below if this needs
// revisiting after a rerun.
const double ratioMin = 0., ratioMax = 10.;

// Variable-width binning: 5 GeV steps from 20 GeV (the analysis jetPtCut, where
// both spectra turn on) up to 100 GeV, then the same coarse edges calculateRAA.C
// uses above 100 GeV, where raw statistics run thin.
const int    NEdge = 29;
double       ptEdge[NEdge] = {
  20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,
  110,120,130,140,150,160,180,200,240,280,350,500
};

// Okabe-Ito, centre-symmetric markers only. pp is the black reference; the four
// PbPb classes each get a distinct colour and shape so nothing depends on
// colour alone.
const char *classHex[NClass]  = { "#D55E00", "#0072B2", "#009E73", "#CC79A7" };
const int   classMark[NClass] = { 20, 21, 33, 34 };
const double classSize[NClass]= { 1.0, 1.0, 1.2, 1.1 };

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

void plotJetPt_minBiasOnly_PbPbOverPP_rawPtCut()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fP = TFile::Open(pbpbFile);
  TFile *fp = TFile::Open(ppFile);
  if(!fP || fP->IsZombie()){ printf("ERROR: cannot open %s\n", pbpbFile); return; }
  if(!fp || fp->IsZombie()){ printf("ERROR: cannot open %s\n", ppFile);   return; }

  TH1D *hpp0 = nullptr;
  fp->GetObject(ppHistName, hpp0);
  if(!hpp0){ printf("ERROR: %s not found in %s\n", ppHistName, ppFile); return; }
  TH1D *hpp = buildShape(hpp0, "hpp");
  if(!hpp) return;

  TH1D *hPbPb[NClass], *hRatio[NClass];
  for(int ci = 0; ci < NClass; ci++){
    TH1D *sum = nullptr;
    bool ok = true;
    for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
      TH1D *h = nullptr;
      fP->GetObject(Form("%s_C%d", pbpbHistBase, si), h);
      if(!h){ printf("ERROR: %s_C%d not found\n", pbpbHistBase, si); ok = false; break; }
      if(!sum){ sum = (TH1D*) h->Clone(Form("rawSum_%d", ci)); sum->SetDirectory(nullptr); }
      else sum->Add(h);
    }
    if(!ok) return;
    hPbPb[ci] = buildShape(sum, Form("hPbPb_C%d", ci));
    if(!hPbPb[ci]) return;

    hRatio[ci] = (TH1D*) hPbPb[ci]->Clone(Form("hRatio_C%d", ci));
    hRatio[ci]->Divide(hpp);
  }

  printf("\n%-12s", "pT [GeV]");
  for(int ci = 0; ci < NClass; ci++) printf(" %14s", classLabel[ci]);
  printf("\n");
  for(int b = 1; b <= hpp->GetNbinsX(); b++){
    printf("%4.0f-%-6.0f", hpp->GetXaxis()->GetBinLowEdge(b), hpp->GetXaxis()->GetBinUpEdge(b));
    for(int ci = 0; ci < NClass; ci++)
      printf(" %8.3f+-%.3f", hRatio[ci]->GetBinContent(b), hRatio[ci]->GetBinError(b));
    printf("\n");
  }

  int colpp = kBlack;
  int colClass[NClass];
  for(int ci = 0; ci < NClass; ci++) colClass[ci] = TColor::GetColor(classHex[ci]);

  double ymax = hPbPb[0]->GetMaximum();
  for(int ci = 1; ci < NClass; ci++) ymax = TMath::Max(ymax, hPbPb[ci]->GetMaximum());
  ymax = TMath::Max(ymax, hpp->GetMaximum());

  // ratioMin/ratioMax are fixed above. Warn rather than silently clip if a
  // point in the plotted window exceeds the fixed ceiling.
  for(int ci = 0; ci < NClass; ci++){
    for(int b = 1; b <= hRatio[ci]->GetNbinsX(); b++){
      double edge = hRatio[ci]->GetXaxis()->GetBinLowEdge(b);
      if(edge < plotPtMin || edge >= 500.) continue;
      double v = hRatio[ci]->GetBinContent(b);
      if(v > ratioMax)
        printf("NOTE: PbPb %s (rawPtCut) / pp at %.0f-%.0f GeV is %.2f, above the ratio ceiling %.1f\n",
               classLabel[ci], edge, hRatio[ci]->GetXaxis()->GetBinUpEdge(b), v, ratioMax);
    }
  }

  const double lm = 0.15, rm = 0.05, tm = 0.07, bm = 0.28, split = 0.35;
  TCanvas *c = new TCanvas("cMinBiasRatioRawPtCut", "", 750, 800);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->SetLogy(); pUp->Draw(); pUp->cd();

  hpp->SetLineColor(colpp); hpp->SetMarkerColor(colpp);
  hpp->SetMarkerStyle(25); hpp->SetMarkerSize(1.1); hpp->SetLineWidth(2);
  hpp->SetTitle("");
  hpp->GetYaxis()->SetTitle("1/N_{jet}(p_{T}>200) dN_{jet}/dp_{T}  [GeV^{-1}]");
  hpp->GetYaxis()->SetTitleSize(0.050);
  hpp->GetYaxis()->SetLabelSize(0.044);
  hpp->GetYaxis()->SetTitleOffset(1.35);
  hpp->GetXaxis()->SetRangeUser(plotPtMin, 500.);
  hpp->SetMaximum(ymax * 3.);
  hpp->SetMinimum(1e-6);
  hpp->Draw("ep");

  TLegend *leg = new TLegend(0.50, 0.60, 0.94, 0.90);
  leg->SetBorderSize(0); leg->SetTextSize(0.038); leg->SetFillStyle(0);
  leg->AddEntry(hpp, "pp MinBias", "lp");
  for(int ci = 0; ci < NClass; ci++){
    hPbPb[ci]->SetLineColor(colClass[ci]); hPbPb[ci]->SetMarkerColor(colClass[ci]);
    hPbPb[ci]->SetMarkerStyle(classMark[ci]); hPbPb[ci]->SetMarkerSize(classSize[ci]);
    hPbPb[ci]->SetLineWidth(2);
    hPbPb[ci]->Draw("ep same");
    leg->AddEntry(hPbPb[ci], Form("PbPb %s, raw p_{T}>40", classLabel[ci]), "lp");
  }
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038);
  lat.DrawLatex(0.19, 0.89, "Shape only -- normalized above 200 GeV. PbPb: raw p_{T}>40 GeV cut, pp: none");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  bool first = true;
  for(int ci = 0; ci < NClass; ci++){
    hRatio[ci]->SetLineColor(colClass[ci]); hRatio[ci]->SetMarkerColor(colClass[ci]);
    hRatio[ci]->SetMarkerStyle(classMark[ci]); hRatio[ci]->SetMarkerSize(classSize[ci]);
    hRatio[ci]->SetLineWidth(2);
    hRatio[ci]->SetTitle("");
    hRatio[ci]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
    hRatio[ci]->GetYaxis()->SetTitle("PbPb / pp");
    hRatio[ci]->GetXaxis()->SetTitleSize(0.050*sc);
    hRatio[ci]->GetXaxis()->SetLabelSize(0.045*sc);
    hRatio[ci]->GetYaxis()->SetTitleSize(0.050*sc);
    hRatio[ci]->GetYaxis()->SetLabelSize(0.045*sc);
    hRatio[ci]->GetYaxis()->SetTitleOffset(1.30/sc);
    hRatio[ci]->GetYaxis()->SetNdivisions(505);
    hRatio[ci]->GetXaxis()->SetRangeUser(plotPtMin, 500.);
    hRatio[ci]->SetMinimum(ratioMin); hRatio[ci]->SetMaximum(ratioMax);
    hRatio[ci]->Draw(first ? "ep" : "ep same");
    first = false;
  }

  TLine *one = new TLine(plotPtMin, 1., 500., 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
