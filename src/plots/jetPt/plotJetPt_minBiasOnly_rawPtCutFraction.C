// Fraction of PbPb MinBias jets surviving the raw-pT>40 GeV cut, as a function
// of the JEC-CORRECTED pT (the x axis of both source histograms) -- i.e.
//   h_inclRecoJetPt_rawPtCut_C* / h_inclRecoJetPt_C*
// per coarse centrality class. No pp curve: pp_scan.C has no rawPtCut variant,
// so there is nothing to divide on the pp side.
//
// This is the direct complement to plotJetPt_minBiasOnly_PbPbOverPP_rawPtCut.C,
// which showed the cut essentially truncates the spectrum below ~45 GeV rather
// than differentially suppressing an excess. This plot is that same conclusion
// read off directly: the survival fraction, not a shape ratio against pp.
//
// Both histograms are filled from the same jet loop with the same weight, and
// h_inclRecoJetPt_rawPtCut's fill is gated by an extra cut on top of
// h_inclRecoJetPt's -- so bin by bin, rawPtCut <= nominal always (verified
// against the current scan output before writing this). That makes each bin a
// genuine efficiency, not an independent ratio, so errors use the binomial
// option (TH1::Divide(...,"B")) rather than the usual independent-Gaussian
// combination.
//
// No area/high-pT normalization anywhere in this file: unlike the PbPb/pp shape
// comparisons, numerator and denominator here already share an absolute scale
// (same events, same jets, same weight), so renormalizing either one would
// destroy the thing being measured.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_rawPtCutFraction.C'
// Run from: src/plots/jetPt/

const char *pbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-8-13_ultraFineCentBins.root";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_rawPtCutFraction.pdf";

// PbPb coarse class -> [first, last] ultra-fine slice index, same convention as
// makeFakeJetFile.C and the sibling PbPbOverPP macros:
//   C1 = 0-10%  : slices  1- 2      C3 = 30-50% : slices  7-10
//   C2 = 10-30% : slices  3- 6      C4 = 50-80% : slices 11-16
const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const double plotPtMin = 20., plotPtMax = 500.;
const double ratioMin  = 0.,  ratioMax  = 1.15;

// Same variable-width binning as the sibling macros, for direct comparability.
const int    NEdge = 29;
double       ptEdge[NEdge] = {
  20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,
  110,120,130,140,150,160,180,200,240,280,350,500
};

// Okabe-Ito, centre-symmetric markers -- same per-class assignment as the
// sibling PbPbOverPP macros, so the class identity reads the same across all
// three plots.
const char *classHex[NClass]  = { "#D55E00", "#0072B2", "#009E73", "#CC79A7" };
const int   classMark[NClass] = { 20, 21, 33, 34 };
const double classSize[NClass]= { 1.0, 1.0, 1.2, 1.1 };

static TH1D* sumSlices(TFile *f, const char *base, int ci, const char *name)
{
  TH1D *sum = nullptr;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *h = nullptr;
    f->GetObject(Form("%s_C%d", base, si), h);
    if(!h){ printf("ERROR: %s_C%d not found\n", base, si); return nullptr; }
    if(!sum){ sum = (TH1D*) h->Clone(name); sum->SetDirectory(nullptr); }
    else sum->Add(h);
  }
  return sum;
}

void plotJetPt_minBiasOnly_rawPtCutFraction()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fP = TFile::Open(pbpbFile);
  if(!fP || fP->IsZombie()){ printf("ERROR: cannot open %s\n", pbpbFile); return; }

  TH1D *hFrac[NClass];
  for(int ci = 0; ci < NClass; ci++){
    TH1D *num = sumSlices(fP, "h_inclRecoJetPt_rawPtCut", ci, Form("num_%d", ci));
    TH1D *den = sumSlices(fP, "h_inclRecoJetPt",          ci, Form("den_%d", ci));
    if(!num || !den) return;

    num = (TH1D*) num->Rebin(NEdge-1, Form("numR_%d", ci), ptEdge);
    den = (TH1D*) den->Rebin(NEdge-1, Form("denR_%d", ci), ptEdge);
    num->SetDirectory(nullptr); den->SetDirectory(nullptr);

    hFrac[ci] = (TH1D*) num->Clone(Form("hFrac_C%d", ci));
    hFrac[ci]->Divide(num, den, 1., 1., "B");
  }

  printf("\n%-12s", "pT [GeV]");
  for(int ci = 0; ci < NClass; ci++) printf(" %14s", classLabel[ci]);
  printf("\n");
  for(int b = 1; b <= hFrac[0]->GetNbinsX(); b++){
    printf("%4.0f-%-6.0f", hFrac[0]->GetXaxis()->GetBinLowEdge(b), hFrac[0]->GetXaxis()->GetBinUpEdge(b));
    for(int ci = 0; ci < NClass; ci++)
      printf(" %8.3f+-%.3f", hFrac[ci]->GetBinContent(b), hFrac[ci]->GetBinError(b));
    printf("\n");
  }

  int colClass[NClass];
  for(int ci = 0; ci < NClass; ci++) colClass[ci] = TColor::GetColor(classHex[ci]);

  TCanvas *c = new TCanvas("cRawPtCutFraction", "", 750, 600);
  c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
  c->SetTopMargin(0.08);  c->SetBottomMargin(0.13);

  bool first = true;
  TLegend *leg = new TLegend(0.55, 0.20, 0.94, 0.45);
  leg->SetBorderSize(0); leg->SetTextSize(0.036); leg->SetFillStyle(0);
  for(int ci = 0; ci < NClass; ci++){
    hFrac[ci]->SetLineColor(colClass[ci]); hFrac[ci]->SetMarkerColor(colClass[ci]);
    hFrac[ci]->SetMarkerStyle(classMark[ci]); hFrac[ci]->SetMarkerSize(classSize[ci]);
    hFrac[ci]->SetLineWidth(2);
    hFrac[ci]->SetTitle("");
    hFrac[ci]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]  (JEC-corrected)");
    hFrac[ci]->GetYaxis()->SetTitle("fraction with raw p_{T}^{jet} > 40 GeV");
    hFrac[ci]->GetXaxis()->SetTitleSize(0.045);
    hFrac[ci]->GetXaxis()->SetLabelSize(0.040);
    hFrac[ci]->GetYaxis()->SetTitleSize(0.045);
    hFrac[ci]->GetYaxis()->SetLabelSize(0.040);
    hFrac[ci]->GetYaxis()->SetTitleOffset(1.45);
    hFrac[ci]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hFrac[ci]->SetMinimum(ratioMin); hFrac[ci]->SetMaximum(ratioMax);
    hFrac[ci]->Draw(first ? "ep" : "ep same");
    leg->AddEntry(hFrac[ci], Form("PbPb %s", classLabel[ci]), "lp");
    first = false;
  }
  leg->Draw();

  TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038);
  lat.DrawLatex(0.17, 0.90, "PbPb MinBias -- raw p_{T}>40 GeV cut survival fraction");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
