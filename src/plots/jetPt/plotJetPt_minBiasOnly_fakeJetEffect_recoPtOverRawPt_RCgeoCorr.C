// Same double ratio as plotJetPt_minBiasOnly_fakeJetEffect_recoPtOverRawPt.C,
// but for the RC_geoCorr background-subtraction method instead of plain RC:
//
//   DR(pT) = [with/without]_recoPt(pT)  /  [with/without]_rawPt(pT)
//
// Each factor is built exactly as in the two single-axis RC_geoCorr macros this
// combines:
//   plotJetPt_minBiasOnly_PbPbOverPP_fakeJetEffect_RCgeoCorr.C   (recoPt / JEC axis)
//   plotJetPt_minBiasOnly_rawPt_fakeJetEffect_RCgeoCorr.C         (rawPt axis)
// See those for the full reasoning on the subtraction itself; this file only
// repeats the with/without computation, then divides the two.
//
// The two [with/without] curves come from DIFFERENT samples (the recoPt side
// from the 2026-8-13 PbPb_scan.C output, the rawPt side from the 2026-8-11
// PbPb_pfCandAnalyzer.C mixed-event output), each already normalized against
// its own N_events, so N_events does not cancel in this double ratio and is
// computed separately and correctly for each side, exactly as the originals do.
//
// DR is plotted against a shared bin edge array; the two axes it is built from
// (corrected pT and raw pT) are physically different quantities for the same
// jet, so DR(pT) answers "how much bigger/smaller is the correction's effect
// at this nominal pT value, measured one way vs the other" -- not a statement
// about the same jet population at two different pT values.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_fakeJetEffect_recoPtOverRawPt_RCgeoCorr.C'
// Run from: src/plots/jetPt/

// --- recoPt (JEC-corrected) side ---------------------------------------
const char *recoPbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-8-13_ultraFineCentBins.root";
const char *recoFakeJetsFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets_RCgeoCorr.root";
const char *recoHistBase = "h_inclRecoJetPt";

// --- rawPt side -----------------------------------------------------------
const char *rawPbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";
const char *rawFakeJetsFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets_RCgeoCorr_rawPt.root";
const char *rawHistBase = "h_inclRawJetPt";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_fakeJetEffect_recoPtOverRawPt_RCgeoCorr.pdf";

// PbPb coarse class -> [first, last] ultra-fine slice index, same convention
// used throughout this figure family.
const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const double plotPtMin = 30., plotPtMax = 500.;
const double drMin = 0.6, drMax = 1.4;

// Shared bin edges for both sides -- starts at 20 (not the rawPt macro's 15,
// which the recoPt side has no equivalent bin for).
const int    NEdge = 29;
double       ptEdge[NEdge] = {
  20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,
  110,120,130,140,150,160,180,200,240,280,350,500
};

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

// Builds the (with subtraction)/(without) ratio for one class, on one axis.
static TH1D* withOverWithoutRatio(TFile *fScan, TFile *fFake, const char *histBase,
                                   int ci, const char *tag)
{
  TH1D *rawNoSub = sumSlices(fScan, histBase, ci, Form("rawNoSub_%s_%d", tag, ci));
  if(!rawNoSub) return nullptr;

  double N_events = 0.;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *hvz = nullptr;
    fScan->GetObject(Form("h_vz_C%d", si), hvz);
    if(!hvz){ printf("ERROR: h_vz_C%d not found (%s)\n", si, tag); return nullptr; }
    N_events += hvz->Integral();
  }

  TH1D *hFake = nullptr;
  fFake->GetObject(Form("h_fakeJets_C%d", ci+1), hFake);
  if(!hFake){ printf("ERROR: h_fakeJets_C%d not found (%s)\n", ci+1, tag); return nullptr; }
  TH1D *fakeScaled = (TH1D*) hFake->Clone(Form("fakeScaled_%s_%d", tag, ci));
  fakeScaled->SetDirectory(nullptr);
  fakeScaled->Scale(N_events);

  TH1D *rawWithSub = (TH1D*) rawNoSub->Clone(Form("rawWithSub_%s_%d", tag, ci));
  rawWithSub->SetDirectory(nullptr);
  rawWithSub->Add(fakeScaled, -1.);
  for(int b = 1; b <= rawWithSub->GetNbinsX(); b++){
    if(rawWithSub->GetBinContent(b) < 0.){ rawWithSub->SetBinContent(b, 0.); rawWithSub->SetBinError(b, 0.); }
  }

  TH1D *hNoSub   = (TH1D*) rawNoSub  ->Rebin(NEdge-1, Form("hNoSub_%s_%d",   tag, ci), ptEdge);
  TH1D *hWithSub = (TH1D*) rawWithSub->Rebin(NEdge-1, Form("hWithSub_%s_%d", tag, ci), ptEdge);
  hNoSub->SetDirectory(nullptr); hWithSub->SetDirectory(nullptr);
  // N_events and bin-width factors are common to numerator and denominator of
  // the ratio below and would cancel exactly -- skipped rather than applied
  // and immediately divided back out.

  TH1D *ratio = (TH1D*) hWithSub->Clone(Form("ratio_%s_%d", tag, ci));
  ratio->Divide(hNoSub);
  printf("%-8s %-5s N_events=%.0f  fakes/evt=%.4f  scaled fake yield=%.1f (%.2f%% of raw)\n",
         classLabel[ci], tag, N_events, hFake->Integral(), fakeScaled->Integral(),
         100.*fakeScaled->Integral()/rawNoSub->Integral());
  return ratio;
}

void plotJetPt_minBiasOnly_fakeJetEffect_recoPtOverRawPt_RCgeoCorr()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fRecoScan = TFile::Open(recoPbpbFile);
  TFile *fRecoFake = TFile::Open(recoFakeJetsFile);
  TFile *fRawScan  = TFile::Open(rawPbpbFile);
  TFile *fRawFake  = TFile::Open(rawFakeJetsFile);
  if(!fRecoScan || fRecoScan->IsZombie()){ printf("ERROR: cannot open %s\n", recoPbpbFile);    return; }
  if(!fRecoFake || fRecoFake->IsZombie()){ printf("ERROR: cannot open %s\n", recoFakeJetsFile); return; }
  if(!fRawScan  || fRawScan ->IsZombie()){ printf("ERROR: cannot open %s\n", rawPbpbFile);      return; }
  if(!fRawFake  || fRawFake ->IsZombie()){ printf("ERROR: cannot open %s\n", rawFakeJetsFile);  return; }

  TH1D *hDR[NClass];
  for(int ci = 0; ci < NClass; ci++){
    TH1D *rReco = withOverWithoutRatio(fRecoScan, fRecoFake, recoHistBase, ci, "reco");
    TH1D *rRaw  = withOverWithoutRatio(fRawScan,  fRawFake,  rawHistBase,  ci, "raw");
    if(!rReco || !rRaw) return;

    hDR[ci] = (TH1D*) rReco->Clone(Form("hDR_%d", ci));
    hDR[ci]->Divide(rRaw);
  }

  printf("\n%-12s", "pT [GeV]");
  for(int ci = 0; ci < NClass; ci++) printf(" %14s", classLabel[ci]);
  printf("\n");
  for(int b = 1; b <= hDR[0]->GetNbinsX(); b++){
    double edge = hDR[0]->GetXaxis()->GetBinLowEdge(b);
    if(edge < plotPtMin) continue;
    printf("%4.0f-%-6.0f", edge, hDR[0]->GetXaxis()->GetBinUpEdge(b));
    for(int ci = 0; ci < NClass; ci++)
      printf(" %8.3f+-%.3f", hDR[ci]->GetBinContent(b), hDR[ci]->GetBinError(b));
    printf("\n");
  }

  for(int ci = 0; ci < NClass; ci++){
    for(int b = 1; b <= hDR[ci]->GetNbinsX(); b++){
      double edge = hDR[ci]->GetXaxis()->GetBinLowEdge(b);
      if(edge < plotPtMin || edge >= plotPtMax) continue;
      double v = hDR[ci]->GetBinContent(b);
      if(v != 0. && (v > drMax || v < drMin))
        printf("NOTE: %s DR at %.0f-%.0f is %.2f, outside display range [%.1f,%.1f]\n",
               classLabel[ci], edge, hDR[ci]->GetXaxis()->GetBinUpEdge(b), v, drMin, drMax);
    }
  }

  int colClass[NClass];
  for(int ci = 0; ci < NClass; ci++) colClass[ci] = TColor::GetColor(classHex[ci]);

  TCanvas *c = new TCanvas("cDRRCgeoCorr", "", 750, 600);
  c->SetLeftMargin(0.15); c->SetRightMargin(0.05);
  c->SetTopMargin(0.08);  c->SetBottomMargin(0.13);

  TLegend *leg = new TLegend(0.55, 0.68, 0.94, 0.90);
  leg->SetBorderSize(0); leg->SetTextSize(0.036); leg->SetFillStyle(0);

  bool first = true;
  for(int ci = 0; ci < NClass; ci++){
    hDR[ci]->SetLineColor(colClass[ci]); hDR[ci]->SetMarkerColor(colClass[ci]);
    hDR[ci]->SetMarkerStyle(classMark[ci]); hDR[ci]->SetMarkerSize(classSize[ci]);
    hDR[ci]->SetLineWidth(2);
    hDR[ci]->SetTitle("");
    hDR[ci]->GetXaxis()->SetTitle("nominal p_{T}^{jet} [GeV]");
    hDR[ci]->GetYaxis()->SetTitle("[with/without]_{recoPt} / [with/without]_{rawPt}");
    hDR[ci]->GetXaxis()->SetTitleSize(0.045);
    hDR[ci]->GetXaxis()->SetLabelSize(0.040);
    hDR[ci]->GetYaxis()->SetTitleSize(0.045);
    hDR[ci]->GetYaxis()->SetLabelSize(0.040);
    hDR[ci]->GetYaxis()->SetTitleOffset(1.45);
    hDR[ci]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hDR[ci]->SetMinimum(drMin); hDR[ci]->SetMaximum(drMax);
    hDR[ci]->Draw(first ? "ep" : "ep same");
    leg->AddEntry(hDR[ci], Form("PbPb %s", classLabel[ci]), "lp");
    first = false;
  }
  leg->Draw();

  TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.036);
  lat.DrawLatex(0.17, 0.90, "Double ratio (RC_geoCorr) of the fake-jet subtraction's (with/without) effect");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
