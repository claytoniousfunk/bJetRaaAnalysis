// Same as plotJetPt_minBiasOnly_PbPbOverPP_fakeJetEffect.C, but the fake-jet
// estimate uses the Delta-pT (PF-PFCs) background subtraction instead of the
// random-cone one: rootFiles/fakeJets/fakeJets_dPT.root, built by
// makeFakeJetFile_dPT.C from h_fastJetPt_PF_JEC_bkgSub_dPT. NOT the file
// calculateRAA.C actually subtracts -- that is still the RC-based fakeJets.root.
// This exists to compare the two background-subtraction methods against each
// other, not to replace the nominal chain.
//
// The dPT method estimates substantially more fakes than RC in central events:
// 8.73/event vs 3.46/event in 0-10% (both measured from the same mixed-event
// file). Worth keeping in mind when reading the size of the correction below.
//
// Top panel: PbPb/pp, without and with the subtraction, overlaid per class
// (open marker = without, filled = with).
// Bottom panel: (with)/(without) of those same PbPb/pp curves -- algebraically
// this is just hPbPb_withSub/hPbPb_noSub (pp cancels), but it is computed as
// the ratio of the two displayed top-panel curves so what is plotted top and
// bottom is always numerically consistent.
//
// --- How the subtraction is applied here ---
// h_fakeJets_C{1..4} in fakeJets_dPT.root is a PER-EVENT fake-jet rate (fakes per
// event per GeV, JEC-corrected pT axis), built from a *different* mixed-event
// scan (PbPb_..._mixedEventPFClustering_..._2026-8-11...). To subtract it from
// this file's raw (un-normalized) jet counts, it is scaled up by THIS file's
// own N_events per class (summed from h_vz_C* over the same ultra-fine slices,
// same convention as makeFakeJetFile_dPT.C) before the bin-by-bin subtraction. This
// is a cross-file rate transfer -- it assumes the per-event fake rate measured
// in the 8-11 mixed-event scan applies to the 8-13 scan's own event sample,
// which is the same assumption calculateRAA.C's stitching makes when it pulls
// in fakeJets_dPT.root. Bins that go negative after subtraction are zeroed, same
// convention as makeFakeJetFile.C.
//
// No pp-side subtraction: the existing fake-jet estimate is PbPb-only
// (calculateRAA.C only ever subtracts it when isPbPb), consistent with pp
// MinBias having no PbPb-like underlying event to build a combinatorial jet
// from.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_PbPbOverPP_fakeJetEffect_dPT.C'
// Run from: src/plots/jetPt/

const char *pbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-8-13_ultraFineCentBins.root";
const char *ppFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root";
const char *ppHistName = "h_inclRecoJetPt";
const char *fakeJetsFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets_dPT.root";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_PbPbOverPP_fakeJetEffect_dPT.pdf";

// PbPb coarse class -> [first, last] ultra-fine slice index, same convention
// used throughout this figure family and by makeFakeJetFile.C:
//   C1 = 0-10%  : slices  1- 2      C3 = 30-50% : slices  7-10
//   C2 = 10-30% : slices  3- 6      C4 = 50-80% : slices 11-16
// fakeJets_dPT.root indexes the same four classes as h_fakeJets_C1..C4.
const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const double ptNormMin = 200.;
const double plotPtMin = 30.;

const double topMin = 0., topMax = 10.;
const double botMin = 0., botMax = 1.3;

const int    NEdge = 29;
double       ptEdge[NEdge] = {
  20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,
  110,120,130,140,150,160,180,200,240,280,350,500
};

// Okabe-Ito, centre-symmetric filled/open pairs per class: filled = with the
// fake-jet subtraction, open = without.
const char *classHex[NClass]   = { "#D55E00", "#0072B2", "#009E73", "#CC79A7" };
const int   markFilled[NClass] = { 20, 21, 33, 34 };
const int   markOpen[NClass]   = { 24, 25, 27, 28 };
const double classSize[NClass] = { 1.0, 1.0, 1.2, 1.1 };

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

static TH1D* buildShape(TH1D *raw, const char *name)
{
  TH1D *h = (TH1D*) raw->Rebin(NEdge-1, name, ptEdge);
  h->SetDirectory(nullptr);
  int b0 = h->FindBin(ptNormMin);
  double norm = h->Integral(b0, h->GetNbinsX()+1);
  if(norm <= 0.){ printf("ERROR: %s has zero integral above %.0f GeV\n", name, ptNormMin); return nullptr; }
  h->Scale(1./norm);
  h->Scale(1., "width");
  return h;
}

void plotJetPt_minBiasOnly_PbPbOverPP_fakeJetEffect_dPT()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fP = TFile::Open(pbpbFile);
  TFile *fp = TFile::Open(ppFile);
  TFile *fF = TFile::Open(fakeJetsFile);
  if(!fP || fP->IsZombie()){ printf("ERROR: cannot open %s\n", pbpbFile);    return; }
  if(!fp || fp->IsZombie()){ printf("ERROR: cannot open %s\n", ppFile);      return; }
  if(!fF || fF->IsZombie()){ printf("ERROR: cannot open %s\n", fakeJetsFile); return; }

  TH1D *hpp0 = nullptr;
  fp->GetObject(ppHistName, hpp0);
  if(!hpp0){ printf("ERROR: %s not found in %s\n", ppHistName, ppFile); return; }
  TH1D *hpp = buildShape(hpp0, "hpp");
  if(!hpp) return;

  TH1D *hRatioNoSub[NClass], *hRatioWithSub[NClass], *hRatioOfRatios[NClass];
  for(int ci = 0; ci < NClass; ci++){

    TH1D *rawNoSub = sumSlices(fP, "h_inclRecoJetPt", ci, Form("rawNoSub_%d", ci));
    if(!rawNoSub) return;

    double N_events = 0.;
    for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
      TH1D *hvz = nullptr;
      fP->GetObject(Form("h_vz_C%d", si), hvz);
      if(!hvz){ printf("ERROR: h_vz_C%d not found\n", si); return; }
      N_events += hvz->Integral();
    }

    TH1D *hFake = nullptr;
    fF->GetObject(Form("h_fakeJets_C%d", ci+1), hFake);
    if(!hFake){ printf("ERROR: h_fakeJets_C%d not found in %s\n", ci+1, fakeJetsFile); return; }
    TH1D *fakeScaled = (TH1D*) hFake->Clone(Form("fakeScaled_%d", ci));
    fakeScaled->SetDirectory(nullptr);
    fakeScaled->Scale(N_events);

    printf("%-8s N_events=%.0f  fakes/evt=%.4f  scaled fake yield=%.1f  (%.2f%% of raw)\n",
           classLabel[ci], N_events, hFake->Integral(), fakeScaled->Integral(),
           100.*fakeScaled->Integral()/rawNoSub->Integral());

    TH1D *rawWithSub = (TH1D*) rawNoSub->Clone(Form("rawWithSub_%d", ci));
    rawWithSub->SetDirectory(nullptr);
    rawWithSub->Add(fakeScaled, -1.);
    int nNeg = 0;
    for(int b = 1; b <= rawWithSub->GetNbinsX(); b++){
      if(rawWithSub->GetBinContent(b) < 0.){
        rawWithSub->SetBinContent(b, 0.);
        rawWithSub->SetBinError(b, 0.);
        nNeg++;
      }
    }
    if(nNeg > 0) printf("  zeroed %d bins that went negative after subtraction\n", nNeg);

    TH1D *hPbPbNoSub   = buildShape(rawNoSub,   Form("hPbPbNoSub_%d", ci));
    TH1D *hPbPbWithSub = buildShape(rawWithSub, Form("hPbPbWithSub_%d", ci));
    if(!hPbPbNoSub || !hPbPbWithSub) return;

    hRatioNoSub[ci]   = (TH1D*) hPbPbNoSub  ->Clone(Form("hRatioNoSub_%d", ci));
    hRatioNoSub[ci]->Divide(hpp);
    hRatioWithSub[ci] = (TH1D*) hPbPbWithSub->Clone(Form("hRatioWithSub_%d", ci));
    hRatioWithSub[ci]->Divide(hpp);

    hRatioOfRatios[ci] = (TH1D*) hRatioWithSub[ci]->Clone(Form("hRatioOfRatios_%d", ci));
    hRatioOfRatios[ci]->Divide(hRatioNoSub[ci]);
  }

  printf("\n%-12s", "pT [GeV]");
  for(int ci = 0; ci < NClass; ci++) printf(" %16s", classLabel[ci]);
  printf("\n");
  for(int b = 1; b <= hpp->GetNbinsX(); b++){
    double edge = hpp->GetXaxis()->GetBinLowEdge(b);
    if(edge < plotPtMin) continue;
    printf("%4.0f-%-6.0f", edge, hpp->GetXaxis()->GetBinUpEdge(b));
    for(int ci = 0; ci < NClass; ci++)
      printf(" %8.3f->%.3f", hRatioNoSub[ci]->GetBinContent(b), hRatioWithSub[ci]->GetBinContent(b));
    printf("\n");
  }

  // Warn rather than silently clip if the fixed display ranges don't cover the
  // actual data in the plotted window.
  for(int ci = 0; ci < NClass; ci++){
    for(int b = 1; b <= hpp->GetNbinsX(); b++){
      double edge = hpp->GetXaxis()->GetBinLowEdge(b);
      if(edge < plotPtMin || edge >= 500.) continue;
      double vNo = hRatioNoSub[ci]->GetBinContent(b), vWith = hRatioWithSub[ci]->GetBinContent(b);
      if(vNo > topMax || vWith > topMax)
        printf("NOTE: %s PbPb/pp at %.0f-%.0f exceeds top-panel ceiling %.1f (noSub=%.2f, withSub=%.2f)\n",
               classLabel[ci], edge, hpp->GetXaxis()->GetBinUpEdge(b), topMax, vNo, vWith);
      double r = hRatioOfRatios[ci]->GetBinContent(b);
      if(r != 0. && (r > botMax || r < botMin))
        printf("NOTE: %s (with/without) at %.0f-%.0f is %.2f, outside bottom-panel range [%.1f,%.1f]\n",
               classLabel[ci], edge, hpp->GetXaxis()->GetBinUpEdge(b), r, botMin, botMax);
    }
  }

  int colClass[NClass];
  for(int ci = 0; ci < NClass; ci++) colClass[ci] = TColor::GetColor(classHex[ci]);

  const double lm = 0.15, rm = 0.05, tm = 0.07, bm = 0.28, split = 0.40;
  TCanvas *c = new TCanvas("cFakeJetEffectDPT", "", 800, 800);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->Draw(); pUp->cd();

  TLegend *leg = new TLegend(0.40, 0.55, 0.94, 0.90);
  leg->SetBorderSize(0); leg->SetTextSize(0.032); leg->SetFillStyle(0); leg->SetNColumns(2);

  bool first = true;
  for(int ci = 0; ci < NClass; ci++){
    hRatioNoSub[ci]->SetLineColor(colClass[ci]); hRatioNoSub[ci]->SetMarkerColor(colClass[ci]);
    hRatioNoSub[ci]->SetMarkerStyle(markOpen[ci]); hRatioNoSub[ci]->SetMarkerSize(classSize[ci]);
    hRatioNoSub[ci]->SetLineWidth(2); hRatioNoSub[ci]->SetLineStyle(2);
    hRatioNoSub[ci]->SetTitle("");
    hRatioNoSub[ci]->GetYaxis()->SetTitle("PbPb / pp");
    hRatioNoSub[ci]->GetYaxis()->SetTitleSize(0.048);
    hRatioNoSub[ci]->GetYaxis()->SetLabelSize(0.042);
    hRatioNoSub[ci]->GetYaxis()->SetTitleOffset(1.35);
    hRatioNoSub[ci]->GetXaxis()->SetRangeUser(plotPtMin, 500.);
    hRatioNoSub[ci]->SetMinimum(topMin); hRatioNoSub[ci]->SetMaximum(topMax);
    hRatioNoSub[ci]->Draw(first ? "ep" : "ep same");
    first = false;

    hRatioWithSub[ci]->SetLineColor(colClass[ci]); hRatioWithSub[ci]->SetMarkerColor(colClass[ci]);
    hRatioWithSub[ci]->SetMarkerStyle(markFilled[ci]); hRatioWithSub[ci]->SetMarkerSize(classSize[ci]);
    hRatioWithSub[ci]->SetLineWidth(2);
    hRatioWithSub[ci]->Draw("ep same");

    leg->AddEntry(hRatioNoSub[ci],   Form("%s, no sub",   classLabel[ci]), "lp");
    leg->AddEntry(hRatioWithSub[ci], Form("%s, with sub", classLabel[ci]), "lp");
  }
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.036);
  lat.DrawLatex(0.19, 0.91, "Open = no fake-jet subtraction, filled = with (dPT-based fakeJets_dPT.root)");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  bool firstB = true;
  for(int ci = 0; ci < NClass; ci++){
    hRatioOfRatios[ci]->SetLineColor(colClass[ci]); hRatioOfRatios[ci]->SetMarkerColor(colClass[ci]);
    hRatioOfRatios[ci]->SetMarkerStyle(markFilled[ci]); hRatioOfRatios[ci]->SetMarkerSize(classSize[ci]);
    hRatioOfRatios[ci]->SetLineWidth(2);
    hRatioOfRatios[ci]->SetTitle("");
    hRatioOfRatios[ci]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
    hRatioOfRatios[ci]->GetYaxis()->SetTitle("with / without");
    hRatioOfRatios[ci]->GetXaxis()->SetTitleSize(0.048*sc);
    hRatioOfRatios[ci]->GetXaxis()->SetLabelSize(0.042*sc);
    hRatioOfRatios[ci]->GetYaxis()->SetTitleSize(0.048*sc);
    hRatioOfRatios[ci]->GetYaxis()->SetLabelSize(0.042*sc);
    hRatioOfRatios[ci]->GetYaxis()->SetTitleOffset(1.35/sc);
    hRatioOfRatios[ci]->GetYaxis()->SetNdivisions(505);
    hRatioOfRatios[ci]->GetXaxis()->SetRangeUser(plotPtMin, 500.);
    hRatioOfRatios[ci]->SetMinimum(botMin); hRatioOfRatios[ci]->SetMaximum(botMax);
    hRatioOfRatios[ci]->Draw(firstB ? "ep" : "ep same");
    firstB = false;
  }

  TLine *one = new TLine(plotPtMin, 1., 500., 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
