// Same as plotJetPt_minBiasOnly_rawPt_fakeJetEffect.C, but the fake-jet
// estimate uses the geometrically-corrected random cone (RC_geoCorr) instead of
// the plain random cone: rootFiles/fakeJets/fakeJets_RCgeoCorr_rawPt.root, built
// by makeFakeJetFile_RCgeoCorr_rawPt.C from h_fastJetPt_PF_bkgSub_RC_geoCorr.
// PbPb only, no pp side, same as the plain-RC version. Not the file
// calculateRAA.C actually subtracts -- exists to compare background-subtraction
// methods.
//
// Everything here comes from ONE file: the PbPb_pfCandAnalyzer.C mixed-event
// scan (the same file makeFakeJetFile.C / makeFakeJetFile_rawPt.C read). That
// file already fills h_inclRawJetPt (genuine raw-pT axis, not JEC-corrected)
// alongside h_inclRecoJetPt, so no rerun was needed for the PbPb side -- unlike
// the JEC-axis version of this plot, which needed h_inclRecoJetPt_rawPtCut
// added to PbPb_scan.C. See that macro for why the axis has to match: the fake
// estimate here is h_fakeJets_C{1..4} from fakeJets_RCgeoCorr_rawPt.root (built
// by makeFakeJetFile_RCgeoCorr_rawPt.C, which points at
// h_fastJetPt_PF_bkgSub_RC_geoCorr -- the RAW, non-JEC, geoCorr-weighted
// FastJet spectrum), so subtracting it from a genuine raw-pT reco-jet
// spectrum keeps both sides in the same pT variable throughout.
//
// Because everything is one file, one sample, this uses genuine per-event
// normalization (divide by N_events, from the same h_vz_C* slices) rather than
// the high-pT shape-matching trick the PbPb/pp macros use to sidestep an
// untrustworthy cross-dataset absolute scale -- there is no second dataset
// here to be untrustworthy against.
//
// Top panel: dN/dpT per event, without and with the subtraction, per class
// (open marker = without, filled = with).
// Bottom panel: with/without, computed from the two displayed top-panel curves
// so what's plotted top and bottom stays numerically consistent.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_rawPt_fakeJetEffect_RCgeoCorr.C'
// Run from: src/plots/jetPt/

const char *pbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";
const char *fakeJetsFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/fakeJets/fakeJets_RCgeoCorr_rawPt.root";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_rawPt_fakeJetEffect_RCgeoCorr.pdf";

// PbPb coarse class -> [first, last] ultra-fine slice index, same convention
// used throughout this figure family:
//   C1 = 0-10%  : slices  1- 2      C3 = 30-50% : slices  7-10
//   C2 = 10-30% : slices  3- 6      C4 = 50-80% : slices 11-16
const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const double plotPtMin = 15., plotPtMax = 500.;
const double botMin = 0., botMax = 1.3;

// Extends 5 GeV lower than the sibling macros' shared ptEdge: h_inclRawJetPt
// has real content in [15,20) (kinematic cuts gate on corrected pT, not raw,
// so a jet just above the corrected-pT threshold can sit below 20 in raw pT).
const int    NEdge = 30;
double       ptEdge[NEdge] = {
  15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,
  110,120,130,140,150,160,180,200,240,280,350,500
};

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

void plotJetPt_minBiasOnly_rawPt_fakeJetEffect_RCgeoCorr()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fP = TFile::Open(pbpbFile);
  TFile *fF = TFile::Open(fakeJetsFile);
  if(!fP || fP->IsZombie()){ printf("ERROR: cannot open %s\n", pbpbFile);    return; }
  if(!fF || fF->IsZombie()){ printf("ERROR: cannot open %s\n", fakeJetsFile); return; }

  TH1D *hNoSub[NClass], *hWithSub[NClass], *hRatio[NClass];
  double ymax = 0.;

  for(int ci = 0; ci < NClass; ci++){

    TH1D *rawNoSub = sumSlices(fP, "h_inclRawJetPt", ci, Form("rawNoSub_%d", ci));
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
        rawWithSub->SetBinContent(b, 0.); rawWithSub->SetBinError(b, 0.); nNeg++;
      }
    }
    if(nNeg > 0) printf("  zeroed %d bins that went negative after subtraction\n", nNeg);

    hNoSub[ci]   = (TH1D*) rawNoSub  ->Rebin(NEdge-1, Form("hNoSub_%d", ci),   ptEdge);
    hWithSub[ci] = (TH1D*) rawWithSub->Rebin(NEdge-1, Form("hWithSub_%d", ci), ptEdge);
    hNoSub[ci]->SetDirectory(nullptr); hWithSub[ci]->SetDirectory(nullptr);

    // genuine per-event dN/dpT: divide by N_events, then by bin width
    hNoSub[ci]  ->Scale(1./N_events); hNoSub[ci]  ->Scale(1., "width");
    hWithSub[ci]->Scale(1./N_events); hWithSub[ci]->Scale(1., "width");

    hRatio[ci] = (TH1D*) hWithSub[ci]->Clone(Form("hRatio_%d", ci));
    hRatio[ci]->Divide(hNoSub[ci]);

    ymax = TMath::Max(ymax, hNoSub[ci]->GetMaximum());
  }

  printf("\n%-12s", "pT [GeV]");
  for(int ci = 0; ci < NClass; ci++) printf(" %16s", classLabel[ci]);
  printf("\n");
  for(int b = 1; b <= hNoSub[0]->GetNbinsX(); b++){
    printf("%4.0f-%-6.0f", hNoSub[0]->GetXaxis()->GetBinLowEdge(b), hNoSub[0]->GetXaxis()->GetBinUpEdge(b));
    for(int ci = 0; ci < NClass; ci++)
      printf(" %8.4g->%.4g", hNoSub[ci]->GetBinContent(b), hWithSub[ci]->GetBinContent(b));
    printf("\n");
  }

  for(int ci = 0; ci < NClass; ci++){
    for(int b = 1; b <= hRatio[ci]->GetNbinsX(); b++){
      double edge = hRatio[ci]->GetXaxis()->GetBinLowEdge(b);
      if(edge < plotPtMin || edge >= plotPtMax) continue;
      double r = hRatio[ci]->GetBinContent(b);
      if(r != 0. && (r > botMax || r < botMin))
        printf("NOTE: %s (with/without) at %.0f-%.0f is %.2f, outside bottom-panel range [%.1f,%.1f]\n",
               classLabel[ci], edge, hRatio[ci]->GetXaxis()->GetBinUpEdge(b), r, botMin, botMax);
    }
  }

  int colClass[NClass];
  for(int ci = 0; ci < NClass; ci++) colClass[ci] = TColor::GetColor(classHex[ci]);

  const double lm = 0.15, rm = 0.05, tm = 0.07, bm = 0.28, split = 0.40;
  TCanvas *c = new TCanvas("cRawPtFakeJetEffectRCgeoCorr", "", 800, 800);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->SetLogy(); pUp->Draw(); pUp->cd();

  TLegend *leg = new TLegend(0.42, 0.55, 0.94, 0.90);
  leg->SetBorderSize(0); leg->SetTextSize(0.032); leg->SetFillStyle(0); leg->SetNColumns(2);

  bool first = true;
  for(int ci = 0; ci < NClass; ci++){
    hNoSub[ci]->SetLineColor(colClass[ci]); hNoSub[ci]->SetMarkerColor(colClass[ci]);
    hNoSub[ci]->SetMarkerStyle(markOpen[ci]); hNoSub[ci]->SetMarkerSize(classSize[ci]);
    hNoSub[ci]->SetLineWidth(2); hNoSub[ci]->SetLineStyle(2);
    hNoSub[ci]->SetTitle("");
    hNoSub[ci]->GetYaxis()->SetTitle("1/N_{evt} dN_{jet}/dp_{T}  [GeV^{-1}]");
    hNoSub[ci]->GetYaxis()->SetTitleSize(0.048);
    hNoSub[ci]->GetYaxis()->SetLabelSize(0.042);
    hNoSub[ci]->GetYaxis()->SetTitleOffset(1.35);
    hNoSub[ci]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hNoSub[ci]->SetMinimum(ymax * 1e-7); hNoSub[ci]->SetMaximum(ymax * 3.);
    hNoSub[ci]->Draw(first ? "ep" : "ep same");
    first = false;

    hWithSub[ci]->SetLineColor(colClass[ci]); hWithSub[ci]->SetMarkerColor(colClass[ci]);
    hWithSub[ci]->SetMarkerStyle(markFilled[ci]); hWithSub[ci]->SetMarkerSize(classSize[ci]);
    hWithSub[ci]->SetLineWidth(2);
    hWithSub[ci]->Draw("ep same");

    leg->AddEntry(hNoSub[ci],   Form("%s, no sub",   classLabel[ci]), "lp");
    leg->AddEntry(hWithSub[ci], Form("%s, with sub", classLabel[ci]), "lp");
  }
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.036);
  lat.DrawLatex(0.19, 0.91, "Raw p_{T} axis, RC_geoCorr background subtraction. Open = no sub, filled = with");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  bool firstB = true;
  for(int ci = 0; ci < NClass; ci++){
    hRatio[ci]->SetLineColor(colClass[ci]); hRatio[ci]->SetMarkerColor(colClass[ci]);
    hRatio[ci]->SetMarkerStyle(markFilled[ci]); hRatio[ci]->SetMarkerSize(classSize[ci]);
    hRatio[ci]->SetLineWidth(2);
    hRatio[ci]->SetTitle("");
    hRatio[ci]->GetXaxis()->SetTitle("raw p_{T}^{jet} [GeV]");
    hRatio[ci]->GetYaxis()->SetTitle("with / without");
    hRatio[ci]->GetXaxis()->SetTitleSize(0.048*sc);
    hRatio[ci]->GetXaxis()->SetLabelSize(0.042*sc);
    hRatio[ci]->GetYaxis()->SetTitleSize(0.048*sc);
    hRatio[ci]->GetYaxis()->SetLabelSize(0.042*sc);
    hRatio[ci]->GetYaxis()->SetTitleOffset(1.35/sc);
    hRatio[ci]->GetYaxis()->SetNdivisions(505);
    hRatio[ci]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hRatio[ci]->SetMinimum(botMin); hRatio[ci]->SetMaximum(botMax);
    hRatio[ci]->Draw(firstB ? "ep" : "ep same");
    firstB = false;
  }

  TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
