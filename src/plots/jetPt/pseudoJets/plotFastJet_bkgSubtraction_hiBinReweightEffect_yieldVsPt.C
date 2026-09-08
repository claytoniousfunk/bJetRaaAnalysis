// hiBin-reweight effect on the RC-subtracted fake-jet spectrum: same event
// sample and background-subtraction method (RC, mixed-event) read from two
// otherwise-identical scans that differ only in whether the hiBin reweight
// (to match the HardProbes Jet80 trigger's hiBin distribution) was applied.
// One panel per coarse centrality class, raw-hiBin vs reweighted-hiBin
// overlaid per panel plus a reweighted/raw ratio underneath -- same style as
// plotFastJet_bkgSubtraction_methodComparison_yieldVsPt.C.
//
// Usage: root -l -b -q 'plotFastJet_bkgSubtraction_hiBinReweightEffect_yieldVsPt.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *rawFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";
const char *reweightedFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_hiBinReweightToHardProbesJet80_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-9_ultraFineCentBins.root";

const char *histBase = "h_fastJetPt_PF_bkgSub_RC";

const char *outDir  = "../../../../figures/jetPt/pseudoJets/";
const char *outName = "fastJet_bkgSub_hiBinReweightEffect_yieldVsPt.pdf";

const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const double plotPtMin = 0., plotPtMax = 160.;

const int    NEdge = 25;
double       ptEdge[NEdge] = {
  0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,150,160
};

struct Condition { const char *path; const char *label; const char *hex; int marker; };
const int NCond = 2;
Condition conds[NCond] = {
  {rawFile,        "raw hiBin",         "#0072B2", 20},
  {reweightedFile, "reweighted hiBin",  "#D55E00", 21},
};

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

static double sumEvents(TFile *f, int ci)
{
  double N = 0.;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *hvz = nullptr;
    f->GetObject(Form("h_vz_C%d", si), hvz);
    N += hvz->Integral();
  }
  return N;
}

void plotFastJet_bkgSubtraction_hiBinReweightEffect_yieldVsPt()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f[NCond];
  for(int mi = 0; mi < NCond; mi++){
    f[mi] = TFile::Open(conds[mi].path);
    if(!f[mi] || f[mi]->IsZombie()){ printf("ERROR: cannot open %s\n", conds[mi].path); return; }
  }

  TH1D *hYield[NClass][NCond];
  int col[NCond];
  double ymax[NClass], ymin_pos[NClass];

  for(int mi = 0; mi < NCond; mi++) col[mi] = TColor::GetColor(conds[mi].hex);

  for(int ci = 0; ci < NClass; ci++){
    ymax[ci] = 0.; ymin_pos[ci] = 1e300;

    for(int mi = 0; mi < NCond; mi++){
      double N_evt = sumEvents(f[mi], ci);

      TH1D *raw = sumSlices(f[mi], histBase, ci, Form("raw_%d_%d", ci, mi));
      if(!raw) return;

      hYield[ci][mi] = (TH1D*) raw->Rebin(NEdge-1, Form("hYield_%d_%d", ci, mi), ptEdge);
      hYield[ci][mi]->SetDirectory(nullptr);
      hYield[ci][mi]->Scale(1./N_evt);
      hYield[ci][mi]->Scale(1., "width");

      ymax[ci] = TMath::Max(ymax[ci], hYield[ci][mi]->GetMaximum());
      for(int b = 1; b <= hYield[ci][mi]->GetNbinsX(); b++){
        double v = hYield[ci][mi]->GetBinContent(b);
        if(v > 0. && v < ymin_pos[ci]) ymin_pos[ci] = v;
      }
    }
  }

  // Ratio: reweighted / raw, per class.
  TH1D *hRatio[NClass];
  for(int ci = 0; ci < NClass; ci++){
    hRatio[ci] = (TH1D*) hYield[ci][1]->Clone(Form("hRatio_%d", ci));
    hRatio[ci]->Divide(hYield[ci][0]);
  }

  // Fixed display range -- warn rather than silently clip if a point falls
  // outside it.
  const double ratioDispMin = 0., ratioDispMax = 2.;
  for(int ci = 0; ci < NClass; ci++){
    for(int b = 1; b <= hRatio[ci]->GetNbinsX(); b++){
      double lo = hRatio[ci]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      double v = hRatio[ci]->GetBinContent(b);
      if(v != 0. && (v > ratioDispMax || v < ratioDispMin))
        printf("NOTE: %s reweighted/raw at %.0f-%.0f is %.2f, outside display range [%.1f,%.1f]\n",
               classLabel[ci], lo, hRatio[ci]->GetXaxis()->GetBinUpEdge(b), v, ratioDispMin, ratioDispMax);
    }
  }

  TCanvas *c = new TCanvas("cFastJetHiBinReweightEffect", "", 1300, 1100);

  const double split = 0.35;   // fraction of each quadrant given to the ratio panel
  const double sc = (1. - split) / split;

  for(int ci = 0; ci < NClass; ci++){
    int row = ci / 2, colIdx = ci % 2;
    double x0 = colIdx * 0.5, x1 = (colIdx + 1) * 0.5;
    double y1 = 1. - row * 0.5, y0 = 1. - (row + 1) * 0.5;

    c->cd();
    TPad *outer = new TPad(Form("outer_%d", ci), "", x0, y0, x1, y1);
    outer->SetMargin(0., 0., 0., 0.);
    outer->Draw();
    outer->cd();

    TPad *pUp = new TPad(Form("pUp_%d", ci), "", 0, split, 1, 1);
    pUp->SetLeftMargin(0.24); pUp->SetRightMargin(0.05);
    pUp->SetTopMargin(0.10);  pUp->SetBottomMargin(0.);
    pUp->SetLogy(); pUp->Draw();

    TPad *pDn = new TPad(Form("pDn_%d", ci), "", 0, 0, 1, split);
    pDn->SetLeftMargin(0.24); pDn->SetRightMargin(0.05);
    pDn->SetTopMargin(0.);    pDn->SetBottomMargin(0.34);
    pDn->Draw();

    pUp->cd();
    bool first = true;
    for(int mi = 0; mi < NCond; mi++){
      TH1D *h = hYield[ci][mi];
      h->SetLineColor(col[mi]); h->SetMarkerColor(col[mi]);
      h->SetMarkerStyle(conds[mi].marker); h->SetMarkerSize(0.9);
      h->SetLineWidth(2);
      h->SetTitle("");
      h->GetXaxis()->SetLabelSize(0.); h->GetXaxis()->SetTitleSize(0.);
      h->GetYaxis()->SetTitle("1/N_{evt} dN_{fakeJet}/dp_{T}  [GeV^{-1}]");
      h->GetYaxis()->SetTitleSize(0.062); h->GetYaxis()->SetLabelSize(0.055);
      h->GetYaxis()->SetTitleOffset(1.25);
      h->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
      h->SetMinimum(ymin_pos[ci] * 0.3); h->SetMaximum(ymax[ci] * 5.);
      h->Draw(first ? "ep" : "ep same");
      first = false;
    }

    // Top-left, inside the frame -- align=13 (left,top) anchors the text's
    // own top-left corner at this point instead of its baseline, so it can't
    // stray above/right of the pad margins.
    TLatex clat; clat.SetNDC(); clat.SetTextSize(0.062); clat.SetTextFont(62);
    clat.SetTextAlign(13);
    clat.DrawLatex(0.24 + 0.03, 1. - 0.10 - 0.03, Form("PbPb %s", classLabel[ci]));

    if(ci == 0){
      TLegend *leg = new TLegend(0.42, 0.72, 0.94, 0.86);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.046);
      for(int mi = 0; mi < NCond; mi++)
        leg->AddEntry(hYield[ci][mi], conds[mi].label, "lp");
      leg->Draw();
    }

    pDn->cd();
    hRatio[ci]->SetLineColor(col[1]); hRatio[ci]->SetMarkerColor(col[1]);
    hRatio[ci]->SetMarkerStyle(conds[1].marker); hRatio[ci]->SetMarkerSize(0.9 * sc / 1.6);
    hRatio[ci]->SetLineWidth(2);
    hRatio[ci]->SetTitle("");
    hRatio[ci]->GetXaxis()->SetTitle("raw p_{T}^{fakeJet} [GeV]");
    hRatio[ci]->GetYaxis()->SetTitle("reweighted / raw");
    hRatio[ci]->GetXaxis()->SetTitleSize(0.062 * sc); hRatio[ci]->GetXaxis()->SetLabelSize(0.055 * sc);
    hRatio[ci]->GetYaxis()->SetTitleSize(0.062 * sc); hRatio[ci]->GetYaxis()->SetLabelSize(0.055 * sc);
    hRatio[ci]->GetYaxis()->SetTitleOffset(1.25 / sc);
    hRatio[ci]->GetYaxis()->SetNdivisions(505);
    hRatio[ci]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hRatio[ci]->SetMinimum(ratioDispMin);
    hRatio[ci]->SetMaximum(ratioDispMax);
    hRatio[ci]->Draw("ep");

    TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();
  }

  c->cd();
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.018);
  lat.DrawLatex(0.06, 0.99, "FastJet fake-jet spectrum, raw p_{T}, RC-subtracted, mixed-event -- hiBin reweight effect");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
