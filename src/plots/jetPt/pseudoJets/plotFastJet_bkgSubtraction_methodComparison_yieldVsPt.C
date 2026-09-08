// Direct method-comparison plot: how does the choice of fake-jet background-
// subtraction method (RC, RC_geoCorr, dPT) change the actual per-event
// fake-jet dN/dpT spectrum, holding the event sample fixed?
//
// Everything same-event, including RC_sigSel (restricted to events passing
// the signal jet-pT selection -- see doSignalSelectedRC in pseudoJets.h,
// only ever computed same-event). A uniform sample means RC_sigSel/RC here
// isolates just the signal-selection effect, unlike comparing it against
// mixed-event RC (which also carries the same/mixed contamination gap --
// see plotFastJet_bkgSubtraction_mixedVsSame.C). One panel per coarse
// centrality class, curves overlaid per panel plus a method/RC ratio
// underneath.
//
// Usage: root -l -b -q 'plotFastJet_bkgSubtraction_methodComparison_yieldVsPt.C'
// Run from: src/plots/jetPt/pseudoJets/

const char *sameFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-18_ultraFineCentBins.root";

const char *outDir  = "../../../../figures/jetPt/pseudoJets/";
const char *outName = "fastJet_bkgSub_methodComparison_yieldVsPt.pdf";

const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const double plotPtMin = 0., plotPtMax = 160.;

const int    NEdge = 25;
double       ptEdge[NEdge] = {
  0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,110,120,130,140,150,160
};

struct Method { const char *suffix; const char *label; const char *hex; int marker; };
const int NMethod = 4;
Method methods[NMethod] = {
  {"RC",                 "RC",                  "#0072B2", 20},
  {"RC_geoCorr",         "RC_{geoCorr}",         "#009E73", 21},
  {"dPT",                "dPT",                  "#D55E00", 33},
  {"RC_sigSel",          "RC (sigSel)",          "#CC79A7", 34},
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

void plotFastJet_bkgSubtraction_methodComparison_yieldVsPt()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fS = TFile::Open(sameFile);
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open %s\n", sameFile); return; }

  TH1D *hYield[NClass][NMethod];
  int col[NMethod];
  double ymax[NClass], ymin_pos[NClass];

  for(int mi = 0; mi < NMethod; mi++) col[mi] = TColor::GetColor(methods[mi].hex);

  for(int ci = 0; ci < NClass; ci++){
    ymax[ci] = 0.; ymin_pos[ci] = 1e300;

    for(int mi = 0; mi < NMethod; mi++){
      double N_evt = sumEvents(fS, ci);

      TString histBase = TString("h_fastJetPt_PF_bkgSub_") + methods[mi].suffix;
      TH1D *raw = sumSlices(fS, histBase.Data(), ci, Form("raw_%d_%s", ci, methods[mi].suffix));
      if(!raw) return;

      hYield[ci][mi] = (TH1D*) raw->Rebin(NEdge-1, Form("hYield_%d_%s", ci, methods[mi].suffix), ptEdge);
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

  // Ratio of each non-reference method to RC (index 0), per class. Built from
  // the same rebinned/normalized hYield histograms, so width/N_evt scaling
  // cancels cleanly in the division.
  TH1D *hRatio[NClass][NMethod];
  double ratioMin[NClass], ratioMax[NClass];
  for(int ci = 0; ci < NClass; ci++){
    ratioMin[ci] = 1e300; ratioMax[ci] = -1e300;
    for(int mi = 1; mi < NMethod; mi++){
      hRatio[ci][mi] = (TH1D*) hYield[ci][mi]->Clone(Form("hRatio_%d_%s", ci, methods[mi].suffix));
      hRatio[ci][mi]->Divide(hYield[ci][0]);
      for(int b = 1; b <= hRatio[ci][mi]->GetNbinsX(); b++){
        double lo = hRatio[ci][mi]->GetXaxis()->GetBinLowEdge(b);
        if(lo < plotPtMin || lo >= plotPtMax) continue;
        double v = hRatio[ci][mi]->GetBinContent(b);
        if(v <= 0.) continue;
        if(v < ratioMin[ci]) ratioMin[ci] = v;
        if(v > ratioMax[ci]) ratioMax[ci] = v;
      }
    }
    if(ratioMax[ci] / ratioMin[ci] > 50.)
      printf("NOTE: %s method/RC ratio spans %.2f-%.2f (>50x) -- likely low-stats tail bins\n",
             classLabel[ci], ratioMin[ci], ratioMax[ci]);
  }

  TCanvas *c = new TCanvas("cFastJetMethodComparison", "", 1300, 1100);

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
    for(int mi = 0; mi < NMethod; mi++){
      TH1D *h = hYield[ci][mi];
      h->SetLineColor(col[mi]); h->SetMarkerColor(col[mi]);
      h->SetMarkerStyle(methods[mi].marker); h->SetMarkerSize(0.9);
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
      TLegend *leg = new TLegend(0.50, 0.63, 0.94, 0.86);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.046);
      for(int mi = 0; mi < NMethod; mi++)
        leg->AddEntry(hYield[ci][mi], methods[mi].label, "lp");
      leg->Draw();
    }

    pDn->cd();
    bool firstR = true;
    for(int mi = 1; mi < NMethod; mi++){
      TH1D *r = hRatio[ci][mi];
      r->SetLineColor(col[mi]); r->SetMarkerColor(col[mi]);
      r->SetMarkerStyle(methods[mi].marker); r->SetMarkerSize(0.9 * sc / 1.6);
      r->SetLineWidth(2);
      r->SetTitle("");
      r->GetXaxis()->SetTitle("raw p_{T}^{fakeJet} [GeV]");
      r->GetYaxis()->SetTitle("method / RC");
      r->GetXaxis()->SetTitleSize(0.062 * sc); r->GetXaxis()->SetLabelSize(0.055 * sc);
      r->GetYaxis()->SetTitleSize(0.062 * sc); r->GetYaxis()->SetLabelSize(0.055 * sc);
      r->GetYaxis()->SetTitleOffset(1.25 / sc);
      r->GetYaxis()->SetNdivisions(505);
      r->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
      r->SetMinimum(TMath::Max(0., ratioMin[ci] * 0.8));
      r->SetMaximum(ratioMax[ci] * 1.2);
      r->Draw(firstR ? "ep" : "ep same");
      firstR = false;
    }
    TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();
  }

  c->cd();
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.018);
  lat.DrawLatex(0.06, 0.99, "FastJet fake-jet spectrum, raw p_{T}, same-event -- background-subtraction method comparison");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
