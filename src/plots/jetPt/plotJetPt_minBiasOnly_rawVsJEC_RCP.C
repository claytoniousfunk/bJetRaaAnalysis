// Task #1 from the 2026-08-18 meeting: does the JEC break the centrality
// hierarchy in 10-30%?
//
// Olga: "make that jetsPerZ plot by rawPt.  Just in case you are messing up
// the corrections in 10-30.  If at rawPt there is no issue with hierarchy,
// then maybe it will point that it is corrections... These two do not
// commute.  In data, your hierarchy is normal.  If you do this plot, with
// JEC, and see the hierarchy is messed up, we will know it's the JEC."
//
// *** WHY THIS IS NOT THE LITERAL PbPb/pp JETS-PER-Z AT RAW pT ***
// h_inclRawJetPt is only produced by the CURRENT pp_scan.C / PbPb_scan.C /
// PbPb_pfCandAnalyzer.C.  Checked every scan output in both repos: it exists
// ONLY in the 16 PbPb MinBias pfCandAnalyzer files.  No pp file has it (0 of
// 58 checked), and no PbPb Jet60/80/100 triggered file has it (0 of 75).  So
// neither the pp reference nor the stitched high-pT region can be built at
// raw pT without rescanning those samples -- see the note at the bottom.
//
// This plot answers the hierarchy question anyway, PbPb-only, by replacing
// the common pp denominator with the most peripheral PbPb class.  That is
// legitimate for THIS question because pp is a single common denominator
// shared by all four centrality classes: dividing every class by the same
// pp(pT) rescales all four curves identically at each pT and therefore
// cannot reorder them.  The centrality ORDERING -- which is the whole
// question -- lives entirely in the PbPb numerators.
//
// So: R_CP-style ratio, each class divided by 50-80%, computed twice from
// the SAME jets under the two pT definitions (identical entry counts,
// 3109553 in C1 for both -- same jets, different pT assignment).
//
// NORMALIZATION: each class is independently scaled so its own integral
// above ptNormMin = 200 GeV equals 1, exactly as the slide-2 PbPb/pp plot
// does (plotJetPt_minBiasOnly_PbPbOverPP.C).  This matters -- a per-event
// normalization instead makes the ratio dominated by the raw jet-yield
// difference between centralities (0-10% sits ~235x above 50-80% at 50 GeV,
// which swamps everything), whereas the shape normalization is what puts the
// "peak height" on the scale Olga reads off the slide: 0-10% ~8, 10-30% ~8
// (should be ~6), 30-50% ~4.
//
// No Ncoll scaling, so this is not R_CP in the suppression-measurement sense.
// Deliberate: we are testing shape hierarchy, and at these pT the spectrum is
// fake-dominated, where fakes scale per event rather than with binary
// collisions.
//
// Usage: root -l -b -q 'plotJetPt_minBiasOnly_rawVsJEC_RCP.C'
// Run from: src/plots/jetPt/

const char *pbpbFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-18_ultraFineCentBins.root";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetPt_minBiasOnly_rawVsJEC_RCP.pdf";

const int   NClass = 4;
const int   sliceLo[NClass] = { 1,  3,  7, 11};
const int   sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };
const char *classHex[NClass]   = { "#D55E00", "#0072B2", "#009E73", "#CC79A7" };
const int   classMark[NClass]  = { 20, 21, 33, 34 };
const double classSize[NClass] = { 1.0, 1.0, 1.2, 1.1 };

// Index of the class used as the denominator (50-80%, most peripheral).
const int REF = NClass - 1;

// Two pT definitions, read from the same file, same events, same jets.
const int   NDef = 2;
const char *defHist[NDef]  = { "h_inclRecoJetPt", "h_inclRawJetPt" };
const char *defLabel[NDef] = { "JEC-corrected p_{T}", "raw p_{T}" };

const double plotPtMin = 30., plotPtMax = 500.;

// Each class scaled so its own integral above this equals 1 -- same scheme as
// plotJetPt_minBiasOnly_PbPbOverPP.C, so peak heights are on the same scale.
const double ptNormMin = 200.;

const int    NEdge = 29;
double       ptEdge[NEdge] = {
  20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,
  110,120,130,140,150,160,180,200,240,280,350,500
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

void plotJetPt_minBiasOnly_rawVsJEC_RCP()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(pbpbFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", pbpbFile); return; }

  TH1D *hY[NDef][NClass], *hR[NDef][NClass];
  int col[NClass];
  for(int ci = 0; ci < NClass; ci++) col[ci] = TColor::GetColor(classHex[ci]);

  double ymax = 0., ymin_pos = 1e300;
  for(int di = 0; di < NDef; di++){
    for(int ci = 0; ci < NClass; ci++){
      TH1D *raw = sumSlices(f, defHist[di], ci, Form("raw_d%d_c%d", di, ci));
      if(!raw) return;

      hY[di][ci] = (TH1D*) raw->Rebin(NEdge-1, Form("hY_d%d_c%d", di, ci), ptEdge);
      hY[di][ci]->SetDirectory(nullptr);
      int b0 = hY[di][ci]->FindBin(ptNormMin);
      double norm = hY[di][ci]->Integral(b0, hY[di][ci]->GetNbinsX()+1);
      if(norm <= 0.){ printf("ERROR: %s class %s has zero integral above %.0f\n", defHist[di], classLabel[ci], ptNormMin); return; }
      hY[di][ci]->Scale(1./norm);
      hY[di][ci]->Scale(1., "width");

      for(int b = 1; b <= hY[di][ci]->GetNbinsX(); b++){
        double lo = hY[di][ci]->GetXaxis()->GetBinLowEdge(b);
        if(lo < plotPtMin || lo >= plotPtMax) continue;
        double v = hY[di][ci]->GetBinContent(b);
        if(v <= 0.) continue;
        if(v > ymax) ymax = v;
        if(v < ymin_pos) ymin_pos = v;
      }
    }
  }

  // R_CP-style: each class over the most peripheral class, same pT definition.
  double rmax = 0.;
  for(int di = 0; di < NDef; di++){
    for(int ci = 0; ci < NClass; ci++){
      hR[di][ci] = (TH1D*) hY[di][ci]->Clone(Form("hR_d%d_c%d", di, ci));
      hR[di][ci]->Divide(hY[di][REF]);
      if(ci == REF) continue;
      for(int b = 1; b <= hR[di][ci]->GetNbinsX(); b++){
        double lo = hR[di][ci]->GetXaxis()->GetBinLowEdge(b);
        if(lo < plotPtMin || lo >= plotPtMax) continue;
        double v = hR[di][ci]->GetBinContent(b);
        if(v > rmax) rmax = v;
      }
    }
  }
  const double rMin = 0., rMax = rmax * 1.15;

  // Peak of each class's ratio -- the number Olga reads off the plot. Search
  // is restricted to the low-pT bump region; over the full range the weak
  // 30-50% bump loses to high-pT statistical noise and the reported "peak"
  // lands on a fluctuation near 200-350 GeV instead.
  const double peakSearchMax = 120.;
  printf("\nPeak of (class / %s) below %.0f GeV, per pT definition:\n", classLabel[REF], peakSearchMax);
  printf("%-22s %12s %12s %12s\n", "pT definition", classLabel[0], classLabel[1], classLabel[2]);
  for(int di = 0; di < NDef; di++){
    printf("%-22s", defLabel[di]);
    for(int ci = 0; ci < REF; ci++){
      double best = 0., bestPt = 0.;
      for(int b = 1; b <= hR[di][ci]->GetNbinsX(); b++){
        double lo = hR[di][ci]->GetXaxis()->GetBinLowEdge(b);
        if(lo < plotPtMin || lo >= peakSearchMax) continue;
        double v = hR[di][ci]->GetBinContent(b);
        if(v > best){ best = v; bestPt = hR[di][ci]->GetXaxis()->GetBinCenter(b); }
      }
      printf("  %6.2f@%-4.0f", best, bestPt);
    }
    printf("\n");
  }

  // Value at 50 GeV specifically -- Olga: "at 50, this is horrendous".
  printf("\n(class / %s) at 50 GeV:\n", classLabel[REF]);
  printf("%-22s %12s %12s %12s\n", "pT definition", classLabel[0], classLabel[1], classLabel[2]);
  for(int di = 0; di < NDef; di++){
    printf("%-22s", defLabel[di]);
    int b50 = hR[di][0]->GetXaxis()->FindBin(52.5);
    for(int ci = 0; ci < REF; ci++) printf(" %12.3f", hR[di][ci]->GetBinContent(b50));
    printf("\n");
  }

  TCanvas *c = new TCanvas("cRawVsJEC_RCP", "", 1300, 1000);

  const double split = 0.42;
  const double sc = (1. - split) / split;

  for(int di = 0; di < NDef; di++){
    double x0 = di * 0.5, x1 = (di + 1) * 0.5;

    c->cd();
    TPad *outer = new TPad(Form("outer_%d", di), "", x0, 0., x1, 1.);
    outer->SetMargin(0., 0., 0., 0.);
    outer->Draw();
    outer->cd();

    TPad *pUp = new TPad(Form("pUp_%d", di), "", 0, split, 1, 1);
    pUp->SetLeftMargin(0.20); pUp->SetRightMargin(0.05);
    pUp->SetTopMargin(0.11);  pUp->SetBottomMargin(0.);
    pUp->SetLogy(); pUp->Draw();

    TPad *pDn = new TPad(Form("pDn_%d", di), "", 0, 0, 1, split);
    pDn->SetLeftMargin(0.20); pDn->SetRightMargin(0.05);
    pDn->SetTopMargin(0.);    pDn->SetBottomMargin(0.26);
    pDn->Draw();

    pUp->cd();
    for(int ci = 0; ci < NClass; ci++){
      TH1D *h = hY[di][ci];
      h->SetLineColor(col[ci]); h->SetMarkerColor(col[ci]);
      h->SetMarkerStyle(classMark[ci]); h->SetMarkerSize(classSize[ci] * 0.8);
      h->SetLineWidth(2);
      h->SetTitle("");
      h->GetXaxis()->SetLabelSize(0.); h->GetXaxis()->SetTitleSize(0.);
      h->GetYaxis()->SetTitle("1/N_{jet}(p_{T}>200) dN_{jet}/dp_{T}  [GeV^{-1}]");
      h->GetYaxis()->SetTitleSize(0.055); h->GetYaxis()->SetLabelSize(0.048);
      h->GetYaxis()->SetTitleOffset(1.55);
      h->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
      h->SetMinimum(ymin_pos * 0.3); h->SetMaximum(ymax * 5.);
      h->Draw(ci == 0 ? "ep" : "ep same");
    }

    TLatex dlat; dlat.SetNDC(); dlat.SetTextSize(0.052); dlat.SetTextFont(62);
    dlat.SetTextAlign(13);
    dlat.DrawLatex(0.23, 0.87, defLabel[di]);

    if(di == 0){
      TLegend *leg = new TLegend(0.52, 0.46, 0.94, 0.76);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.044);
      for(int ci = 0; ci < NClass; ci++)
        leg->AddEntry(hY[di][ci], Form("PbPb %s", classLabel[ci]), "lp");
      leg->Draw();
    }

    pDn->cd();
    bool first = true;
    for(int ci = 0; ci < REF; ci++){
      TH1D *r = hR[di][ci];
      r->SetLineColor(col[ci]); r->SetMarkerColor(col[ci]);
      r->SetMarkerStyle(classMark[ci]); r->SetMarkerSize(classSize[ci] * 0.8 * sc / 1.4);
      r->SetLineWidth(2);
      r->SetTitle("");
      r->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
      r->GetYaxis()->SetTitle(Form("class / %s", classLabel[REF]));
      r->GetXaxis()->SetTitleSize(0.055*sc); r->GetXaxis()->SetLabelSize(0.048*sc);
      r->GetYaxis()->SetTitleSize(0.055*sc); r->GetYaxis()->SetLabelSize(0.048*sc);
      r->GetYaxis()->SetTitleOffset(1.55/sc);
      r->GetYaxis()->SetNdivisions(505);
      r->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
      r->SetMinimum(rMin); r->SetMaximum(rMax);
      r->Draw(first ? "ep" : "ep same");
      first = false;
    }
    TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();
  }

  c->cd();
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.021);
  lat.DrawLatex(0.05, 0.985, "PbPb MinBias, no subtraction, each class shape-normalized above 200 GeV -- centrality hierarchy at raw vs JEC-corrected p_{T}");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
