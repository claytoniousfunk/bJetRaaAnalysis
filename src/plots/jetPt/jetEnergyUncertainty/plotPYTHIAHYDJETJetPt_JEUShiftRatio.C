// PYTHIA+HYDJET reco jet pT spectrum with the JEU shifted up and down, as
// shifted / nominal. The PbPb counterpart of plotPYTHIAJetPt_JEUShiftRatio.C;
// reco level only, nothing unfolded.
//
// Reads h_inclRecoJetPt_flavor_{nominal,JEUShiftUp,JEUShiftDown}_C<slice> from a
// PYTHIAHYDJET_scan ultraFine output (JEU variants booked 2026-10-02), flavors
// summed -- inclusive reco jets, fakes and all.
//
// Outputs (figures/jetKinematics/jetEnergyUncertainty/)
//   PYTHIAHYDJETJetPt_JEUShiftRatio_<0to10pct|10to30pct|30to50pct|50to80pct>.pdf
//   PYTHIAHYDJETJetPt_JEUShiftRatio_slice<01..16>_<lo>to<hi>pct.pdf
//       spectra (nominal, up, down) in 5 GeV bins on top, shifted / nominal below
//       in wide bins, from 50 GeV; the first set merges the 5% slices of
//       coarseCent.h into each class, the second is one 5% slice per figure.
//   PYTHIAHYDJETJetPt_JEUShift_vsCentrality.pdf
//       shifted / nominal in jet-pT windows against centrality, one point per 5%
//       slice (0-80%); up filled, down open.
//
// The shift moves the jet pT only: the jetPtCut and every filter follow the
// NOMINAL pT, so the lowest bins of the down ratio are not a clean shift of a
// complete spectrum. Uncertainty file: Autumn18_HI_V8_MC_Uncertainty_AK4Calo
// (~5% at |eta|<0.09, 100-250 GeV; the pp file is ~2.3%). Ratio errors are the
// shifted spectrum's statistics only (same jets as nominal, RatioErr::kNumerator).
//
// Usage: root -l -b -q 'plotPYTHIAHYDJETJetPt_JEUShiftRatio.C("<scan.root>")'
// Run from: src/plots/jetPt/jetEnergyUncertainty/

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "../../../../headers/functions/divideByBinwidth.h"
#include "../../../../headers/plotting/plotStyle.h"
#include "../../../../headers/plotting/ratioPanel.h"
#include "../../../../headers/plotting/coarseCent.h"

namespace {

const char *outDir = "../../../../figures/jetKinematics/jetEnergyUncertainty/";

const double ptLo = 50, ptHi = 400;

const char *var[3] = {"nominal", "JEUShiftUp", "JEUShiftDown"};

// 1D pT spectrum (flavors summed) of one summed 2D histogram
// bOnly: |flavor| = 5, i.e. the b (+5) and bbar (-5) bins of the signed flavor axis
TH1D* spectrum(TH2D *H, const char *name, bool bOnly = false)
{
  TH1D *h;
  if(bOnly){
    int bm = H->GetYaxis()->FindBin(-4.5), bp = H->GetYaxis()->FindBin(5.5);
    h = H->ProjectionX(name, bm, bm);
    TH1D *hp = H->ProjectionX(Form("%s_bbar", name), bp, bp);
    h->Add(hp); delete hp;
  }
  else h = H->ProjectionX(name, 1, H->GetNbinsY());
  h->SetDirectory(nullptr);
  return h;
}

// every edge on the 5 GeV grid; the ratio is quoted in these wider bins so
// sparse low-statistics bins do not drive it
const double wideEdges[] = {50, 60, 70, 80, 100, 120, 150, 200, 250, 300, 400};
const int    nWide = sizeof(wideEdges)/sizeof(double) - 1;

TH1D* wide(TH1D *h, const char *name)
{
  TH1D *r = (TH1D*) h->Rebin(nWide, name, wideEdges);   // new histogram
  r->SetDirectory(nullptr);
  return r;
}

// One figure: nominal / up / down spectra on top (5 GeV bins), shifted / nominal
// below in wide bins. H[v] is the summed (pT, flavor) histogram of variation v.
void oneFigure(TH2D *H[3], const char *label, const char *tag, bool bOnly = false, const char *sel = "calo jets, reco")
{
  TH1D *s[3], *w[3];
  for(int v = 0; v < 3; v++){
    s[v] = spectrum(H[v], Form("sp_%s_%s", var[v], tag), bOnly);
    w[v] = wide(s[v], Form("wd_%s_%s", var[v], tag));
  }
  TH1D *rUp = makeRatio(w[1], w[0], Form("rUp_%s", tag));
  TH1D *rDn = makeRatio(w[2], w[0], Form("rDn_%s", tag));
  for(int v = 0; v < 3; v++) divideByBinwidth(s[v]);

  styleH(s[0], okabeHex[0], markFilledCircle, 0.8);
  styleH(s[1], okabeHex[6], markFilledSquare, 0.8);
  styleH(s[2], okabeHex[5], markFilledDiamond, 1.0);
  styleH(rUp,  okabeHex[6], markFilledSquare, 1.0);
  styleH(rDn,  okabeHex[5], markFilledDiamond, 1.2);

  TCanvas *c = new TCanvas(Form("c_%s", tag), "", 700, 800);
  TPad *top, *bot; splitPads(top, bot);

  top->cd(); top->SetLogy();   // steeply falling spectrum; no negative bins to hide here
  s[0]->GetXaxis()->SetRangeUser(ptLo, ptHi);
  s[0]->SetTitle("");
  s[0]->GetYaxis()->SetTitle("dN/dp_{T}^{reco} (arb.)");
  s[0]->GetYaxis()->SetTitleSize(0.055); s[0]->GetYaxis()->SetLabelSize(0.050);
  s[0]->GetXaxis()->SetLabelSize(0);
  s[0]->Draw("E1"); s[1]->Draw("E1 SAME"); s[2]->Draw("E1 SAME");
  TLegend *leg = makeLegend(0.55, 0.66, 0.90, 0.88, 0.048);
  leg->AddEntry(s[0], "Nominal", "lp");
  leg->AddEntry(s[1], "JEU up", "lp");
  leg->AddEntry(s[2], "JEU down", "lp");
  leg->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.045);
  t.DrawLatex(0.20, 0.17, Form("PYTHIA+HYDJET %s", label));
  t.SetTextSize(0.038);
  t.DrawLatex(0.20, 0.11, sel);

  bot->cd();
  rUp->GetXaxis()->SetRangeUser(ptLo, ptHi);
  styleRatioAxes(rUp, "Reco jet p_{T} (GeV)", "Shifted / nominal");
  rUp->GetYaxis()->SetRangeUser(0.5, 1.8);   // central slices leave 0.6-1.5 at low pT
  rUp->Draw("E1"); rDn->Draw("E1 SAME");
  unityLine(ptLo, ptHi)->Draw();

  gSystem->mkdir(outDir, kTRUE);
  savePdfTight(c, Form("%sPYTHIAHYDJETJetPt_JEUShiftRatio_%s.pdf", outDir, tag));
  delete c;
  for(int v = 0; v < 3; v++){ delete w[v]; }
  delete rUp; delete rDn;
}

void perClass(TFile *f, int ci, const char *base = "h_inclRecoJetPt_flavor", bool bOnly = false,
              const char *sel = "calo jets, reco", const char *selTag = "")
{
  TH2D *H[3];
  for(int v = 0; v < 3; v++){
    H[v] = coarseSum2(f, Form("%s_%s", base, var[v]), ci, Form("pc%d", v));
    if(!H[v]){ printf("missing %s_%s slices for %s -- scan predates the JEU histograms?\n", base, var[v], coarseLabel[ci]); return; }
  }
  oneFigure(H, coarseLabel[ci], Form("%s%s", selTag, coarseTag[ci]), bOnly, sel);
  for(int v = 0; v < 3; v++) delete H[v];
}

// one 5% slice (1..16 -> 0-5% ... 75-80%)
void perSlice(TFile *f, int si)
{
  TH2D *H[3];
  for(int v = 0; v < 3; v++){
    f->GetObject(Form("h_inclRecoJetPt_flavor_%s_C%d", var[v], si), H[v]);
    if(!H[v]){ printf("missing h_inclRecoJetPt_flavor_%s_C%d -- scan predates the JEU histograms?\n", var[v], si); return; }
  }
  const int lo = 5*(si - 1), hi = 5*si;
  oneFigure(H, Form("%d-%d%%", lo, hi), Form("slice%02d_%dto%dpct", si, lo, hi));
}

// counts in [lo, hi) on the 5 GeV grid; FindBin(hi) starts at hi, so use hi - eps
double window(TH1D *h, double lo, double hi, double *err = nullptr)
{
  int b1 = h->GetXaxis()->FindBin(lo + 1e-6), b2 = h->GetXaxis()->FindBin(hi - 1e-6);
  double e = 0.;
  double v = h->IntegralAndError(b1, b2, e);
  if(err) *err = e;
  return v;
}

void vsCentrality(TFile *f, const char *base = "h_inclRecoJetPt_flavor", bool bOnly = false,
                  const char *sel = "reco, inclusive", const char *selTag = "")
{
  const int nWin = 4;
  const double winLo[nWin] = { 50, 100, 150, 250};
  const double winHi[nWin] = {100, 150, 250, 400};
  const char  *hex[nWin]   = {okabeHex[0], okabeHex[1], okabeHex[5], okabeHex[3]};   // black, orange, blue, green
  const int    mkF[nWin]   = {markFilledCircle, markFilledSquare, markFilledDiamond, markCross};
  const int    mkO[nWin]   = {markOpenCircle,   markOpenSquare,   markOpenDiamond,   markCross};
  TGraphErrors *g[nWin][2];
  for(int w = 0; w < nWin; w++) for(int s = 0; s < 2; s++) g[w][s] = new TGraphErrors();

  printf("\nslice  cent  | window GeV | up     down   (reco, shifted / nominal)\n");
  for(int si = 1; si <= 16; si++){
    TH1D *h[3];
    for(int v = 0; v < 3; v++){
      TH2D *H = nullptr; f->GetObject(Form("%s_%s_C%d", base, var[v], si), H);
      if(!H){ printf("missing %s_%s_C%d -- scan predates the JEU histograms?\n", base, var[v], si); return; }
      h[v] = spectrum(H, Form("vc_%s_%d", var[v], si), bOnly);
    }
    const double cent = 5.*(si - 1) + 2.5;
    for(int w = 0; w < nWin; w++){
      double nom = window(h[0], winLo[w], winHi[w]), r[2];
      for(int s = 0; s < 2; s++){
        double e, sh = window(h[s+1], winLo[w], winHi[w], &e);
        r[s] = nom > 0 ? sh / nom : 0.;
        int np = g[w][s]->GetN();
        g[w][s]->SetPoint(np, cent, r[s]);
        g[w][s]->SetPointError(np, 0., nom > 0 ? e / nom : 0.);
      }
      printf("%4d  %4.1f%% | %3.0f-%-3.0f    | %.3f  %.3f\n", si, cent, winLo[w], winHi[w], r[0], r[1]);
    }
    for(int v = 0; v < 3; v++) delete h[v];
  }

  TCanvas *c = new TCanvas("c_vc", "", 800, 700);
  c->SetLeftMargin(0.14); c->SetBottomMargin(0.13); c->SetTopMargin(0.06); c->SetRightMargin(0.04);
  TH1F *fr0 = c->DrawFrame(0, 0.6, 80, 1.8);
  fr0->GetXaxis()->SetTitle("Centrality (%)"); fr0->GetYaxis()->SetTitle("Shifted / nominal, reco");
  fr0->GetXaxis()->SetTitleSize(0.05); fr0->GetYaxis()->SetTitleSize(0.05);
  fr0->GetXaxis()->SetLabelSize(0.045); fr0->GetYaxis()->SetLabelSize(0.045);
  fr0->GetYaxis()->SetTitleOffset(1.35);
  TLine one; one.SetLineStyle(7); one.SetLineColor(kGray+1); one.DrawLine(0, 1., 80, 1.);
  TLegend *leg = makeLegend(0.22, 0.405, 0.95, 0.53, 0.037);
  leg->SetNColumns(2);
  for(int w = 0; w < nWin; w++) for(int s = 0; s < 2; s++){
    int col = TColor::GetColor(hex[w]);
    g[w][s]->SetMarkerColor(col); g[w][s]->SetLineColor(col);
    g[w][s]->SetMarkerStyle(s == 0 ? mkF[w] : mkO[w]); g[w][s]->SetMarkerSize(1.2); g[w][s]->SetLineWidth(2);
    g[w][s]->Draw("PE same");
    if(s == 0) leg->AddEntry(g[w][s], Form("%.0f-%.0f GeV", winLo[w], winHi[w]), "p");
  }
  leg->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.045);
  t.DrawLatex(0.18, 0.88, "PYTHIA+HYDJET calo jets");
  t.SetTextSize(0.038);
  t.DrawLatex(0.18, 0.83, Form("%s; filled: JEU up, open: JEU down", sel));

  gSystem->mkdir(outDir, kTRUE);
  savePdfTight(c, Form("%sPYTHIAHYDJETJetPt_JEUShift_vsCentrality%s.pdf", outDir, selTag));
  delete c;
}

} // namespace

void plotPYTHIAHYDJETJetPt_JEUShiftRatio(const char *scanFile)
{
  initPlotStyle();
  TFile *f = TFile::Open(scanFile);
  if(!f || f->IsZombie()){ printf("cannot open %s\n", scanFile); return; }
  for(int ci = 0; ci < NCoarse; ci++) perClass(f, ci);
  for(int si = 1; si <= 16; si++) perSlice(f, si);
  vsCentrality(f);

  // muon-tagged jets (incl. reco muon tag, no trigger gating), all flavors and b only
  const char *mu = "h_inclRecoJetPt_inclRecoMuonTag_flavor";
  for(int ci = 0; ci < NCoarse; ci++){
    perClass(f, ci, mu, false, "calo jets, #mu-tagged", "muTagged_allJets_");
    perClass(f, ci, mu, true,  "calo jets, #mu-tagged b", "muTagged_bJets_");
  }
  // no vs-centrality version for the muon-tagged jets: 5% slices are too sparse
}
