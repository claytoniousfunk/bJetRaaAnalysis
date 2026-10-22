// PbPb data all-track ptRel against PYTHIA+HYDJET light-jet muon ptRel: where the
// data-driven light template differs from the light muons it stands in for, and
// whether the central excess follows the share of combinatorial jets. 2026-10-20,
// user. Includes studyDataDrivenLightAllTrk_PF_PbPb.C (its inputs and helpers).
// Per coarse class (data and P+H summed over the class's 5% slices):
//   allTrkVsMu_meanPtRel_PbPb.pdf
//     <ptRel> over 0-3 GeV vs jet pT: data all tracks reweighted to the 0-80% P+H
//     light-muon pT (the study's template), P+H light-jet muons of the class (reco
//     level, no extra smearing), and the nominal light template (pp light smeared by
//     the class's fitted axis sigma). Errors: sqrt(var / Neff) over 0-3 GeV. Bottom:
//     data - P+H light muon.
//   allTrkVsMu_shapes_PbPb.pdf
//     the ptRel shapes over jet pT 80-150 GeV merged (P+H light muons have 6-77
//     effective entries per central window alone): data all tracks reweighted to the
//     0-80% target in that range, P+H light muons of the class; fit binning, unit
//     area. Bottom: data / P+H.
//   allTrkVsMu_excessVsTrackFrac_PbPb.pdf
//     left: f15, the fraction of jets (|eta| < 1.6) with a selected track above
//     15 GeV within dR < 0.4, vs jet pT, per class -- low where most jets are
//     upward fluctuations or combinatorial (central MinBias). Right: the data - P+H
//     <ptRel> excess of each class and window against its f15. An anticorrelation
//     SUPPORTS (does not prove) that the central excess comes from tracks in
//     combinatorial / fluctuated jets, which the P+H light muons (gen-matched jets,
//     removeUnmatchedJets) do not contain. The proof needs the MC closure.
// Output: figures/bPurity/PFJets_PbPb/dataDrivenLight/
// Usage, from src/plots/bPurity/:  root -l -b -q plotAllTrkVsLightMuon_PF_PbPb.C
#include "studyDataDrivenLightAllTrk_PF_PbPb.C"

namespace avm {
// mean and its error over 0 - hi GeV, from contents and errors (Neff = sum^2 / sum err^2)
void meanErr(TH1 *h, double &m, double &e, double hi = 3.){
  double s = 0, sx = 0, sxx = 0, se2 = 0;
  for(int i = 1; i <= h->FindBin(hi - 1e-6); i++){ const double x = h->GetBinCenter(i), w = h->GetBinContent(i); s += w; sx += w*x; sxx += w*x*x; se2 += pow(h->GetBinError(i), 2); }
  if(s <= 0){ m = 0; e = 0; return; }
  m = sx/s; const double var = sxx/s - m*m, ne = se2 > 0 ? s*s/se2 : 0.; e = ne > 0 ? sqrt(TMath::Max(var, 0.)/ne) : 0.;
}
TGraphErrors* graph(int n, const double *x, const double *ex, const double *y, const double *ey, const char *hex, int mark, double size = 1.2){
  TGraphErrors *g = new TGraphErrors(n, x, y, ex, ey); const int col = TColor::GetColor(hex);
  g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mark); g->SetMarkerSize(size); g->SetLineWidth(2); return g;
}
TPad* pane(TCanvas *cv, int c, int n, bool top, double split, double yTop = 0.89){
  const double x0 = 0.05, w = (1. - x0 - 0.01)/n, xl = x0 + c*w - (c == 0 ? x0 : 0.), xh = x0 + (c+1)*w;
  cv->cd(); TPad *p = new TPad(Form("avm%d%d_%d", c, top, dla::uid++), "", xl, top ? split : 0., xh, top ? yTop : split);
  p->SetLeftMargin(c == 0 ? x0/(xh - xl) + 0.12 : 0.02); p->SetRightMargin(0.02);
  if(top){ p->SetTopMargin(0.05); p->SetBottomMargin(0.015); } else { p->SetTopMargin(0.02); p->SetBottomMargin(0.30); }
  p->Draw(); p->cd(); return p;
}
void text(double x, double y, const char *s, double size){ TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size); t.SetTextAlign(13); t.DrawLatex(x, y, s); }
}
using namespace avm;

void plotAllTrkVsLightMuon_PF_PbPb(){
  setJetCollection("PFJets"); RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL); initPlotStyle(); gStyle->SetEndErrorSize(0);
  TFile *fT = TFile::Open(trkFile), *fL = TFile::Open(phLightFile); PHFiles ph;
  if(!fT || !fL || !openPH(ph)){ printf("ERROR: inputs\n"); return; }
  gSystem->mkdir(outDir, kTRUE);
  auto lightMu = [&](int s0, int s1){ TH3 *m = nullptr; for(const char *fl : {"u", "d", "s", "g"}){ TH3 *h = sumSlices(fL, Form("h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_%sJets", fl), s0, s1); if(!h) return (TH3*) nullptr; if(!m) m = h; else { m->Add(h); delete h; } } return m; };
  TH3 *m80 = lightMu(1, 16);
  if(!m80){ printf("ERROR: P+H light-muon 3D missing\n"); return; }

  double x[nWin], ex[nWin];
  for(int k = 0; k < nWin; k++){ x[k] = 0.5*(win[k][0] + win[k][1]); ex[k] = 0.5*(win[k][1] - win[k][0]); }
  double mD[NCls][nWin], eD[NCls][nWin], mM[NCls][nWin], eM[NCls][nWin], mN[NCls][nWin], f15[NCls][nWin], ef15[NCls][nWin], dif[NCls][nWin], edif[NCls][nWin];
  TH1D *shD[NCls], *shM[NCls];
  printf("class   window   data rw        P+H light mu   nominal   data-P+H        f15\n");
  for(int c = 0; c < nCls; c++){
    const int s0 = clsSlice[c][0], s1 = clsSlice[c][1];
    TH3 *tA = sumSlices(fT, "h_jetTrkPtRel_trkPt_recoJetPt", s0, s1), *tL = sumSlices(fT, "h_jetTrkMaxPtRel_trkPt_recoJetPt", s0, s1), *mC = lightMu(s0, s1);
    TH1 *jets = nullptr; for(int i = s0; i <= s1; i++){ TH1 *h = (TH1*) fT->Get(Form("h_inclRecoJetPt_C%d", i)); if(!jets){ jets = (TH1*) h->Clone(Form("jets%d", c)); jets->SetDirectory(nullptr); } else jets->Add(h); }
    if(!tA || !tL || !mC || !jets){ printf("ERROR: %s inputs\n", clsName[c]); return; }
    for(int k = 0; k < nWin; k++){
      const double lo = win[k][0], hi = win[k][1];
      TH2D *T = proj2(tA, lo, hi), *M = proj2(mC, lo, hi), *M80 = proj2(m80, lo, hi), *TL = proj2(tL, lo, hi);
      TH1D *tg = M80->ProjectionY(Form("tg%d", uid++)), *La = reweighted(T, tg), *Mx = M->ProjectionX(Form("mx%d", uid++));
      Tpl t = buildTemplatesPbPb(ph, c, lo, hi);
      meanErr(La, mD[c][k], eD[c][k]); meanErr(Mx, mM[c][k], eM[c][k]); double en; meanErr(t.l, mN[c][k], en);
      dif[c][k] = mD[c][k] - mM[c][k]; edif[c][k] = hypot(eD[c][k], eM[c][k]);
      // f15: jets with a leading selected track above 15 GeV (all ptRel, overflow included) / all jets
      const double nLead = TL->Integral(0, TL->GetNbinsX()+1, TL->GetYaxis()->FindBin(15.001), TL->GetNbinsY()+1);
      const double nJet = jets->Integral(jets->FindBin(lo + 1e-6), jets->FindBin(hi - 1e-6));
      f15[c][k] = nJet > 0 ? nLead/nJet : 0.; ef15[c][k] = nJet > 0 ? sqrt(f15[c][k]*(1. - f15[c][k])/nJet) : 0.;
      printf("%-7s %3.0f-%-3.0f  %.3f+-%.3f   %.3f+-%.3f    %.3f    %+.3f+-%.3f   %.3f\n", clsName[c], lo, hi, mD[c][k], eD[c][k], mM[c][k], eM[c][k], mN[c][k], dif[c][k], edif[c][k], f15[c][k]);
      delete T; delete M; delete M80; delete TL; delete tg; delete La; delete Mx;
    }
    // shapes over 80-150 GeV merged
    { TH2D *T = proj2(tA, 80, 150), *M = proj2(mC, 80, 150), *M80 = proj2(m80, 80, 150);
      TH1D *tg = M80->ProjectionY(Form("tg%d", uid++)), *La = reweighted(T, tg), *Mx = M->ProjectionX(Form("mx%d", uid++));
      shD[c] = range05(La, Form("shD%d", c)); shM[c] = range05(Mx, Form("shM%d", c));
      printf("%-7s 80-150 shapes: P+H light-muon Neff %.0f, data tracks %.0f\n", clsName[c], neff(Mx), T->Integral());
      for(TH1D *h : {shD[c], shM[c]}){ h->Scale(1./h->Integral()); h->Scale(1., "width"); }
      delete T; delete M; delete M80; delete tg; delete La; delete Mx; }
    delete tA; delete tL; delete mC; delete jets;
  }

  const char *hdr = "PbPb 5.02 TeV, anti-k_{T} R = 0.4 PF jets; data tracks: charged PF candidates, p_{T} > 15 GeV, #DeltaR < 0.4";
  const double split = 0.36;
  // ---------------- figure 1: <ptRel> vs jet pT
  {
    TCanvas *cv = new TCanvas("cvM", "", 1600, 600);
    for(int c = 0; c < nCls; c++){
      pane(cv, c, nCls, true, split);
      TH1D *fr = new TH1D(Form("frM%d", c), "", 1, 80, 300); fr->SetStats(0); fr->SetMinimum(0.6); fr->SetMaximum(1.45);
      fr->GetYaxis()->SetTitle(c == 0 ? "#LTp_{T}^{rel}#GT_{0-3 GeV} [GeV]" : ""); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetTitleOffset(1.1);
      fr->GetYaxis()->SetLabelSize(c == 0 ? 0.065 : 0.); fr->GetXaxis()->SetLabelSize(0.); fr->Draw("axis");
      double z[nWin] = {0};
      graph(nWin, x, ex, mN[c], z, hexBkg2, 1, 0.)->Draw("z same");
      graph(nWin, x, ex, mM[c], eM[c], hexMC, markOpenSquare)->Draw("pz same");
      graph(nWin, x, ex, mD[c], eD[c], hexCorrected, markFilledCircle)->Draw("pz same");
      text(c == 0 ? 0.45 : 0.08, 0.92, clsName[c], 0.07);
      pane(cv, c, nCls, false, split);
      TH1D *fb = new TH1D(Form("fbM%d", c), "", 1, 80, 300); fb->SetStats(0); fb->SetMinimum(-0.3); fb->SetMaximum(0.45);
      styleRatioAxes(fb, "p_{T}^{jet} [GeV]", c == 0 ? "data #minus P+H" : ""); fb->GetYaxis()->SetLabelSize(c == 0 ? 0.085 : 0.); fb->Draw("axis");
      TLine *l0 = new TLine(80, 0, 300, 0); l0->SetLineStyle(2); l0->SetLineColor(kGray+2); l0->Draw();
      graph(nWin, x, ex, dif[c], edif[c], hexCorrected, markFilledCircle)->Draw("pz same");
    }
    cv->cd(); TLegend *lg = makeLegend(0.30, 0.90, 0.99, 0.99, 0.030); lg->SetNColumns(3);
    TH1D *a = new TH1D("la", "", 1, 0, 1), *b = new TH1D("lb", "", 1, 0, 1), *n = new TH1D("ln", "", 1, 0, 1);
    styleH(a, hexCorrected, markFilledCircle); styleH(b, hexMC, markOpenSquare); styleLine(n, hexBkg2, 2);
    lg->AddEntry(a, "data all tracks, reweighted", "pe"); lg->AddEntry(b, "P+H light-jet muons", "pe"); lg->AddEntry(n, "nominal light template", "l"); lg->Draw();
    text(0.01, 0.985, hdr, 0.026);
    savePdfTight(cv, Form("%s/allTrkVsMu_meanPtRel_PbPb.pdf", outDir));
  }
  // ---------------- figure 2: shapes, 80-150 GeV
  {
    TCanvas *cv = new TCanvas("cvS", "", 1600, 600);
    for(int c = 0; c < nCls; c++){
      pane(cv, c, nCls, true, split);
      styleH(shD[c], hexCorrected, markFilledCircle, 0.9); styleH(shM[c], hexMC, markOpenSquare, 0.9);
      TH1D *fr = (TH1D*) shM[c]->Clone(Form("frS%d", c)); fr->Reset(); fr->SetDirectory(nullptr); fr->SetMinimum(0.); fr->SetMaximum(1.4*TMath::Max(shD[c]->GetMaximum(), shM[c]->GetMaximum()));
      fr->GetYaxis()->SetTitle(c == 0 ? "(1/N) dN/dp_{T}^{rel} [GeV^{-1}]" : ""); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetTitleOffset(1.1);
      fr->GetYaxis()->SetLabelSize(c == 0 ? 0.065 : 0.); fr->GetXaxis()->SetLabelSize(0.); fr->Draw("axis");
      shM[c]->Draw("same e1"); shD[c]->Draw("same e1");
      text(c == 0 ? 0.45 : 0.08, 0.92, Form("%s, 80 < p_{T}^{jet} < 150 GeV", clsName[c]), 0.062);
      pane(cv, c, nCls, false, split);
      TH1D *r = makeRatio(shD[c], shM[c], Form("rS%d", c), RatioErr::kBoth); styleH(r, hexCorrected, markFilledCircle, 0.9);
      r->SetMinimum(0.); r->SetMaximum(3.); styleRatioAxes(r, "p_{T}^{rel} [GeV]", c == 0 ? "data / P+H" : ""); r->GetYaxis()->SetLabelSize(c == 0 ? 0.085 : 0.);
      r->Draw("e1"); unityLine(0, fitHi)->Draw();
    }
    cv->cd(); TLegend *lg = makeLegend(0.45, 0.90, 0.99, 0.99, 0.030); lg->SetNColumns(2);
    lg->AddEntry(shD[0], "data all tracks, reweighted", "pe"); lg->AddEntry(shM[0], "P+H light-jet muons", "pe"); lg->Draw();
    text(0.01, 0.985, hdr, 0.026);
    savePdfTight(cv, Form("%s/allTrkVsMu_shapes_PbPb.pdf", outDir));
  }
  // ---------------- figure 3: excess vs the fraction of jets with a 15 GeV track
  {
    TCanvas *cv = new TCanvas("cvF", "", 1300, 560);
    const char *hex[NCls] = {okabeHex[6], okabeHex[1], okabeHex[3], okabeHex[5]}; const int mk[NCls] = {markFilledCircle, markFilledSquare, markFilledDiamond, markOpenCircle};
    cv->cd(); TPad *pl = new TPad("pl", "", 0., 0., 0.5, 0.92), *pr = new TPad("pr", "", 0.5, 0., 1., 0.92);
    for(TPad *p : {pl, pr}){ p->SetLeftMargin(0.15); p->SetBottomMargin(0.14); p->SetTopMargin(0.04); p->SetRightMargin(0.03); p->Draw(); }
    pl->cd(); TH1D *fl = new TH1D("fl", "", 1, 80, 300); fl->SetStats(0); fl->SetMinimum(0.); fl->SetMaximum(1.);
    fl->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]"); fl->GetYaxis()->SetTitle("f_{15}: jets with a track > 15 GeV"); fl->GetXaxis()->SetTitleSize(0.05); fl->GetYaxis()->SetTitleSize(0.05); fl->Draw("axis");
    TLegend *lg = makeLegend(0.62, 0.18, 0.95, 0.42, 0.045);
    for(int c = 0; c < nCls; c++){ TGraphErrors *g = graph(nWin, x, ex, f15[c], ef15[c], hex[c], mk[c]); g->Draw("pz same"); lg->AddEntry(g, clsName[c], "p"); }
    lg->Draw();
    pr->cd(); TH1D *fx = new TH1D("fx", "", 1, 0., 1.); fx->SetStats(0); fx->SetMinimum(-0.3); fx->SetMaximum(0.45);
    fx->GetXaxis()->SetTitle("f_{15}"); fx->GetYaxis()->SetTitle("#LTp_{T}^{rel}#GT_{0-3}: data #minus P+H light muon [GeV]"); fx->GetXaxis()->SetTitleSize(0.05); fx->GetYaxis()->SetTitleSize(0.045); fx->GetYaxis()->SetTitleOffset(1.4); fx->Draw("axis");
    TLine *l0 = new TLine(0, 0, 1, 0); l0->SetLineStyle(2); l0->SetLineColor(kGray+2); l0->Draw();
    for(int c = 0; c < nCls; c++){ TGraphErrors *g = new TGraphErrors(nWin, f15[c], dif[c], ef15[c], edif[c]); const int col = TColor::GetColor(hex[c]);
      g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mk[c]); g->SetMarkerSize(1.2); g->Draw("pz same"); }
    text(0.18, 0.93, "each point: one class and jet-p_{T} window (80-300 GeV)", 0.040);
    cv->cd(); text(0.01, 0.995, hdr, 0.026);
    savePdfTight(cv, Form("%s/allTrkVsMu_excessVsTrackFrac_PbPb.pdf", outDir));
  }
}
