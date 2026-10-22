// MC closure of the data-driven light template's premise: do light-jet TRACKS and
// light-jet MUONS share a ptRel distribution at fixed pT? 2026-10-21, user. Inside one
// PYTHIA+HYDJET scan with doTrackPtRelHists (config_PYTHIAHYDJET.h): the same events,
// jet selection, weights and truth flavor give both
//   tracks  h_jetTrkPtRel_trkPt_recoJetPt_<fl>_C<i>    every track (pT > 15, dR < 0.4)
//   muons   h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_<fl>_C<i>
// light = u + d + s + g, summed over each coarse class's 5% slices. The tracks are
// reweighted in track pT to the light-muon pT of the same class and window
// (studyDataDrivenLightAllTrk_PF_PbPb.C's reweighted(), as done to the data) and
// compared with the muons' own ptRel. If they agree, the data-driven template stands
// in for light muons and its difference from P+H in data is data/MC; if not, the gap
// here is a correction the method needs.
// Figures (figures/bPurity/PFJets_PbPb/dataDrivenLight/), <tag> = trkMuClosure or
// trkMuClosure_selftest:
//   <tag>_meanPtRel.pdf      <ptRel>_0-3 vs jet pT per class: muons, tracks
//                            reweighted, tracks unweighted; bottom: tracks rw - muons
//   <tag>_shapes.pdf         ptRel shapes, jet pT 80-150 GeV merged: muons vs tracks
//                            reweighted; bottom: tracks / muons
//   <tag>_fixedPt.pdf        the premise itself, no reweighting: <ptRel>_0-3 in
//                            slices of track / muon pT, jet pT 80-150 GeV, per class
// "selftest": the file's light MUONS stand in for the tracks, so every comparison must
// close exactly -- checks this macro on any template scan, closure file or not.
// Usage, from src/plots/bPurity/:
//   root -l -b -q 'plotTrackMuonClosure_PF_PbPb.C("<P+H _trkPtRel scan>.root")'
//   root -l -b -q 'plotTrackMuonClosure_PF_PbPb.C("selftest")'
#include "plotAllTrkVsLightMuon_PF_PbPb.C"

void plotTrackMuonClosure_PF_PbPb(const char *phFile = "selftest"){
  setJetCollection("PFJets"); initPlotStyle(); gStyle->SetEndErrorSize(0);
  const bool self = TString(phFile) == "selftest";
  TFile *f = TFile::Open(self ? phLightFile : phFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", self ? phLightFile : phFile); return; }
  const char *muName = "h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn", *trkName = self ? muName : "h_jetTrkPtRel_trkPt_recoJetPt";
  const TString tag = self ? "trkMuClosure_selftest" : "trkMuClosure";
  gSystem->mkdir(outDir, kTRUE);
  auto light = [&](const char *base, int s0, int s1){ TH3 *m = nullptr; for(const char *fl : {"u", "d", "s", "g"}){ TH3 *h = sumSlices(f, Form("%s_%sJets", base, fl), s0, s1); if(!h) return (TH3*) nullptr; if(!m) m = h; else { m->Add(h); delete h; } } return m; };

  double x[nWin], ex[nWin];
  for(int k = 0; k < nWin; k++){ x[k] = 0.5*(win[k][0] + win[k][1]); ex[k] = 0.5*(win[k][1] - win[k][0]); }
  const int nPS = 5; const double pS[nPS+1] = {15, 20, 25, 35, 50, 100};   // track / muon pT slices (bin edges of ptRelMuPt3D.h)
  double xp[nPS], exp_[nPS];
  for(int j = 0; j < nPS; j++){ xp[j] = 0.5*(pS[j] + pS[j+1]); exp_[j] = 0.5*(pS[j+1] - pS[j]); }
  double mM[NCls][nWin], eM[NCls][nWin], mR[NCls][nWin], eR[NCls][nWin], mU[NCls][nWin], eU[NCls][nWin], dRM[NCls][nWin], edRM[NCls][nWin];
  double fM[NCls][nPS], efM[NCls][nPS], fT[NCls][nPS], efT[NCls][nPS];
  TH1D *shM[NCls], *shR[NCls];
  printf("%s: %s\nclass   window   muons          tracks rw      tracks unw     rw - mu          mu Neff\n", tag.Data(), self ? phLightFile : phFile);
  for(int c = 0; c < nCls; c++){
    TH3 *tr = light(trkName, clsSlice[c][0], clsSlice[c][1]), *mu = light(muName, clsSlice[c][0], clsSlice[c][1]);
    if(!tr || !mu){ printf("ERROR: %s light %s missing -- %s\n", clsName[c], !tr ? "track 3D" : "muon 3D", !tr ? "a scan with doTrackPtRelHists is needed" : ""); return; }
    for(int k = 0; k < nWin; k++){
      TH2D *T = proj2(tr, win[k][0], win[k][1]), *M = proj2(mu, win[k][0], win[k][1]);
      TH1D *tg = M->ProjectionY(Form("tg%d", uid++)), *Mx = M->ProjectionX(Form("mx%d", uid++)), *Tx = T->ProjectionX(Form("tx%d", uid++)), *R = reweighted(T, tg);
      meanErr(Mx, mM[c][k], eM[c][k]); meanErr(R, mR[c][k], eR[c][k]); meanErr(Tx, mU[c][k], eU[c][k]);
      dRM[c][k] = mR[c][k] - mM[c][k]; edRM[c][k] = self ? 0. : hypot(eR[c][k], eM[c][k]);   // selftest: identical inputs, fully correlated
      printf("%-7s %3.0f-%-3.0f  %.3f+-%.3f   %.3f+-%.3f   %.3f+-%.3f   %+.3f+-%.3f   %5.0f\n", clsName[c], win[k][0], win[k][1], mM[c][k], eM[c][k], mR[c][k], eR[c][k], mU[c][k], eU[c][k], dRM[c][k], edRM[c][k], neff(Mx));
      delete T; delete M; delete tg; delete Mx; delete Tx; delete R;
    }
    TH2D *T = proj2(tr, 80, 150), *M = proj2(mu, 80, 150);
    TH1D *tg = M->ProjectionY(Form("tg%d", uid++)), *R = reweighted(T, tg), *Mx = M->ProjectionX(Form("mx%d", uid++));
    shR[c] = range05(R, Form("shR%d", c)); shM[c] = range05(Mx, Form("shM%d", c));
    for(TH1D *h : {shR[c], shM[c]}){ h->Scale(1./h->Integral()); h->Scale(1., "width"); }
    for(int j = 0; j < nPS; j++){
      const int ylo = T->GetYaxis()->FindBin(pS[j] + 1e-6), yhi = T->GetYaxis()->FindBin(pS[j+1] - 1e-6);
      TH1D *a = T->ProjectionX(Form("pa%d", uid++), ylo, yhi), *b = M->ProjectionX(Form("pb%d", uid++), ylo, yhi);
      meanErr(a, fT[c][j], efT[c][j]); meanErr(b, fM[c][j], efM[c][j]); delete a; delete b;
    }
    delete T; delete M; delete tg; delete R; delete Mx; delete tr; delete mu;
  }

  const char *hdr = self ? "PYTHIA+HYDJET, light jets; SELF-TEST: muons used as the tracks (must close exactly)" : "PYTHIA+HYDJET, light (u, d, s, g) jets: every track (p_{T} > 15 GeV, #DeltaR < 0.4) vs the tagging muon, same events and weights";
  const double split = 0.36;
  { // figure 1
    TCanvas *cv = new TCanvas("cvC1", "", 1600, 600);
    for(int c = 0; c < nCls; c++){
      pane(cv, c, nCls, true, split);
      TH1D *fr = new TH1D(Form("frC1%d", c), "", 1, 80, 300); fr->SetStats(0); fr->SetMinimum(0.5); fr->SetMaximum(1.5);
      fr->GetYaxis()->SetTitle(c == 0 ? "#LTp_{T}^{rel}#GT_{0-3 GeV} [GeV]" : ""); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetTitleOffset(1.1);
      fr->GetYaxis()->SetLabelSize(c == 0 ? 0.065 : 0.); fr->GetXaxis()->SetLabelSize(0.); fr->Draw("axis");
      graph(nWin, x, ex, mU[c], eU[c], hexData, markOpenCircle)->Draw("pz same");
      graph(nWin, x, ex, mM[c], eM[c], hexMC, markOpenSquare)->Draw("pz same");
      graph(nWin, x, ex, mR[c], eR[c], hexCorrected, markFilledCircle)->Draw("pz same");
      text(c == 0 ? 0.45 : 0.08, 0.92, clsName[c], 0.07);
      pane(cv, c, nCls, false, split);
      TH1D *fb = new TH1D(Form("fbC1%d", c), "", 1, 80, 300); fb->SetStats(0); fb->SetMinimum(-0.3); fb->SetMaximum(0.3);
      styleRatioAxes(fb, "p_{T}^{jet} [GeV]", c == 0 ? "tracks #minus muons" : ""); fb->GetYaxis()->SetLabelSize(c == 0 ? 0.085 : 0.); fb->Draw("axis");
      TLine *l0 = new TLine(80, 0, 300, 0); l0->SetLineStyle(2); l0->SetLineColor(kGray+2); l0->Draw();
      graph(nWin, x, ex, dRM[c], edRM[c], hexCorrected, markFilledCircle)->Draw("pz same");
    }
    cv->cd(); TLegend *lg = makeLegend(0.30, 0.90, 0.99, 0.99, 0.030); lg->SetNColumns(3);
    TH1D *a = new TH1D("cla", "", 1, 0, 1), *b = new TH1D("clb", "", 1, 0, 1), *u = new TH1D("clu", "", 1, 0, 1);
    styleH(a, hexCorrected, markFilledCircle); styleH(b, hexMC, markOpenSquare); styleH(u, hexData, markOpenCircle);
    lg->AddEntry(b, "light-jet muons", "pe"); lg->AddEntry(a, "tracks, reweighted to muon p_{T}", "pe"); lg->AddEntry(u, "tracks, unweighted", "pe"); lg->Draw();
    text(0.01, 0.985, hdr, 0.026);
    savePdfTight(cv, Form("%s/%s_meanPtRel.pdf", outDir, tag.Data()));
  }
  { // figure 2
    TCanvas *cv = new TCanvas("cvC2", "", 1600, 600);
    for(int c = 0; c < nCls; c++){
      pane(cv, c, nCls, true, split);
      styleH(shR[c], hexCorrected, markFilledCircle, 0.9); styleH(shM[c], hexMC, markOpenSquare, 0.9);
      TH1D *fr = (TH1D*) shM[c]->Clone(Form("frC2%d", c)); fr->Reset(); fr->SetDirectory(nullptr); fr->SetMinimum(0.); fr->SetMaximum(1.4*TMath::Max(shR[c]->GetMaximum(), shM[c]->GetMaximum()));
      fr->GetYaxis()->SetTitle(c == 0 ? "(1/N) dN/dp_{T}^{rel} [GeV^{-1}]" : ""); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetTitleOffset(1.1);
      fr->GetYaxis()->SetLabelSize(c == 0 ? 0.065 : 0.); fr->GetXaxis()->SetLabelSize(0.); fr->Draw("axis");
      shM[c]->Draw("same e1"); shR[c]->Draw("same e1");
      text(c == 0 ? 0.45 : 0.08, 0.92, Form("%s, 80 < p_{T}^{jet} < 150 GeV", clsName[c]), 0.062);
      pane(cv, c, nCls, false, split);
      TH1D *r = makeRatio(shR[c], shM[c], Form("rC2%d", c), self ? RatioErr::kNumerator : RatioErr::kBoth); styleH(r, hexCorrected, markFilledCircle, 0.9);
      r->SetMinimum(0.); r->SetMaximum(2.5); styleRatioAxes(r, "p_{T}^{rel} [GeV]", c == 0 ? "tracks / muons" : ""); r->GetYaxis()->SetLabelSize(c == 0 ? 0.085 : 0.);
      r->Draw("e1"); unityLine(0, fitHi)->Draw();
    }
    cv->cd(); TLegend *lg = makeLegend(0.45, 0.90, 0.99, 0.99, 0.030); lg->SetNColumns(2);
    lg->AddEntry(shM[0], "light-jet muons", "pe"); lg->AddEntry(shR[0], "tracks, reweighted to muon p_{T}", "pe"); lg->Draw();
    text(0.01, 0.985, hdr, 0.026);
    savePdfTight(cv, Form("%s/%s_shapes.pdf", outDir, tag.Data()));
  }
  { // figure 3: fixed pT
    TCanvas *cv = new TCanvas("cvC3", "", 1600, 600);
    for(int c = 0; c < nCls; c++){
      pane(cv, c, nCls, true, split);
      TH1D *fr = new TH1D(Form("frC3%d", c), "", 1, 15, 100); fr->SetStats(0); fr->SetMinimum(0.3); fr->SetMaximum(1.6);
      fr->GetYaxis()->SetTitle(c == 0 ? "#LTp_{T}^{rel}#GT_{0-3 GeV} [GeV]" : ""); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetTitleOffset(1.1);
      fr->GetYaxis()->SetLabelSize(c == 0 ? 0.065 : 0.); fr->GetXaxis()->SetLabelSize(0.); fr->Draw("axis");
      graph(nPS, xp, exp_, fM[c], efM[c], hexMC, markOpenSquare)->Draw("pz same");
      graph(nPS, xp, exp_, fT[c], efT[c], hexCorrected, markFilledCircle)->Draw("pz same");
      text(c == 0 ? 0.45 : 0.08, 0.92, Form("%s, 80 < p_{T}^{jet} < 150 GeV", clsName[c]), 0.062);
      pane(cv, c, nCls, false, split);
      double d[nPS], ed[nPS]; for(int j = 0; j < nPS; j++){ d[j] = fT[c][j] - fM[c][j]; ed[j] = self ? 0. : hypot(efT[c][j], efM[c][j]); }
      TH1D *fb = new TH1D(Form("fbC3%d", c), "", 1, 15, 100); fb->SetStats(0); fb->SetMinimum(-0.4); fb->SetMaximum(0.4);
      styleRatioAxes(fb, "track / muon p_{T} [GeV]", c == 0 ? "tracks #minus muons" : ""); fb->GetYaxis()->SetLabelSize(c == 0 ? 0.085 : 0.); fb->Draw("axis");
      TLine *l0 = new TLine(15, 0, 100, 0); l0->SetLineStyle(2); l0->SetLineColor(kGray+2); l0->Draw();
      graph(nPS, xp, exp_, d, ed, hexCorrected, markFilledCircle)->Draw("pz same");
    }
    cv->cd(); TLegend *lg = makeLegend(0.45, 0.90, 0.99, 0.99, 0.030); lg->SetNColumns(2);
    TH1D *a = new TH1D("c3a", "", 1, 0, 1), *b = new TH1D("c3b", "", 1, 0, 1); styleH(a, hexCorrected, markFilledCircle); styleH(b, hexMC, markOpenSquare);
    lg->AddEntry(b, "light-jet muons", "pe"); lg->AddEntry(a, "tracks (no reweighting)", "pe"); lg->Draw();
    text(0.01, 0.985, hdr, 0.026);
    savePdfTight(cv, Form("%s/%s_fixedPt.pdf", outDir, tag.Data()));
  }
}
