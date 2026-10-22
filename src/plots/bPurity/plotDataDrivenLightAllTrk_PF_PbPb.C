// Figures for the construction of the PbPb data-driven light templates from ALL
// tracks (studyDataDrivenLightAllTrk_PF_PbPb.C, included below so the drawn
// templates are exactly the fitted ones). 2026-10-19, user. The PbPb counterpart of
// plotDataDrivenLightReweighting_pp_PF.C. One page per coarse class, jet-pT windows
// 80-100 ... 200-300 GeV across:
//   allTrkLight_trackPtWeight_PbPb.pdf
//     top: the pT spectrum T(p) of every charged track (pT > 15 GeV, dR < 0.4) of the
//     class's jets in PbPb MinBias data, and the target M_L(p), the PYTHIA+HYDJET light
//     (u, d, s, g) tagged-jet muon pT of 0-80% (the study's default); unit area, per
//     GeV (log y: steeply falling, all positive). Bottom: w(p) = M_L(p) / T(p).
//   allTrkLight_ptRelTemplates_PbPb.pdf
//     top: the data all-track ptRel, unweighted and reweighted (sum over p of
//     T(ptRel, p) w(p)), against the nominal light template the fit would otherwise
//     use (pp PYTHIA light smeared by the class's fitted axis sigma, kLightSmearedPP);
//     in the fit binning, unit area, per GeV. Bottom: ratio to the nominal.
// Output: figures/bPurity/PFJets_PbPb/dataDrivenLight/
// Usage, from src/plots/bPurity/:  root -l -b -q plotDataDrivenLightAllTrk_PF_PbPb.C
#include "studyDataDrivenLightAllTrk_PF_PbPb.C"

namespace pda {
const int nW = 5;
// one panel of a 1 x nW row: top above split, bottom below
TPad* cell(TCanvas *c, int k, bool top, double split){
  const double x0 = 0.04, w = (1. - x0 - 0.01)/nW;
  const double xl = x0 + k*w - (k == 0 ? x0 : 0.), xh = x0 + (k+1)*w;
  TPad *p = new TPad(Form("pd%d%d_%d", k, top, dla::uid++), "", xl, top ? split : 0., xh, top ? 0.93 : split);
  p->SetLeftMargin(k == 0 ? x0/(xh - xl) + 0.10 : 0.02); p->SetRightMargin(0.02);
  if(top){ p->SetTopMargin(0.04); p->SetBottomMargin(0.015); } else { p->SetTopMargin(0.02); p->SetBottomMargin(0.30); }
  c->cd(); p->Draw(); p->cd(); return p;
}
void label(double x, double y, const char *s, double size = 0.065, int align = 13){
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size); t.SetTextAlign(align); t.DrawLatex(x, y, s);
}
void density(TH1D *h){ const double s = h->Integral(1, h->GetNbinsX()); if(s > 0) h->Scale(1./s); h->Scale(1., "width"); }
}
using namespace pda;

void plotDataDrivenLightAllTrk_PF_PbPb(){
  setJetCollection("PFJets"); RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL); initPlotStyle(); gStyle->SetEndErrorSize(0);
  TFile *fT = TFile::Open(trkFile), *fL = TFile::Open(phLightFile); PHFiles ph;
  if(!fT || !fL || !openPH(ph)){ printf("ERROR: inputs\n"); return; }
  gSystem->mkdir(outDir, kTRUE);
  // the study's default target: PYTHIA+HYDJET light muons of 0-80%, every class
  TH3 *mL = nullptr;
  for(const char *fl : {"u", "d", "s", "g"}){ TH3 *h = sumSlices(fL, Form("h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_%sJets", fl), 1, 16); if(!h){ printf("ERROR: light-muon 3D missing\n"); return; } if(!mL) mL = h; else { mL->Add(h); delete h; } }

  const TString fPt = Form("%s/allTrkLight_trackPtWeight_PbPb.pdf", outDir), fTp = Form("%s/allTrkLight_ptRelTemplates_PbPb.pdf", outDir);
  TCanvas *c1 = new TCanvas("c1", "", 1700, 640), *c2 = new TCanvas("c2", "", 1700, 640);
  const double split = 0.36;
  for(int c = 0; c < nCls; c++){
    TH3 *tAll = sumSlices(fT, "h_jetTrkPtRel_trkPt_recoJetPt", clsSlice[c][0], clsSlice[c][1]);
    if(!tAll){ printf("ERROR: %s: all-track 3D missing\n", clsName[c]); return; }
    c1->Clear(); c2->Clear();
    printf("%s  window   tracks(all, pT>15)  <pT> tracks / mu   <ptRel 0-3> unweighted / reweighted / nominal light\n", clsName[c]);
    for(int k = 0; k < nW; k++){
      const double lo = win[k][0], hi = win[k][1];
      TH2D *T = proj2(tAll, lo, hi), *Mmu = proj2(mL, lo, hi);
      TH1D *target = Mmu->ProjectionY(Form("tg%d", uid++));
      TH1D *Lrw = reweighted(T, target), *Tx = T->ProjectionX(Form("tx%d", uid++));
      Tpl t = buildTemplatesPbPb(ph, c, lo, hi);
      TH1D *lRw = range05(Lrw, Form("lRw%d", uid++)), *lUn = range05(Tx, Form("lUn%d", uid++)), *lNom = (TH1D*) t.l->Clone(Form("lNom%d", uid++)); lNom->SetDirectory(nullptr);
      // track pT spectra (means from the counts) and the weight
      TH1D *tPt = T->ProjectionY(Form("tPt%d", uid++)), *mPt = (TH1D*) target->Clone(Form("mPt%d", uid++));
      const double nTrk = tPt->Integral(1, tPt->GetNbinsX()), tMean = tPt->GetMean(), mMean = mPt->GetMean();
      TH1D *tn = (TH1D*) tPt->Clone(Form("tn%d", uid++)), *mn = (TH1D*) mPt->Clone(Form("mn%d", uid++));
      tn->Scale(1./tn->Integral(1, tn->GetNbinsX())); mn->Scale(1./mn->Integral(1, mn->GetNbinsX()));
      TH1D *wPt = makeRatio(mn, tn, Form("w%d", uid++), RatioErr::kBoth);   // the weight reweighted() applies
      density(tPt); density(mPt);
      for(TH1D *h : {lRw, lUn, lNom}){ h->Scale(1./h->Integral()); h->Scale(1., "width"); }
      printf("        %3.0f-%-3.0f  %8.0f            %5.1f / %5.1f        %.3f / %.3f / %.3f\n", lo, hi, nTrk, tMean, mMean, mean3(lUn), mean3(lRw), mean3(lNom));

      // figure 1
      TPad *p = cell(c1, k, true, split); p->SetLogy(); p->SetLogx();
      styleH(tPt, hexData, markOpenCircle, 0.8); styleH(mPt, hexMC, markFilledSquare, 0.8);
      TH1D *fr = (TH1D*) tPt->Clone(Form("fr%d", uid++)); fr->Reset(); fr->SetDirectory(nullptr);
      fr->GetXaxis()->SetRangeUser(15, 150); fr->SetMinimum(2e-5); fr->SetMaximum(1.);
      fr->GetYaxis()->SetTitle(k == 0 ? "(1/N) dN/dp_{T} [GeV^{-1}]" : ""); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetTitleOffset(1.25);
      fr->GetYaxis()->SetLabelSize(k == 0 ? 0.065 : 0.); fr->GetXaxis()->SetLabelSize(0.);
      fr->Draw("axis"); mPt->Draw("same e1"); tPt->Draw("same e1");
      label(k == 0 ? 0.42 : 0.06, 0.94, Form("%.0f < p_{T}^{jet} < %.0f GeV", lo, hi), 0.068);
      label(k == 0 ? 0.42 : 0.06, 0.85, Form("#LTp_{T}#GT: %.1f / %.1f GeV", tMean, mMean), 0.060);
      if(k == nW - 1){ TLegend *lg = makeLegend(0.06, 0.05, 0.98, 0.25, 0.056); lg->AddEntry(tPt, "data, all tracks", "pe"); lg->AddEntry(mPt, "P+H light-jet muon, 0-80%", "pe"); lg->Draw(); }
      TPad *q = cell(c1, k, false, split); q->SetLogx();
      styleH(wPt, hexCorrected, markFilledCircle, 0.8); wPt->GetXaxis()->SetRangeUser(15, 150); wPt->SetMinimum(0.); wPt->SetMaximum(4.);
      styleRatioAxes(wPt, "p_{T} [GeV]", k == 0 ? "w = M_{L} / T" : "");
      wPt->GetXaxis()->SetTitleSize(0.10); wPt->GetXaxis()->SetLabelSize(0.085); wPt->GetXaxis()->SetMoreLogLabels(); wPt->GetXaxis()->SetNoExponent();
      wPt->GetYaxis()->SetTitleSize(0.10); wPt->GetYaxis()->SetTitleOffset(0.75); wPt->GetYaxis()->SetLabelSize(k == 0 ? 0.085 : 0.);
      wPt->Draw("e1"); unityLine(15, 150)->Draw();

      // figure 2
      p = cell(c2, k, true, split);
      styleH(lUn, hexData, markOpenCircle, 0.8); styleH(lRw, hexCorrected, markFilledCircle, 0.8); styleLine(lNom, hexMC, 3);
      TH1D *f2 = (TH1D*) lNom->Clone(Form("f2%d", uid++)); f2->Reset(); f2->SetDirectory(nullptr);
      f2->SetMinimum(0.); f2->SetMaximum(1.6*TMath::Max(lNom->GetMaximum(), TMath::Max(lRw->GetMaximum(), lUn->GetMaximum())));
      f2->GetYaxis()->SetTitle(k == 0 ? "(1/N) dN/dp_{T}^{rel} [GeV^{-1}]" : ""); f2->GetYaxis()->SetTitleSize(0.075); f2->GetYaxis()->SetTitleOffset(1.25);
      f2->GetYaxis()->SetLabelSize(k == 0 ? 0.065 : 0.); f2->GetXaxis()->SetLabelSize(0.);
      f2->Draw("axis"); lNom->Draw("same hist"); lUn->Draw("same e1"); lRw->Draw("same e1");
      label(k == 0 ? 0.42 : 0.06, 0.94, Form("%.0f < p_{T}^{jet} < %.0f GeV", lo, hi), 0.068);
      label(k == 0 ? 0.42 : 0.06, 0.85, Form("#LTp_{T}^{rel}#GT_{0-3}: %.2f / %.2f / %.2f", mean3(lUn), mean3(lRw), mean3(lNom)), 0.052);
      if(k == nW - 1){ TLegend *lg = makeLegend(0.12, 0.50, 0.98, 0.78, 0.052);
        lg->AddEntry(lUn, "data, unweighted", "pe"); lg->AddEntry(lRw, "data, reweighted", "pe"); lg->AddEntry(lNom, "nominal (smeared pp)", "l"); lg->Draw(); }
      q = cell(c2, k, false, split);
      TH1D *rU = makeRatio(lUn, lNom, Form("rU%d", uid++)), *rR = makeRatio(lRw, lNom, Form("rR%d", uid++));
      styleH(rU, hexData, markOpenCircle, 0.8); styleH(rR, hexCorrected, markFilledCircle, 0.8);
      rR->SetMinimum(0.); rR->SetMaximum(3.); styleRatioAxes(rR, "p_{T}^{rel} [GeV]", k == 0 ? "/ nominal" : "");
      rR->GetXaxis()->SetTitleSize(0.10); rR->GetXaxis()->SetLabelSize(0.085);
      rR->GetYaxis()->SetTitleSize(0.10); rR->GetYaxis()->SetTitleOffset(0.75); rR->GetYaxis()->SetLabelSize(k == 0 ? 0.085 : 0.);
      rR->Draw("e1"); rU->Draw("same e1"); rR->Draw("same e1"); unityLine(0, fitHi)->Draw();

      delete T; delete Mmu; delete target; delete Lrw; delete Tx; delete tn; delete mn;
    }
    for(TCanvas *cv : {c1, c2}){ cv->cd(); label(0.01, 0.995, Form("PbPb 5.02 TeV, %s, anti-k_{T} R = 0.4 PF jets; tracks: charged PF candidates, p_{T} > 15 GeV, #DeltaR < 0.4", clsName[c]), 0.036, 13); }
    const char *open = c == 0 ? "(" : (c == nCls - 1 ? ")" : "");
    c1->Print(fPt + open, Form("Title:%s", clsName[c])); c2->Print(fTp + open, Form("Title:%s", clsName[c]));
    delete tAll;
  }
  for(const TString &f : {fPt, fTp}){
    const TString tmp = f + ".crop.pdf";
    if(gSystem->Exec(Form("pdfcrop --margins 8 %s %s > /dev/null 2>&1", f.Data(), tmp.Data())) == 0) gSystem->Rename(tmp, f);
  }
}
