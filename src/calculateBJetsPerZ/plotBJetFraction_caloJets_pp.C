// pp CALO jets: inclusive b-jet spectrum and inclusive-jet spectrum on one
// plot, and their ratio, the inclusive b-jet fraction.
//
//   b jets     correctMuonTagFrequency_caloJets_pp.C: purity x muon-tagged
//              spectrum (SingleMuon, mu12 trigger-efficiency corrected),
//              unfolded N = 2, muon-reco and tag-frequency corrected. dN/dpT.
//   inclusive  inclusiveSpectrumSys_caloJets_pp.C, unfolded N = 2 (as the b jets), b-jet
//              binning. Jets per Z per bin, in Jet100-sample units.
//
// NORMALIZATION. The inclusive spectrum is per Z, with N_Z counted in the
// 2026-05-04 SingleMuon scan; the b-jet spectrum is a raw yield in the
// 2026-09-29 calo SingleMuon scan. Both are the same dataset (N_Z 112760 vs
// 112956, h_vz entries within 0.2%), and the inclusive per-Z normalization
// already assumes the Jet100 and SingleMuon samples share a luminosity. Each
// spectrum is divided by the N_Z of its own file, so the fraction is
// (b / N_Z^b) / (incl / N_Z^incl), and both are drawn x N_Z^incl, i.e. in the
// Jet100-sample-equivalent yield of the inclusive spectrum figure. The two N_Z
// count the same events, so their statistical errors cancel and are not added.
//
// UNCERTAINTIES ON THE FRACTION.
//   spectrum JEU   both spectra's JEU source is the same pair of shifted
//                  PYTHIA calo reco spectra (rootFiles/JEU/unfold_PYTHIA_caloJets_JEUShift.root),
//                  unfolded with each chain's own response and N. Fully
//                  correlated: the fraction's JEU systematic is the difference
//                  of the two spectra's signed JEU shifts,
//                    max_k |d_b,k - d_incl,k|,   k = up, down.
//   everything else uncorrelated, quadrature: b purity, iterations, response
//                  MC stat., 70-80 GeV turn-on and muon corrections on the b side; iterations and
//                  response MC stat. on the inclusive side (different responses,
//                  different MC samples).
//   stat.          quadrature (b jets from SingleMuon, inclusive from the jet-
//                  triggered samples).
// The b-purity JEU source (jet energy scale moving the purity-fit windows) is
// also correlated with the spectrum JEU but enters the b purity, not the
// spectrum; it stays in the b-purity total, uncorrelated, as in the b chain.
//
// MC TRUTH. PYTHIA gen-level b-jet fraction from the calo PYTHIA scan (mcPath,
// caloJetBPurityFit_pp.h): h_inclGenJetPt_flavor, b = |flavor| 5 + gluon
// splitting 17 (as correctMuonTagFrequency_caloJets_pp.C) over all flavours, vs
// gen jet pT, |eta| < 1.6, no trigger. This gen-jet histogram is filled outside
// the muon-tag code and is not affected by the pre-2026-10-04 gen-muon-tag bug.
//
// Output  rootFiles/CorrectedBJetSpectra/Data/bJetFraction_caloJets_pp.root
//         figures/bJetSpectra/bAndInclusiveJetSpectra_caloJets_pp.pdf
//         figures/bJetSpectra/bJetFraction_caloJets_pp.pdf
//
// Usage: root -l -b -q plotBJetFraction_caloJets_pp.C
// Run from: src/calculateBJetsPerZ/

#include "TFile.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TPad.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"
#include "TMath.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/functions/caloJetBPurityFit_pp.h"   // ptEdges, NPt, ptReportMin, dataPath

const char *bPath     = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_inclusive_caloJets_pp.root";
const char *bJEUPath  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_caloJets_pp.root";
const char *inclPath  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/inclJetSpectrumSys_caloJets_pp.root";
const char *nzInclPath= "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateJetsPerZ/rootFiles/JetsPerZ/jetsPerZSpectra_caloJets_coarseMuEff.root";
const char *outRoot   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/bJetFraction_caloJets_pp.root";
const char *figDir    = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/bJetSpectra";
const char *inclCfg   = "N2";           // inclusive: N = 2, as the b-jet unfolding
const double etaWidth = 3.2;            // |eta| < 1.6, as calculateJetsPerZ.cc
const double muEffZ   = 0.9708;         // same as calculateJetsPerZ.cc's muEffPP
const double zLo = 60., zHi = 120.;

// b-jet systematic sources other than the spectrum JEU (uncorrelated with the inclusive side)
const int NBSrc = 6;
const char *bSrc[NBSrc]    = {"bPurity", "iterations", "responseStat", "turnOn", "muonRecoEff", "muonTagFrequency"};
const int NISrc = 2;
const char *iSrc[NISrc]    = {"iterations", "responseStat"};

void plotBJetFraction_caloJets_pp(){
  initPlotStyle();
  gSystem->mkdir(figDir, kTRUE);

  TFile *fB = TFile::Open(bPath), *fBJ = TFile::Open(bJEUPath), *fI = TFile::Open(inclPath), *fZ = TFile::Open(nzInclPath), *fD = TFile::Open(dataPath);
  for(TFile *f : {fB, fBJ, fI, fZ, fD}) if(!f || f->IsZombie()){ printf("ERROR: cannot open an input file\n"); return; }

  // ---- inputs
  TH1D *bSpec = nullptr, *bStat = nullptr, *bJu = nullptr, *bJd = nullptr, *bJ = nullptr, *bS[NBSrc];
  fB->GetObject("h_bJetPt_inclusive", bSpec); fB->GetObject("h_incl_statRel", bStat); fB->GetObject("h_incl_sysRel_spectrumJEU", bJ);
  fBJ->GetObject("h_bJetPt_sysRel_spectrumJEU_up", bJu); fBJ->GetObject("h_bJetPt_sysRel_spectrumJEU_down", bJd);
  for(int k = 0; k < NBSrc; k++){ bS[k] = nullptr; fB->GetObject(Form("h_incl_sysRel_%s", bSrc[k]), bS[k]); }
  TH1D *iSpec = nullptr, *iStat = nullptr, *iJu = nullptr, *iJd = nullptr, *iS[NISrc];
  fI->GetObject(Form("h_incl_unfolded_%s_bjet", inclCfg), iSpec); fI->GetObject(Form("h_incl_statRel_%s_bjet", inclCfg), iStat);
  fI->GetObject(Form("h_incl_sysRel_spectrumJEU_up_%s_bjet", inclCfg), iJu); fI->GetObject(Form("h_incl_sysRel_spectrumJEU_down_%s_bjet", inclCfg), iJd);
  for(int k = 0; k < NISrc; k++){ iS[k] = nullptr; fI->GetObject(Form("h_incl_sysRel_%s_%s_bjet", iSrc[k], inclCfg), iS[k]); }
  bool ok = bSpec && bStat && bJu && bJd && bJ && iSpec && iStat && iJu && iJd;
  for(int k = 0; k < NBSrc; k++) ok = ok && bS[k];
  for(int k = 0; k < NISrc; k++) ok = ok && iS[k];
  if(!ok){ printf("ERROR: missing input histograms (rerun inclusiveSpectrumSys_caloJets_pp.C for the signed JEU shifts)\n"); return; }

  // the b chain carries the gen-level JEU unchanged: its |JEU| must equal max(|up|,|down|)
  for(int i = 1; i <= NPt; i++){
    double m = TMath::Max(fabs(bJu->GetBinContent(i)), fabs(bJd->GetBinContent(i)));
    if(fabs(m - bJ->GetBinContent(i)) > 1e-6){ printf("ERROR: b-jet JEU in %s does not match the signed shifts in %s (bin %d)\n", bPath, bJEUPath, i); return; }
  }

  // ---- Z counts
  TVectorD *nzv = nullptr; fZ->GetObject("NZ", nzv);
  TH1D *hm = nullptr; fD->GetObject("h_dimuonMass", hm);
  if(!nzv || !hm){ printf("ERROR: missing NZ / h_dimuonMass\n"); return; }
  const double nzIncl = (*nzv)[0];
  const double nzB    = hm->Integral(hm->FindBin(zLo), hm->FindBin(zHi - 1e-6))/(muEffZ*muEffZ);
  // ---- MC truth fraction
  TFile *fMC = TFile::Open(mcPath);
  TH2D *hFl = nullptr; if(fMC && !fMC->IsZombie()) fMC->GetObject("h_inclGenJetPt_flavor", hFl);
  if(!hFl){ printf("ERROR: missing h_inclGenJetPt_flavor in %s\n", mcPath); return; }
  TH1D *mcAll = hFl->ProjectionX("mcAll", 0, hFl->GetNbinsY()+1), *mcB = nullptr;
  for(int fl : {-5, 5, 17}){
    int yb = hFl->GetYaxis()->FindBin(fl + 0.1);   // flavour bins are [v, v+1)
    TH1D *p = hFl->ProjectionX(Form("mcB_f%d", fl), yb, yb);
    if(!mcB) mcB = (TH1D*) p->Clone("mcB"); else mcB->Add(p);
  }
  TH1D *mcAllR = (TH1D*) mcAll->Rebin(NPt, "mcAllR", ptEdges), *mcBR = (TH1D*) mcB->Rebin(NPt, "mcBR", ptEdges);
  TH1D *fMCt = (TH1D*) mcBR->Clone("h_bJetFraction_MCtruth"); fMCt->SetDirectory(nullptr);
  fMCt->Divide(mcBR, mcAllR, 1., 1., "B");
  fMCt->SetTitle("PYTHIA gen-level b-jet fraction (|flavor| 5 + bGS 17) / all;gen #it{p}_{T}^{jet} [GeV];b jets / inclusive jets");

  printf("N_Z: inclusive file %.1f, b-jet (SingleMuon calo 2026-09-29) %.1f, ratio %.4f\n", nzIncl, nzB, nzB/nzIncl);

  // ---- spectra in d2N/dpT deta, Jet100-sample-equivalent yield (x N_Z^incl)
  TH1D *hB = (TH1D*) bSpec->Clone("h_bJetSpectrum"); hB->SetDirectory(nullptr);
  hB->Scale(nzIncl/nzB/etaWidth);
  hB->SetTitle("b jets (inclusive);#it{p}_{T}^{jet} [GeV];d^{2}#it{N}/d#it{p}_{T}d#eta [GeV^{-1}]");
  TH1D *hI = (TH1D*) iSpec->Clone("h_inclJetSpectrum"); hI->SetDirectory(nullptr);
  for(int i = 1; i <= NPt; i++){
    double w = hI->GetBinWidth(i)*etaWidth;
    hI->SetBinContent(i, iSpec->GetBinContent(i)*nzIncl/w); hI->SetBinError(i, iSpec->GetBinError(i)*nzIncl/w);
  }
  hI->SetTitle("inclusive jets;#it{p}_{T}^{jet} [GeV];d^{2}#it{N}/d#it{p}_{T}d#eta [GeV^{-1}]");

  // ---- relative systematics
  auto book = [&](const char *n, const char *t){ TH1D *h = new TH1D(n, Form("%s;#it{p}_{T}^{jet} [GeV];relative uncertainty", t), NPt, ptEdges); h->SetDirectory(nullptr); return h; };
  TH1D *bSys  = book("h_bJet_sysRel",  "b jets: total systematic");
  TH1D *iSys  = book("h_incl_sysRel",  "inclusive jets: total systematic");
  TH1D *fF    = book("h_bJetFraction", "inclusive b-jet fraction");
  fF->GetYaxis()->SetTitle("b jets / inclusive jets");
  TH1D *fStat = book("h_bJetFraction_statRel",              "fraction: stat.");
  TH1D *fSys  = book("h_bJetFraction_sysRel",               "fraction: total systematic (JEU correlated)");
  TH1D *fJEU  = book("h_bJetFraction_sysRel_spectrumJEU",   "fraction: spectrum JEU, correlated (b and inclusive shifted together)");
  TH1D *fJEUu = book("h_bJetFraction_sysRel_spectrumJEU_uncorr", "fraction: spectrum JEU if treated as uncorrelated (comparison only)");
  TH1D *fOth  = book("h_bJetFraction_sysRel_other",         "fraction: all other sources, quadrature");
  TH1D *fNoC  = book("h_bJetFraction_sysRel_noCancel",      "fraction: total systematic, JEU uncorrelated (comparison only)");

  printf("\n%-9s %10s %10s %9s %8s | %7s %7s | %7s %7s %7s %7s | %7s\n", "jet pT", "b", "incl", "fraction", "PYTHIA",
         "b sys", "inc sys", "JEU cor", "JEU unc", "other", "TOTAL", "stat");
  for(int i = 1; i <= NPt; i++){
    double vb = hB->GetBinContent(i), vi = hI->GetBinContent(i);
    if(vb <= 0. || vi <= 0.) continue;
    double qb = 0., qi = 0., qo = 0.;
    for(int k = 0; k < NBSrc; k++){ double s = bS[k]->GetBinContent(i); qb += s*s; qo += s*s; }
    for(int k = 0; k < NISrc; k++){ double s = iS[k]->GetBinContent(i); qi += s*s; qo += s*s; }
    double jb = bJ->GetBinContent(i), ji = TMath::Max(fabs(iJu->GetBinContent(i)), fabs(iJd->GetBinContent(i)));
    double bs = sqrt(qb + jb*jb), is = sqrt(qi + ji*ji);
    double ru = bJu->GetBinContent(i) - iJu->GetBinContent(i);
    double rd = bJd->GetBinContent(i) - iJd->GetBinContent(i);
    double jc = TMath::Max(fabs(ru), fabs(rd)), ju = sqrt(jb*jb + ji*ji);
    double tot = sqrt(qo + jc*jc), st = sqrt(pow(bStat->GetBinContent(i), 2) + pow(iStat->GetBinContent(i), 2));
    double fr = vb/vi;
    bSys->SetBinContent(i, bs); iSys->SetBinContent(i, is);
    fF->SetBinContent(i, fr); fF->SetBinError(i, st*fr);
    fStat->SetBinContent(i, st); fSys->SetBinContent(i, tot); fJEU->SetBinContent(i, jc); fJEUu->SetBinContent(i, ju);
    fOth->SetBinContent(i, sqrt(qo)); fNoC->SetBinContent(i, sqrt(qo + ju*ju));
    printf("%3.0f-%-5.0f %10.4g %10.4g %9.5f %8.5f | %6.2f%% %6.2f%% | %6.2f%% %6.2f%% %6.2f%% %6.2f%% | %6.2f%%%s\n", ptEdges[i-1], ptEdges[i], vb, vi, fr, fMCt->GetBinContent(i),
           100*bs, 100*is, 100*jc, 100*ju, 100*sqrt(qo), 100*tot, 100*st, ptEdges[i-1] < ptReportMin - 1e-6 ? "  (underflow, not reported)" : "");
  }
  // the b-jet spectrum's own stat error is already its bin error; the inclusive one's too

  TFile *fo = TFile::Open(outRoot, "recreate");
  hB->Write(); bSys->Write(); hI->Write(); iSys->Write();
  fF->Write(); fMCt->Write(); fStat->Write(); fSys->Write(); fJEU->Write(); fJEUu->Write(); fOth->Write(); fNoC->Write();
  TNamed info("info", TString(Form("pp calo jets: inclusive b-jet fraction = b jets (correctMuonTagFrequency_caloJets_pp.C, unfolded N = 2) / "
                                   "inclusive jets (inclusiveSpectrumSys_caloJets_pp.C, unfolded N = 2), each per its own N_Z "
                                   "(incl %.1f, b %.1f); spectra drawn x N_Z^incl / %.1f in eta. Spectrum JEU fully correlated "
                                   "(same shifted PYTHIA reco spectra): fraction JEU = max over up/down of |d_b - d_incl| (signed shifts); all other sources and stat. in quadrature. "
                                   "Bins below %.0f GeV are unfolding underflow, not reported.", nzIncl, nzB, etaWidth, ptReportMin)).Data());
  info.Write(); fo->Close();

  auto shown = [&](int i){ return ptEdges[i-1] >= ptReportMin - 1e-6; };
  const double xlo = ptReportMin, xhi = ptEdges[NPt];
  auto drawBoxes = [&](TH1D *h, TH1D *rel, const char *hex, double alpha){
    for(int i = 1; i <= NPt; i++){
      double v = h->GetBinContent(i); if(v <= 0. || !shown(i)) continue;
      double e = rel->GetBinContent(i)*v;
      TBox *bx = new TBox(h->GetBinLowEdge(i), v - e, h->GetBinLowEdge(i+1), v + e);
      bx->SetFillColorAlpha(TColor::GetColor(hex), alpha); bx->SetLineColor(TColor::GetColor(hex)); bx->SetLineWidth(1);
      bx->Draw("l"); bx->Draw();
    }
  };
  auto blankHidden = [&](TH1D *h){ TH1D *c = (TH1D*) h->Clone(Form("%s_draw", h->GetName())); c->SetDirectory(nullptr);
    for(int i = 1; i <= NPt; i++) if(!shown(i)){ c->SetBinContent(i, -999.); c->SetBinError(i, 0.); } return c; };
  const char *hexB = "#D55E00", *hexI = "#000000", *hexIBox = "#999999";

  // ================================================== figure 1: both spectra
  {
    TCanvas *c = new TCanvas("cSpec", "", 700, 850);
    TPad *pU = new TPad("pU", "", 0, 0.30, 1, 1), *pL = new TPad("pL", "", 0, 0, 1, 0.30);
    pU->SetLeftMargin(0.16); pU->SetRightMargin(0.04); pU->SetTopMargin(0.05); pU->SetBottomMargin(0.015); pU->SetLogy();
    pL->SetLeftMargin(0.16); pL->SetRightMargin(0.04); pL->SetTopMargin(0.03); pL->SetBottomMargin(0.33);
    pU->Draw(); pL->Draw();

    pU->cd();
    double ymax = 0., ymin = 1e30;
    for(int i = 1; i <= NPt; i++) if(shown(i)) for(TH1D *h : {hB, hI}){ ymax = TMath::Max(ymax, h->GetBinContent(i)); ymin = TMath::Min(ymin, h->GetBinContent(i)); }
    TH1D *fr = (TH1D*) hI->Clone("frSpec"); fr->Reset(); fr->SetTitle("");
    fr->GetXaxis()->SetRangeUser(xlo, xhi); fr->SetMinimum(0.2*ymin); fr->SetMaximum(ymax*1e3);   // headroom for header + legend
    fr->GetXaxis()->SetLabelSize(0);
    fr->GetYaxis()->SetTitle("d^{2}#it{N}_{jet} / d#it{p}_{T} d#eta [GeV^{-1}]");
    fr->GetYaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetTitleOffset(1.45);
    fr->Draw("axis");
    drawBoxes(hI, iSys, hexIBox, 0.45);
    drawBoxes(hB, bSys, hexB, 0.30);
    TH1D *dI = blankHidden(hI), *dB = blankHidden(hB);
    styleH(dI, hexI, markFilledCircle, 1.2); styleH(dB, hexB, markFilledSquare, 1.2);
    dI->Draw("E1 X0 same"); dB->Draw("E1 X0 same");
    fr->Draw("axis same");

    TLatex t; t.SetNDC(); t.SetTextFont(42);
    t.SetTextSize(0.046); t.DrawLatex(0.20, 0.885, "pp 5.02 TeV, calo jets");
    t.SetTextSize(0.036); t.DrawLatex(0.20, 0.835, "anti-#it{k}_{T} #it{R} = 0.4, |#eta| < 1.6, unfolded");
    t.DrawLatex(0.20, 0.790, "Jet100-sample-equivalent yield");
    TLegend *leg = makeLegend(0.56, 0.56, 0.94, 0.75, 0.034);
    TH1D *legIb = (TH1D*) dI->Clone("legIb"); legIb->SetFillColorAlpha(TColor::GetColor(hexIBox), 0.45); legIb->SetLineColor(TColor::GetColor(hexIBox));
    TH1D *legBb = (TH1D*) dB->Clone("legBb"); legBb->SetFillColorAlpha(TColor::GetColor(hexB), 0.30);
    leg->AddEntry(dI, "inclusive jets (#it{N} = 2)", "pe");
    leg->AddEntry(legIb, "inclusive, syst.", "f");
    leg->AddEntry(dB, "b jets (#it{N} = 2)", "pe");
    leg->AddEntry(legBb, "b jets, syst.", "f");
    leg->Draw();

    pL->cd();
    TH1D *lBs = (TH1D*) bSys->Clone("lBs"), *lIs = (TH1D*) iSys->Clone("lIs"), *lBt = (TH1D*) bStat->Clone("lBt"), *lIt = (TH1D*) iStat->Clone("lIt");
    double m = 0.;
    for(TH1D *h : {lBs, lIs, lBt, lIt}){
      h->SetDirectory(nullptr); h->Scale(100.);
      for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.);
      for(int i = 1; i <= NPt; i++) if(shown(i)) m = TMath::Max(m, h->GetBinContent(i));
      h->GetXaxis()->SetRangeUser(xlo, xhi);
    }
    TH1D *rF = (TH1D*) lBs->Clone("rFspec2"); rF->Reset(); rF->SetTitle("");
    rF->SetMinimum(0.); rF->SetMaximum(1.55*m);
    rF->GetXaxis()->SetRangeUser(xlo, xhi);
    rF->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); rF->GetYaxis()->SetTitle("rel. unc. [%]");
    rF->GetXaxis()->SetTitleSize(0.12); rF->GetXaxis()->SetLabelSize(0.10); rF->GetXaxis()->SetTitleOffset(1.15);
    rF->GetYaxis()->SetTitleSize(0.105); rF->GetYaxis()->SetLabelSize(0.095); rF->GetYaxis()->SetTitleOffset(0.66); rF->GetYaxis()->SetNdivisions(505);
    rF->Draw("axis");
    styleLine(lIs, hexI, 3); styleLine(lBs, hexB, 3);
    styleLine(lIt, hexI, 2); lIt->SetLineStyle(2); styleLine(lBt, hexB, 2); lBt->SetLineStyle(2);
    lIs->Draw("hist same"); lBs->Draw("hist same"); lIt->Draw("hist same"); lBt->Draw("hist same");
    TLegend *l2 = makeLegend(0.19, 0.76, 0.94, 0.95, 0.075); l2->SetNColumns(4);
    l2->AddEntry(lIs, "incl. syst.", "l"); l2->AddEntry(lIt, "incl. stat.", "l");
    l2->AddEntry(lBs, "b syst.", "l");      l2->AddEntry(lBt, "b stat.", "l");
    l2->Draw();

    c->SaveAs(Form("%s/bAndInclusiveJetSpectra_caloJets_pp.pdf", figDir));
  }

  // ================================================== figure 2: b-jet fraction
  {
    // drawn in percent; the output file keeps the fraction
    TH1D *fPct = (TH1D*) fF->Clone("fPct"), *mcPct = (TH1D*) fMCt->Clone("mcPct");
    for(TH1D *h : {fPct, mcPct}){ h->SetDirectory(nullptr); h->Scale(100.); }
    TCanvas *c = new TCanvas("cFrac", "", 700, 850);
    TPad *pU = new TPad("pU2", "", 0, 0.30, 1, 1), *pL = new TPad("pL2", "", 0, 0, 1, 0.30);
    pU->SetLeftMargin(0.16); pU->SetRightMargin(0.04); pU->SetTopMargin(0.05); pU->SetBottomMargin(0.015);
    pL->SetLeftMargin(0.16); pL->SetRightMargin(0.04); pL->SetTopMargin(0.03); pL->SetBottomMargin(0.33);
    pU->Draw(); pL->Draw();

    pU->cd();
    double ymax = 0.;
    for(int i = 1; i <= NPt; i++) if(shown(i)) ymax = TMath::Max(ymax, TMath::Max(fPct->GetBinContent(i)*(1. + fSys->GetBinContent(i)), mcPct->GetBinContent(i)));
    TH1D *fr = (TH1D*) fPct->Clone("frFrac"); fr->Reset(); fr->SetTitle("");
    fr->GetXaxis()->SetRangeUser(xlo, xhi); fr->SetMinimum(0.); fr->SetMaximum(1.6*ymax);   // headroom for header + legend
    fr->GetXaxis()->SetLabelSize(0);
    fr->GetYaxis()->SetTitle("b jets / inclusive jets [%]");
    fr->GetYaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetLabelSize(0.042); fr->GetYaxis()->SetTitleOffset(1.45);
    fr->Draw("axis");
    drawBoxes(fPct, fSys, hexB, 0.30);
    TH1D *dF = blankHidden(fPct); styleH(dF, "#000000", markFilledCircle, 1.2);
    TH1D *dMC = (TH1D*) mcPct->Clone("dMC"); dMC->SetDirectory(nullptr);
    for(int i = 0; i <= NPt+1; i++){ dMC->SetBinError(i, 0.); if(i < 1 || i > NPt || !shown(i)) dMC->SetBinContent(i, -999.); }
    dMC->GetXaxis()->SetRangeUser(xlo, xhi);
    styleLine(dMC, "#0072B2", 3); dMC->SetLineStyle(3);
    dMC->Draw("hist same");
    dF->Draw("E1 X0 same");
    fr->Draw("axis same");

    TLatex t; t.SetNDC(); t.SetTextFont(42);
    t.SetTextSize(0.046); t.DrawLatex(0.20, 0.885, "pp 5.02 TeV, calo jets");
    t.SetTextSize(0.036); t.DrawLatex(0.20, 0.835, "anti-#it{k}_{T} #it{R} = 0.4, |#eta| < 1.6, unfolded");
    t.DrawLatex(0.20, 0.790, "b jets / inclusive jets, both #it{N} = 2");
    TLegend *leg = makeLegend(0.20, 0.60, 0.60, 0.76, 0.034);
    TH1D *legBx = (TH1D*) dF->Clone("legFbx"); legBx->SetFillColorAlpha(TColor::GetColor(hexB), 0.30); legBx->SetLineColor(TColor::GetColor(hexB));
    leg->AddEntry(dF, "inclusive b-jet fraction (stat.)", "pe");
    leg->AddEntry(legBx, "syst. (JEU correlated)", "f");
    leg->AddEntry(dMC, "PYTHIA, gen-level truth", "l");
    leg->Draw();

    pL->cd();
    TH1D *lT = (TH1D*) fSys->Clone("lTf"), *lN = (TH1D*) fNoC->Clone("lNf"), *lJ = (TH1D*) fJEU->Clone("lJf"), *lS = (TH1D*) fStat->Clone("lSf");
    double m = 0.;
    for(TH1D *h : {lT, lN, lJ, lS}){
      h->SetDirectory(nullptr); h->Scale(100.);
      for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.);
      for(int i = 1; i <= NPt; i++) if(shown(i)) m = TMath::Max(m, h->GetBinContent(i));
      h->GetXaxis()->SetRangeUser(xlo, xhi);
    }
    TH1D *rF = (TH1D*) lT->Clone("rFfrac"); rF->Reset(); rF->SetTitle("");
    rF->SetMinimum(0.); rF->SetMaximum(1.55*m);
    rF->GetXaxis()->SetRangeUser(xlo, xhi);
    rF->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); rF->GetYaxis()->SetTitle("rel. unc. [%]");
    rF->GetXaxis()->SetTitleSize(0.12); rF->GetXaxis()->SetLabelSize(0.10); rF->GetXaxis()->SetTitleOffset(1.15);
    rF->GetYaxis()->SetTitleSize(0.105); rF->GetYaxis()->SetLabelSize(0.095); rF->GetYaxis()->SetTitleOffset(0.66); rF->GetYaxis()->SetNdivisions(505);
    rF->Draw("axis");
    styleLine(lT, hexB, 3); styleLine(lN, "#999999", 3); lN->SetLineStyle(7);
    styleLine(lJ, "#0072B2", 3); styleLine(lS, "#000000", 2); lS->SetLineStyle(2);
    lN->Draw("hist same"); lT->Draw("hist same"); lJ->Draw("hist same"); lS->Draw("hist same");
    TLegend *l2 = makeLegend(0.19, 0.74, 0.94, 0.95, 0.072); l2->SetNColumns(2);
    l2->AddEntry(lT, "total syst.", "l");          l2->AddEntry(lN, "total, JEU not cancelled", "l");
    l2->AddEntry(lJ, "JEU (correlated)", "l");     l2->AddEntry(lS, "stat.", "l");
    l2->Draw();

    c->SaveAs(Form("%s/bJetFraction_caloJets_pp.pdf", figDir));
  }
}
