// Correct the unfolded pp CALO-jet b-jet spectrum for the tagging muon's
// reconstruction + tight-ID efficiency.
//
// The efficiency is the Z -> mu mu tag-and-probe measurement of
// src/muonReconstructionEfficiencyCalculator/muonTagAndProbeEfficiency.C on
// rootFiles/scanningOuput/pp/pp_HighEGJet_tnp.root (not muon-triggered, so no
// trigger bias), step "all" (any muon object -> analysis tight ID), |eta| < 2.0,
// in probe-pT bins. That macro prints its results and writes no ROOT file, so
// the table below is copied from its output (run 2026-10-04); re-run it and
// update the table if the measurement changes.
//
// The correction uses the AVERAGE of the probe-pT bins (unweighted mean):
//   corrected = unfolded / eps_avg
// Its uncertainty is a new systematic source, constant in jet pT:
//   pT dependence   max |eps_bin - eps_avg| / eps_avg   (using the average
//                   instead of the pT-dependent value)
//   stat            error on the mean, sqrt(sum sigma_bin^2) / N / eps_avg
//   method          background-subtraction systematic of the integrated value
//                   (largest shift of the sideband / no-subtraction variations)
// added in quadrature, and in quadrature with the unfolded spectrum's existing
// systematics.
//
// CAVEATS: tag-and-probe measures isolated muons from Z decays; the tagging
// muons here sit inside jets, where the reconstruction and ID can behave
// differently. The analysis MC value 0.9708 (tight | gen) also includes a step
// tag-and-probe cannot see (muons that leave no muon object), so the two are
// not like-for-like (see the tag-and-probe macro's header).
//
// Input   rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_caloJets_pp.root
// Output  rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_muRecoCorr_caloJets_pp.root
//           h_bJetPt_corrected          unfolded / eps_avg, stat error
//           h_bJetPt_corrected_sysAbs   same, error = total systematic
//           h_corr_statRel, h_corr_sysRel
//           h_corr_sysRel_bPurity, _spectrumJEU, _iterations, _muonRecoEff
//           h_muonRecoEff_vsProbePt     the tag-and-probe table, for reference
//         figures/bJetSpectra/bJetSpectrum_muRecoCorr_caloJets_pp.pdf
//         figures/systematics/bJetSpectrum_muRecoCorr_sysBreakdown_caloJets_pp.pdf
//
// Usage: root -l -b -q correctMuonRecoEff_caloJets_pp.C
// Run from: src/calculateBJetsPerZ/

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TColor.h"
#include "TStyle.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/ratioPanel.h"
#include "../../headers/functions/caloJetBPurityFit_pp.h"   // ptEdges, NPt, ptReportMin

// ---- tag-and-probe table (muonTagAndProbeEfficiency.C, pp_HighEGJet_tnp.root, step "all") ----
const int NTnP = 7;
const double tnpEdges[NTnP+1] = {15, 20, 25, 30, 40, 50, 70, 100};
const double tnpEff[NTnP]     = {0.964, 0.970, 0.965, 0.962, 0.958, 0.973, 0.972};
const double tnpErrUp[NTnP]   = {0.008, 0.005, 0.003, 0.002, 0.002, 0.002, 0.004};
const double tnpErrDn[NTnP]   = {0.008, 0.005, 0.004, 0.002, 0.002, 0.002, 0.004};
const double tnpMethodSys     = 0.0030;   // integrated "all": background-method systematic
const double tnpIntegrated    = 0.9638;   // integrated "all", for reference only

const char *inPath  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_caloJets_pp.root";
const char *outRoot = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/bJetSpectrum_muRecoCorr_caloJets_pp.root";
const char *outFig  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/bJetSpectra/bJetSpectrum_muRecoCorr_caloJets_pp.pdf";
const char *outFigS = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics/bJetSpectrum_muRecoCorr_sysBreakdown_caloJets_pp.pdf";

bool reported(int i){ return ptEdges[i-1] >= ptReportMin - 1e-6; }

TH1D* bookRel(const char *name, const char *title)
{
  TH1D *h = new TH1D(name, Form("%s;#it{p}_{T}^{jet} [GeV];relative uncertainty", title), NPt, ptEdges);
  h->SetDirectory(nullptr); return h;
}

void correctMuonRecoEff_caloJets_pp()
{
  initPlotStyle();
  gSystem->mkdir(gSystem->DirName(outRoot), kTRUE);
  gSystem->mkdir(gSystem->DirName(outFig), kTRUE);
  gSystem->mkdir(gSystem->DirName(outFigS), kTRUE);

  // ---- average efficiency and its uncertainty --------------------------------
  double sum = 0., s2 = 0.;
  for(int b = 0; b < NTnP; b++){ sum += tnpEff[b]; s2 += pow(0.5*(tnpErrUp[b] + tnpErrDn[b]), 2); }
  const double effAvg = sum/NTnP;
  double spread = 0.; for(int b = 0; b < NTnP; b++) spread = TMath::Max(spread, fabs(tnpEff[b] - effAvg));
  const double relSpread = spread/effAvg, relStat = sqrt(s2)/NTnP/effAvg, relMethod = tnpMethodSys/effAvg;
  const double relEff = sqrt(relSpread*relSpread + relStat*relStat + relMethod*relMethod);
  printf("\nmuon reco+ID efficiency (tag-and-probe, average of %d probe-pT bins): %.4f\n", NTnP, effAvg);
  printf("  uncertainty: pT dependence %.2f%%, stat %.2f%%, method %.2f%%  ->  %.2f%%\n",
         100*relSpread, 100*relStat, 100*relMethod, 100*relEff);
  printf("  (integrated value %.4f, for reference)\n", tnpIntegrated);

  // ---- inputs -----------------------------------------------------------------
  TFile *fI = TFile::Open(inPath);
  if(!fI || fI->IsZombie()){ printf("ERROR: cannot open %s -- run unfoldBJetSpectrum_caloJets_pp.C\n", inPath); return; }
  TH1D *hU = nullptr, *rPur = nullptr, *rJEU = nullptr, *rIter = nullptr, *rStat = nullptr;
  fI->GetObject("h_bJetPt_unfolded", hU);
  fI->GetObject("h_unf_sysRel_bPurity", rPur);
  fI->GetObject("h_unf_sysRel_spectrumJEU", rJEU);
  fI->GetObject("h_unf_sysRel_iterations", rIter);
  fI->GetObject("h_unf_statRel", rStat);
  if(!hU || !rPur || !rJEU || !rIter || !rStat){ printf("ERROR: missing histograms in %s\n", inPath); return; }

  // ---- correct ----------------------------------------------------------------
  TH1D *hC = (TH1D*) hU->Clone("h_bJetPt_corrected"); hC->SetDirectory(nullptr);
  hC->Scale(1./effAvg);   // relative stat errors unchanged
  hC->SetTitle("b jets, unfolded, corrected for muon reco efficiency, stat.;#it{p}_{T}^{jet} [GeV];d#it{N}/d#it{p}_{T} [GeV^{-1}]");
  TH1D *hCsys = (TH1D*) hC->Clone("h_bJetPt_corrected_sysAbs"); hCsys->SetDirectory(nullptr);
  hCsys->SetTitle("b jets, corrected, error = total systematic;#it{p}_{T}^{jet} [GeV];d#it{N}/d#it{p}_{T} [GeV^{-1}]");

  TH1D *cPur  = (TH1D*) rPur->Clone("h_corr_sysRel_bPurity");      cPur->SetDirectory(nullptr);
  TH1D *cJEU  = (TH1D*) rJEU->Clone("h_corr_sysRel_spectrumJEU");  cJEU->SetDirectory(nullptr);
  TH1D *cIter = (TH1D*) rIter->Clone("h_corr_sysRel_iterations");  cIter->SetDirectory(nullptr);
  TH1D *cStat = (TH1D*) rStat->Clone("h_corr_statRel");            cStat->SetDirectory(nullptr);
  TH1D *cEff  = bookRel("h_corr_sysRel_muonRecoEff", "muon reco efficiency (tag-and-probe average)");
  TH1D *cTot  = bookRel("h_corr_sysRel", "total systematic");

  printf("\n%-9s %11s %11s | %8s %8s %8s %8s | %8s %8s\n", "jet pT", "unfolded", "corrected",
         "purity", "specJEU", "iter", "muEff", "TOTAL", "stat");
  for(int i = 1; i <= NPt; i++){
    double t = sqrt(pow(rPur->GetBinContent(i), 2) + pow(rJEU->GetBinContent(i), 2) + pow(rIter->GetBinContent(i), 2) + relEff*relEff);
    if(hU->GetBinContent(i) <= 0.) t = 0.;
    cEff->SetBinContent(i, hU->GetBinContent(i) > 0. ? relEff : 0.);
    cTot->SetBinContent(i, t);
    hCsys->SetBinError(i, t*hC->GetBinContent(i));
    printf("%3.0f-%-5.0f %11.4g %11.4g | %7.2f%% %7.2f%% %7.2f%% %7.2f%% | %7.2f%% %7.2f%%%s\n", ptEdges[i-1], ptEdges[i],
           hU->GetBinContent(i), hC->GetBinContent(i), 100*rPur->GetBinContent(i), 100*rJEU->GetBinContent(i),
           100*rIter->GetBinContent(i), 100*cEff->GetBinContent(i), 100*t, 100*rStat->GetBinContent(i),
           reported(i) ? "" : "  (underflow)");
  }

  TH1D *hT = new TH1D("h_muonRecoEff_vsProbePt", "muon reco+ID efficiency, Z tag-and-probe, step all;probe #it{p}_{T} [GeV];efficiency", NTnP, tnpEdges);
  hT->SetDirectory(nullptr);
  for(int b = 0; b < NTnP; b++){ hT->SetBinContent(b+1, tnpEff[b]); hT->SetBinError(b+1, 0.5*(tnpErrUp[b] + tnpErrDn[b])); }

  TFile *fo = TFile::Open(outRoot, "recreate");
  hC->Write(); hCsys->Write(); cStat->Write(); cTot->Write(); cPur->Write(); cJEU->Write(); cIter->Write(); cEff->Write(); hT->Write();
  TNamed info("info", Form("pp calo b jets: unfolded spectrum / muon reco+ID efficiency %.4f (average of %d tag-and-probe probe-pT "
                            "bins, pp_HighEGJet_tnp.root, step all, |eta|<2). Efficiency uncertainty %.2f%% (pT dependence %.2f%%, "
                            "stat %.2f%%, method %.2f%%), added in quadrature to purity, spectrum JEU and iteration systematics.",
                            effAvg, NTnP, 100*relEff, 100*relSpread, 100*relStat, 100*relMethod));
  info.Write();
  fo->Close();
  printf("\nwritten %s\n", outRoot);

  int colSys = TColor::GetColor("#D55E00");

  // ---- figure: corrected spectrum ---------------------------------------------
  {
    TCanvas *c = new TCanvas("c", "", 700, 800);
    TPad *pTop, *pBot; splitPads(pTop, pBot);
    pTop->cd(); pTop->SetLogy();   // log: falls by ~3 orders of magnitude, no negative bins
    TH1D *fr = (TH1D*) hC->Clone("fr"); fr->Reset(); fr->SetTitle("");
    double ymax = 0., ymin = 1e30;
    for(int i = 1; i <= NPt; i++) if(reported(i)){ ymax = TMath::Max(ymax, hC->GetBinContent(i)); ymin = TMath::Min(ymin, hC->GetBinContent(i)); }
    fr->SetMinimum(0.3*ymin); fr->SetMaximum(30.*ymax);
    fr->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); fr->GetXaxis()->SetLabelSize(0);
    fr->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} [GeV^{-1}]"); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
    fr->Draw("AXIS");
    for(int i = 1; i <= NPt; i++){
      double v = hC->GetBinContent(i); if(v <= 0. || !reported(i)) continue;
      TBox *bx = new TBox(ptEdges[i-1], v - hCsys->GetBinError(i), ptEdges[i], v + hCsys->GetBinError(i));
      bx->SetFillColorAlpha(colSys, 0.30); bx->SetLineColor(colSys); bx->Draw("l"); bx->Draw();
    }
    TH1D *d = (TH1D*) hC->Clone("d");
    for(int i = 1; i <= NPt; i++) if(!reported(i)){ d->SetBinContent(i, -999.); d->SetBinError(i, 0.); }
    styleH(d, "#000000", markFilledCircle, 1.2);
    d->Draw("E1 X0 same");
    TLegend *leg = makeLegend(0.25, 0.05, 0.60, 0.20, 0.036);
    leg->AddEntry(d, "b jets, stat.", "lp");
    TBox *lb = new TBox(); lb->SetFillColorAlpha(colSys, 0.30); lb->SetLineColor(colSys);
    leg->AddEntry(lb, "total systematic", "f");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
    la.DrawLatex(0.21, 0.84, "pp 5.02 TeV, calo jets: b jets");
    la.SetTextSize(0.034);
    la.DrawLatex(0.21, 0.785, "unfolded, corrected for muon reco efficiency");
    la.DrawLatex(0.21, 0.74, Form("(#varepsilon_{#mu} = %.3f, tag-and-probe average)", effAvg));

    pBot->cd();
    TH1D *rS = (TH1D*) cTot->Clone("rS"), *rT = (TH1D*) cStat->Clone("rT");
    for(TH1D *h : {rS, rT}){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); h->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); }
    styleLine(rS, "#D55E00", 3); styleLine(rT, "#0072B2", 3); rT->SetLineStyle(2);
    double m = 0.; for(int i = 1; i <= NPt; i++) if(reported(i)) m = TMath::Max(m, TMath::Max(rS->GetBinContent(i), rT->GetBinContent(i)));
    TH1D *rF = (TH1D*) rS->Clone("rF"); rF->Reset();
    rF->SetMinimum(0.); rF->SetMaximum(1.7*m);   // headroom for the legend row
    rF->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "rel. unc. [%]");
    rF->Draw("AXIS");
    rS->Draw("HIST same"); rT->Draw("HIST same");
    TLegend *legB = makeLegend(0.19, 0.78, 0.95, 0.97, 0.075);
    legB->SetNColumns(2);
    legB->AddEntry(rS, "total syst.", "l"); legB->AddEntry(rT, "stat.", "l");
    legB->Draw();
    c->SaveAs(outFig);
    delete c;
  }

  // ---- figure: systematic breakdown -------------------------------------------
  {
    TCanvas *c = new TCanvas("cS", "", 800, 700);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    TH1D *lT = (TH1D*) cTot->Clone("lT"), *lP = (TH1D*) cPur->Clone("lP"), *lJ = (TH1D*) cJEU->Clone("lJ");
    TH1D *lI = (TH1D*) cIter->Clone("lI"), *lE = (TH1D*) cEff->Clone("lE");
    for(TH1D *h : {lT, lP, lJ, lI, lE}){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); h->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); }
    double ymax = 0.; for(int i = 1; i <= NPt; i++) if(reported(i)) ymax = TMath::Max(ymax, lT->GetBinContent(i));
    TH1F *fr = c->DrawFrame(ptReportMin, 0., ptEdges[NPt], 1.9*ymax);   // headroom for header + legend
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle("relative systematic uncertainty [%]");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
    styleLine(lP, "#009E73", 3); lP->SetLineStyle(2);
    styleLine(lJ, "#0072B2", 3); lJ->SetLineStyle(3);
    styleLine(lI, "#E69F00", 3); lI->SetLineStyle(7);
    styleLine(lE, "#CC79A7", 3); lE->SetLineStyle(9);
    styleLine(lT, "#000000", 4);
    lP->Draw("HIST same"); lJ->Draw("HIST same"); lI->Draw("HIST same"); lE->Draw("HIST same"); lT->Draw("HIST same");
    TLegend *leg = makeLegend(0.40, 0.60, 0.95, 0.86, 0.032);
    leg->AddEntry(lP, "b purity (propagated through unfolding)", "l");
    leg->AddEntry(lJ, "spectrum JEU (gen level)", "l");
    leg->AddEntry(lI, "iterations", "l");
    leg->AddEntry(lE, "muon reco efficiency", "l");
    leg->AddEntry(lT, "total (quadrature)", "l");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
    la.DrawLatex(0.17, 0.885, "pp 5.02 TeV, calo jets: corrected b-jet spectrum systematics");
    savePdfTight(c, outFigS);
    delete c;
  }
  printf("written %s\nwritten %s\n", outFig, outFigS);
}
