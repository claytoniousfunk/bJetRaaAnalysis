// Unfold the pp CALO-jet b-jet spectrum (muon-tagged x b purity) with the
// muon-tagged b-jet response, and propagate its uncertainties.
//
// INPUT  rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_caloJets_pp.root
//        (constructBJetSpectrum_caloJets_pp.C): h_bJetPt in dN/dpT with stat
//        errors, and the relative systematics per source.
//
// RESPONSE  h_matchedRecoJetPt_genJetPt_bJets_muTagged (X reco, Y gen) from the
//        PYTHIA calo response scans with PF-matched flavour, even + odd halves
//        added (2026-10-01). Muon-tagged = a reco muon in the jet, no trigger
//        requirement (the data spectrum is muon-triggered; the trigger
//        efficiency is a separate correction). Rebinned from the native 5 GeV
//        bins to ptEdges. Reco bins below recoFloor are zeroed and the truth is
//        taken from the untruncated matrix, so jets reconstructed below the
//        floor become misses -- the data have no jets below 70 GeV, and this is
//        the treatment unfoldClosureTest.C applied when measuring the optimal
//        iteration count. Matched jets only: the truth is the gen projection of
//        the matched matrix (no b-jet miss histogram exists), so the unfolded
//        spectrum is per matched muon-tagged b jet; reconstruction and tagging
//        efficiencies are separate corrections.
//
// ITERATIONS  nIterOpt = 2, the minimum-MSE iteration of
//        unfoldClosureTest.C("pp",15,"weighted","tiltScan:0.5",80.,500.,true,1,
//          <PFflavor even>,<PFflavor odd>,"_bMuTagPFflavor",-1,
//          "h_matchedRecoJetPt_genJetPt_bJets_muTagged",
//          "h_matchedRecoJetPt_genJetPt_bJets_muTagged",true,70.,
//          "50,60,70,80,100,120,150,200,300,500")
//        (2026-10-04: MSE 0.0110 at N = 2; flat to within ~15% over N = 2-15,
//        N = 1 under-converged, max|r - 1| 10.7%). The iteration systematic is
//        max(|U(N-1) - U(N)|, |U(N+1) - U(N)|) / U(N).
//
// UNCERTAINTIES on the unfolded spectrum
//   stat           the measured spectrum's stat errors (counts + purity fit)
//                  propagated through the unfolding covariance (RooUnfold).
//   b purity       each purity source propagated through the unfolding: the
//                  measured spectrum is scaled by (1 +/- rel_s) in every bin,
//                  unfolded, and the larger relative change is that source's
//                  systematic. Sources added in quadrature.
//   spectrum JEU   already a GEN-level uncertainty (unfoldPYTHIAJetPt_JEUShift.C
//                  unfolds the shifted spectra), so applied directly, not
//                  re-unfolded.
//   iterations     as above.
//   total          quadrature sum of the three.
//
// Underflow bins (50-60, 60-70, 70-80 GeV) are unfolded with the rest but not
// reported; figures show 80-500 GeV only (the matrix figure shows everything).
//
// Output  rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_caloJets_pp.root
//           h_bJetPt_measured            measured b-jet dN/dpT (input), stat
//           h_bJetPt_unfolded            unfolded dN/dpT at nIterOpt, stat
//           h_bJetPt_unfolded_N<n>       unfolded at N = nIterOpt-1, +1
//           h_bJetPt_unfolded_sysAbs     unfolded, error = total systematic
//           h_unf_statRel, h_unf_sysRel  relative stat / total systematic
//           h_unf_sysRel_bPurity, h_unf_sysRel_spectrumJEU, h_unf_sysRel_iterations
//           h_unf_sysRel_<purity source> each purity source after unfolding
//           h_response                   the rebinned, floor-truncated response
//         figures/bJetSpectra/unfold_caloJets_pp_response.pdf
//         figures/bJetSpectra/unfold_caloJets_pp_measuredVsUnfolded.pdf
//         figures/bJetSpectra/unfoldedBJetSpectrum_caloJets_pp.pdf
//         figures/systematics/unfoldedBJetSpectrum_sysBreakdown_caloJets_pp.pdf
//
// Usage: root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");' \
//                         unfoldBJetSpectrum_caloJets_pp.C
// Run from: src/calculateBJetsPerZ/

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TColor.h"
#include "TStyle.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/ratioPanel.h"
#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/functions/caloJetBPurityFit_pp.h"   // ptEdges, NPt, ptReportMin

const char *specPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_caloJets_pp.root";
const char *respDir  = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIA/";
const char *respE    = "PYTHIA_DiJet_response_caloJets_PFflavor_manualJEC_evenEvents_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-10-1.root";
const char *respO    = "PYTHIA_DiJet_response_caloJets_PFflavor_manualJEC_oddEvents_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-10-1.root";
const char *respName = "h_matchedRecoJetPt_genJetPt_bJets_muTagged";

const double recoFloor = 70.;
const int    nIterOpt  = 2;

const char *outRoot   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_caloJets_pp.root";
const char *figDir    = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/bJetSpectra";
const char *figSysDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics";

const int NSrc = 6;
const char *srcKey[NSrc]   = {"bGS", "cMult", "fitRangeLow", "fitRangeHigh", "JEU", "JER"};

// native fixed-bin 2D -> ptEdges x ptEdges, by bin centre (as unfoldClosureTest's rebinToVarBins)
TH2D* rebin2D(TH2D *h, const char *name)
{
  TH2D *o = new TH2D(name, h->GetTitle(), NPt, ptEdges, NPt, ptEdges); o->SetDirectory(nullptr); o->Sumw2();
  for(int ix = 1; ix <= h->GetNbinsX(); ix++){
    int jx = o->GetXaxis()->FindBin(h->GetXaxis()->GetBinCenter(ix)); if(jx < 1 || jx > NPt) continue;
    for(int iy = 1; iy <= h->GetNbinsY(); iy++){
      int jy = o->GetYaxis()->FindBin(h->GetYaxis()->GetBinCenter(iy)); if(jy < 1 || jy > NPt) continue;
      o->SetBinContent(jx, jy, o->GetBinContent(jx, jy) + h->GetBinContent(ix, iy));
      o->SetBinError(jx, jy, sqrt(pow(o->GetBinError(jx, jy), 2) + pow(h->GetBinError(ix, iy), 2)));
    }
  }
  return o;
}

TH1D* bookRel(const char *name, const char *title)
{
  TH1D *h = new TH1D(name, Form("%s;#it{p}_{T}^{jet} [GeV];relative uncertainty", title), NPt, ptEdges);
  h->SetDirectory(nullptr); return h;
}

bool reported(int i){ return ptEdges[i-1] >= ptReportMin - 1e-6; }

// dN/dpT <-> counts per bin
TH1D* toCounts(TH1D *h, const char *name){ TH1D *c = (TH1D*) h->Clone(name); c->SetDirectory(nullptr);
  for(int i = 1; i <= c->GetNbinsX(); i++){ double w = c->GetBinWidth(i); c->SetBinContent(i, c->GetBinContent(i)*w); c->SetBinError(i, c->GetBinError(i)*w); } return c; }
TH1D* toDensity(TH1D *h, const char *name){ TH1D *c = (TH1D*) h->Clone(name); c->SetDirectory(nullptr); divideByBinwidth(c); return c; }

void unfoldBJetSpectrum_caloJets_pp()
{
  initPlotStyle();
  gSystem->mkdir(gSystem->DirName(outRoot), kTRUE);
  gSystem->mkdir(figDir, kTRUE); gSystem->mkdir(figSysDir, kTRUE);

  // ---- inputs -----------------------------------------------------------------
  TFile *fS = TFile::Open(specPath);
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open %s -- run constructBJetSpectrum_caloJets_pp.C\n", specPath); return; }
  TH1D *hMeasD = nullptr, *hRelJEU = nullptr, *hRelSrc[NSrc];
  fS->GetObject("h_bJetPt", hMeasD);
  fS->GetObject("h_bJetPt_sysRel_spectrumJEU", hRelJEU);
  for(int s = 0; s < NSrc; s++) fS->GetObject(Form("h_bJetPt_sysRel_%s", srcKey[s]), hRelSrc[s]);
  if(!hMeasD || !hRelJEU){ printf("ERROR: missing input histograms in %s\n", specPath); return; }
  for(int s = 0; s < NSrc; s++) if(!hRelSrc[s]){ printf("ERROR: missing h_bJetPt_sysRel_%s\n", srcKey[s]); return; }
  if(hMeasD->GetNbinsX() != NPt){ printf("ERROR: spectrum has %d bins, expected %d\n", hMeasD->GetNbinsX(), NPt); return; }

  TFile *fE = TFile::Open(Form("%s%s", respDir, respE)), *fO = TFile::Open(Form("%s%s", respDir, respO));
  TH2D *rE = nullptr, *rO = nullptr;
  if(fE) fE->GetObject(respName, rE);
  if(fO) fO->GetObject(respName, rO);
  if(!rE || !rO){ printf("ERROR: %s missing in the response files\n", respName); return; }
  TH2D *rNative = (TH2D*) rE->Clone("rNative"); rNative->SetDirectory(nullptr); rNative->Add(rO);

  // ---- response ---------------------------------------------------------------
  TH2D *resp2D = rebin2D(rNative, "h_response_full");
  TH1D *truthMC = resp2D->ProjectionY("truthMC", 1, NPt); truthMC->SetDirectory(nullptr);   // before the floor
  TH2D *respT = (TH2D*) resp2D->Clone("h_response"); respT->SetDirectory(nullptr);
  for(int ix = 1; ix <= NPt; ix++)
    if(respT->GetXaxis()->GetBinUpEdge(ix) <= recoFloor + 1e-6)
      for(int iy = 0; iy <= NPt+1; iy++){ respT->SetBinContent(ix, iy, 0.); respT->SetBinError(ix, iy, 0.); }
  TH1D *measMC = respT->ProjectionX("measMC", 1, NPt); measMC->SetDirectory(nullptr);
  RooUnfoldResponse resp(measMC, truthMC, respT, "resp_bMuTag_pp_calo", "PYTHIA calo muon-tagged b response");

  // ---- unfold -----------------------------------------------------------------
  TH1D *hMeasC = toCounts(hMeasD, "hMeasC");
  auto unfoldN = [&](TH1D *measCounts, int n, const char *name){
    RooUnfoldBayes u(&resp, measCounts, n); u.SetVerbose(0);
    TH1D *h = (TH1D*) u.Hunfold(); h->SetName(name); h->SetDirectory(nullptr);
    return toDensity(h, name);
  };
  TH1D *hU  = unfoldN(hMeasC, nIterOpt,     "h_bJetPt_unfolded");
  TH1D *hUm = unfoldN(hMeasC, nIterOpt - 1, Form("h_bJetPt_unfolded_N%d", nIterOpt - 1));
  TH1D *hUp = unfoldN(hMeasC, nIterOpt + 1, Form("h_bJetPt_unfolded_N%d", nIterOpt + 1));

  // ---- systematics on the unfolded spectrum -----------------------------------
  TH1D *uSrc[NSrc];
  for(int s = 0; s < NSrc; s++){
    uSrc[s] = bookRel(Form("h_unf_sysRel_%s", srcKey[s]), Form("unfolded, b-purity source %s (propagated)", srcKey[s]));
    TH1D *mUp = (TH1D*) hMeasC->Clone(Form("mUp%d", s)), *mDn = (TH1D*) hMeasC->Clone(Form("mDn%d", s));
    for(int i = 1; i <= NPt; i++){
      double r = hRelSrc[s]->GetBinContent(i);
      mUp->SetBinContent(i, hMeasC->GetBinContent(i)*(1. + r));
      mDn->SetBinContent(i, hMeasC->GetBinContent(i)*(1. - r));
    }
    TH1D *uUp = unfoldN(mUp, nIterOpt, Form("uUp%d", s)), *uDn = unfoldN(mDn, nIterOpt, Form("uDn%d", s));
    for(int i = 1; i <= NPt; i++){
      double n = hU->GetBinContent(i);
      if(n > 0.) uSrc[s]->SetBinContent(i, TMath::Max(fabs(uUp->GetBinContent(i)/n - 1.), fabs(uDn->GetBinContent(i)/n - 1.)));
    }
  }
  TH1D *uPur  = bookRel("h_unf_sysRel_bPurity",      "unfolded, b purity (all sources, propagated)");
  TH1D *uJEU  = bookRel("h_unf_sysRel_spectrumJEU",  "unfolded, spectrum JEU (gen level)");
  TH1D *uIter = bookRel("h_unf_sysRel_iterations",   Form("unfolded, iterations N = %d +/- 1", nIterOpt));
  TH1D *uTot  = bookRel("h_unf_sysRel",              "unfolded, total systematic");
  TH1D *uStat = bookRel("h_unf_statRel",             "unfolded, stat. (propagated through the unfolding)");
  TH1D *hUsys = (TH1D*) hU->Clone("h_bJetPt_unfolded_sysAbs"); hUsys->SetDirectory(nullptr);
  hUsys->SetTitle("unfolded b jets, error = total systematic;#it{p}_{T}^{jet} [GeV];d#it{N}/d#it{p}_{T} [GeV^{-1}]");

  printf("\n%-9s %10s %10s %8s | %8s %8s %8s %8s %8s %8s | %8s %8s %8s %8s\n", "jet pT", "measured", "unfolded", "U/M",
         "bGS", "cMult", "fitLow", "fitHigh", "JEU(p)", "JER", "purity", "specJEU", "iter", "TOTAL");
  for(int i = 1; i <= NPt; i++){
    double n = hU->GetBinContent(i);
    double q = 0.; for(int s = 0; s < NSrc; s++) q += pow(uSrc[s]->GetBinContent(i), 2);
    double rp = sqrt(q), rj = hRelJEU->GetBinContent(i);
    double ri = n > 0. ? TMath::Max(fabs(hUm->GetBinContent(i)/n - 1.), fabs(hUp->GetBinContent(i)/n - 1.)) : 0.;
    double rt = sqrt(rp*rp + rj*rj + ri*ri);
    uPur->SetBinContent(i, rp); uJEU->SetBinContent(i, rj); uIter->SetBinContent(i, ri); uTot->SetBinContent(i, rt);
    uStat->SetBinContent(i, n > 0. ? hU->GetBinError(i)/n : 0.);
    hUsys->SetBinError(i, rt*n);
    double m = hMeasD->GetBinContent(i);
    printf("%3.0f-%-5.0f %10.4g %10.4g %8.3f |", ptEdges[i-1], ptEdges[i], m, n, m > 0. ? n/m : 0.);
    for(int s = 0; s < NSrc; s++) printf(" %7.2f%%", 100*uSrc[s]->GetBinContent(i));
    printf(" | %7.2f%% %7.2f%% %7.2f%% %7.2f%%  stat %.2f%%%s\n", 100*rp, 100*rj, 100*ri, 100*rt, 100*uStat->GetBinContent(i),
           reported(i) ? "" : "  (underflow)");
  }

  // ---- write ------------------------------------------------------------------
  TFile *fo = TFile::Open(outRoot, "recreate");
  TH1D *hMeasOut = (TH1D*) hMeasD->Clone("h_bJetPt_measured"); hMeasOut->Write();
  hU->Write(); hUm->Write(); hUp->Write(); hUsys->Write();
  uStat->Write(); uTot->Write(); uPur->Write(); uJEU->Write(); uIter->Write();
  for(int s = 0; s < NSrc; s++) uSrc[s]->Write();
  respT->Write();
  TNamed info("info", Form("pp calo jets: b-jet spectrum (muon-tagged x purity) unfolded with %s (PYTHIA calo, PF flavour, "
                            "even+odd 2026-10-01), Bayes N = %d (unfoldClosureTest min MSE), reco floor %.0f GeV, matched jets "
                            "only. Syst: purity sources propagated through the unfolding (+/- shift), spectrum JEU at gen level, "
                            "iterations N+/-1; quadrature. Stat propagated by RooUnfold. Underflow 50-80 GeV not reported.",
                            respName, nIterOpt, recoFloor));
  info.Write();
  fo->Close();
  printf("\nwritten %s\n", outRoot);

  int colSys = TColor::GetColor("#D55E00");

  // ---- figure 1: response matrix ----------------------------------------------
  {
    // row-normalised: P(reco | gen) for each gen bin, so the colour reads as a
    // migration probability rather than a cross-section. Drawn with one
    // equal-size cell per bin and the bin ranges as labels: the bins run from
    // 10 to 200 GeV wide, and on a pT axis the low ones collapse.
    TH2D *pr = new TH2D("prob", "", NPt, 0, NPt, NPt, 0, NPt); pr->SetDirectory(nullptr);
    for(int iy = 1; iy <= NPt; iy++){
      double tot = truthMC->GetBinContent(iy);   // includes jets lost below the floor
      for(int ix = 1; ix <= NPt; ix++) pr->SetBinContent(ix, iy, tot > 0. ? respT->GetBinContent(ix, iy)/tot : 0.);
    }
    for(int i = 1; i <= NPt; i++){
      TString lab = Form("%.0f-%.0f", ptEdges[i-1], ptEdges[i]);
      pr->GetXaxis()->SetBinLabel(i, lab); pr->GetYaxis()->SetBinLabel(i, lab);
    }
    gStyle->SetPalette(kViridis);
    TCanvas *c = new TCanvas("cResp", "", 850, 750);
    c->SetLeftMargin(0.17); c->SetRightMargin(0.15); c->SetBottomMargin(0.16); c->SetTopMargin(0.12);
    pr->GetXaxis()->SetTitle("reco #it{p}_{T}^{jet} [GeV]"); pr->GetYaxis()->SetTitle("gen #it{p}_{T}^{jet} [GeV]");
    pr->GetZaxis()->SetTitle("P(reco | gen)");
    pr->GetXaxis()->SetLabelSize(0.036); pr->GetYaxis()->SetLabelSize(0.036);
    pr->GetXaxis()->SetTitleSize(0.045); pr->GetYaxis()->SetTitleSize(0.045); pr->GetZaxis()->SetTitleSize(0.042);
    pr->GetXaxis()->SetTitleOffset(1.55); pr->GetYaxis()->SetTitleOffset(1.85); pr->GetZaxis()->SetTitleOffset(1.1);
    pr->SetMinimum(0.); pr->SetMaximum(1.);
    pr->SetMarkerSize(1.2); gStyle->SetPaintTextFormat(".2f");
    pr->Draw("COLZ TEXT");
    // reco floor sits on a bin edge: index of the first reco bin at or above it
    double xFloor = 0.; for(int i = 0; i <= NPt; i++) if(ptEdges[i] <= recoFloor + 1e-6) xFloor = i;
    TLine *lf = new TLine(xFloor, 0, xFloor, NPt); lf->SetLineColor(kRed+1); lf->SetLineWidth(3); lf->SetLineStyle(2); lf->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.17, 0.945, "PYTHIA, pp 5.02 TeV, calo jets: muon-tagged b-jet response");
    la.SetTextSize(0.030); la.SetTextColor(kRed+1);
    la.DrawLatex(0.17, 0.90, Form("dashed: reco floor %.0f GeV (jets reconstructed below it are misses)", recoFloor));
    savePdfTight(c, Form("%s/unfold_caloJets_pp_response.pdf", figDir));
    delete c;
  }

  // ---- figure 2: measured vs unfolded (N-1, N, N+1) ---------------------------
  {
    TCanvas *c = new TCanvas("cMU", "", 700, 800);
    TPad *pTop, *pBot; splitPads(pTop, pBot);
    pTop->cd(); pTop->SetLogy();   // log: falls by ~3 orders of magnitude, no negative bins
    TH1D *fr = (TH1D*) hU->Clone("frMU"); fr->Reset(); fr->SetTitle("");
    double ymax = 0., ymin = 1e30;
    for(int i = 1; i <= NPt; i++) if(reported(i)){ ymax = TMath::Max(ymax, TMath::Max(hU->GetBinContent(i), hMeasD->GetBinContent(i)));
                                                     ymin = TMath::Min(ymin, TMath::Min(hU->GetBinContent(i), hMeasD->GetBinContent(i))); }
    fr->SetMinimum(0.3*ymin); fr->SetMaximum(30.*ymax);
    fr->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); fr->GetXaxis()->SetLabelSize(0);
    fr->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} [GeV^{-1}]"); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
    fr->Draw("AXIS");
    TH1D *dM = (TH1D*) hMeasD->Clone("dM"), *dU = (TH1D*) hU->Clone("dU"), *dUm = (TH1D*) hUm->Clone("dUm"), *dUp = (TH1D*) hUp->Clone("dUp");
    for(TH1D *h : {dM, dU, dUm, dUp}) for(int i = 1; i <= NPt; i++) if(!reported(i) || h->GetBinContent(i) <= 0.){ h->SetBinContent(i, -999.); h->SetBinError(i, 0.); }
    styleH(dM,  "#0072B2", markOpenCircle,   1.2);
    styleH(dU,  "#000000", markFilledCircle, 1.2);
    styleH(dUm, "#E69F00", markOpenSquare,   1.2);
    styleH(dUp, "#CC79A7", markOpenDiamond,  1.5);
    for(TH1D *h : {dUm, dUp}) for(int i = 1; i <= NPt; i++) h->SetBinError(i, 0.);
    dM->Draw("E1 X0 same"); dUm->Draw("P same"); dUp->Draw("P same"); dU->Draw("E1 X0 same");
    TLegend *leg = makeLegend(0.25, 0.05, 0.62, 0.30, 0.036);   // lower left is empty: the spectrum falls
    leg->AddEntry(dM,  "measured (reco), stat.", "lp");
    leg->AddEntry(dU,  Form("unfolded, #it{N} = %d, stat.", nIterOpt), "lp");
    leg->AddEntry(dUm, Form("unfolded, #it{N} = %d", nIterOpt - 1), "p");
    leg->AddEntry(dUp, Form("unfolded, #it{N} = %d", nIterOpt + 1), "p");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
    la.DrawLatex(0.21, 0.84, "pp 5.02 TeV, calo jets: b jets");
    la.SetTextSize(0.036);
    la.DrawLatex(0.21, 0.78, "muon-tagged #times b purity, no efficiency correction");

    pBot->cd();
    TH1D *r0 = makeRatio(hU, hMeasD, "r0"), *rm = makeRatio(hUm, hMeasD, "rm"), *rp = makeRatio(hUp, hMeasD, "rp");
    for(TH1D *h : {r0, rm, rp}){ for(int i = 1; i <= NPt; i++){ h->SetBinError(i, 0.); if(!reported(i)) h->SetBinContent(i, -999.); } }
    styleH(r0, "#000000", markFilledCircle, 1.2); styleH(rm, "#E69F00", markOpenSquare, 1.2); styleH(rp, "#CC79A7", markOpenDiamond, 1.5);
    double lo = 1e9, hi = -1e9;
    for(TH1D *h : {r0, rm, rp}) for(int i = 1; i <= NPt; i++) if(reported(i)){ lo = TMath::Min(lo, h->GetBinContent(i)); hi = TMath::Max(hi, h->GetBinContent(i)); }
    TH1D *rF = (TH1D*) r0->Clone("rFMU"); rF->Reset();
    double pad = 0.15*(hi - lo + 0.05);
    rF->SetMinimum(lo - pad); rF->SetMaximum(hi + 2.5*pad);
    rF->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "unfolded / meas.");
    rF->Draw("AXIS");
    unityLine(ptReportMin, ptEdges[NPt])->Draw();
    rm->Draw("P same"); rp->Draw("P same"); r0->Draw("P same");
    c->SaveAs(Form("%s/unfold_caloJets_pp_measuredVsUnfolded.pdf", figDir));
    delete c;
  }

  // ---- figure 3: the unfolded spectrum with stat + total syst -----------------
  {
    TCanvas *c = new TCanvas("cU", "", 700, 800);
    TPad *pTop, *pBot; splitPads(pTop, pBot);
    pTop->cd(); pTop->SetLogy();   // log: falls by ~3 orders of magnitude, no negative bins
    TH1D *fr = (TH1D*) hU->Clone("frU"); fr->Reset(); fr->SetTitle("");
    double ymax = 0., ymin = 1e30;
    for(int i = 1; i <= NPt; i++) if(reported(i)){ ymax = TMath::Max(ymax, hU->GetBinContent(i)); ymin = TMath::Min(ymin, hU->GetBinContent(i)); }
    fr->SetMinimum(0.3*ymin); fr->SetMaximum(30.*ymax);
    fr->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); fr->GetXaxis()->SetLabelSize(0);
    fr->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} [GeV^{-1}]"); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
    fr->Draw("AXIS");
    for(int i = 1; i <= NPt; i++){
      double v = hU->GetBinContent(i); if(v <= 0. || !reported(i)) continue;
      TBox *bx = new TBox(ptEdges[i-1], v - hUsys->GetBinError(i), ptEdges[i], v + hUsys->GetBinError(i));
      bx->SetFillColorAlpha(colSys, 0.30); bx->SetLineColor(colSys); bx->Draw("l"); bx->Draw();
    }
    TH1D *dU = (TH1D*) hU->Clone("dU3");
    for(int i = 1; i <= NPt; i++) if(!reported(i)){ dU->SetBinContent(i, -999.); dU->SetBinError(i, 0.); }
    styleH(dU, "#000000", markFilledCircle, 1.2);
    dU->Draw("E1 X0 same");
    TLegend *leg = makeLegend(0.25, 0.05, 0.60, 0.20, 0.036);
    leg->AddEntry(dU, "unfolded b jets, stat.", "lp");
    TBox *lb = new TBox(); lb->SetFillColorAlpha(colSys, 0.30); lb->SetLineColor(colSys);
    leg->AddEntry(lb, "total systematic", "f");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
    la.DrawLatex(0.21, 0.84, "pp 5.02 TeV, calo jets: b jets");
    la.SetTextSize(0.036);
    la.DrawLatex(0.21, 0.78, Form("unfolded (Bayes, #it{N} = %d), no efficiency correction", nIterOpt));

    pBot->cd();
    TH1D *rS = (TH1D*) uTot->Clone("rS3"), *rT = (TH1D*) uStat->Clone("rT3");
    for(TH1D *h : {rS, rT}){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); }
    styleLine(rS, "#D55E00", 3); styleLine(rT, "#0072B2", 3); rT->SetLineStyle(2);
    double m = 0.; for(int i = 1; i <= NPt; i++) if(reported(i)) m = TMath::Max(m, TMath::Max(rS->GetBinContent(i), rT->GetBinContent(i)));
    TH1D *rF = (TH1D*) rS->Clone("rF3"); rF->Reset();
    rF->SetMinimum(0.); rF->SetMaximum(1.7*m);   // headroom for the legend row
    rF->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "rel. unc. [%]");
    rF->Draw("AXIS");
    rS->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); rT->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    rS->Draw("HIST same"); rT->Draw("HIST same");
    TLegend *legB = makeLegend(0.19, 0.78, 0.95, 0.97, 0.075);
    legB->SetNColumns(2);
    legB->AddEntry(rS, "total syst.", "l"); legB->AddEntry(rT, "stat.", "l");
    legB->Draw();
    c->SaveAs(Form("%s/unfoldedBJetSpectrum_caloJets_pp.pdf", figDir));
    delete c;
  }

  // ---- figure 4: systematic breakdown of the unfolded spectrum ----------------
  {
    TCanvas *c = new TCanvas("cUS", "", 800, 700);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    TH1D *lT = (TH1D*) uTot->Clone("lT4"), *lP = (TH1D*) uPur->Clone("lP4"), *lJ = (TH1D*) uJEU->Clone("lJ4"), *lI = (TH1D*) uIter->Clone("lI4");
    for(TH1D *h : {lT, lP, lJ, lI}){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); h->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); }
    double ymax = 0.; for(int i = 1; i <= NPt; i++) if(reported(i)) ymax = TMath::Max(ymax, lT->GetBinContent(i));
    TH1F *fr = c->DrawFrame(ptReportMin, 0., ptEdges[NPt], 1.8*ymax);   // headroom for header + legend
    fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
    fr->GetYaxis()->SetTitle("relative systematic uncertainty [%]");
    fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
    fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
    styleLine(lP, "#009E73", 3); lP->SetLineStyle(2);
    styleLine(lJ, "#0072B2", 3); lJ->SetLineStyle(3);
    styleLine(lI, "#E69F00", 3); lI->SetLineStyle(7);
    styleLine(lT, "#000000", 4);
    lP->Draw("HIST same"); lJ->Draw("HIST same"); lI->Draw("HIST same"); lT->Draw("HIST same");
    TLegend *leg = makeLegend(0.40, 0.62, 0.95, 0.85, 0.032);
    leg->AddEntry(lP, "b purity (propagated through unfolding)", "l");
    leg->AddEntry(lJ, "spectrum JEU (gen level)", "l");
    leg->AddEntry(lI, Form("iterations (#it{N} = %d #pm 1)", nIterOpt), "l");
    leg->AddEntry(lT, "total (quadrature)", "l");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
    la.DrawLatex(0.17, 0.875, "pp 5.02 TeV, calo jets: unfolded b-jet spectrum systematics");
    savePdfTight(c, Form("%s/unfoldedBJetSpectrum_sysBreakdown_caloJets_pp.pdf", figSysDir));
    delete c;
  }
  printf("figures written to %s and %s\n", figDir, figSysDir);
}
