// "R_AA" of inclusive calo jets in MC (no quenching), before and after
// unfolding, per coarse centrality class. The point is to see what the
// unfolding does to an R_AA built exactly as the data one is, in a sample where
// the right answer is known.
//
//   R_AA(pT) = (1/N_evt) dN_AA/dpT / ( T_AA dsigma_pp/dpT )
//
// PbPb = PYTHIA+HYDJET, pp = PYTHIA, both calo jets + manual JEC.
//
// T_AA NORMALIZATION IN MC. Each PYTHIA+HYDJET event carries exactly ONE
// embedded PYTHIA hard scatter, whatever its centrality, so the per-event jet
// yield of the embedded sample does not grow with Ncoll (checked: gen jets
// above 100 GeV per weighted event are 5.21e-4 in PYTHIA and 5.19-5.24e-4 in
// every 5% PYTHIA+HYDJET slice 0-75%). To make the MC stand in for an
// unquenched PbPb class, each embedded event is promoted to the
//   Ncoll * sigma_gen / sigma_NN = T_AA * sigma_gen
// hard scatters a real event of that class would contain, with sigma_gen the
// PYTHIA cross section the sample represents. The pp side is the cross section
// dsigma/dpT = sigma_gen * (per-event yield). Hence
//
//   R_AA = [T_AA sigma_gen Y_AA] / [T_AA sigma_gen Y_pp] = Y_AA / Y_pp
//
// i.e. T_AA cancels EXACTLY, as it must for an unquenched model, and sigma_gen
// with it (the forest weights are relative, so sigma_gen is not known in
// absolute units here and is set to 1). The T_AA factors are applied explicitly
// below anyway, so the macro mirrors the data chain and the values used are on
// record. Scaling the per-event MC yield by 1/T_AA WITHOUT the sigma_gen T_AA
// promotion would give R_AA ~ sigma_NN/(Ncoll sigma_gen), i.e. 1/Ncoll -- not a
// meaningful number for an embedded sample.
//
// Glauber <Ncoll> per coarse class (src/plots/glauber/plotNcollFromNZ.C):
// 1558, 783.9, 267.0, 50.37, sigma_NN = 68.6 mb -> T_AA = 22.71, 11.43, 3.892,
// 0.734 mb^-1. The Glauber uncertainty (1.8-5%) is a common scale on each class
// and is not drawn.
//
// INPUTS
//   spectra   2026-10-02 PYTHIA_scan and PYTHIAHYDJET_scan (ultraFine, 5% slices
//             summed into the coarse classes with coarseCent.h)
//             reco  h_inclRecoJetPt_flavor_nominal[_C<s>], flavors summed
//             gen   h_inclGenJetPt_flavor[_C<s>],          flavors summed
//             N_evt h_vz[_C<s>] (all events; the inclusive jets are ungated)
//   response  the ones the data chain uses (calculateJetsPerZ.cc):
//             pp    PYTHIA_scan_response calo manualJEC 2026-09-22
//             PbPb  PYTHIAHYDJET_scan_response calo manualJEC 2026-09-25, C1..C4
//             h_matchedRecoJetPt_genJetPt_allJets[_C<i>], misses from
//             h_unmatchedGenJetPt[_C<i>] added to the truth.
//   Reco floor 60 GeV and Bayes 1 iteration, as the data chain; N = 1..4 are
//   printed.
//
// WHAT TO EXPECT. Gen-level R_AA is 1 by construction (drawn as the reference).
// At reco level the PbPb spectrum is smeared harder than pp (underlying-event
// fluctuations) and carries fake jets, both of which push R_AA above 1 at low pT
// in central classes. The response holds NO fakes, so unfolding removes the
// smearing but not the fakes: whatever is left above 1 at low pT after
// unfolding is the fake contribution (the data chain's doFakeJetSubtraction is
// not applied here). The responses are scanned from the same MC samples as the
// spectra, so this is not a statistically independent closure test.
//
// ERRORS. Statistical, PbPb and pp propagated (independent samples). Unfolded
// errors are RooUnfold's diagonal ones; bin-to-bin correlations ignored.
//
// EXCLUDING UNMATCHED (argument true). Drops reco jets with no gen match
// (flavor 18) from the measured spectra, i.e. removes the fakes the response
// does not model. Only possible in MC; it shows what the unfolding does when
// the fakes are handled perfectly. Outputs gain _matchedRecoOnly.
//
// Output  figures/unfoldTest/MCRAA_unfoldingEffect_caloJets.pdf (2x2 classes)
//         rootFiles/unfoldTest/MCRAA_unfoldingEffect_caloJets.root
//           RAA_{reco,unfold,gen}_C<i>, i = 1..4 -> 0-10 ... 50-80%
//
// Usage: root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");
//                          gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");'
//                     'plotMCRAA_unfoldingEffect.C+(false)'   (true: unmatched excluded)
// Run from: src/unfoldTest/
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/ratioPanel.h"
#include "../../headers/plotting/coarseCent.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <cstdio>

namespace {

const char *dir = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/";
const char *scanPP = "PYTHIA/PYTHIA_DiJet_caloJets_manualJEC_pThat-15_mu12_pTmu-15to999_tight_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_2026-10-2.root";
const char *scanAA = "PYTHIAHYDJET/PYTHIAHYDJET_manualJEC_DiJet_caloJets_pThat-15_mu12_pTmu-15to999_tight_vzReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_WDecayFilter_2026-10-2_ultraFineCentBins.root";
const char *respPP = "PYTHIA/PYTHIA_DiJet_response_caloJets_manualJEC_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_2026-9-22.root";
const char *respAA = "PYTHIAHYDJET/PYTHIAHYDJET_response_DiJet_caloJets_manualJEC_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-9-25.root";
TString outPdf  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/unfoldTest/MCRAA_unfoldingEffect_caloJets.pdf";
TString outRoot = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/unfoldTest/MCRAA_unfoldingEffect_caloJets.root";

const double recoFloor = 60.;
const int    nIterNominal = 1, nIterMax = 4;
const double yLo = 0.6, yHi = 2.0;

// data chain's axis (calculateJetsPerZ.cc newAxis); every edge on the 5 GeV grid
const double edges[] = {60,70,80,90,100,110,130,150,180,200,240,280,350,500};
const int    nEdges  = sizeof(edges)/sizeof(double) - 1;

const double sigmaNN = 68.6;                                   // mb
const double ncoll[NCoarse] = {1558., 783.9, 267.0, 50.37};    // Glauber, coarse classes
const double sigmaGen = 1.;                                    // cancels; see header

TH1* getH(TFile *f, const char *name)
{
  TH1 *h = nullptr; f->GetObject(name, h);
  if(!h) printf("ERROR: %s missing in %s\n", name, f->GetName());
  return h;
}

// Flavors summed. excludeFlavor >= 0 leaves that flavor out: the scans tag a reco
// jet with no gen match as flavor 18 (PYTHIAHYDJET_scan.C; PYTHIA_scan.C never
// fills 18, so it is a no-op on the pp side).
const int flavorUnmatched = 18;
TH1D* flavorSum(TH2D *H, const char *name, int excludeFlavor = -999)
{
  TH1D *h = H->ProjectionX(name, 1, H->GetNbinsY());
  h->SetDirectory(nullptr);
  if(excludeFlavor != -999){
    const int yb = H->GetYaxis()->FindBin(excludeFlavor + 1e-3);
    TH1D *x = H->ProjectionX(Form("%s_excl", name), yb, yb);
    h->Add(x, -1.);
    delete x;
  }
  return h;
}

void zeroBelowFloor(TH1D *h)
{
  for(int b = 0; b <= h->GetNbinsX() + 1; b++)
    if(h->GetXaxis()->GetBinUpEdge(b) < recoFloor + 1e-6){ h->SetBinContent(b, 0.); h->SetBinError(b, 0.); }
}

// Response with misses in the truth and the reco floor applied, as
// calculateJetsPerZ.cc: truth from the untruncated matrix, so matched jets
// reconstructed below the floor count as inefficiency.
RooUnfoldResponse* makeResponse(TFile *f, const char *mName, const char *missName, const char *tag)
{
  TH2D *m = (TH2D*) getH(f, mName);
  TH1D *miss = (TH1D*) getH(f, missName);
  if(!m || !miss) return nullptr;
  m = (TH2D*) m->Clone(Form("resp2D_%s", tag)); m->SetDirectory(nullptr);
  TH1D *truth = m->ProjectionY(Form("truth_%s", tag), 1, m->GetNbinsX());
  truth->SetDirectory(nullptr); truth->Add(miss);
  for(int ix = 0; ix <= m->GetNbinsX() + 1; ix++)
    if(m->GetXaxis()->GetBinUpEdge(ix) < recoFloor + 1e-6)
      for(int iy = 0; iy <= m->GetNbinsY() + 1; iy++){ m->SetBinContent(ix, iy, 0.); m->SetBinError(ix, iy, 0.); }
  TH1D *meas = m->ProjectionX(Form("meas_%s", tag), 1, m->GetNbinsY());
  meas->SetDirectory(nullptr);
  return new RooUnfoldResponse(meas, truth, m, Form("resp_%s", tag), tag);
}

// unfolded spectra for N = 1..nIterMax
void unfoldAll(RooUnfoldResponse *r, TH1D *reco, TH1D *out[], const char *tag)
{
  RooUnfoldBayes u(r, reco, 1);
  u.SetVerbose(0);
  for(int n = 1; n <= nIterMax; n++){
    u.SetIterations(n);
    out[n] = (TH1D*) u.Hunfold();
    out[n]->SetName(Form("unf_%s_N%d", tag, n)); out[n]->SetDirectory(nullptr);
  }
}

// counts on the 5 GeV axis -> per-event (scale) yield per GeV on the wide axis
TH1D* yieldOf(TH1D *h, double scale, const char *name)
{
  TH1D *r = (TH1D*) h->Rebin(nEdges, name, edges);   // new histogram
  r->SetDirectory(nullptr);
  r->Scale(scale);
  r->Scale(1., "width");
  return r;
}

} // namespace

void plotMCRAA_unfoldingEffect(bool excludeUnmatched = false)
{
  const int excl = excludeUnmatched ? flavorUnmatched : -999;
  if(excludeUnmatched){ outPdf.ReplaceAll(".pdf", "_matchedRecoOnly.pdf"); outRoot.ReplaceAll(".root", "_matchedRecoOnly.root"); }
  initPlotStyle();
  TFile *fsPP = TFile::Open(Form("%s%s", dir, scanPP)), *fsAA = TFile::Open(Form("%s%s", dir, scanAA));
  TFile *frPP = TFile::Open(Form("%s%s", dir, respPP)), *frAA = TFile::Open(Form("%s%s", dir, respAA));
  for(TFile *f : {fsPP, fsAA, frPP, frAA})
    if(!f || f->IsZombie()){ printf("ERROR: cannot open an input\n"); return; }

  // ---- pp: dsigma/dpT = sigma_gen * per-event yield
  const double nEvtPP = getH(fsPP, "h_vz")->Integral();
  TH1D *recoPP = flavorSum((TH2D*) getH(fsPP, "h_inclRecoJetPt_flavor_nominal"), "recoPP", excl);
  TH1D *genPP  = flavorSum((TH2D*) getH(fsPP, "h_inclGenJetPt_flavor"), "genPP");
  zeroBelowFloor(recoPP);
  RooUnfoldResponse *rPP = makeResponse(frPP, "h_matchedRecoJetPt_genJetPt_allJets", "h_unmatchedGenJetPt", "pp");
  if(!rPP) return;
  TH1D *unfPP[nIterMax + 1];
  unfoldAll(rPP, recoPP, unfPP, "pp");

  const double sPP = sigmaGen/nEvtPP;
  TH1D *xsRecoPP = yieldOf(recoPP, sPP, "xsRecoPP"), *xsGenPP = yieldOf(genPP, sPP, "xsGenPP");
  TH1D *xsUnfPP[nIterMax + 1];
  for(int n = 1; n <= nIterMax; n++) xsUnfPP[n] = yieldOf(unfPP[n], sPP, Form("xsUnfPP_N%d", n));

  gSystem->mkdir(gSystem->DirName(outPdf.Data()), kTRUE);
  gSystem->mkdir(gSystem->DirName(outRoot.Data()), kTRUE);
  TFile *fout = TFile::Open(outRoot, "recreate");
  TNamed("inputs", Form("pp scan %s; PbPb scan %s; pp response %s; PbPb response %s; floor %.0f GeV; Bayes N=%d; unmatched reco jets %s",
                        scanPP, scanAA, respPP, respAA, recoFloor, nIterNominal, excludeUnmatched ? "excluded" : "included")).Write();

  TCanvas *c = new TCanvas("cMCRAA", "", 1100, 900);
  c->Divide(2, 2, 0.001, 0.001);
  const double xLo = edges[0], xHi = edges[nEdges];

  for(int ci = 0; ci < NCoarse; ci++){
    const double taa = ncoll[ci]/sigmaNN;                       // mb^-1
    const double nEvt = coarseEvents(fsAA, ci, "h_vz");
    TH2D *R = coarseSum2(fsAA, "h_inclRecoJetPt_flavor_nominal", ci, "reco");
    TH2D *G = coarseSum2(fsAA, "h_inclGenJetPt_flavor", ci, "gen");
    if(nEvt <= 0 || !R || !G){ printf("ERROR: class %s incomplete in the PbPb scan\n", coarseLabel[ci]); return; }
    TH1D *recoAA = flavorSum(R, Form("recoAA_C%d", ci + 1), excl), *genAA = flavorSum(G, Form("genAA_C%d", ci + 1));
    delete R; delete G;
    zeroBelowFloor(recoAA);

    RooUnfoldResponse *rAA = makeResponse(frAA, Form("h_matchedRecoJetPt_genJetPt_allJets_C%d", ci + 1),
                                          Form("h_unmatchedGenJetPt_C%d", ci + 1), Form("C%d", ci + 1));
    if(!rAA) return;
    TH1D *unfAA[nIterMax + 1];
    unfoldAll(rAA, recoAA, unfAA, Form("C%d", ci + 1));

    // PbPb per-event yield, each embedded event promoted to T_AA * sigma_gen hard scatters
    const double sAA = taa*sigmaGen/nEvt;
    TH1D *yReco = yieldOf(recoAA, sAA, Form("yReco_C%d", ci + 1));
    TH1D *yGen  = yieldOf(genAA,  sAA, Form("yGen_C%d", ci + 1));
    TH1D *yUnf[nIterMax + 1];
    for(int n = 1; n <= nIterMax; n++) yUnf[n] = yieldOf(unfAA[n], sAA, Form("yUnf_C%d_N%d", ci + 1, n));

    // R_AA = yield / (T_AA dsigma/dpT)
    TH1D *tReco = (TH1D*) xsRecoPP->Clone("tReco"), *tGen = (TH1D*) xsGenPP->Clone("tGen");
    tReco->Scale(taa); tGen->Scale(taa);
    TH1D *raaReco = makeRatio(yReco, tReco, Form("RAA_reco_C%d", ci + 1), RatioErr::kBoth);
    TH1D *raaGen  = makeRatio(yGen,  tGen,  Form("RAA_gen_C%d",  ci + 1), RatioErr::kBoth);
    TH1D *raaUnf[nIterMax + 1];
    for(int n = 1; n <= nIterMax; n++){
      TH1D *t = (TH1D*) xsUnfPP[n]->Clone("tUnf"); t->Scale(taa);
      raaUnf[n] = makeRatio(yUnf[n], t, Form("RAA_unfold_N%d_C%d", n, ci + 1), RatioErr::kBoth);
      delete t;
    }
    delete tReco; delete tGen;

    printf("\n%s  T_AA = %.3f mb^-1  N_evt(w) = %.4g\n  pT bin      reco    gen   | unfolded N=1    N=2    N=3    N=4\n",
           coarseLabel[ci], taa, nEvt);
    for(int b = 1; b <= nEdges; b++){
      printf("  %3.0f-%-4.0f  %.3f  %.3f  |", edges[b-1], edges[b], raaReco->GetBinContent(b), raaGen->GetBinContent(b));
      for(int n = 1; n <= nIterMax; n++) printf("  %.3f", raaUnf[n]->GetBinContent(b));
      printf("\n");
    }

    fout->cd();
    raaReco->Write(Form("RAA_reco_C%d", ci + 1));
    raaGen->Write(Form("RAA_gen_C%d", ci + 1));
    raaUnf[nIterNominal]->Write(Form("RAA_unfold_C%d", ci + 1));
    for(int n = 1; n <= nIterMax; n++) raaUnf[n]->Write();

    // ---- panel
    c->cd(ci + 1);
    gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.06); gPad->SetBottomMargin(0.13);
    TH1D *fr = new TH1D(Form("fr_C%d", ci + 1), "", 1, xLo, xHi);
    fr->SetDirectory(nullptr);
    fr->SetMinimum(yLo); fr->SetMaximum(yHi);
    fr->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
    fr->GetYaxis()->SetTitle("MC \"R_{AA}\"");
    fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05);
    fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045);
    fr->GetYaxis()->SetTitleOffset(1.3);
    fr->Draw("AXIS");
    unityLine(xLo, xHi)->Draw();

    styleH(raaGen,  hexMC,          markFilledDiamond, 1.3);
    styleH(raaReco, hexData,        markFilledCircle,  1.0);
    styleH(raaUnf[nIterNominal], hexCorrected, markFilledSquare, 1.0);
    raaGen->Draw("E1 SAME"); raaReco->Draw("E1 SAME"); raaUnf[nIterNominal]->Draw("E1 SAME");

    TLatex t; t.SetNDC(); t.SetTextSize(0.05);
    t.DrawLatex(0.38, 0.87, Form("PYTHIA+HYDJET %s / PYTHIA", coarseLabel[ci]));
    t.SetTextSize(0.04);
    t.DrawLatex(0.38, 0.81, Form("T_{AA} = %.3g mb^{-1}, calo jets, |#eta| < 1.6", taa));
    if(excludeUnmatched) t.DrawLatex(0.38, 0.19, "unmatched reco jets excluded");
    // points above the frame (0-10% reco at low pT is ~20): print them rather than
    // let a log axis or a huge range flatten the region of interest
    TString off;
    for(TH1D *h : {raaReco, raaUnf[nIterNominal]})
      for(int b = 1; b <= nEdges; b++)
        if(h->GetBinContent(b) > yHi)
          off += Form("%s %.0f-%.0f: %.2f;  ", h == raaReco ? "reco" : "unf.", edges[b-1], edges[b], h->GetBinContent(b));
    if(off.Length()){
      t.SetTextSize(0.032);
      t.DrawLatex(0.38, 0.75, "Off scale:");
      TObjArray *parts = off.Tokenize(";");
      for(int k = 0; k < parts->GetEntries(); k++)
        t.DrawLatex(0.40 + 0.27*(k/3), 0.71 - 0.045*(k%3), TString(((TObjString*)parts->At(k))->GetString().Strip(TString::kBoth)).Data());
      delete parts;
    }
    if(ci == NCoarse - 1){
      TLegend *lg = makeLegend(0.45, 0.58, 0.94, 0.76, 0.045);
      lg->AddEntry(raaReco, "Reco, not unfolded", "pe");
      lg->AddEntry(raaUnf[nIterNominal], Form("Unfolded (Bayes, %d it.)", nIterNominal), "pe");
      lg->AddEntry(raaGen, "Gen (truth)", "pe");
      lg->Draw();
    }
  }

  savePdfTight(c, outPdf);
  fout->Close();
  printf("\nwrote %s\n      %s\n", outPdf.Data(), outRoot.Data());
}
