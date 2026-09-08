// Effect of the ptRel background subtraction on the b-jet template fit purity.
//
// UPDATED 2026-09-08. Two changes from the version that produced the figures in
// the 2026-09-01 talk:
//
//   1. BOTH background pairings are now subtracted, matching
//      S = D - T2 - T3 in plotPtRelTemplateDecomposition.C:
//        T2  h_mixedMuonPtRel_recoJetPt       fake mu + real reco jet
//        T3  h_fastJetMuonPtRel_..._bkgSub_RC fake mu + mixed fastJet
//      Only T3 came off before; T2 did not exist in the scan until 2026-09-08.
//      Both describe a muon that is not from the jet it was tagged to.
//
//   2. The source is the 100-resample scan, ~95x the mixed-event statistics of
//      the file used before, so the subtracted templates are far less noisy.
//
//   Also relevant: T3's jet pT axis moved from raw rcSub to JEC-corrected on
//   2026-09-08, so it now shares the data's pT scale. In the earlier figures a
//   given jet pT window selected a systematically different jet population on
//   the template side than on the data side.
//
//   3. The fit range is 0-4 GeV rather than 0-3. T2 pairs a muon
//      against a jet axis computed without that muon (see the geometry caveat
//      in plotPtRelTemplateDecomposition.C), which puts its mean ptRel near
//      5 GeV in central events -- a 0-3 window excluded nearly all of it, so
//      the template could not influence the fit whatever its size. See the
//      note on low_x/high_x: the range choice moves the baseline purity by
//      more than the subtraction does, so it is not a neutral knob.
//
// Reproduces the data-preparation and RooFit template-fit machinery from
// templateFitter() in calculateBJetsPerZ.cc (that file's PbPb, isData=1,
// 2-template branch only -- see that function, lines ~196-513, for the
// original), factored into a standalone function so it can be called twice
// per (class, jet pT) point: once with the fake-muon subtraction applied,
// once without. calculateBJetsPerZ.cc itself has no such toggle -- the
// subtraction there is unconditional -- so this is a new comparison, not a
// reproduction of an existing plot.
//
// *** FAKE-MUON SOURCE FILE: DELIBERATELY DIFFERENT FROM calculateBJetsPerZ.cc ***
// templateFitter() hardcodes
//   PbPb_SingleMuon_..._mixedEventPFClustering_..._2026-8-27_coarseBins_partial.root
// which was built (via makeCoarseBins_muonPtRel.C) from
//   PbPb_SingleMuon_..._mixedEventPFClustering_..._2026-8-27_ultraFineCentBins_partial.root
// That file has only 29539 events in slice C1 against the same-event sample's
// 731408 -- a ~4% subsample. It is almost certainly the mixed-event job that
// crashed on the null mixedEventPFCandidates_pt/eta/phi/id pointers (see that
// debugging session): the "_partial" suffix in the parent file's name was
// dropped somewhere in the chain, so nothing about the coarse-bins filename
// flags it as thin. Checked directly: the 2026-9-1 mixed-event file has
// 732005 events in C1, matching same-event to within 0.1%, and postdates the
// crash fix. This macro merges its fine slices into the same 4 coarse classes
// on the fly (see mergeFakeCoarse below) rather than pointing at either
// pre-built coarse-bins file. If calculateBJetsPerZ.cc's own fake-muon
// subtraction is ever revisited, its hardcoded path should be updated to
// match.
//
// DATA AND MC TEMPLATE FILES are the exact final values openTemplateFitterFiles()
// leaves in file_PbPb_SingleMuon / file_PH_DiJet / file_PH_BJet / file_PH_MuJet
// -- note file_PbPb_SingleMuon is assigned TWICE in that function (once to a
// 2026-5-4 file with jetTrkMaxFilter, then overwritten by a 2026-2-12 file
// without it); the second assignment is what templateFitter() actually reads,
// so that is what is used here too, to keep this an apples-to-apples
// before/after comparison against the SAME data callateBJetsPerZ.cc currently
// fits. That shadowed first assignment looks like an unintentional bug in the
// production script, worth a separate look.
//
// FIT PROCEDURE, replicated exactly from templateFitter(): project each 2D
// (ptRel, jetPt) template onto ptRel in the jetPt window; merge muJet+bJet
// samples into the DiJet b/bGS templates and muJet into the c template
// (mergeBtemplates/mergeCtemplates = true, matching calculateBPurity()'s
// defaults); compute c_truth/l_truth from the PRE-merge DiJet-only sample
// (order matters -- this must happen before the muJet/bJet merge, exactly as
// in the source); apply the bGS_multiplier = 1.2 gluon-splitting
// enhancement; fold the c template into the light template at its MC-truth
// fraction (2-template fit: b vs. light+c, not a free c fraction); rebin data
// and templates onto muRelPtAxis (0-5 GeV, 0.1 GeV steps) and convert to
// density; fit fb (b fraction) with RooFit (RooHistPdf b + light templates,
// SumW2Error, range 0-5 GeV). fb.getVal()/getError() is the b-purity and its
// fit error, matching templateFitter()'s returnValueIndex 1/2 exactly.
//
// The subtraction itself follows templateFitter() lines 313-319: per-event
// normalise data and each template, subtract, then scale back to the data's raw
// count so the fit sees the data's original normalisation. The only change is
// that two templates come off instead of one; they share fakeFile's event count
// so a single division serves both.
//
// RATIO PANEL ERRORS treat the two fits (with/without subtraction) as
// independent when propagating to the ratio. They are not independent -- both
// are fit to the same underlying data histogram, just before vs. after one
// histogram subtraction -- so the true ratio uncertainty is smaller than
// quadrature-of-relative-errors suggests. Read the ratio panel's error bars
// as an upper bound, not an exact uncertainty.
//
// Usage: root -l -b -q 'plotBPurity_fakeMuonSubtractionEffect.C'
// Run from: src/calculateBJetsPerZ/

#include "../../headers/functions/divideByBinwidth.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"

const char *dataFile =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_WDecayFilter_2026-2-12.root";
const char *diJetFile =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/PYTHIAHYDJET_DiJet_pThat-30_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_vzReweight_hiBinReweight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_weightCut_2026-2-12.root";
const char *bJetFile =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/PYTHIAHYDJET_BJet_pThat-30_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_vzReweight_hiBinReweight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_weightCut_2026-2-12.root";
const char *muJetFile =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/PYTHIAHYDJET_MuJet_pThat-30_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_vzReweight_hiBinReweight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_weightCut_2026-2-12.root";
const char *fakeFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/"
  "PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_"
  "mixedEventPFClustering_fastJetResamples-100_2026-9-8_ultraFineCentBins.root";

// The two background templates of the ptRel decomposition
// (plotPtRelTemplateDecomposition.C). Both describe a muon that is NOT from the
// jet it was tagged to; they differ only in whether the jet is real.
const char *bkgT2 = "h_mixedMuonPtRel_recoJetPt";                 // fake mu + real reco jet
const char *bkgT3 = "h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC";  // fake mu + mixed fastJet

const char *outDir = "../../figures/bPurity/";

const int NClass = 4;
const int sliceLo[NClass] = { 1,  3,  7, 11};   // fine-slice range in fakeFile per coarse class
const int sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };
const char *classTag[NClass]   = { "C1", "C2", "C3", "C4" };

const int    M = 51;
double muRelPtAxis[M] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0};

const int    N_jetPtAxisEdges = 8;
double       jetPtAxisEdges[N_jetPtAxisEdges] = {80,90,100,120,150,200,300,500};

// RooFit fit range, the named "fit_region" below. The observable itself spans
// the full muRelPtAxis (0-5 GeV); this restricts where the likelihood is
// evaluated. This is NOT a neutral choice -- the uncorrected purity in 0-10%,
// 80-90 GeV reads 0.767 / 0.843 / 0.828 at 0-3 / 0-5 / 0-4, so the range moves
// the baseline by more than the background subtraction does. Set deliberately.
//
// History: 0-3 originally; widened to 0-5 on 2026-09-08 because T2's weight
// sits near 5 GeV (mean ~5.1 in central events) and a 0-3 window excluded
// almost all of it, leaving the template unable to affect the fit whatever its
// size; then set to 0-4.
const double low_x = 0.0, high_x = 4.0;
const double c_multiplier = 1.0;
const double bGS_multiplier = 1.2;   // gluon-splitting enhancement factor

const char *noSubHex   = "#0072B2";   // Okabe-Ito blue
const char *withSubHex = "#D55E00";   // Okabe-Ito vermillion

const double ratioMin = 0.8, ratioMax = 1.2;

TFile *fData = nullptr, *fDiJet = nullptr, *fBJet = nullptr, *fMuJet = nullptr, *fFake = nullptr;

// Merge fine slices sliceLo[classIdx]..sliceHi[classIdx] of the fake-muon
// ptRel-vs-jetPt map into one coarse-class histogram, and its event count.
static TH2D* mergeFakeCoarse(int classIdx, const char *base, double &Nevents, int tag){
  TH2D *sum = nullptr;
  Nevents = 0.;
  for(int si = sliceLo[classIdx]; si <= sliceHi[classIdx]; si++){
    TH2D *h = nullptr; TH1D *v = nullptr;
    fFake->GetObject(Form("%s_C%d", base, si), h);
    fFake->GetObject(Form("h_vz_C%d", si), v);
    if(!h || !v){ printf("WARNING: %s slice C%d missing\n", base, si); continue; }
    if(!sum){ sum = (TH2D*) h->Clone(Form("bkgSum_%s_%d_%d", base, classIdx, tag)); sum->SetDirectory(nullptr); }
    else sum->Add(h);
    Nevents += v->Integral();
  }
  return sum;
}

// Replicates templateFitter()'s PbPb / isData=1 / 2-template-fit path.
// Returns the RooFit b-fraction (purity) and its fit error via out-params.
static void fitBPurity(int classIdx, double lowPt, double highPt, bool doFakeSub,
                        double &purity, double &purityErr, int tag){

  int centBin = classIdx + 1;

  TH2D *H_data = nullptr; TH1D *h_vz_data = nullptr;
  fData->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_C%d", centBin), H_data);
  fData->GetObject(Form("h_vz_C%d", centBin), h_vz_data);
  if(!H_data || !h_vz_data){ printf("ERROR: data histograms missing for C%d\n", centBin); purity = purityErr = -1.; return; }

  TH2D *H_b = nullptr, *H_bGS = nullptr, *H_c = nullptr, *H_u = nullptr, *H_d = nullptr, *H_s = nullptr, *H_g = nullptr;
  fDiJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bJets_C%dT0", centBin), H_b);
  fDiJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bGSJets_C%dT0", centBin), H_bGS);
  fDiJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_cJets_C%dT0", centBin), H_c);
  fDiJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_uJets_C%dT0", centBin), H_u);
  fDiJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_dJets_C%dT0", centBin), H_d);
  fDiJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_sJets_C%dT0", centBin), H_s);
  fDiJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_gJets_C%dT0", centBin), H_g);

  TH2D *H_bJet_b = nullptr, *H_bJet_bGS = nullptr, *H_muJet_b = nullptr, *H_muJet_bGS = nullptr, *H_muJet_c = nullptr;
  fBJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bJets_C%dT0", centBin), H_bJet_b);
  fBJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bGSJets_C%dT0", centBin), H_bJet_bGS);
  fMuJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bJets_C%dT0", centBin), H_muJet_b);
  fMuJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bGSJets_C%dT0", centBin), H_muJet_bGS);
  fMuJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_cJets_C%dT0", centBin), H_muJet_c);

  if(!H_b||!H_bGS||!H_c||!H_u||!H_d||!H_s||!H_g||!H_bJet_b||!H_bJet_bGS||!H_muJet_b||!H_muJet_bGS||!H_muJet_c){
    printf("ERROR: template histograms missing for C%d\n", centBin); purity = purityErr = -1.; return;
  }

  double smallShift = 0.01;
  int b1 = H_data->GetYaxis()->FindBin(lowPt + smallShift);
  int b2 = H_data->GetYaxis()->FindBin(highPt - smallShift);

  TH1D *h_data = H_data->ProjectionX(Form("h_data_%d", tag), b1, b2);
  h_data->SetDirectory(nullptr);

  if(doFakeSub){
    // Subtract BOTH background pairings, matching S = D - T2 - T3 in
    // plotPtRelTemplateDecomposition.C. Previously only T3 came off.
    //
    // Both templates and the data are put on a per-event footing before the
    // subtraction and the result scaled back to the data's raw count, so the
    // fit sees a histogram with the data's original normalisation. That is
    // templateFitter()'s convention, kept so this stays comparable to it; the
    // two templates share fakeFile's event count so one division serves both.
    const char *base[2] = {bkgT2, bkgT3};
    TH1D *h_bkg[2] = {nullptr, nullptr};
    double Nbkg[2] = {0., 0.};
    for(int t = 0; t < 2; t++){
      TH2D *H = mergeFakeCoarse(classIdx, base[t], Nbkg[t], tag*10 + t);
      if(!H || Nbkg[t] <= 0.){
        printf("ERROR: %s missing/empty for C%d\n", base[t], centBin);
        purity = purityErr = -1.; return;
      }
      h_bkg[t] = H->ProjectionX(Form("h_bkg%d_%d", t, tag), b1, b2);
      h_bkg[t]->SetDirectory(nullptr);
      h_bkg[t]->Scale(1./Nbkg[t]);
      delete H;
    }

    h_data->Scale(1./h_vz_data->Integral());
    h_data->Add(h_bkg[0], -1.);
    h_data->Add(h_bkg[1], -1.);
    h_data->Scale(h_vz_data->Integral());

    delete h_bkg[0]; delete h_bkg[1];
  }

  TH1D *h_b   = H_b  ->ProjectionX(Form("h_b_%d", tag),   b1, b2); h_b  ->SetDirectory(nullptr);
  TH1D *h_bGS = H_bGS->ProjectionX(Form("h_bGS_%d", tag), b1, b2); h_bGS->SetDirectory(nullptr);
  TH1D *h_c   = H_c  ->ProjectionX(Form("h_c_%d", tag),   b1, b2); h_c  ->SetDirectory(nullptr);
  TH1D *h_u   = H_u  ->ProjectionX(Form("h_u_%d", tag),   b1, b2); h_u  ->SetDirectory(nullptr);
  TH1D *h_d   = H_d  ->ProjectionX(Form("h_d_%d", tag),   b1, b2); h_d  ->SetDirectory(nullptr);
  TH1D *h_s   = H_s  ->ProjectionX(Form("h_s_%d", tag),   b1, b2); h_s  ->SetDirectory(nullptr);
  TH1D *h_g   = H_g  ->ProjectionX(Form("h_g_%d", tag),   b1, b2); h_g  ->SetDirectory(nullptr);

  TH1D *h_l = (TH1D*) h_u->Clone(Form("h_l_%d", tag));
  h_l->SetDirectory(nullptr);
  h_l->Add(h_d); h_l->Add(h_s); h_l->Add(h_g);

  // c_truth / l_truth from the PRE-merge DiJet-only sample -- must be
  // computed before muJet/bJet are folded into h_b/h_bGS/h_c below, exactly
  // as in templateFitter() (that ordering is load-bearing, not incidental)
  double c_truth = h_c->Integral() / (h_l->Integral() + h_b->Integral() + h_c->Integral() + h_bGS->Integral());
  double l_truth = h_l->Integral() / (h_l->Integral() + h_b->Integral() + h_c->Integral() + h_bGS->Integral());

  TH1D *h_bJet_b   = H_bJet_b  ->ProjectionX(Form("h_bJet_b_%d", tag),   b1, b2); h_bJet_b  ->SetDirectory(nullptr);
  TH1D *h_bJet_bGS = H_bJet_bGS->ProjectionX(Form("h_bJet_bGS_%d", tag), b1, b2); h_bJet_bGS->SetDirectory(nullptr);
  TH1D *h_muJet_b  = H_muJet_b ->ProjectionX(Form("h_muJet_b_%d", tag),  b1, b2); h_muJet_b ->SetDirectory(nullptr);
  TH1D *h_muJet_bGS= H_muJet_bGS->ProjectionX(Form("h_muJet_bGS_%d", tag), b1, b2); h_muJet_bGS->SetDirectory(nullptr);
  TH1D *h_muJet_c  = H_muJet_c ->ProjectionX(Form("h_muJet_c_%d", tag),  b1, b2); h_muJet_c ->SetDirectory(nullptr);

  // mergeCtemplates / mergeBtemplates = true (calculateBPurity() defaults)
  h_c->Add(h_muJet_c);
  h_b->Add(h_muJet_b); h_b->Add(h_bJet_b);
  h_bGS->Add(h_muJet_bGS); h_bGS->Add(h_bJet_bGS);

  double N_h_b = h_b->Integral();
  double N_h_bGS = h_bGS->Integral();
  double f_bGS_truth = N_h_bGS / (N_h_b + N_h_bGS);
  double GS_enhancement_shift = bGS_multiplier - 1.;

  h_bGS->Scale(N_h_b*(f_bGS_truth + GS_enhancement_shift) / (N_h_bGS*(1-(f_bGS_truth + GS_enhancement_shift))));
  h_b->Add(h_bGS);

  // 2-template fit: fold c into light at its MC-truth fraction rather than
  // fitting it freely
  h_c->Scale(c_multiplier * c_truth / h_c->Integral());
  h_l->Scale(l_truth / h_l->Integral());
  h_l->Add(h_c);

  h_l->Scale(1./h_l->Integral());
  h_b->Scale(1./h_b->Integral());

  TH1D *h_roo = (TH1D*) h_data->Clone(Form("h_roo_%d", tag));
  h_roo->SetDirectory(nullptr);
  h_data->Scale(1./h_data->Integral());

  h_roo  = (TH1D*) h_roo ->Rebin(M-1, Form("h_roo_r_%d", tag), muRelPtAxis);
  h_l    = (TH1D*) h_l   ->Rebin(M-1, Form("h_l_r_%d", tag), muRelPtAxis);
  h_b    = (TH1D*) h_b   ->Rebin(M-1, Form("h_b_r_%d", tag), muRelPtAxis);
  divideByBinwidth(h_roo);
  divideByBinwidth(h_l);
  divideByBinwidth(h_b);

  RooRealVar x(Form("x_%d", tag), "Observable", muRelPtAxis[0], muRelPtAxis[M-1]);
  x.setRange("fit_region", low_x, high_x);
  RooDataHist dataHist(Form("dataHist_%d", tag), "data histogram", x, RooFit::Import(*h_roo));
  RooDataHist bHist(Form("bHist_%d", tag), "MC b-jet histogram", x, RooFit::Import(*h_b));
  RooDataHist lHist(Form("lHist_%d", tag), "MC light+c-jet histogram", x, RooFit::Import(*h_l));
  RooHistPdf bPdf(Form("bPdf_%d", tag), "MC b-jet PDF", x, bHist, 0);
  RooHistPdf lPdf(Form("lPdf_%d", tag), "MC light+c-jet PDF", x, lHist, 0);
  RooRealVar fb(Form("fb_%d", tag), "b fraction", 0.5, 0.0, 1.0);
  RooAddPdf model(Form("model_%d", tag), "b+light templates", RooArgList(bPdf, lPdf), RooArgList(fb));
  model.fitTo(dataHist,
              RooFit::Minimizer("Minuit","migrad"),
              RooFit::SumW2Error(true),
              RooFit::Range("fit_region"),
              RooFit::PrintLevel(-1));

  purity = fb.getVal();
  purityErr = fb.getError();
}

void plotBPurity_fakeMuonSubtractionEffect(){
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);   // 56 fits' worth of INFO chatter otherwise

  fData  = TFile::Open(dataFile);
  fDiJet = TFile::Open(diJetFile);
  fBJet  = TFile::Open(bJetFile);
  fMuJet = TFile::Open(muJetFile);
  fFake  = TFile::Open(fakeFile);
  if(!fData||fData->IsZombie()){ printf("ERROR: cannot open %s\n", dataFile); return; }
  if(!fDiJet||fDiJet->IsZombie()){ printf("ERROR: cannot open %s\n", diJetFile); return; }
  if(!fBJet||fBJet->IsZombie()){ printf("ERROR: cannot open %s\n", bJetFile); return; }
  if(!fMuJet||fMuJet->IsZombie()){ printf("ERROR: cannot open %s\n", muJetFile); return; }
  if(!fFake||fFake->IsZombie()){ printf("ERROR: cannot open %s\n", fakeFile); return; }

  const int NPt = N_jetPtAxisEdges - 1;
  double ptCenter[NPt], ptWidth[NPt];
  for(int i = 0; i < NPt; i++){
    ptCenter[i] = 0.5*(jetPtAxisEdges[i] + jetPtAxisEdges[i+1]);
    ptWidth[i]  = 0.5*(jetPtAxisEdges[i+1] - jetPtAxisEdges[i]);
  }

  int colNo = TColor::GetColor(noSubHex), colWith = TColor::GetColor(withSubHex);

  for(int ci = 0; ci < NClass; ci++){

    double pNo[NPt], eNo[NPt], pWith[NPt], eWith[NPt];
    double ratio[NPt], eRatio[NPt];
    int tag = ci;

    printf("\n===== %s =====\n", classLabel[ci]);
    printf("%-14s %12s %12s %12s\n", "jetPt", "no-sub", "with-sub", "with/no");
    for(int i = 0; i < NPt; i++){
      fitBPurity(ci, jetPtAxisEdges[i], jetPtAxisEdges[i+1], false, pNo[i], eNo[i], tag*100 + i);
      fitBPurity(ci, jetPtAxisEdges[i], jetPtAxisEdges[i+1], true,  pWith[i], eWith[i], tag*100 + 50 + i);

      ratio[i] = (pNo[i] > 0.) ? pWith[i]/pNo[i] : 0.;
      eRatio[i] = (pNo[i] > 0. && pWith[i] > 0.) ?
        ratio[i]*TMath::Sqrt(TMath::Power(eWith[i]/pWith[i],2) + TMath::Power(eNo[i]/pNo[i],2)) : 0.;

      printf("%4.0f-%-9.0f %9.4f+-%.4f %9.4f+-%.4f %9.4f\n",
             jetPtAxisEdges[i], jetPtAxisEdges[i+1], pNo[i], eNo[i], pWith[i], eWith[i], ratio[i]);
    }

    TGraphErrors *gNo   = new TGraphErrors(NPt, ptCenter, pNo,   ptWidth, eNo);
    TGraphErrors *gWith = new TGraphErrors(NPt, ptCenter, pWith, ptWidth, eWith);
    TGraphErrors *gRatio= new TGraphErrors(NPt, ptCenter, ratio, ptWidth, eRatio);

    gNo->SetLineColor(colNo); gNo->SetMarkerColor(colNo);
    gNo->SetMarkerStyle(20); gNo->SetMarkerSize(1.1); gNo->SetLineWidth(2);
    gWith->SetLineColor(colWith); gWith->SetMarkerColor(colWith);
    gWith->SetMarkerStyle(21); gWith->SetMarkerSize(1.1); gWith->SetLineWidth(2);
    gRatio->SetLineColor(kBlack); gRatio->SetMarkerColor(kBlack);
    gRatio->SetMarkerStyle(20); gRatio->SetMarkerSize(1.0); gRatio->SetLineWidth(2);

    const double lm = 0.15, rm = 0.05, tm = 0.10, bm = 0.28, split = 0.35;
    TCanvas *c = new TCanvas(Form("c_%d", ci), "", 700, 800);

    TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
    pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
    pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
    pUp->Draw(); pUp->cd();

    TH1F *fr = pUp->DrawFrame(jetPtAxisEdges[0], 0., jetPtAxisEdges[NPt], 1.0);
    fr->GetYaxis()->SetTitle("b-jet purity");
    fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetLabelSize(0.045);
    fr->GetYaxis()->SetTitleOffset(1.25);
    fr->GetXaxis()->SetLabelSize(0.);

    gNo->Draw("P same");
    gWith->Draw("P same");

    TLegend *leg = new TLegend(0.55, 0.70, 0.90, 0.86);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);
    leg->AddEntry(gNo, "no background subtraction", "lp");
    leg->AddEntry(gWith, "with T2 + T3 subtraction", "lp");
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.19, 0.90, Form("PbPb SingleMuon (5.02 TeV), %s", classLabel[ci]));
    la.DrawLatex(0.19, 0.855, Form("2-template fit (b vs. light+c), %.0f < p_{T}^{rel} < %.0f GeV", low_x, high_x));

    c->cd();
    TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
    pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
    pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
    pDn->Draw(); pDn->cd();

    const double sc = (1.-split)/split;
    // fixed window -- some points' error bars run well outside [0.8,1.2]
    // (see the printed table), so they will visibly clip off the top/bottom
    // rather than being auto-fit; that clipping is deliberate zoom, not a bug
    TH1F *fr2 = pDn->DrawFrame(jetPtAxisEdges[0], ratioMin, jetPtAxisEdges[NPt], ratioMax);
    fr2->GetXaxis()->SetTitle("p_{T}^{recoJet} [GeV]");
    fr2->GetYaxis()->SetTitle("with / no sub");
    fr2->GetXaxis()->SetTitleSize(0.055*sc); fr2->GetXaxis()->SetLabelSize(0.045*sc);
    fr2->GetYaxis()->SetTitleSize(0.055*sc); fr2->GetYaxis()->SetLabelSize(0.045*sc);
    fr2->GetYaxis()->SetTitleOffset(1.25/sc);
    fr2->GetYaxis()->SetNdivisions(505);

    TLine *one = new TLine(jetPtAxisEdges[0], 1., jetPtAxisEdges[NPt], 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    gRatio->Draw("P same");

    TString out = TString(outDir) + Form("bPurity_fakeMuonSubtractionEffect_%s.pdf", classTag[ci]);
    c->SaveAs(out);
    printf("Saved %s\n", out.Data());
  }
}
