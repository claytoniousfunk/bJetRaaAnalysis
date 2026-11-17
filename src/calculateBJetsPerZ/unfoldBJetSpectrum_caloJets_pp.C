// Unfold the pp CALO-jet b-jet spectrum (muon-tagged x b purity) with the
// muon-tagged b-jet response, and propagate its uncertainties.
//
// INPUT  rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_caloJets_pp.root
//        (constructBJetSpectrum_caloJets_pp.C): h_bJetPt in dN/dpT with stat
//        errors, and the relative systematics per source.
//
// RESPONSE  h_matchedRecoJetPt_genJetPt_bJets_muTagged (X reco, Y gen) from the
//        PYTHIA calo response scans with PF-matched flavor, even + odd halves
//        added (2026-10-01). Muon-tagged = a reco muon in the jet, no trigger
//        requirement (the data spectrum is muon-triggered; the trigger
//        efficiency is a separate correction). Rebinned from the native 5 GeV
//        bins to ptEdges. Reco bins below recoFloor are zeroed and the truth is
//        taken from the untruncated matrix, so jets reconstructed below the
//        floor become misses. FLOOR = 70 GeV: the SingleMuon calo forest has no
//        jets below 70. The 70-80 GeV bin is a turn-on and arrives here already
//        corrected for it (constructBJetSpectrum_caloJets_pp.C). It must be:
//        muon-tagged calo b jets reconstruct at ~0.7 of gen pT, so reco 70-80
//        feeds gen 80-120 directly, and the uncorrected half-empty bin gave an
//        unfolding factor at gen 80-100 of 1.32 against a PYTHIA truth of 2.25.
//        Raising the floor to 80 instead removes that information: only 12% of
//        gen 80-100 tagged b jets reconstruct above 80, and N = 2 then misses the
//        tilted closure by 10% (min MSE moves to N = 9). Matched jets only: the truth is the gen projection of
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
//        N = 1 under-converged, max|r - 1| 10.7%). Measured with the floor at
//        70, the floor used here. The iteration systematic is
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
//   response stat  the response matrix's finite MC statistics, by toys: every
//                  bin of the (untruncated, rebinned) response is fluctuated
//                  within its MC error -- Poisson in its effective entries
//                  (c/sigma)^2, scaled back by sigma^2/c, so weighted bins stay
//                  positive and thin bins fluctuate correctly -- and the truth
//                  (prior), reco-floor truncation, measured projection and
//                  RooUnfoldResponse are rebuilt from it, so prior and
//                  migrations move together. The data are unfolded with each of
//                  nToys toys at nIterOpt; the RMS over toys / nominal is the
//                  systematic. The toys' mean shift is stored as a check.
//                  Treated as a systematic: it comes from the MC sample size,
//                  is independent of the data statistics and does not shrink
//                  with more data.
//   total          quadrature sum of the four.
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
//           h_unf_sysRel_responseStat    response-matrix MC statistics (toy RMS)
//           h_unf_respToys_meanShift     toys' mean / nominal - 1 (check)
//           h_respToys_values            every toy's unfolded dN/dpT (X = toy, Y = jet pT), for
//                                        src/plots/responseMatrix/plotResponseStatToys_caloJets_pp.C
//           h_response_neff, h_response_relErr  per-cell effective entries (c/sigma)^2 and
//                                        sigma/c of the rebinned, untruncated response
//           h_unf_sysRel_<purity source> each purity source after unfolding
//           h_response                   the rebinned, floor-truncated response
//         figures/bJetSpectra/unfold_caloJets_pp_response.pdf
//         figures/bJetSpectra/unfold_caloJets_pp_measuredVsUnfolded.pdf
//         figures/bJetSpectra/unfold_caloJets_pp_measuredVsUnfolded_simple.pdf   measured and nominal only, for talks
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
#include "../../headers/functions/bJetUnfoldExt.h"
#include "../../headers/functions/clipResponseOutliers.h"
#include "../../headers/functions/jerSmear.h"   // JER scheme (PFJets): nominal response x1.10, variation x1.10
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
#include "TRandom3.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/functions/dataDrivenPrior.h"
#include "../../headers/functions/jeuShift.h"
#include "TVectorD.h"
#include "../../headers/plotting/ratioPanel.h"
#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/functions/caloJetBPurityFit_pp.h"   // ptEdges, NPt, ptReportMin

// specPath / outRoot follow the jet collection (collTag() in the names, set in
// the macro); the response is the sum of ppRespFiles (jetCollectionInputs.h).
const char *specPath = nullptr;
// response histogram: ppRespName (jetCollectionInputs.h)

// calo: data have no jets below 70; 70-80 turn-on corrected upstream.
// caloMatchedPF: 50, the turn-on correction starts there (every bin kept).
// PFJets: 60, see CALO-MATCH ID below.
double recoFloor = 70.;

// CALO-MATCH ID (only when PFJets data meet a CALO-MATCHED response -- the
// 2026-10-07 first PF run; PFJets now unfolds with the plain-PF sibling-repo
// muon-tagged response and none of this applies). Then the plain-PF data were
// unfolded with the calo-matched muon-tagged b response,
// whose reco axis lacks the muon-tagged b jets without a calo partner: 9% are
// kept at 50-60 GeV, 56% at 70-80, 86% at 90-100, > 97% above 120
// (measureCaloMatchEff.C). Unfolding plain-PF data with it made the result
// swing by 35% at 80-120 GeV between N = 1 and 3. Each reco row of the native
// response is therefore scaled by w = 1 + idCorrFrac (1/eps - 1), eps =
// h_eff_muTagB_pp at the row's reco pT, nominal idCorrFrac = 1; the truth is the
// gen projection of the scaled matrix. Rows with eps < idEpsMin are not scaled
// (weights > 10 on a few surviving jets), which is why the floor is 60, not 50.
// Systematic: idCorrFrac = 0.5 and 1.5, envelope, carried in the turnOn slot
// (h_unf_sysRel_turnOn), which is otherwise zero for PFJets: it is a turn-on of
// the response, and every downstream macro already propagates that slot.
const char  *idEffPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/correctionFactors/caloMatchEff.root";
const double idEpsMin  = 0.1;
double idCorrFrac = 1.;
int          nIterOpt  = 2;   // PFJets with its plain-PF response: 3 (below)
const int    nToys     = 1000;   // response-matrix MC-statistics toys
const int    toySeed   = 20261004;

const char *outRoot   = nullptr;
const char *figDir    = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/bJetSpectra";
const char *figSysDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics";

const int NSrc = 6;   // + tplStat (template MC statistics, 2026-10-08)
const char *srcKey[NSrc]   = {"bGS", "cMult", "fitRangeHigh", "JEU", "JER", "tplStat"};

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

// JER variation of the b response (PFJets): jerSmear.h -- since 2026-11-05 the nominal response is
// itself smeared to a 10% worse resolution (HIN-18-014) and this is a further 10% on it (until then:
// unsmeared nominal, +20%, extra width 0.663). BJET_JER_SCHEME=old restores that.
TH2D* smearRespJER(TH2D *r, const char *name){ return jerSmear::variation(r, name); }
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

// coll: jet collection (jetCollectionInputs.h; "" = BJET_JETCOLL env or caloJets)
void unfoldBJetSpectrum_caloJets_pp(const char *coll = "")
{
  initPlotStyle();
  setJetCollection(coll);
  if(!ptBinGuard(purVar())) return;   // 7-bin jet-pT axis variant (caloJetBPurityFitCore.h)
  specPath = strdup(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_%s_pp%s.root", collTag(), purVar()));
  outRoot  = strdup(Form("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/CorrectedBJetSpectra/Data/unfoldedBJetSpectrum_%s_pp%s.root", collTag(), purVar()));
  if(jetColl != kCaloJets) recoFloor = 50.;
  bool respCaloMatched = false;
  for(const std::string &fn : ppRespFiles) if(fn.find("caloJetMatched") != std::string::npos) respCaloMatched = true;
  const bool idCorr = jetColl == kPFJets && respCaloMatched;
  if(idCorr) recoFloor = 60.;
  // PFJets, plain-PF 2026-02-18 muon-tagged response: N = 3, min MSE of unfoldClosureTest.C("pp",15,
  // "weighted","tilt:0.3",80.,500.,false,1,<resp>,<resp>,...,"h_matchedRecoJetPt_genJetPt_bJets","",true,50.,
  // "50,60,70,80,100,120,150,200,300,500"), 2026-10-07 (train = test: one file)
  if(jetColl == kPFJets && !idCorr) nIterOpt = bNIterPP;   // bJetUnfoldExt.h
  // PFJets with its plain-PF response, 2026-10-08: both axes extended to 20 GeV
  // (bJetUnfoldExt.h), reco floor 70, weight outliers of the response clipped
  // (clipResponseOutliers.h) -- as bJetChain_caloJets_PbPb.C, so pp and PbPb
  // unfold alike. The floor systematic becomes 70 -> 80 (was 50 -> 100).
  const bool extB = jetColl == kPFJets && !idCorr;
  if(extB) recoFloor = bRecoFloorPF();   // 70; 50 with BJET_UNFOLD_LOWFLOOR=1 (bJetUnfoldExt.h)
  const std::vector<double> extE = bExtEdges(NPt + 1, ptEdges);
  gSystem->mkdir(gSystem->DirName(outRoot), kTRUE);
  gSystem->mkdir(figDir, kTRUE); gSystem->mkdir(figSysDir, kTRUE);

  // ---- inputs -----------------------------------------------------------------
  TFile *fS = TFile::Open(specPath);
  if(!fS || fS->IsZombie()){ printf("ERROR: cannot open %s -- run constructBJetSpectrum_caloJets_pp.C\n", specPath); return; }
  TH1D *hMeasD = nullptr, *hRelJEU = nullptr, *hRelTurn = nullptr, *hRelSrc[NSrc];
  fS->GetObject("h_bJetPt", hMeasD);
  fS->GetObject("h_bJetPt_sysRel_spectrumJEU", hRelJEU);
  fS->GetObject("h_bJetPt_sysRel_turnOn", hRelTurn);
  for(int s = 0; s < NSrc; s++) fS->GetObject(Form("h_bJetPt_sysRel_%s", srcKey[s]), hRelSrc[s]);
  if(!hMeasD || !hRelJEU || !hRelTurn){ printf("ERROR: missing input histograms in %s\n", specPath); return; }
  for(int s = 0; s < NSrc; s++) if(!hRelSrc[s]){ printf("ERROR: missing h_bJetPt_sysRel_%s\n", srcKey[s]); return; }
  if(hMeasD->GetNbinsX() != NPt){ printf("ERROR: spectrum has %d bins, expected %d\n", hMeasD->GetNbinsX(), NPt); return; }

  TH2D *rNative = nullptr;
  TH1D *jerResSim = nullptr, *jerResNom = nullptr, *jerResVar = nullptr, *jerDataOut = nullptr;   // JER figure (2026-11-17)
  for(const std::string &fn : ppRespFiles){
    TFile *f = TFile::Open(Form("%s%s", ppRespDir, fn.c_str())); TH2D *r = nullptr;
    if(f && !f->IsZombie()) f->GetObject(ppRespName, r);
    if(!r){ printf("ERROR: %s missing in %s\n", ppRespName, fn.c_str()); return; }
    if(!rNative){ rNative = (TH2D*) r->Clone("rNative"); rNative->SetDirectory(nullptr); } else rNative->Add(r);
  }
  if(extB) clipResponseOutliers(rNative, 100., true, "pp");
  if(extB){   // JER nominal (2026-11-05, user): the response smeared to a 10% worse resolution (jerSmear.h)
    jerResSim = jerSmear::resolution(rNative, "h_bJER_resSim");   // simulation width (JER figure, 2026-11-17)
    TH2D *sm = jerSmear::nominal(rNative, "rNativeJERnom"); delete rNative; rNative = sm; rNative->SetName("rNative");
    jerSmear::announce("pp b response");
  }
  // calo-match ID correction of the response's reco axis (PFJets; see CALO-MATCH ID)
  TH1D *idEff = nullptr;
  TH2D *rNativeRaw = (TH2D*) rNative->Clone("rNativeRaw"); rNativeRaw->SetDirectory(nullptr);
  auto idScale = [&](TH2D *src, double frac, const char *name){
    TH2D *o = (TH2D*) src->Clone(name); o->SetDirectory(nullptr);
    for(int ix = 1; ix <= o->GetNbinsX(); ix++){
      double e = idEff->GetBinContent(idEff->FindBin(o->GetXaxis()->GetBinCenter(ix)));
      if(e < idEpsMin) continue;
      double w = 1. + frac*(1./e - 1.);
      for(int iy = 0; iy <= o->GetNbinsY()+1; iy++){ o->SetBinContent(ix, iy, o->GetBinContent(ix, iy)*w); o->SetBinError(ix, iy, o->GetBinError(ix, iy)*w); }
    }
    return o;
  };
  if(idCorr){
    TFile *fE = TFile::Open(idEffPath);
    if(!fE || fE->IsZombie()){ printf("ERROR: %s missing -- run measureCaloMatchEff.C\n", idEffPath); return; }
    fE->GetObject("h_eff_muTagB_pp", idEff);
    if(!idEff){ printf("ERROR: h_eff_muTagB_pp missing\n"); return; }
    idEff->SetDirectory(nullptr);
    delete rNative; rNative = idScale(rNativeRaw, idCorrFrac, "rNative");
    printf("calo-match ID correction of the response: reco rows x 1/eps (eps >= %.2f), floor %.0f GeV\n", idEpsMin, recoFloor);
  }

  // ---- data-driven prior (PFJets with its plain-PF response; dataDrivenPrior.h) --
  // the gen rows of the fine response are weighted by (g/100)^-k, k fitted so the
  // prior follows the unfolded data; the MC-prior result and a reco floor of 100
  // GeV are the systematics (h_unf_sysRel_prior, h_unf_sysRel_recoFloor), carried
  // in the turnOn slot
  const bool ddPrior = jetColl == kPFJets && !idCorr && bUseDataDrivenPrior;   // MC prior since 2026-10-21 (bJetUnfoldExt.h)
  const double floorAlt = extB ? 70. : 100.;   // extB: 60 -> 70 since 2026-11-05 (70 -> 80 before)
  TH1D *hMeasCounts = toCounts(hMeasD, "hMeasCountsPrior");
  auto unfoldFine = [&](TH2D *fine, double fl, int n, const char *name, TH1D **truthOut) -> TH1D* {
    if(extB){ TH2D *rx = rebin2DTo(fine, extE, Form("%s_rx", name)); TH1D *U = unfoldOnExt<RooUnfoldResponse, RooUnfoldBayes>(rx, hMeasCounts, NPt + 1, ptEdges, fl, n, name, truthOut); delete rx; return U; }
    TH2D *r2 = rebin2D(fine, Form("%s_r2", name));
    TH1D *tr = r2->ProjectionY(Form("%s_tr", name), 1, NPt); tr->SetDirectory(nullptr);
    for(int ix = 1; ix <= NPt; ix++) if(r2->GetXaxis()->GetBinUpEdge(ix) <= fl + 1e-6) for(int iy = 0; iy <= NPt+1; iy++){ r2->SetBinContent(ix, iy, 0.); r2->SetBinError(ix, iy, 0.); }
    TH1D *mm = r2->ProjectionX(Form("%s_mm", name), 1, NPt); mm->SetDirectory(nullptr);
    TH1D *in = (TH1D*) hMeasCounts->Clone(Form("%s_in", name)); in->SetDirectory(nullptr);
    for(int i = 1; i <= NPt; i++) if(ptEdges[i] <= fl + 1e-6){ in->SetBinContent(i, 0.); in->SetBinError(i, 0.); }
    RooUnfoldResponse rr(mm, tr, r2, Form("%s_rr", name), "");
    RooUnfoldBayes u(&rr, in, n); u.SetVerbose(0);
    TH1D *U = (TH1D*) u.Hunfold()->Clone(name); U->SetDirectory(nullptr);
    if(truthOut) *truthOut = tr; else delete tr;
    delete r2; delete mm; delete in;
    return U;
  };
  double priorTilt = 0.;
  TH2D *rNativeMCPrior = (TH2D*) rNative->Clone("rNativeMCPrior"); rNativeMCPrior->SetDirectory(nullptr);
  if(ddPrior){
    priorTilt = fitPriorTilt(rNative, [&](TH2D *f, TH1D *&t){ return unfoldFine(f, recoFloor, nIterOpt, "priorFit", &t); }, ptReportMin, extB ? bPriorFitHi() : 300., 8, true);
    TH2D *t = tiltGenRows(rNative, priorTilt, "rNativeTilted"); delete rNative; rNative = t;
    printf("data-driven prior: response gen rows x (g/100)^-%.3f\n", priorTilt);
  }
  // ---- response ---------------------------------------------------------------
  TH2D *resp2D = extB ? rebin2DTo(rNative, extE, "h_response_full") : rebin2D(rNative, "h_response_full");
  const int nR = resp2D->GetNbinsX();   // NPt, or NPt + 3 on the extended axes
  TH1D *truthMC = resp2D->ProjectionY("truthMC", 1, nR); truthMC->SetDirectory(nullptr);   // before the floor
  TH2D *respT = (TH2D*) resp2D->Clone("h_response"); respT->SetDirectory(nullptr);
  for(int ix = 1; ix <= nR; ix++)
    if(respT->GetXaxis()->GetBinUpEdge(ix) <= recoFloor + 1e-6)
      for(int iy = 0; iy <= nR+1; iy++){ respT->SetBinContent(ix, iy, 0.); respT->SetBinError(ix, iy, 0.); }
  TH1D *measMC = respT->ProjectionX("measMC", 1, nR); measMC->SetDirectory(nullptr);
  // analysis-axis counts in, analysis-axis counts out, whatever the response axes
  auto bayesA = [&](RooUnfoldResponse &rr, TH1D *measCounts, int n, const char *name) -> TH1D* {
    if(!extB){ RooUnfoldBayes u(&rr, measCounts, n); u.SetVerbose(0); TH1D *h = (TH1D*) u.Hunfold(); h->SetName(name); h->SetDirectory(nullptr); return h; }
    TH1D *in = bToExtAxis(measCounts, extE, recoFloor, Form("%s_inX", name));
    RooUnfoldBayes u(&rr, in, n); u.SetVerbose(0); TH1D *hx = (TH1D*) u.Hunfold(); hx->SetDirectory(nullptr);
    TH1D *h = bToAnaAxis(hx, NPt + 1, ptEdges, name); delete in; delete hx; return h; };
  RooUnfoldResponse resp(measMC, truthMC, respT, Form("resp_bMuTag_pp_%s", collTag()), Form("PYTHIA %s muon-tagged b response", collLabel()));

  // ---- unfold -----------------------------------------------------------------
  TH1D *hMeasC = toCounts(hMeasD, "hMeasC");
  auto unfoldN = [&](TH1D *measCounts, int n, const char *name){
    TH1D *h = bayesA(resp, measCounts, n, name);
    return toDensity(h, name);
  };
  TH1D *hU  = unfoldN(hMeasC, nIterOpt,     "h_bJetPt_unfolded");
  TH1D *hUm = unfoldN(hMeasC, nIterOpt - 1, Form("h_bJetPt_unfolded_N%d", nIterOpt - 1));
  TH1D *hUp = unfoldN(hMeasC, nIterOpt + 1, Form("h_bJetPt_unfolded_N%d", nIterOpt + 1));

  // JEU on the b-jet spectrum itself (PFJets, 2026-10-21, user): the measured spectrum shifted by
  // +/- the AK4PF JEC uncertainty (pp, |eta| < 1.6 average; headers/functions/jeuShift.h) and unfolded
  // with this response and N; signed up/down kept (h_unf_varRel_JEUshift_up/_down) for the correlated
  // ratios of bJetsPerZSys_PF_PbPb.C. Replaces the all-flavour gen-level ratio as spectrumJEU.
  const bool jeuOnB = jetColl == kPFJets && !idCorr;
  TH1D *jeuUp = nullptr, *jeuDn = nullptr;
  if(jeuOnB){
    const std::vector<double> dPP = jeuDelta("/home/clayton/Analysis/code/bJetRaaAnalysis/JetEnergyCorrections/Spring18_ppRef5TeV_V6_MC_Uncertainty_AK4PF.txt");
    TH1D *mU = jeuShiftMeasured(hMeasC, dPP, +1, "hMeasC_JEUup"), *mD = jeuShiftMeasured(hMeasC, dPP, -1, "hMeasC_JEUdown");   // fit-based shape since 2026-11-16 (jeuShift.h; BJET_JEU_SHAPE=measured = before)
    TH1D *uU = unfoldN(mU, nIterOpt, "h_bJetPt_unfolded_JEUup"), *uD = unfoldN(mD, nIterOpt, "h_bJetPt_unfolded_JEUdown");
    jeuUp = (TH1D*) hU->Clone("h_unf_varRel_JEUshift_up"); jeuUp->Reset(); jeuUp->SetDirectory(nullptr);
    jeuDn = (TH1D*) hU->Clone("h_unf_varRel_JEUshift_down"); jeuDn->Reset(); jeuDn->SetDirectory(nullptr);
    for(int i = 1; i <= NPt; i++){ const double n = hU->GetBinContent(i);
      jeuUp->SetBinContent(i, n > 0. ? uU->GetBinContent(i)/n - 1. : 0.); jeuDn->SetBinContent(i, n > 0. ? uD->GetBinContent(i)/n - 1. : 0.); }
    jeuUp->SetTitle("unfolded, JEU up (measured b spectrum shifted), signed;#it{p}_{T}^{jet} [GeV];var / nom - 1");
    jeuDn->SetTitle("unfolded, JEU down (measured b spectrum shifted), signed;#it{p}_{T}^{jet} [GeV];var / nom - 1");
    delete mU; delete mD; delete uU; delete uD;
  }

  // JER of the b response (PFJets, 2026-10-21, user: the pp b spectrum had none -- only the PbPb
  // response was smeared): the fine (nominal-smeared) response smeared further (jerSmear.h), unfolded with this N;
  // signed h_unf_varRel_JERresp, read by bJetsPerZSys_PF_PbPb.C as the pp part of JER (pp-PbPb uncorr.)
  TH1D *jerResp = (TH1D*) hU->Clone("h_unf_varRel_JERresp"); jerResp->Reset(); jerResp->SetDirectory(nullptr);
  jerResp->SetTitle(Form("unfolded, b response JER variation (%s), signed;#it{p}_{T}^{jet} [GeV];var / nom - 1", jerSmear::tagText()));
  if(jeuOnB){
    TH2D *sm = smearRespJER(rNative, "rNativeJER");
    jerResNom = jerSmear::resolution(rNative, "h_bJER_resNom"); jerResVar = jerSmear::resolution(sm, "h_bJER_resVar");   // JER figure, 2026-11-17
    TH2D *r2 = extB ? rebin2DTo(sm, extE, "respJER") : rebin2D(sm, "respJER");
    TH1D *tT = r2->ProjectionY("respJER_t", 1, nR); tT->SetDirectory(nullptr);
    for(int ix = 1; ix <= nR; ix++) if(r2->GetXaxis()->GetBinUpEdge(ix) <= recoFloor + 1e-6) for(int iy = 0; iy <= nR+1; iy++){ r2->SetBinContent(ix, iy, 0.); r2->SetBinError(ix, iy, 0.); }
    TH1D *tM = r2->ProjectionX("respJER_m", 1, nR); tM->SetDirectory(nullptr);
    RooUnfoldResponse rj(tM, tT, r2, "respJER_rr", "");
    TH1D *uC = bayesA(rj, hMeasC, nIterOpt, "uJERresp"), *uJ = toDensity(uC, "uJERrespD");
    // 2026-11-05 (user): measured on the MC (jerSmear.h) -- the nominal response's reco spectrum unfolded
    // with the nominal and the varied response; BJET_JER_INPUT=data = the data (uJ / hU above). Both printed.
    TH1D *mcA = extB ? bToAnaAxis(measMC, NPt + 1, ptEdges, "jerMCa") : (TH1D*) measMC->Clone("jerMCa"); mcA->SetDirectory(nullptr);
    TH1D *m0 = bayesA(resp, mcA, nIterOpt, "uJERmc0"), *m1 = bayesA(rj, mcA, nIterOpt, "uJERmc1");
    TH1D *jerData = (TH1D*) jerResp->Clone("jerData"); jerData->SetDirectory(nullptr);
    for(int i = 1; i <= NPt; i++){
      const bool rep = bReportedBin(ptEdges[i-1], ptEdges[i]);
      const double n = hU->GetBinContent(i), a = m0->GetBinContent(i);
      jerData->SetBinContent(i, (n > 0. && rep) ? uJ->GetBinContent(i)/n - 1. : 0.);
      const double mcV = (a > 0. && rep) ? m1->GetBinContent(i)/a - 1. : 0.;
      jerResp->SetBinContent(i, jerSmear::onData() ? jerData->GetBinContent(i) : mcV);
    }
    printf("b-response JER (pp, %s), MC | data:", jerSmear::inputText());
    for(int i = 1; i <= NPt; i++) if(bReportedBin(ptEdges[i-1], ptEdges[i])){ const double a = m0->GetBinContent(i);
      printf("  %.0f: %+.3f | %+.3f", ptEdges[i-1], a > 0. ? m1->GetBinContent(i)/a - 1. : 0., jerData->GetBinContent(i)); }
    printf("\n");
    jerDataOut = (TH1D*) jerData->Clone("h_unf_varRel_JERresp_data"); jerDataOut->SetDirectory(nullptr);   // the data-based shift, for the JER figure
    delete sm; delete r2; delete tT; delete tM; delete uJ; delete mcA; delete m0; delete m1; delete jerData;
  }

  // unfolding bias (PFJets, 2026-10-21, user): same-sample closure at the DATA-LIKE tilt k (fitted to the
  // unfolded data, fitPriorTilt; measured only, the prior stays the MC one), |unfolded / true - 1|;
  // read by bJetsPerZSys_PF_PbPb.C as the pp part of its unfoldingBias source
  TH1D *uBias = bookRel("h_unf_sysRel_unfoldingBias", "unfolded, non-closure at the data-like truth (same-sample closure)");
  if(jeuOnB){
    const double kD = bClosureTilt > -90. ? bClosureTilt : fitPriorTilt(rNative, [&](TH2D *f, TH1D *&t){ return unfoldFine(f, recoFloor, nIterOpt, "ubFit", &t); }, ptReportMin, extB ? bPriorFitHi() : 300., 8, false);   // 0.3 since 2026-11-05
    TH2D *P = bClosureTilt > -90. ? closureTiltRows(rNative, bClosureTilt, "ubP") : tiltGenRows(rNative, kD, "ubP");
    TH1D *pd = (TH1D*) hMeasC->Clone("ubPd"); pd->Reset(); pd->SetDirectory(nullptr);
    TH1D *tt = (TH1D*) hMeasC->Clone("ubTt"); tt->Reset(); tt->SetDirectory(nullptr);
    for(int ix = 1; ix <= P->GetNbinsX(); ix++) for(int iy = 1; iy <= P->GetNbinsY(); iy++){
      const double c = P->GetBinContent(ix, iy); if(c == 0.) continue;
      pd->Fill(P->GetXaxis()->GetBinCenter(ix), c); tt->Fill(P->GetYaxis()->GetBinCenter(iy), c);
    }
    for(int i = 1; i <= NPt; i++){ const double mc = hMeasC->GetBinContent(i); pd->SetBinError(i, mc > 0. ? pd->GetBinContent(i)*hMeasC->GetBinError(i)/mc : 0.); }
    TH1D *U = unfoldN(pd, nIterOpt, "ubUn");   // density
    printf("unfolding bias: closure tilt k = %+.3f, unfolded / true:", kD);
    for(int i = 1; i <= NPt; i++){ const double t = tt->GetBinContent(i)/tt->GetBinWidth(i); const double r = t > 0. ? U->GetBinContent(i)/t : 0.;
      uBias->SetBinContent(i, (t > 0. && bReportedBin(ptEdges[i-1], ptEdges[i])) ? fabs(r - 1.) : 0.); if(ptEdges[i-1] >= 80. && ptEdges[i] <= 200.) printf("  %.3f", r); }
    printf("\n");
    delete P; delete pd; delete tt; delete U;
  }

  // calo-match ID correction at idCorrFrac 0.5 and 1.5 (PFJets): envelope -> turnOn slot
  TH1D *uID = bookRel("h_unf_sysRel_caloMatchID", "unfolded, calo-match ID correction of the response (+/-50%)");
  if(idCorr) for(double fr : {0.5, 1.5}){
    TH2D *rn = idScale(rNativeRaw, fr, Form("rNativeID%.1f", fr));
    TH2D *r2 = rebin2D(rn, Form("respID%.1f", fr));
    TH1D *tr = r2->ProjectionY(Form("truthID%.1f", fr), 1, NPt); tr->SetDirectory(nullptr);
    for(int ix = 1; ix <= NPt; ix++)
      if(r2->GetXaxis()->GetBinUpEdge(ix) <= recoFloor + 1e-6)
        for(int iy = 0; iy <= NPt+1; iy++){ r2->SetBinContent(ix, iy, 0.); r2->SetBinError(ix, iy, 0.); }
    TH1D *mm = r2->ProjectionX(Form("measID%.1f", fr), 1, NPt); mm->SetDirectory(nullptr);
    RooUnfoldResponse rv(mm, tr, r2, Form("respIDv%.1f", fr), "");
    RooUnfoldBayes u(&rv, hMeasC, nIterOpt); u.SetVerbose(0);
    TH1D *h = toDensity((TH1D*) u.Hunfold(), Form("uID%.1f", fr));
    for(int i = 1; i <= NPt; i++){
      double n = hU->GetBinContent(i);
      if(n > 0.) uID->SetBinContent(i, TMath::Max(uID->GetBinContent(i), fabs(h->GetBinContent(i)/n - 1.)));
    }
  }

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
  // ---- turn-on correction (70-80 GeV measured bin), propagated the same way ---
  TH1D *uTurn = bookRel("h_unf_sysRel_turnOn", "unfolded, 70-80 GeV turn-on correction (propagated)");
  {
    TH1D *mUp = (TH1D*) hMeasC->Clone("mUpT"), *mDn = (TH1D*) hMeasC->Clone("mDnT");
    for(int i = 1; i <= NPt; i++){
      double r = hRelTurn->GetBinContent(i);
      mUp->SetBinContent(i, hMeasC->GetBinContent(i)*(1. + r));
      mDn->SetBinContent(i, hMeasC->GetBinContent(i)*(1. - r));
    }
    TH1D *uUp = unfoldN(mUp, nIterOpt, "uUpT"), *uDn = unfoldN(mDn, nIterOpt, "uDnT");
    for(int i = 1; i <= NPt; i++){
      double n = hU->GetBinContent(i);
      if(n > 0.) uTurn->SetBinContent(i, TMath::Max(fabs(uUp->GetBinContent(i)/n - 1.), fabs(uDn->GetBinContent(i)/n - 1.)));
    }
  }
  // data-driven prior and reco floor (PFJets): relative to the nominal, into the turn-on slot
  TH1D *uPrior = bookRel("h_unf_sysRel_prior", "unfolded, MC prior instead of the data-driven one");
  TH1D *uFloor = bookRel("h_unf_sysRel_recoFloor", Form("unfolded, reco floor %.0f instead of %.0f GeV", floorAlt, recoFloor));
  if(ddPrior){
    TH1D *a = toDensity(unfoldFine(rNativeMCPrior, recoFloor, nIterOpt, "uMCPrior", nullptr), "uMCPriorD");
    TH1D *b = toDensity(unfoldFine(rNative, floorAlt, nIterOpt, "uFloorAlt", nullptr), "uFloorAltD");
    for(int i = 1; i <= NPt; i++){ double n = hU->GetBinContent(i); if(n <= 0.) continue;
      uPrior->SetBinContent(i, fabs(a->GetBinContent(i)/n - 1.));
      if(ptEdges[i-1] >= floorAlt - 1e-6) uFloor->SetBinContent(i, fabs(b->GetBinContent(i)/n - 1.)); }
    for(int i = 1; i <= NPt; i++) uTurn->SetBinContent(i, sqrt(pow(uTurn->GetBinContent(i), 2) + pow(uPrior->GetBinContent(i), 2) + pow(uFloor->GetBinContent(i), 2)));
  }
  // the calo-match ID correction (PFJets) rides in the turn-on slot
  for(int i = 1; i <= NPt; i++) uTurn->SetBinContent(i, std::hypot(uTurn->GetBinContent(i), uID->GetBinContent(i)));
  // ---- response-matrix MC statistics: toys -------------------------------------
  TH1D *uResp  = bookRel("h_unf_sysRel_responseStat",  Form("unfolded, response-matrix MC stat. (RMS of %d toys)", nToys));
  TH1D *uRespM = bookRel("h_unf_respToys_meanShift",   "unfolded, response toys: mean / nominal - 1 (check)");
  TH2D *hToyVals = new TH2D("h_respToys_values", "unfolded dN/dpT per response toy;toy;#it{p}_{T}^{jet} [GeV]",
                            nToys, 0, nToys, NPt, ptEdges);
  hToyVals->SetDirectory(nullptr);
  TH2D *hNeff = (TH2D*) resp2D->Clone("h_response_neff"), *hRelE = (TH2D*) resp2D->Clone("h_response_relErr");
  hNeff->SetDirectory(nullptr); hRelE->SetDirectory(nullptr); hNeff->Reset(); hRelE->Reset();
  hNeff->SetTitle("response effective entries (c/#sigma)^{2};reco #it{p}_{T}^{jet} [GeV];gen #it{p}_{T}^{jet} [GeV]");
  hRelE->SetTitle("response relative MC error #sigma/c;reco #it{p}_{T}^{jet} [GeV];gen #it{p}_{T}^{jet} [GeV]");
  for(int ix = 1; ix <= nR; ix++) for(int iy = 1; iy <= nR; iy++){
    double c = resp2D->GetBinContent(ix, iy), e = resp2D->GetBinError(ix, iy);
    if(c > 0. && e > 0.){ hNeff->SetBinContent(ix, iy, (c/e)*(c/e)); hRelE->SetBinContent(ix, iy, e/c); }
  }
  {
    TRandom3 rnd(toySeed);
    std::vector<double> s1(NPt+1, 0.), s2(NPt+1, 0.);
    for(int t = 0; t < nToys; t++){
      TH2D *r2 = (TH2D*) resp2D->Clone(Form("toyResp%d", t)); r2->SetDirectory(nullptr);
      for(int ix = 1; ix <= nR; ix++) for(int iy = 1; iy <= nR; iy++){
        double c = resp2D->GetBinContent(ix, iy), e = resp2D->GetBinError(ix, iy);
        if(c <= 0. || e <= 0.) continue;
        double neff = (c/e)*(c/e);
        r2->SetBinContent(ix, iy, rnd.Poisson(neff)*e*e/c);
      }
      TH1D *tTruth = r2->ProjectionY(Form("toyTruth%d", t), 1, nR); tTruth->SetDirectory(nullptr);   // prior, before the floor
      for(int ix = 1; ix <= nR; ix++)
        if(r2->GetXaxis()->GetBinUpEdge(ix) <= recoFloor + 1e-6)
          for(int iy = 0; iy <= nR+1; iy++){ r2->SetBinContent(ix, iy, 0.); r2->SetBinError(ix, iy, 0.); }
      TH1D *tMeas = r2->ProjectionX(Form("toyMeas%d", t), 1, nR); tMeas->SetDirectory(nullptr);
      RooUnfoldResponse rt(tMeas, tTruth, r2, Form("toyRooResp%d", t), "");
      TH1D *h = bayesA(rt, hMeasC, nIterOpt, Form("toyU%d", t));
      for(int i = 1; i <= NPt; i++){ double v = h->GetBinContent(i)/h->GetBinWidth(i); s1[i] += v; s2[i] += v*v; hToyVals->SetBinContent(t+1, i, v); }
      delete h; delete tMeas; delete tTruth; delete r2;
    }
    for(int i = 1; i <= NPt; i++){
      double n = hU->GetBinContent(i); if(n <= 0.) continue;
      double mean = s1[i]/nToys, rms = sqrt(TMath::Max(0., s2[i]/nToys - mean*mean));
      uResp->SetBinContent(i, rms/n); uRespM->SetBinContent(i, mean/n - 1.);
    }
  }

  TH1D *uPur  = bookRel("h_unf_sysRel_bPurity",      "unfolded, b purity (all sources, propagated)");
  TH1D *uJEU  = bookRel("h_unf_sysRel_spectrumJEU",  "unfolded, spectrum JEU (gen level)");
  TH1D *uIter = bookRel("h_unf_sysRel_iterations",   Form("unfolded, iterations N = %d +/- 1", nIterOpt));
  TH1D *uTot  = bookRel("h_unf_sysRel",              "unfolded, total systematic");
  TH1D *uStat = bookRel("h_unf_statRel",             "unfolded, stat. (propagated through the unfolding)");
  TH1D *hUsys = (TH1D*) hU->Clone("h_bJetPt_unfolded_sysAbs"); hUsys->SetDirectory(nullptr);
  hUsys->SetTitle("unfolded b jets, error = total systematic;#it{p}_{T}^{jet} [GeV];d#it{N}/d#it{p}_{T} [GeV^{-1}]");

  printf("\n%-9s %10s %10s %8s | %8s %8s %8s %8s %8s | %8s %8s %8s %8s %8s %8s\n", "jet pT", "measured", "unfolded", "U/M",
         "bGS", "cMult", "fitHi", "JEU(p)", "JER", "purity", "specJEU", "iter", "respMC", "turnOn", "TOTAL");
  for(int i = 1; i <= NPt; i++){
    double n = hU->GetBinContent(i);
    double q = 0.; for(int s = 0; s < NSrc; s++) q += pow(uSrc[s]->GetBinContent(i), 2);
    double rp = sqrt(q), rj = jeuOnB ? TMath::Max(fabs(jeuUp->GetBinContent(i)), fabs(jeuDn->GetBinContent(i))) : hRelJEU->GetBinContent(i);
    double ri = n > 0. ? TMath::Max(fabs(hUm->GetBinContent(i)/n - 1.), fabs(hUp->GetBinContent(i)/n - 1.)) : 0.;
    double rr = uResp->GetBinContent(i);
    double ro = uTurn->GetBinContent(i);
    double rt = sqrt(rp*rp + rj*rj + ri*ri + rr*rr + ro*ro);
    uPur->SetBinContent(i, rp); uJEU->SetBinContent(i, rj); uIter->SetBinContent(i, ri); uTot->SetBinContent(i, rt);
    uStat->SetBinContent(i, n > 0. ? hU->GetBinError(i)/n : 0.);
    hUsys->SetBinError(i, rt*n);
    double m = hMeasD->GetBinContent(i);
    printf("%3.0f-%-5.0f %10.4g %10.4g %8.3f |", ptEdges[i-1], ptEdges[i], m, n, m > 0. ? n/m : 0.);
    for(int s = 0; s < NSrc; s++) printf(" %7.2f%%", 100*uSrc[s]->GetBinContent(i));
    printf(" | %7.2f%% %7.2f%% %7.2f%% %7.2f%% %7.2f%% %7.2f%%  stat %.2f%%  (toy mean shift %+.2f%%)%s\n", 100*rp, 100*rj, 100*ri, 100*rr, 100*ro, 100*rt,
           100*uStat->GetBinContent(i), 100*uRespM->GetBinContent(i), reported(i) ? "" : "  (underflow)");
  }

  // ---- write ------------------------------------------------------------------
  TFile *fo = TFile::Open(outRoot, "recreate");
  TH1D *hMeasOut = (TH1D*) hMeasD->Clone("h_bJetPt_measured"); hMeasOut->Write();
  hU->Write(); hUm->Write(); hUp->Write(); hUsys->Write();
  if(jeuUp){ jeuUp->Write(); jeuDn->Write(); }
  uBias->Write(); jerResp->Write();
  for(TH1D *h : {jerResSim, jerResNom, jerResVar, jerDataOut}) if(h) h->Write();
  uStat->Write(); uTot->Write(); uPur->Write(); uJEU->Write(); uIter->Write(); uTurn->Write(); uID->Write(); uPrior->Write(); uFloor->Write();
  { TVectorD kv(1); kv[0] = priorTilt; kv.Write("priorTilt"); } uResp->Write(); uRespM->Write(); hToyVals->Write(); hNeff->Write(); hRelE->Write();
  for(int s = 0; s < NSrc; s++) uSrc[s]->Write();
  respT->Write();
  TNamed info("info", Form("pp calo jets: b-jet spectrum (muon-tagged x purity) unfolded with %s (PYTHIA calo, PF flavor, "
                            "even+odd 2026-10-01), Bayes N = %d (unfoldClosureTest min MSE), reco floor %.0f GeV, matched jets "
                            "only. Syst: purity sources propagated through the unfolding (+/- shift), spectrum JEU at gen level, "
                            "iterations N+/-1, 70-80 GeV turn-on correction (propagated), response-matrix MC stat (RMS of %d Poisson toys of the response); quadrature. Stat "
                            "propagated by RooUnfold. Underflow 50-80 GeV not reported.",
                            ppRespName, nIterOpt, recoFloor, nToys));
  info.Write();
  fo->Close();
  { TFile *fc = TFile::Open(outRoot, "update"); TNamed collN("jetCollection", collTag()); collN.Write(); fc->Close(); }
  printf("\nwritten %s\n", outRoot);

  int colSys = TColor::GetColor("#D55E00");

  // ---- figure 1: response matrix ----------------------------------------------
  {
    // row-normalised: P(reco | gen) for each gen bin, so the colour reads as a
    // migration probability rather than a cross-section. Drawn with one
    // equal-size cell per bin and the bin ranges as labels: the bins run from
    // 10 to 200 GeV wide, and on a pT axis the low ones collapse.
    TH2D *pr = new TH2D("prob", "", nR, 0, nR, nR, 0, nR); pr->SetDirectory(nullptr);
    for(int iy = 1; iy <= nR; iy++){
      double tot = truthMC->GetBinContent(iy);   // includes jets lost below the floor
      for(int ix = 1; ix <= nR; ix++) pr->SetBinContent(ix, iy, tot > 0. ? respT->GetBinContent(ix, iy)/tot : 0.);
    }
    for(int i = 1; i <= nR; i++){
      TString lab = Form("%.0f-%.0f", respT->GetXaxis()->GetBinLowEdge(i), respT->GetXaxis()->GetBinUpEdge(i));
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
    double xFloor = 0.; for(int i = 0; i <= nR; i++) if(respT->GetXaxis()->GetBinLowEdge(i + 1) <= recoFloor + 1e-6) xFloor = i;
    TLine *lf = new TLine(xFloor, 0, xFloor, nR); lf->SetLineColor(kRed+1); lf->SetLineWidth(3); lf->SetLineStyle(2); lf->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.036);
    la.DrawLatex(0.17, 0.945, Form("PYTHIA, pp 5.02 TeV, %s: muon-tagged b-jet response", collLabel()));
    la.SetTextSize(0.030); la.SetTextColor(kRed+1);
    la.DrawLatex(0.17, 0.90, Form("dashed: reco floor %.0f GeV (jets reconstructed below it are misses)", recoFloor));
    savePdfTight(c, Form("%s/unfold_%s_pp_response%s.pdf", figDir, collTag(), purVar()));
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
    la.DrawLatex(0.21, 0.84, Form("pp 5.02 TeV, %s: b jets", collLabel()));
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
    savePdfTight(c, Form("%s/unfold_%s_pp_measuredVsUnfolded%s.pdf", figDir, collTag(), purVar()));
    delete c;
  }

  // ---- figure 2b: the same for talks: measured and the nominal unfolding only
  {
    TCanvas *c = new TCanvas("cMUs", "", 700, 800);
    TPad *pTop, *pBot; splitPads(pTop, pBot);
    pTop->cd(); pTop->SetLogy();
    TH1D *fr = (TH1D*) hU->Clone("frMUs"); fr->Reset(); fr->SetTitle("");
    double ymax = 0., ymin = 1e30;
    for(int i = 1; i <= NPt; i++) if(reported(i)){ ymax = TMath::Max(ymax, TMath::Max(hU->GetBinContent(i), hMeasD->GetBinContent(i)));
                                                     ymin = TMath::Min(ymin, TMath::Min(hU->GetBinContent(i), hMeasD->GetBinContent(i))); }
    fr->SetMinimum(0.3*ymin); fr->SetMaximum(30.*ymax);
    fr->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); fr->GetXaxis()->SetLabelSize(0);
    fr->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} [GeV^{-1}]"); fr->GetYaxis()->SetTitleSize(0.055);
    fr->GetYaxis()->SetTitleOffset(1.35); fr->GetYaxis()->SetLabelSize(0.045);
    fr->Draw("AXIS");
    TH1D *dM = (TH1D*) hMeasD->Clone("dMs"), *dU = (TH1D*) hU->Clone("dUs");
    for(TH1D *h : {dM, dU}) for(int i = 1; i <= NPt; i++) if(!reported(i) || h->GetBinContent(i) <= 0.){ h->SetBinContent(i, -999.); h->SetBinError(i, 0.); }
    styleH(dM, "#0072B2", markOpenCircle, 1.3); styleH(dU, "#000000", markFilledCircle, 1.3);
    dM->Draw("E1 X0 same"); dU->Draw("E1 X0 same");
    TLegend *leg = makeLegend(0.25, 0.08, 0.62, 0.24, 0.045);
    leg->AddEntry(dM, "measured", "lp"); leg->AddEntry(dU, "unfolded", "lp"); leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.050);
    la.DrawLatex(0.21, 0.84, Form("pp, #it{b} jets, %s", collLabel()));
    pBot->cd();
    TH1D *r0 = makeRatio(hU, hMeasD, "r0s");
    for(int i = 1; i <= NPt; i++){ r0->SetBinError(i, 0.); if(!reported(i)) r0->SetBinContent(i, -999.); }
    styleH(r0, "#000000", markFilledCircle, 1.3);
    double lo = 1e9, hi = -1e9;
    for(int i = 1; i <= NPt; i++) if(reported(i)){ lo = TMath::Min(lo, r0->GetBinContent(i)); hi = TMath::Max(hi, r0->GetBinContent(i)); }
    TH1D *rF = (TH1D*) r0->Clone("rFMUs"); rF->Reset();
    double pad = 0.15*(hi - lo + 0.05);
    rF->SetMinimum(lo - pad); rF->SetMaximum(hi + 2.5*pad);
    rF->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]);
    styleRatioAxes(rF, "#it{p}_{T}^{jet} [GeV]", "unfolded / meas.");
    rF->Draw("AXIS");
    unityLine(ptReportMin, ptEdges[NPt])->Draw();
    r0->Draw("P same");
    savePdfTight(c, Form("%s/unfold_%s_pp_measuredVsUnfolded_simple%s.pdf", figDir, collTag(), purVar()));
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
    la.DrawLatex(0.21, 0.84, Form("pp 5.02 TeV, %s: b jets", collLabel()));
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
    savePdfTight(c, Form("%s/unfoldedBJetSpectrum_%s_pp%s.pdf", figDir, collTag(), purVar()));
    delete c;
  }

  // ---- figure 4: systematic breakdown of the unfolded spectrum ----------------
  {
    TCanvas *c = new TCanvas("cUS", "", 800, 700);
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.06);
    TH1D *lT = (TH1D*) uTot->Clone("lT4"), *lP = (TH1D*) uPur->Clone("lP4"), *lJ = (TH1D*) uJEU->Clone("lJ4"), *lI = (TH1D*) uIter->Clone("lI4");
    TH1D *lR = (TH1D*) uResp->Clone("lR4");
    for(TH1D *h : {lT, lP, lJ, lI, lR}){ h->Scale(100.); for(int i = 0; i <= NPt+1; i++) h->SetBinError(i, 0.); h->GetXaxis()->SetRangeUser(ptReportMin, ptEdges[NPt]); }
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
    styleLine(lR, "#CC79A7", 3); lR->SetLineStyle(9);
    lP->Draw("HIST same"); lJ->Draw("HIST same"); lI->Draw("HIST same"); lR->Draw("HIST same"); lT->Draw("HIST same");
    TLegend *leg = makeLegend(0.40, 0.58, 0.95, 0.86, 0.032);
    leg->AddEntry(lP, "b purity (propagated through unfolding)", "l");
    leg->AddEntry(lJ, "spectrum JEU (gen level)", "l");
    leg->AddEntry(lI, Form("iterations (#it{N} = %d #pm 1)", nIterOpt), "l");
    leg->AddEntry(lR, "response-matrix MC stat.", "l");
    leg->AddEntry(lT, "total (quadrature)", "l");
    leg->Draw();
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
    la.DrawLatex(0.17, 0.875, Form("pp 5.02 TeV, %s: unfolded b-jet spectrum systematics", collLabel()));
    savePdfTight(c, Form("%s/unfoldedBJetSpectrum_sysBreakdown_%s_pp%s.pdf", figSysDir, collTag(), purVar()));
    delete c;
  }
  printf("figures written to %s and %s\n", figDir, figSysDir);
}
