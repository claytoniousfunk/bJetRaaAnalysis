#pragma once
#include <cstdlib>
#include <cstdio>
#include "TString.h"
// PURITY VARIANT (2026-10-21, user): environment BJET_PURITY_VARIANT=gsAbs0p175_fit5 runs the
// purity fits with the bGS share at natural + 0.175 (absolute, calo convention, also for PFJets;
// the bGS systematic then scans +/-0.025 absolute) and the fixed upper fit edge at 5 GeV (0-5 GeV
// fit). Every purity / pp b-spectrum file and figure of the chain then carries purVar()
// ("_gsAbs0p175fit5") in its name, so the nominal outputs are never touched. Unset: nominal.
inline bool purityVariantOn(){
  static int v = -1;
  if(v < 0){ const char *e = std::getenv("BJET_PURITY_VARIANT"); v = (e && TString(e) == "gsAbs0p175_fit5") ? 1 : 0;
             if(v) printf("PURITY VARIANT gsAbs0p175_fit5: bGS share natural + 0.175 (absolute), fit 0-5 GeV, outputs suffixed _gsAbs0p175fit5\n"); }
  return v == 1;
}
// LIGHT SIGMA (2026-10-22, user): the PbPb light template is smeared with the MC-closure sigma*
// (studyLightSigmaClosure_PF_PbPb.C; table headers/fitParameters/lightSigmaClosure_PF_PbPb.h) --
// NOMINAL since 2026-10-22 -- and the lightTpl systematic shifts it by the closure's toy RMS.
// BJET_PURITY_VARIANT=bcSigma restores the previous nominal (the sigma of the b + c shape fit, its fit
// error as the lightTpl shift), outputs suffixed "_bcSigma". The gsAbs0p175_fit5 variant keeps the
// b + c sigma too: the table was derived at the nominal 3 GeV edge and natural g->bb share.
// (BJET_PURITY_VARIANT=lightSigStar was the 2026-10-21 variant of the same thing; now identical to the
// nominal apart from its "_lightSigStar" suffix.)
inline bool bcSigmaVariantOn(){
  static int v = -1;
  if(v < 0){ const char *e = std::getenv("BJET_PURITY_VARIANT"); v = (e && TString(e) == "bcSigma") ? 1 : 0;
             if(v) printf("PURITY VARIANT bcSigma: PbPb light template smeared with the b + c shape-fit sigma (pre-2026-10-22 nominal), outputs suffixed _bcSigma\n"); }
  return v == 1;
}
inline bool lightSigStarOn(){
  static int v = -1;
  if(v < 0){ const char *e = std::getenv("BJET_PURITY_VARIANT"); v = (e && TString(e) == "lightSigStar") ? 1 : 0; }
  return v == 1;
}
// CSVv2 PURITY (2026-10-24, user: "unfold and compare to gen"): BJET_PURITY_VARIANT=csvV2_<wp>
// (wp = loose, medium, shapes, tight) feeds the pp b-jet chain (construct -> unfold -> muon reco
// eff -> tag frequency) the purity of ALL muon-tagged jets implied by a CSVv2 cut,
// P_tag f_pass / eps_b(PYTHIA) from src/plots/bPurity/fitPtRelCSVv2Cut_pp_PF.C, packaged as
// rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2<wp>.root by
// src/plots/bPurity/makeCSVv2PurityFile_pp_PF.C. Outputs suffixed "_csvV2<wp>". pp only: no
// fit or template setting changes, and no PbPb input carries the suffix.
// IP3D PURITY (2026-10-27, user: "what is the inclusive b-jet fraction measurement when we use ip3d"): the same
// mechanism with BJET_PURITY_VARIANT=ip3d_<tag>, suffix "_ip3d<tag>", reading the all-jet purity of the CSVv2-split
// ln|IP3D| fits (plotCSVv2TagEffFromIP3DFits_pp_PF.C) packaged by src/plots/bPurity/makeIP3DPurityFile_pp_PF.C.
inline const char* csvV2Var(){
  static TString v = "?";
  if(v == "?"){ const char *e = std::getenv("BJET_PURITY_VARIANT"); v = (e && TString(e).BeginsWith("csvV2_")) ? TString::Format("_csvV2%s", e + 6)
                                                                    : (e && TString(e).BeginsWith("ip3d_")) ? TString::Format("_ip3d%s", e + 5) : TString("");
                if(v.Length()) printf("PURITY VARIANT %s: pp purity from the CSVv2 cut fits, outputs suffixed %s\n", e, v.Data()); }
  return v.Data();
}
// LOW-FLOOR b-jet UNFOLDING (2026-10-26, user: "show the bJet RZ all the way down to this low pt"):
// BJET_UNFOLD_LOWFLOOR=1 moves the PF b-jet unfolding's reco floor from 70 to 50 GeV (the pp data start
// at 50) and reports from 60 GeV (bJetUnfoldExt.h: bUfSplit 50, bUfHi 60). Every file of the pp chain and
// every purity file it reads then carries "_lowFloor" after the purity suffix (copy the purity files to
// the suffixed names first); the PbPb chain takes it through its outTag. Unset: nominal.
inline bool bLowFloorOn(){
  static int v = -1;
  if(v < 0){ const char *e = std::getenv("BJET_UNFOLD_LOWFLOOR"); v = (e && TString(e) == "1") ? 1 : 0;
             if(v) printf("LOW-FLOOR b-jet UNFOLDING: reco floor 50 GeV, truth [20,50) [50,60) | 60 70 80 ..., reported from 60 GeV, outputs suffixed _lowFloor\n"); }
  return v == 1;
}
inline const char* lowFloorVar(){ return bLowFloorOn() ? "_lowFloor" : ""; }
// MU5 pp CHAIN (2026-10-29, user: "the b-jet inclusive fraction, derived from mu5"): BJET_MU5=1 runs the pp
// PF b-jet chain on the mu5 scans (HLT_HIL3Mu5_NHitQ10, hybridSoft, muon pT 7-15 GeV): the prescale-weighted
// pp SingleMuon data, the mu5 PYTHIA templates (turn-on shape, tag-frequency comparison) and the mu5 PYTHIA
// response (its muon-tagged b response, and the tag frequency from it), all set in setJetCollection
// (jetCollectionInputs.h); the muon data / MC scale factor is 1 (no 7-15 GeV hybridSoft tag and probe).
// Every pp chain file carries "_mu5" after the purity suffix, purity files included (makeCSVv2PurityFile_pp_PF.C
// with inTag "_mu5"). PFJets only.
inline bool mu5On(){
  static int v = -1;
  if(v < 0){ const char *e = std::getenv("BJET_MU5"); v = (e && TString(e) == "1") ? 1 : 0;
             if(v) printf("MU5 pp CHAIN: mu5 data, templates and response, muon SF 1, outputs suffixed _mu5\n"); }
  return v == 1;
}
inline const char* mu5Var(){ return mu5On() ? "_mu5" : ""; }
// b-jet BAYES-ITERATION VARIANT (2026-11-11, user: "still use the 0.3 tilt, use these new N_opt values"):
// BJET_NITER=bias sets the b-jet N (bNIterPP / bNIterPbPb, bJetUnfoldExt.h) to the 5% plateau on the
// BIAS^2 of the +0.3 tilted-truth closure instead of the MSE: pp 6, PbPb 10 / 8 / 6 / 6 (0-10 ... 50-80%);
// BJET_NITER=ref keeps the nominal N under its own suffix, a same-code reference for the comparison.
// Every pp chain file and every purity file it reads carries "_nIterBias" / "_nIterRef" after the other
// suffixes (link the purity files to the suffixed names first); the PbPb chain takes it through its outTag.
inline int bNIterVariant(){
  static int v = -1;
  if(v < 0){ const char *e = std::getenv("BJET_NITER"); v = !e ? 0 : TString(e) == "bias" ? 1 : TString(e) == "ref" ? 2 : TString(e) == "opt" ? 3 : TString(e) == "bias5" ? 4 : 0;
             if(v) printf("b-jet N VARIANT %s: outputs suffixed %s\n", e, v == 1 ? "_nIterBias (pp 6, PbPb 10/8/6/6)" : v == 3 ? "_nIterOpt (pp 2, PbPb 3/3/2/1)" : v == 4 ? "_nIterBias5 (pp 2, PbPb 6/3/2/1)" : "_nIterRef (nominal N)"); }
  return v;
}
// BJET_NITER=opt (2026-11-16): the MSE 5%-plateau N re-derived with the CURRENT inputs (2D-fit purity errors),
// optimalIterations_bJets_PF_MCprior.C(0., "_csvV2_2D") and (0., "_csvV2_2D7bins"): pp 2, PbPb 3 / 3 / 2 / 1 on both
// the 4-bin and the 7-bin axis (the nominal 3, 2 / 2 / 1 / 1 dates from the 2026-11-05 inputs); suffix "_nIterOpt".
inline const char* nIterVar(){ return bNIterVariant() == 1 ? "_nIterBias" : bNIterVariant() == 2 ? "_nIterRef" : bNIterVariant() == 3 ? "_nIterOpt" : bNIterVariant() == 4 ? "_nIterBias5" : ""; }
// BJET_NITER=bias5 (2026-11-16, user: "Yes, 5% criterion"): the MSE optimum of BJET_NITER=opt with the added requirement that
// every reported bin of the tilt-0.3 closure is within 5% of truth, on the 7-BIN axis (optimalIterations_bJets_PF_MCprior.C,
// "within x% of truth" printout): only 0-10% moves, 3 -> 6 (pp 2, PbPb 6 / 3 / 2 / 1); suffix "_nIterBias5". On the 4-bin
// axis the same rule gives 4 in 0-10%, which this switch does NOT set.
// 7-BIN JET-pT AXIS (2026-11-16; caloJetBPurityFitCore.h): with BJET_PT_7BINS defined a chain macro's histograms are
// on another axis than the nominal's, so its output tag must say so; without it, a tag naming 7 bins is refused too.
inline bool ptBinGuard(const TString &tag){
#ifdef BJET_PT_7BINS
  if(!tag.Contains("7bins")){ printf("ERROR: BJET_PT_7BINS is defined but the output tag \"%s\" does not contain \"7bins\": refusing to write nominal-named files on the 7-bin axis\n", tag.Data()); return false; }
  printf("JET-pT AXIS: 7-bin variant, 80 90 100 110 120 135 150 200 GeV (tag %s)\n", tag.Data());
#else
  if(tag.Contains("7bins")){ printf("ERROR: tag \"%s\" names the 7-bin axis but BJET_PT_7BINS is not defined (root -e '#define BJET_PT_7BINS 1' ...)\n", tag.Data()); return false; }
#endif
  return true;
}
// RECO-FLOOR SOURCE (2026-11-16, user: "Just include the non-closure while using the 60 GeV floor. No extra source needed"):
// the reco-floor systematic (unfolding with the reco floor at 70 instead of 60 GeV) is no longer computed or written by
// bJetsPerZSys_PF_PbPb.C (recoFloor, incl_recoFloor) and inclusiveJetsPerZSys_PF.C (recoFloor). The unfolding non-closure
// (unfoldingBias, tilted-truth closure) is evaluated at the nominal 60 GeV floor and stands for the unfolding's
// uncertainty. The results never carried it since 2026-10-21 (first tools/dropSysSources.C pass, "_noRecoFloor" in the
// file names, kept so every downstream name still resolves; that pass is now a plain copy). BJET_SYS_RECOFLOOR=1 computes
// and writes the source again.
inline bool recoFloorSourceOn(){
  static int v = -1;
  if(v < 0){ const char *e = std::getenv("BJET_SYS_RECOFLOOR"); v = (e && TString(e) == "1") ? 1 : 0;
             if(v) printf("RECO-FLOOR SOURCE ON (BJET_SYS_RECOFLOOR=1): floor 60 -> 70 GeV variation computed and written\n"); }
  return v == 1;
}
// the closure sigma* is used unless a variant asks for the b + c sigma
inline bool useClosureLightSigma(){ return !bcSigmaVariantOn() && !purityVariantOn(); }
inline const char* purVarBase(){ return csvV2Var()[0] ? csvV2Var() : purityVariantOn() ? "_gsAbs0p175fit5" : bcSigmaVariantOn() ? "_bcSigma" : (lightSigStarOn() ? "_lightSigStar" : ""); }
inline const char* purVar(){ static TString s; s = TString(purVarBase()) + mu5Var() + lowFloorVar() + nIterVar(); return s.Data(); }
// suffix of the pp b-spectrum / purity inputs: only the variants that change the pp fits
inline const char* ppPurVar(){ static TString s; s = TString(csvV2Var()[0] ? csvV2Var() : purityVariantOn() ? "_gsAbs0p175fit5" : "") + mu5Var() + lowFloorVar() + nIterVar(); return s.Data(); }
