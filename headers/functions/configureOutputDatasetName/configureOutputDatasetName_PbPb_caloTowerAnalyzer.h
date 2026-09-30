// Output naming for PbPb_caloTowerAnalyzer.C.
//
// Copied from configureOutputDatasetName_PbPb_pfCandAnalyzer.h with one
// argument repurposed: the last string is towerTag() from caloTowers.h rather
// than constituentTag() from pseudoJets.h. Everything else is identical, so
// the two scans' filenames stay directly comparable field by field.
#include "TDatime.h"

TString configureOutputDatasetName(bool doSingleMuonSample,
				   bool doMinBiasSample,
				   bool doHardProbesSample,
				   bool doNoRhoModificationSample,
				   bool doWithRhoModificationSample,
				   bool applyMinBiasTrigger,
				   bool applyJet60Trigger,
				   bool applyJet80Trigger,
				   bool applyJet100Trigger,
				   bool applyMu12TriggerEfficiencyCorrection,
				   bool doJetTrkMaxFilter,
				   bool doEtaPhiMask,
				   bool doWDecayFilter,
				   bool doBJetNeutrinoEnergyShift,
				   bool doJERCorrection,
				   bool apply_JER_smear,
				   bool apply_JEU_shift_up,
				   bool apply_JEU_shift_down,
				   double muPtCut,
				   double muPtMaxCut,
				   bool fillMu5,
				   bool fillMu7,
				   bool fillMu12,
				   bool doEventMixing,
				   bool skipSingleConstituentJets,
				   bool doHiBinReweightToHardProbesJet80,
				   bool useCaloJetsOverride,
				   bool useFlowJetsOverride,
				   int N_fastJetMixedEventResamples,
				   double pseudoJetCandPt_min,
				   const char *bkgMapTag,
				   const char *towerTagStr,
				   bool doTowerPUSub)
{

  TString result = "output";
  result.Append("_PbPb");

  TString datasetIndicator = "";
  if(doSingleMuonSample) datasetIndicator = "_SingleMuon";
  else if(doMinBiasSample) datasetIndicator = "_MinBias_Part1";
  else if(doHardProbesSample) datasetIndicator = "_HardProbes";
  else if(doNoRhoModificationSample) datasetIndicator = "_noRhoModification";
  else if(doWithRhoModificationSample) datasetIndicator = "_withRhoModification";
  else{};
  result.Append(datasetIndicator);

  if(useCaloJetsOverride) result.Append("_caloJets");
  if(useFlowJetsOverride) result.Append("_flowJets");

  // general information
  if(applyMinBiasTrigger) result.Append("_MinBiasHLT");
  if(applyJet60Trigger) result.Append("_Jet60HLT");
  if(applyJet80Trigger) result.Append("_Jet80HLT");
  if(applyJet100Trigger) result.Append("_Jet100HLT");

  if(doHiBinReweightToHardProbesJet80) result.Append("_hiBinReweightToHardProbesJet80");

  if(fillMu5) result.Append(Form("_mu5_pTmu-%1.0fto%1.0f_hybridSoft",muPtCut,muPtMaxCut));
  else if(fillMu7) result.Append(Form("_mu7_pTmu-%1.0fto%2.0f_hybridSoft",muPtCut,muPtMaxCut));
  else if(fillMu12){
    if(muPtMaxCut < 100) result.Append(Form("_mu12_pTmu-%2.0fto%2.0f_tight",muPtCut,muPtMaxCut));
    else result.Append(Form("_mu12_pTmu-%2.0fto%3.0f_tight",muPtCut,muPtMaxCut));
  }
  else{};
  if(applyMu12TriggerEfficiencyCorrection) result.Append("_mu12TriggerEfficiencyCorrection");

  // jet-based filters
  if(doJetTrkMaxFilter) result.Append("_jetTrkMaxFilter");
  if(doEtaPhiMask) result.Append("_etaPhiMask");
  if(doWDecayFilter) result.Append("_WDecayFilter");
  if(doBJetNeutrinoEnergyShift) result.Append("_BJetNeutrinoEnergyShift");
  if(doJERCorrection) result.Append("_JERCorrection");
  if(apply_JER_smear) result.Append("_applyJERSmear");
  if(apply_JEU_shift_up) result.Append("_applyJEUShiftUp");
  if(apply_JEU_shift_down) result.Append("_applyJEUShiftDown");

  // pfCand / event-mixing options
  if(doEventMixing){
    result.Append("_mixedEventPFClustering");
    // only meaningful for doEventMixing -- same-event clustering has no
    // randomness to resample, so this is omitted there rather than printing
    // a misleading "_fastJetResamples-N" on a file it did not affect
    result.Append(Form("_fastJetResamples-%i",N_fastJetMixedEventResamples));
  }
  else result.Append("_sameEventPFClustering");
  // PF-candidate pT floor for the random cones and FastJet inputs. Dropped from
  // the name on 2026-09-07 and restored: the RC and dPT maps depend on it
  // directly, and without it a 2 GeV scan is indistinguishable from a 0 GeV one.
  result.Append(Form("_pseudoJetCandPtMin-%.1f", pseudoJetCandPt_min));
  // Tower selection: ET threshold, clustering eta acceptance and any em/had-only
  // variant (towerTag() in caloTowers.h). This is not cosmetic -- 76% of towers
  // sit above the default 0.3 GeV threshold, so two scans at different
  // thresholds cluster substantially different input and must not share a name.
  result.Append(towerTagStr);
  // forest-style tower-level PU subtraction (doTowerPUSub in caloTowers.h). It
  // only adds histograms, but a file without it has them booked and empty, so
  // the name says which kind this is. Tagged only when it actually ran: the
  // scan skips it in mixed-event mode.
  if(doTowerPUSub && !doEventMixing) result.Append("_towerPUSub");
  // which background map the subtracted spectra used (bkgMapTag() in pseudoJets.h)
  result.Append(bkgMapTag);
  if(skipSingleConstituentJets) result.Append("_skipSingleConstituentJets");

  TDatime dt;
  result.Append(Form("_%i-%i-%i",dt.GetYear(),dt.GetMonth(),dt.GetDay()));

  return result;

}
