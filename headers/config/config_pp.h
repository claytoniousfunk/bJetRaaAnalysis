int isMC_status = 0;
int AASetup_status = 0;
TString jetTreeString = "jetTree";
TString muonTreeString = "muonTree";
TString hltString = "hltTree";
// event filters
const int NeventFilters_SingleMuon = 3;
std::string eventFilters_SingleMuon[NeventFilters_SingleMuon] = {"pPAprimaryVertexFilter","HBHENoiseFilterResultRun2Loose","pBeamScrapingFilter"};
const int NeventFilters_MinBias = 1;
std::string eventFilters_MinBias[NeventFilters_MinBias] = {"pprimaryVertexFilter"};
// jet options
bool doJetTrkMaxFilter = true;
bool doEtaPhiMask = false;
bool doWDecayFilter = true;
// triggers
bool applyJet15Trigger = false;
bool applyJet30Trigger = false;
bool applyJet40Trigger = false;
bool applyJet60Trigger = false;
bool applyJet80Trigger = false;
bool applyJet100Trigger = false;
bool applyAntiMu5Jet30Trigger = false;
bool applyAntiMu5Jet40Trigger = false;
bool applyAntiMu5Jet60Trigger = false;
bool fillMu5 = false;
bool fillMu7 = false;
bool fillMu12 = true;    // 2026-11-16: back to mu12 (HLT_HIL3Mu12, tight ID, muon pT > 15 GeV), as every other scan config
// mu5 only: weight every trigger-gated fill (h_vz_triggerOn included) by the event's
// HLT_HIL3Mu5_NHitQ10_v1 prescale (3 or 5 in the 2017 pp run), so the mu5 spectra are
// luminosity-normalized. Adds _mu5PrescaleWeighted to the output name. 2026-10-29.
bool weightMu5Prescale = true;
// spectra alterations
bool applyMu12TriggerEfficiencyCorrection = true;  // mu12 HLT efficiency (per-muon weight); set false for a mu5 scan
bool doJESCorrection = false;
bool doBJetNeutrinoEnergyShift = false;
bool doJERCorrection = false;
// options for systematic studies
bool apply_JER_smear = false;
bool apply_JEU_shift_up = false;
bool apply_JEU_shift_down = false;
// select dataset
bool doSingleMuonSample = true;
bool doMinBiasSample = false;
bool doHighEGJetSample = false;
// ptRel x jet pT x CSVv2 map of the triggered muon-tagged jets (headers/AnalysisSetup/ptRelDeepCSV3D.h,
// as PYTHIA_scan.C; eventMap::loadCSVv2 reads discr_csvV2 from the jet tree). Additive: no
// existing histogram changes, so the output name is unchanged.
// CSVv2 replaced DeepCSV on 2026-10-22 (the tagger of the CMS b-jet shapes analyses).
bool doCSVv2Map = true;
// The same map with the DeepCSV discriminant (discr_deepCSV, eventMap::loadDeepCSV), a cross-check of the CSVv2 purity
// fits, booked NEXT TO the CSVv2 map so pp and PbPb carry both taggers in one scan generation (2026-11-16, user; the
// 2026-10-22 scans had DeepCSV instead of CSVv2). Same jets, weights and axes (first z bin [-2, 0) = jets the tagger did
// not evaluate, stored as -1), names h_muptrel_recoJetPt_deepCSV_inclRecoMuonTag_triggerOn[_<flavor>]. Additive: the output
// name is unchanged. ~0.84 MB. Without a discr_deepCSV branch the scan prints so and every jet lands in the underflow.
bool doDeepCSVMap = true;
// System8 maps for a data measurement of the CSVv2 b-tag efficiency (headers/AnalysisSetup/system8.h:
// away jet = leading other jet, |eta| < 1.6, pT > 30 GeV, |dphi| > 2.5, tagged if CSVv2 > 0.84; maps in CSVv2
// working-point categories). Additive: no existing histogram changes, output name unchanged. 2026-10-26.
bool doSystem8Maps = true;
// Muon IP3D maps (headers/AnalysisSetup/muonIP3D.h: ln|IP3D| x jet pT x CSVv2 category and ptRel x ln|IP3D| x jet
// pT, for the IP3D template fit of CMS BTV-12-001). Additive: no existing histogram changes, output name unchanged.
// Needs muIP3D in the muon tree (ggHiNtuplizer); the scan aborts if it is missing. 2026-10-27.
bool doMuonIP3DMaps = true;
// jet mods
bool useCaloJetsOverride = false;
bool useManualJEC = true;
// Drop the data-only L2L3Residual from the calo-jet JEC (useCaloJetsOverride &&
// useManualJEC only). The AK4Calo residual files in JetEnergyCorrections/ are
// byte-identical copies of the AK4PF ones, so calo jets get a PF-derived
// data/MC correction: +2.3-2.9% in pp (Spring18_ppRef5TeV_V6) against +0.8-1.0%
// in PbPb (Autumn18_HI_V8) at 200-450 GeV, |eta| < 1.6. The MC responses get
// no residual, so that offset survives unfolding and raises pp relative to
// PbPb, by roughly (1.029/1.010)^5 ~ 1.10 in R_AA terms. Flip both pp and PbPb
// together, or the comparison mixes conventions. Adds _noL2L3Residual to the
// output name.
bool skipCaloL2L3Residual = false;
// PF jets only: drop every reco PF jet with no ak4Calo jet within dR < 0.2
// (src/scanning/scan_calo_jet_match.h). Applies to ALL reco-jet histograms and
// the muon tag; event counts and gen-level histograms are unaffected. The scan
// aborts if the forest has no calo tree. Adds _caloJetMatched to the output name.
bool requireCaloJetMatch = false;
