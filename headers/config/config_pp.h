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
bool applyJet80Trigger = true;
bool applyJet100Trigger = false;
bool applyAntiMu5Jet30Trigger = false;
bool applyAntiMu5Jet40Trigger = false;
bool applyAntiMu5Jet60Trigger = false;
bool fillMu5 = false;
bool fillMu7 = false;
bool fillMu12 = true;
// spectra alterations
bool applyMu12TriggerEfficiencyCorrection = false;
bool doJESCorrection = false;
bool doBJetNeutrinoEnergyShift = false;
bool doJERCorrection = false;
// options for systematic studies
bool apply_JER_smear = false;
bool apply_JEU_shift_up = false;
bool apply_JEU_shift_down = false;
// select dataset
bool doSingleMuonSample = false;
bool doMinBiasSample = false;
bool doHighEGJetSample = true;
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
bool requireCaloJetMatch = true;
