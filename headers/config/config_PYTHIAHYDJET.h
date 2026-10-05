// general analysis options
int isMC_status = 1;
int AASetup_status = 1;
TString jetTreeString = "jetTree";
TString muonTreeString = "muonTree";
TString hltString = "hltTree";
const int NeventFilters = 5;
std::string eventFilters[NeventFilters] = {"pprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "collisionEventSelectionAODv2", "phfCoincFilter2Th4", "pclusterCompatibilityFilter"};
// options for systematic studies
bool apply_JER_smear = false;
bool apply_JEU_shift_up = false;
bool apply_JEU_shift_down = false;
// options to skip genParticle loops (for forests with no genParticle info)
// TRUE only for forests without HiGenParticleAna/hi (the old MuJet / BJet
// lists). Leave false for the fix2 forests: the gen-jet muon tags and the
// tag frequency need the gen particles.
bool skipGenParticles = false;
// reweighting functions
bool doPThatWeight = true;
bool doHiBinReweight = false;
bool doHiBinReweightToHardProbesJet80 = false;
bool doVzReweight = true;
bool doJetPtReweight = false;
// jet-based filters
bool doGenJetPthatFilter = false;
bool doLeadingXjetDumpFilter = false;
bool doXdumpReweight = false;
bool doJetTrkMaxFilter = true;
bool doRemoveHYDJETjet = false;
bool doEtaPhiMask = false;
bool doBJetEnergyShift = false;
bool doJERCorrection = false;
bool doJESCorrection = false;
bool doHadronPtRelReweight = false;
bool doHadronPtRelReweightToMuon = false;
bool doDRReweight = false;
bool doWeightCut = false;
bool doJetAxisSmearing = false;
bool doWDecayFilter = true;
bool doBJetSpectraReweightToData = false;
bool doPThatCorrelationFilter = true;
// drop reco jets with no matched gen jet (refpt < 0): the HYDJET background fakes.
// Their count is pThat-independent but they carry the pThat weight of the event,
// so a few low-pThat events dominate the error at low pT; the fake background is
// subtracted from data instead. Only PYTHIAHYDJET_scan.C reads this.
bool doRemoveUnmatchedJets = true;
// shifting the hiBin distribution by this amount
int hiBinShift = 10;
// jet-alterations for closure
bool doBJetNeutrinoEnergyShift = false;
// experimental filters
bool onlyOneMuonTaggedJetPerEvent = false;
// jet-axis smearing parameters
double mu_phi = 0.0;
double sigma_phi = 0.005;
double mu_eta = 0.0;
double sigma_eta = 0.005;
// stuff for dataset naming
TString generator = "PYTHIAHYDJET";
// data set
bool doDiJetSample = false;
bool doMuJetSample = false;
bool doBJetSample = true;
bool doDiJetSample_batch1 = false;
bool doDiJetSample_batch2 = false;
bool doDiJetSample_batch3 = false;
bool doDiJetSample_batch4 = false;
bool doDiJetSample_batch5 = false;
bool doDiJetSample_batch6 = false;
bool doDiJetSample_batch7 = false;
bool doDiJetSample_batch8 = false;
bool doDiJetSample_batch9 = false;
bool doDiJetSample_batch10 = false;
bool doDiJetSample_batch11 = false;
bool doDiJetSample_batch12 = false;
bool doDiJetSample_batch13 = false;
bool doDiJetSample_batch14 = false;
bool doDiJetSample_batch15 = false;
// jet triggers
bool applyJet60Trigger = false;
bool applyJet80Trigger = false;
bool applyMu12TriggerEfficiencyCorrection = false;
bool fillMu5 = false;
bool fillMu7 = false;
bool fillMu12 = true;
bool useCaloJetsOverride = true;
bool useManualJEC = true;
// for response scan. The PbPb b-jet unfolding needs BOTH halves (even + odd
// summed for the nominal response, split for the closure test): run once with
// onlyEvenEvents, then again with onlyOddEvents.
bool onlyEvenEvents = true;
bool onlyOddEvents = false;
// response scan: fill only jets whose jet-branch muon (mupt/mueta, no
// muon-tree matching) passes muPtCut and |eta| < 2
bool onlyMuTaggedJets = false;
