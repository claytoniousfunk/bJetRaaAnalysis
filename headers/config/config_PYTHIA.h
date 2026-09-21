// general analysis options
int isMC_status = 1;
int AASetup_status = 0;
TString jetTreeString = "jetTree";
TString muonTreeString = "muonTree";
TString hltString = "hltTree";
const int NeventFilters = 3;
std::string eventFilters[NeventFilters] = {"pPAprimaryVertexFilter","HBHENoiseFilterResultRun2Loose","pBeamScrapingFilter"};
// options for systematic studies
bool apply_JER_smear = false;
bool apply_JEU_shift_up = false;
bool apply_JEU_shift_down = false;
// options to skip genParticle loops (for forests with no genParticle info)
bool skipGenParticles = false;
// reweighting functions
bool doVzReweight = true;
bool doJetPtReweight = false;
// jet-based filers
bool doGenJetPthatFilter = false;
bool doLeadingXjetDumpFilter = false;
bool doXdumpReweight = false;
bool doJetTrkMaxFilter = true;
bool doRemoveHYDJETjet = true;
bool doEtaPhiMask = false;
bool doPThatCorrelationFilter = false;
// jet-energy alterations
bool doBJetEnergyShift = false;
bool doBJetNeutrinoEnergyShift = false;
bool doNeutrinoEnergyAddition = false;
bool doInclNeutrinoEnergySubtraction = false;
bool apply_genJetPt_smear = false;
bool doJERCorrection = false;
bool doJESCorrection = false;
// muptrel reweight
bool doMuPtRelReweight = false;
bool doHadronPtRelReweight = false;
bool doHadronPtRelReweightToMuon = false;
bool doDRReweight = false;
bool doWeightCut = false;
bool doJetAxisSmearing = false;
bool doWDecayFilter = true;
// select dataset
bool doDiJetSample = true;
// stuff for dataset naming
TString generator = "PYTHIA";
// jet-axis smearing parameters
double mu_phi = 0.0;
double sigma_phi = 0.005;
double mu_eta = 0.0;
double sigma_eta = 0.005;
// triggers
bool applyJet60Trigger = false;
bool applyJet80Trigger = false;
bool applyJet100Trigger = false;
bool applyMu12TriggerEfficiencyCorrection = false;
bool applyMu5Jet30Trigger = false;
bool applyMu5Jet40Trigger = false;
bool applyMu5Jet60Trigger = false;
bool fillMu5 = false;
bool fillMu7 = false;
bool fillMu12 = true;
// jet mods
bool useCaloJetsOverride = true; // overrides typical PF jet collection with caloJets
// Correct rawpt with the JEC text files (calo or PF, matching the collection)
// instead of taking the forest's jtpt. The forest corrects calo jets with the
// AK4PF payload (ak4CaloJetSequence_pp_mc_cff.py: payload = "AK4PF"), leaving
// them 22% low against gen; AK4Calo L2Relative on rawpt closes at 1.014. For PF
// the two agree exactly. Adds "_manualJEC" to the output name.
bool useManualJEC = true;


