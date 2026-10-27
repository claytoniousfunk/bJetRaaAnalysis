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
// "reverse unfolding" (2026-10-08), PYTHIAHYDJET_scan_response.C only: an EVENT
// weight at the leading muon-tagged reco jet's pT, the PbPb b-jet data / MC reco
// ratio (headers/fitParameters/leadMuTagEventWeight_PYTHIAHYDJET_PF.h, generated
// by src/calculateBJetsPerZ/reverseUnfoldReweight_PF_PbPb.C). Adds
// _leadMuTagEventReweight to the response file name.
bool doLeadMuTagEventReweight = false;   // off: template and ptRel x muon pT scans must not carry it; the reverse unfolding (nominal since late 2026-10-09) is done offline
// Track-ptRel closure (2026-10-21, user): with this on, every selected track (pT > 15
// GeV, dR < 0.4; PbPbTracks/trackTree or the charged PF candidates,
// headers/functions/jetTrackPtRel.h) of every jet passing the template selection
// fills h_jetTrkPtRel_trkPt_recoJetPt_<flavor>_C<i>, same weight, jet pT and truth
// flavor as the muon 3D h_muptrel_mupt_recoJetPt_inclRecoMuonTag_triggerOn_<flavor>_C<i>
// (no trigger requirement on the tracks). Tests the data-driven light template's
// premise: light-jet tracks and light-jet muons share a ptRel distribution at fixed
// pT. Adds _trkPtRel to the output name. Off for the template scans.
bool doTrackPtRelHists = false;
// ptRel x jet pT x CSVv2 maps of the muon-tagged jets, per centrality index and flavor,
// nominal jet pT (headers/AnalysisSetup/ptRelDeepCSV3D.h, as PYTHIA_scan.C; eventMap::loadCSVv2
// reads discr_csvV2 from akCs4PF / akPu4Calo, which carry it in these forests). Additive: no
// existing histogram changes, so the output name is unchanged. ~0.84 MB per map in memory,
// NCentralityIndices x kNJetFlavors maps (~140 MB with the ultraFine scheme).
// CSVv2 replaced DeepCSV on 2026-10-24 (the tagger of the CMS b-jet shapes analyses).
bool doCSVv2Map = true;
// System8 maps for a data measurement of the CSVv2 b-tag efficiency (headers/AnalysisSetup/system8.h:
// away jet = leading other jet, |eta| < 1.6, pT > 30 GeV, |dphi| > 2.5, tagged if CSVv2 > 0.84; maps in CSVv2
// working-point categories). Additive: no existing histogram changes, output name unchanged. 2026-10-26.
bool doSystem8Maps = true;
// Muon IP3D maps (headers/AnalysisSetup/muonIP3D.h: ln|IP3D| x jet pT x CSVv2 category and ptRel x ln|IP3D| x jet
// pT, for the IP3D template fit of CMS BTV-12-001). Additive: no existing histogram changes, output name unchanged.
// Needs muIP3D in the muon tree (ggHiNtuplizer); the scan aborts if it is missing. 2026-10-27.
bool doMuonIP3DMaps = true;
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
bool doDiJetSample = true;
bool doMuJetSample = false;
bool doBJetSample = false;
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
bool useCaloJetsOverride = false;
bool useManualJEC = true;
// PF jets only: drop every reco PF jet with no akPu4Calo jet within dR < 0.2
// (src/scanning/scan_calo_jet_match.h). Applies to ALL reco-jet histograms and
// the muon tag; event counts and gen-level histograms are unaffected. The scan
// aborts if the forest has no calo tree. Adds _caloJetMatched to the output name.
bool requireCaloJetMatch = false;
// for response scan. The PbPb b-jet unfolding needs BOTH halves (even + odd
// summed for the nominal response, split for the closure test): run once with
// onlyEvenEvents, then again with onlyOddEvents.
bool onlyEvenEvents = false;
bool onlyOddEvents = false;
// response scan: fill only jets whose jet-branch muon (mupt/mueta, no
// muon-tree matching) passes muPtCut and |eta| < 2
bool onlyMuTaggedJets = false;
