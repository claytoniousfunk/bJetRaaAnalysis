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
// pThat cross-section weight (em->weight); false fills every event with weight
// 1, as in PYTHIAHYDJET. Output name: _pThat-<cut> or _pThat-unweighted.
bool doPThatWeight = true;
bool doVzReweight = true;
bool doJetPtReweight = false;
// jet-based filers
bool doGenJetPthatFilter = false;
bool doLeadingXjetDumpFilter = false;
bool doXdumpReweight = false;
bool doJetTrkMaxFilter = true;
bool doRemoveHYDJETjet = true;
bool doEtaPhiMask = false;
bool doPThatCorrelationFilter = true;
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
// Multi-sigma jet-axis smearing scan (PYTHIA_scan.C). On top of the nominal
// templates, the muon ptRel templates are refilled with the reco jet axis
// smeared by each sigma below, the SAME sigma in eta and in phi, nominal jet pT
// only. One standard-normal (eta, phi) offset is drawn per jet and scaled by
// each sigma, so neighbouring sigmas differ only by the smearing, not by noise.
// Purpose: pp templates with the axis smeared to PbPb-like resolution, the
// sigma fixed per centrality class and jet-pT window by matching the
// PYTHIA+HYDJET b and c templates and then applied to light (whose PbPb MC
// statistics are too thin). Output names ..._<flavor>_axisSmear<k>; the sigmas
// are stored as h_axisSmearSigma. Tagged _axisSmearScan.
bool doAxisSmearScan = true;
const int NAxisSmear = 9;
const double axisSmearSigma[NAxisSmear] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10};
// ptRel x jet pT x CSVv2 maps of the muon-tagged jets, per flavor, nominal jet pT
// (headers/AnalysisSetup/ptRelDeepCSV3D.h; eventMap::loadCSVv2 reads discr_csvV2 from the
// jet tree). Additive: no existing histogram changes, so the output name is unchanged.
// CSVv2 replaced DeepCSV on 2026-10-22 (the tagger of the CMS b-jet shapes analyses).
bool doCSVv2Map = true;
// The same ptRel x jet pT x CSVv2 map for the b-purity JES / JER systematics of the 2D fit (2026-11-17, user: "JER for the purity
// should only change which jets fill the templates in MC"): filled with the T1 (JER-smeared), T2 (JES up), T3 (JES down) jet pT
// of the 1D templates, so the systematic is a refit with templates from the varied jets instead of the 1D ratio T_var / T0
// applied to the 2D templates. Names ..._csvV2_inclRecoMuonTag_triggerOn_<flavor>[_C<i>]_T<t>. Additive, output name unchanged.
bool doCSVv2JetVarMaps = true;
// Even-indexed events only (evi % 2 == 0) of the T0 and T1 CSVv2 maps, names ..._T0_even, ..._T1_even; the odd half is the
// full map minus the even one. Statistical test of the JER refit: T1 - T0 between independent halves against within one.
// ~0.7 MB per map, 9 flavors: ~19 MB for the three T maps, ~13 MB for the two even maps.
bool doCSVv2EvenOdd = true;
// The same map with the DeepCSV discriminant (discr_deepCSV, eventMap::loadDeepCSV), a cross-check of the CSVv2 purity
// fits, booked NEXT TO the CSVv2 map so pp and PbPb carry both taggers in one scan generation (2026-11-16, user; the
// 2026-10-22 scans had DeepCSV instead of CSVv2). Same jets, weights and axes (first z bin [-2, 0) = jets the tagger did
// not evaluate, stored as -1), names h_muptrel_recoJetPt_deepCSV_inclRecoMuonTag_triggerOn[_<flavor>]. Additive: the output
// name is unchanged. ~0.84 MB per flavor. Without a discr_deepCSV branch the scan prints so and every jet lands in the underflow.
bool doDeepCSVMap = true;
// System8 maps for a data measurement of the CSVv2 b-tag efficiency (headers/AnalysisSetup/system8.h:
// away jet = leading other jet, |eta| < 1.6, pT > 30 GeV, |dphi| > 2.5, tagged if CSVv2 > 0.84; maps in CSVv2
// working-point categories). Additive: no existing histogram changes, output name unchanged. 2026-10-26.
bool doSystem8Maps = true;
// Muon IP3D maps (headers/AnalysisSetup/muonIP3D.h: ln|IP3D| x jet pT x CSVv2 category and ptRel x ln|IP3D| x jet
// pT, for the IP3D template fit of CMS BTV-12-001). Additive: no existing histogram changes, output name unchanged.
// Needs muIP3D in the muon tree (ggHiNtuplizer); the scan aborts if it is missing. 2026-10-27.
bool doMuonIP3DMaps = true;
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
bool fillMu12 = true;    // 2026-11-06: back to mu12 (the nominal templates) for the even / odd template scans
// jet mods
bool useCaloJetsOverride = false;
bool useManualJEC = true;
// PF jets only: drop every reco PF jet with no ak4Calo jet within dR < 0.2
// (src/scanning/scan_calo_jet_match.h). Applies to ALL reco-jet histograms and
// the muon tag; event counts and gen-level histograms are unaffected. The scan
// aborts if the forest has no calo tree. Adds _caloJetMatched to the output name.
bool requireCaloJetMatch = false;
// "reverse unfolding" (2026-10-08), PYTHIA_scan_response.C only: an EVENT weight
// at the leading muon-tagged reco jet's pT, the pp b-jet data / PYTHIA reco ratio
// (headers/fitParameters/leadMuTagEventWeight_PYTHIA_PF.h, generated by
// src/calculateBJetsPerZ/reverseUnfoldReweight_PF_pp.C). Adds
// _leadMuTagEventReweight to the response file name.
bool doLeadMuTagEventReweight = false;   // off: template and ptRel x muon pT scans must not carry it; the reverse unfolding (nominal since late 2026-10-09) is done offline
// keep only even- or odd-indexed events, giving two statistically independent halves:
// the response scan (unfoldClosureTest.C) and, since 2026-11-06, PYTHIA_scan.C (b-purity
// closure; set per job by PYTHIA_scan's parity argument, condor_PYTHIA_scan.py even|odd)
bool onlyEvenEvents = false;
bool onlyOddEvents = false;
// response scan: fill only jets whose jet-branch muon (mupt/mueta, no
// muon-tree matching) passes muPtCut and |eta| < 2
bool onlyMuTaggedJets = false;


