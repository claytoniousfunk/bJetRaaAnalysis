// general ROOT/C includes
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TLatex.h>
#include <TCut.h>
#include "TDatime.h"
#include <vector>
#include "TCanvas.h"
#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <random>
#include <memory>
#include "TSystem.h"
#include "TRandom2.h"

// FastJet
#ifdef DO_FASTJET
#include "fastjet/ClusterSequence.hh"
#endif
#include "../../../headers/fastJetCandInfo.h"

// event map
#include "../../../eventMap/eventMap.h"
// jet corrector
#include "../../../JetEnergyCorrections/JetCorrector.h"
// jet uncertainty
#include "../../../JetEnergyCorrections/JetUncertainty.h"
// general analysis variables
#include "../../../headers/config/config_centrality.h"
// JERCorrection params
#include "../../../headers/fitParameters/JERCorrectionParams_PYTHIA_mu12.h"
TF1 *fitFxn_PYTHIA_JERCorrection;
// JER-correction function
#include "../../../headers/fitFunctions/fitFxn_PYTHIA_JERCorrection.h"

// HLT fit params/fxn
#include "../../../headers/fitParameters/HLTFitParams_PbPb.h"
TF1 *fitFxn_PbPb_HLT_C4, *fitFxn_PbPb_HLT_C3, *fitFxn_PbPb_HLT_C2, *fitFxn_PbPb_HLT_C1;
#include "../../../headers/fitFunctions/fitFxn_PbPb_HLT.h"

// eta-phi mask function
#include "../../../headers/functions/etaPhiMask.h"
// getDr function
#include "../../../headers/functions/getDr.h"
// getJetPtBin function
#include "../../../headers/functions/getJetPtBin.h"
// getPtRel function
#include "../../../headers/functions/getPtRel.h"
// isQualityMuon_hybridSoft function
#include "../../../headers/functions/isQualityMuon_hybridSoft.h"
// isQualityMuon_tight function
#include "../../../headers/functions/isQualityMuon_tight.h"
// isWDecayMuon function
#include "../../../headers/functions/isWDecayMuon.h"
// isWDecayMuon_raw function (input is rawJetPt)
#include "../../../headers/functions/isWDecayMuon_raw.h"
// triggerIsOn function
#include "../../../headers/functions/triggerIsOn.h"
// pthat filter function
#include "../../../headers/functions/passesLeadingGenJetPthatFilter.h"
// JetTrkMax filter function
#include "../../../headers/functions/jet_filter/passesJetTrkMaxFilter.h"
// print introduction
#include "../../../headers/introductions/printIntroduction_PbPb_scan_V3p7.h"
// analysis config
//#include "../../../headers/config/config_PbPb_SingleMuon.h"
//#include "../../../headers/config/config_PbPb_MinBias.h"
//#include "../../../headers/config/config_PbPb_diJet.h"
#include "../../../headers/config/config_PbPb.h"
// read config
#include "../../../headers/config/readConfig.h"
// dataset naming functions
#include "../../../headers/functions/getDatasetName/getDatasetName_PbPb.h"
#include "../../../headers/functions/getInputFileName/getInputFileName_PbPb.h"
#include "../../../headers/functions/configureOutputDatasetName/configureOutputDatasetName_PbPb_pfCandAnalyzer.h"
// dimuon mass calculation
#include "../../../headers/functions/calculateDimuonMass.h"
// shared scanning helpers
#include "../scan_jet_corrections.h"
#include "../scan_muon_tag.h"
// pfCand analysis variables
#include "../../../headers/AnalysisSetup/pseudoJets.h"
// hibin fit parameters / functions
#include "../../../headers/fitParameters/hiBinFitParams_PYTHIAHYDJET.h"
TF1 *fitFxn_hiBin;
#include "../../../headers/fitFunctions/fitFxn_hiBin.h"


// initialize histograms
// ~~~~~~~~~ event variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH1D *h_eventsBeforeSelection;
TH1D *h_eventsAfterSelection;
TH1D *h_NEvents;
// -----------------------------------------inclusive events ------- -------
TH1D *h_vz[NCentralityIndices];
TH1D *h_hiBin;
// -----------------------------------------events w/ muTrigger ------- -------
TH1D *h_vz_triggerOn[NCentralityIndices];
TH1D *h_hiBin_triggerOn;
// -----------------------------------------events w/ jetTrigger ------- -------
TH1D *h_vz_jet60[NCentralityIndices];
TH1D *h_hiBin_jet60;
TH1D *h_vz_jet80[NCentralityIndices];
TH1D *h_hiBin_jet80;
TH1D *h_vz_jet100[NCentralityIndices];
TH1D *h_hiBin_jet100;
TH1D *h_vz_jet100_withJetAbove130[NCentralityIndices];
TH1D *h_hiBin_jet100_withJetAbove130;
// -----------------------------------------events w/ incl. reco jet -------
TH1D *h_vz_jet[NCentralityIndices];
TH1D *h_hiBin_jet;
// ----------------------- events w/ incl. reco jet + incl. reco muon tag --
TH1D *h_vz_inclRecoMuonTag[NCentralityIndices];
TH1D *h_hiBin_inclRecoMuonTag;
// ---------- events w/ incl. reco jet + incl. reco muon tag + trigger on --
TH1D *h_vz_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH1D *h_hiBin_inclRecoMuonTag_triggerOn;
// ~~~~~~~~~ jet variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ----------------------------------------- incl. reco jets --------------
TH1D *h_inclRecoJetPt[NCentralityIndices];
TH1D *h_inclRawJetPt[NCentralityIndices];
TH1D *h_inclRecoJetEta[NCentralityIndices];
TH1D *h_inclRecoJetPhi[NCentralityIndices];
TH2D *h_inclRecoJetPt_inclRecoJetEta[NCentralityIndices];
TH2D *h_inclRecoJetPt_inclRecoJetPhi[NCentralityIndices];
TH2D *h_inclRecoJetEta_inclRecoJetPhi[NCentralityIndices][NJetPtIndices];
// ------------------------------ incl. reco jets + incl. reco muon tag ----
TH1D *h_inclRecoJetPt_inclRecoMuonTag[NCentralityIndices];
TH1D *h_inclRecoJetEta_inclRecoMuonTag[NCentralityIndices];
TH1D *h_inclRecoJetPhi_inclRecoMuonTag[NCentralityIndices];
TH2D *h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag[NCentralityIndices];
TH2D *h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag[NCentralityIndices];
TH2D *h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[NCentralityIndices][NJetPtIndices];
// ------------------ incl. reco jets + incl. reco muon tag + trigger on ----
TH1D *h_inclRecoJetPt_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH1D *h_inclRecoJetEta_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH1D *h_inclRecoJetPhi_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH2D *h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH2D *h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH2D *h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[NCentralityIndices][NJetPtIndices];
// ~~~~~~~~~ muon variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH1D *h_inclMuPt;
TH1D *h_nMu_triggerOn;
TH2D *h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH2D *h_mupt_recoJetPt_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH2D *h_mueta_recoJetPt_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH2D *h_muphi_recoJetPt_inclRecoMuonTag_triggerOn[NCentralityIndices];
TH2D *h_muJetDr_recoJetPt[NCentralityIndices];
// ~~~~~~~~~ analysis histograms ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH1D *h_muPtOverJetPt[NCentralityIndices][NJetPtIndices];
TH2D *h_mupt_jetpt[NCentralityIndices];
TH2D *h_mupt_muptrel[NCentralityIndices][NJetPtIndices];
TH2D *h_muptrel_jetpt[NCentralityIndices];
TH1D *h_NJetPerEvent[NCentralityIndices];
TH1D *h_NMuTaggedJetPerEvent[NCentralityIndices];
TH2D *h_muptrel_hiBin[NJetPtIndices];
TH1D *h_dimuonMass[NCentralityIndices];
TH1D *h_dimuonMass_sameSign[NCentralityIndices];
TH1D *h_Jet60HLT[NCentralityIndices];
TH1D *h_Jet60HLT_Prescale[NCentralityIndices];
TH1D *h_Jet80HLT[NCentralityIndices];
TH1D *h_Jet80HLT_Prescale[NCentralityIndices];
TH1D *h_Jet100HLT[NCentralityIndices];
TH1D *h_Jet100HLT_Prescale[NCentralityIndices];
// ~~~~~~~~~ pfCand / FastJet histograms ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH1D     *h_pfPt[NCentralityIndices];
TH1D     *h_pseudoJetPt[NCentralityIndices];
TH1D     *h_pseudoJetPt_geoCorr[NCentralityIndices];
TH1D     *h_fastJetPt_PF[NCentralityIndices];
TH1D     *h_fastJetPt_PF_JEC[NCentralityIndices];
TH1D     *h_fastJetPt_PFCs[NCentralityIndices];
TH1D     *h_fastJetPt_PFCs_JEC[NCentralityIndices];
TH1D     *h_nPFcand[NCentralityIndices];
TH1D     *h_nPFcandCS[NCentralityIndices];
TH1D     *h_nPFcandFastJet[NCentralityIndices];  // candidates actually handed to FastJet
TProfile2D *h_dPTEtaPhi_PF_PFCs[NCentralityIndices];
// Same as h_dPTEtaPhi_PF_PFCs, but the per-jet-pair dPT is weighted by
// cos(dR_min) between the matched PF and PFCs jets, mirroring the RC map's
// geoCorr treatment (pseudoJetPt_geoCorr_k = pfPt_l*cos(dR)).
TProfile2D *h_dPTEtaPhi_PF_PFCs_geoCorr[NCentralityIndices];
TProfile2D *h_dPTEtaPhi_PF_PFCs_dPTAbove0[NCentralityIndices];
TProfile2D *h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60[NCentralityIndices];
TH1D       *h_dRmin_PF_PFCs[NCentralityIndices];
TH2D     *h_fastJetPtPF_dRmin[NCentralityIndices][NJetPtIndices];
TH2D     *h_fastJetPtPF_etaPFCs[NCentralityIndices][NJetPtIndices];
// Random-cone eta/phi map: mean cone pT at each (eta,phi) cell.
// Cell size ~0.1x0.1 — sub-cone resolution, enough stats per cell.
// Used to build a local UE background map for FastJet jet subtraction.
const int    NRC_EtaBins = 32;   // -1.6 to 1.6, width 0.1
const int    NRC_PhiBins = 64;   // -pi to pi,   width ~0.098
TProfile2D  *h_randConeEtaPhi[NCentralityIndices];
TProfile2D  *h_randConeEtaPhi_geoCorr[NCentralityIndices];
// Same-event RC, restricted to events passing the signal jet-pT selection
// (see doSignalSelectedRC in pseudoJets.h). Only ever filled when
// !doEventMixing; stay empty (and unwritten) otherwise.
TH1D        *h_pseudoJetPt_sigSel[NCentralityIndices];
TProfile2D  *h_randConeEtaPhi_sigSel[NCentralityIndices];

TH1D        *h_fastJetPt_PF_bkgSub_RC[NCentralityIndices];
TH1D        *h_fastJetPt_PF_JEC_bkgSub_RC[NCentralityIndices];

// Same as h_fastJetPt_PF_bkgSub_RC, restricted to events passing the signal
// jet-pT selection (see doSignalSelectedRC in pseudoJets.h). Only ever
// filled when !doEventMixing; stay empty (and unwritten) otherwise -- same
// caveat as h_pseudoJetPt_sigSel.
TH1D        *h_fastJetPt_PF_bkgSub_RC_sigSel[NCentralityIndices];
TH1D        *h_fastJetPt_PF_JEC_bkgSub_RC_sigSel[NCentralityIndices];

TH1D        *h_fastJetPt_PF_bkgSub_RC_geoCorr[NCentralityIndices];
TH1D        *h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr[NCentralityIndices];

TH1D        *h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect[NCentralityIndices];
TH1D        *h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect[NCentralityIndices];

TH1D        *h_fastJetPt_PF_bkgSub_dPT[NCentralityIndices];
TH1D        *h_fastJetPt_PF_JEC_bkgSub_dPT[NCentralityIndices];

TH1D        *h_fastJetPt_PF_bkgSub_dPT_geoCorr[NCentralityIndices];
TH1D        *h_fastJetPt_PF_JEC_bkgSub_dPT_geoCorr[NCentralityIndices];

TH1D        *h_fastJetPt_PF_bkgSub_dPT_dPTAbove0[NCentralityIndices];
TH1D        *h_fastJetPt_PF_JEC_bkgSub_dPT_dPTAbove0[NCentralityIndices];

TH1D        *h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60[NCentralityIndices];
TH1D        *h_fastJetPt_PF_JEC_bkgSub_dPT_PFCsPTAbove60[NCentralityIndices];

TH2D *h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC[NCentralityIndices];
TH2D *h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC[NCentralityIndices];
TH2D *h_fastJetMuonDR_inclusiveClosestFastJet[NCentralityIndices];
TH2D *h_muonDR_inclusiveClosestJet[NCentralityIndices];
// Trigger-gated twin of the above. Identical fill, plus evtTriggerDecision.
// The ungated original is left as it is because other macros already read it on
// that convention; this exists so the dR distribution can be compared against
// h_mixedMuonPtRel_recoJetPt, whose fill IS gated. Without a matching pair the
// two are normalised over different event samples -- the trigger fraction runs
// from 0.12 in 0-10% to 0.025 in 50-80%, and the fraction among muon-near-jet
// PAIRS is different again (0.16 central to 0.93 peripheral, because a 20 GeV
// jet near a muon is common in central events and needs a hard scattering in
// peripheral ones), so the two cannot be reconciled after the fact.
// Divide this one by h_vz_triggerOn, not h_vz.
TH2D *h_muonDR_inclusiveClosestJet_triggerOn[NCentralityIndices];

// --- muon-injection study (Olga, 2026-09-08) ---------------------------------
// A mixed-event muon matched to donor-event fastJets two ways:
//   noInject       donor candidates clustered ALONE, muon matched afterwards.
//                  This is what the existing templates do, and the jet axis and
//                  pT know nothing about the muon.
//   inject         the muon is added to the donor candidates and everything is
//                  reclustered, so the muon pulls the axis and adds its pT --
//                  what would happen if that muon really were in the event.
// The difference between them is the bias currently carried by
// h_mixedMuonPtRel_recoJetPt: a muon-tagged jet is under-represented at high
// jet pT because the tagging muon's pT was never added to it.
//
// Two tagging conventions are filled for the injected case: nearest jet within
// dR < epsilon_mm (matching how data tags, and like-for-like with noInject),
// and the jet that anti-kT actually assigned the muon to (Constit). The gap
// between them measures how often the nearest jet is not the owning jet.
//
// Inclusive donor jet pT for both collections, filled per injection so the two
// share a denominator and any difference is the injection alone.
TH1D *h_donorJetPt_noInject[NCentralityIndices];
TH1D *h_donorJetPt_inject[NCentralityIndices];
TH2D *h_injMuonDR_donorJetPt_noInject[NCentralityIndices];
TH2D *h_injMuonPtRel_donorJetPt_noInject[NCentralityIndices];
TH2D *h_injMuonDR_donorJetPt_inject[NCentralityIndices];
TH2D *h_injMuonPtRel_donorJetPt_inject[NCentralityIndices];
TH2D *h_injMuonDR_donorJetPt_injectConstit[NCentralityIndices];
TH2D *h_injMuonPtRel_donorJetPt_injectConstit[NCentralityIndices];


// --- ptRel background templates for the muon-tagged-jet decomposition ---
// The measured data ptRel (h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn) is a
// sum of four pairings:
//   (real mu, real jet)   the signal we want to isolate
//   (real mu, fake jet)   -> h_realMuonPtRel_mixedFastJetPt   (mixed-event jet)
//   (fake mu, real jet)   -> h_mixedMuonPtRel_recoJetPt       (mixed-event muon)
//   (fake mu, fake jet)   -> h_fastJetMuonPtRel_..._bkgSub_RC (already existed)
// The two new ones require doEventMixing and are filled ONLY inside the
// mixed-event branch; they stay empty in a same-event scan.
//
// Every term of the decomposition -- these two, the doubly-fake RC histogram,
// and the data -- carries JEC-corrected, background-subtracted jet pT. The RC
// histogram used to be filled on the raw rcSub scale; that was corrected so a
// given jet pT window selects the same jets on both sides of the subtraction.
// Output files produced before that change hold raw-scale contents under the
// same histogram name, so they cannot be mixed with newer ones.
//
// All apply the SAME dR < epsilon_mm tagging criterion the real analysis uses
// (findRecoMuonTag), so the per-event rates are directly comparable to data.
// Note h_muonDR_inclusiveClosestJet, by contrast, has NO dR cut -- it records
// the closest jet at any distance -- so it cannot be reused for this.
TH2D *h_mixedMuonPtRel_recoJetPt[NCentralityIndices];
TH2D *h_realMuonPtRel_mixedFastJetPt[NCentralityIndices];


// RC eta/phi maps loaded from external file at run time
TProfile2D  *h_RC_map[NCentralityIndices];
TProfile2D  *h_RC_geoCorr_map[NCentralityIndices];
TProfile2D  *h_dPT_map[NCentralityIndices];
TProfile2D  *h_dPT_geoCorr_map[NCentralityIndices];
TProfile2D  *h_dPT_dPTAbove0_map[NCentralityIndices];
TProfile2D  *h_dPT_PFCsPTAbove60_map[NCentralityIndices];

///////////////////////  start the program
void PbPb_pfCandAnalyzer(int group = 1){

  if(fillMu5){
    muPtCut = 7.0;
    muPtMaxCut = 9.0;
  }
  else if(fillMu7){
    muPtCut = 9.0;
    muPtMaxCut = 15.0;
  }
  else if(fillMu12){
    muPtCut = 15.0;
    //muPtCut = 20.0;
    //muPtMaxCut = 60.0;
    muPtMaxCut = 999.0;
  }
  else{};



  std::string inputFileList = "";
  if(doSingleMuonSample) inputFileList = "../../../fileNames/fileNames_HISingleMuon_withPFCandidates_partial.txt";
  else if(doMinBiasSample) inputFileList = "../../../fileNames/fileNames_HIMinimumBias0_Part1_withTracksAndPFCandidates.txt";
  else if(doHardProbesSample) {
    inputFileList = "";
    std::cout << "no withPFCandidates forest for HardProbes.  Exiting...\n";
    return;
  }
    
  else if(doNoRhoModificationSample) inputFileList = "../../../fileNames/fileNames_HIMinimumBias0_Part1_noRhoModulation.txt";
  else if(doWithRhoModificationSample) inputFileList = "../../../fileNames/fileNames_HIMinimumBias0_Part1_withRhoModulation.txt";
  else{};

  std::ifstream instr(inputFileList.c_str(), std::ifstream::in);
  if(!instr.is_open()){
    cout << "filelist not found!! Exiting..." << endl;
    return;
  }
  std::string filename;
  Int_t ifile = 0;

  while(instr>>filename){

    ifile++;

    if(ifile != group) continue;

    std::string input = filename.c_str();

  
    // TString inputDataset = "";
    // TString inputFileName = "";

    // inputDataset = getDatasetName(doSingleMuonSample,
    // 				  doMinBiasSample,
    // 				  doHardProbesSample);

    // inputFileName = getInputFileName(doSingleMuonSample,
    // 				     doMinBiasSample,
    // 				     doHardProbesSample);

    // TString input = Form("%s%s_%i.root",inputDataset.Data(),inputFileName.Data(),group);
    // std::cout << "input dataset = " << input << std::endl;


    TString outputBaseDir = "/eos/cms/store/group/phys_heavyions/cbennett/scanningOutput/";
    TString outputDatasetName = "";
    outputDatasetName = configureOutputDatasetName(doSingleMuonSample,
						   doMinBiasSample,
						   doHardProbesSample,
						   doNoRhoModificationSample,
						   doWithRhoModificationSample,
						   applyMinBiasTrigger,
						   applyJet60Trigger,
						   applyJet80Trigger,
						   applyJet100Trigger,
						   applyMu12TriggerEfficiencyCorrection,
						   doJetTrkMaxFilter,
						   doEtaPhiMask,
						   doWDecayFilter,
						   doBJetNeutrinoEnergyShift,
						   doJERCorrection,
						   apply_JER_smear,
						   apply_JEU_shift_up,
						   apply_JEU_shift_down,
						   muPtCut,
						   muPtMaxCut,
						   fillMu5,
						   fillMu7,
						   fillMu12,
						   doEventMixing,
						   skipSingleConstituentJets,
						   doHiBinReweightToHardProbesJet80,
						   useCaloJetsOverride,
						   useFlowJetsOverride,
						   N_fastJetMixedEventResamples);


    TString suffixEdit = CENT_SCHEME_SUFFIX;

    TString outputDir = Form("%s%s%s",outputBaseDir.Data(),outputDatasetName.Data(),suffixEdit.Data());

    TString output = Form("%s/PbPb_pfCandAnalyzer_output_%i.root",outputDir.Data(),group);

    std::cout << "output dataset = " << output << std::endl;

    // Create the output directory if it is not already there. With one condor
    // job per input file many jobs race to create the same directory, so a
    // failed mkdir is not itself an error -- another job winning the race looks
    // identical to a real failure. Only give up if the directory is still
    // missing after the attempt.
    if(gSystem->AccessPathName(outputDir)){
      std::cout << " Output directory not found, creating: " << outputDir << std::endl;
      gSystem->mkdir(outputDir, kTRUE);
    }

    if(gSystem->AccessPathName(outputDir)){
      std::cout << "\033[1;31m Could not create output directory: \033[0m " << outputDir << std::endl;
      return;
    }
  

  
    // JET ENERGY CORRECTIONS
    vector<string> Files_PF;
    Files_PF.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_DATA_L2Relative_AK4PF.txt"); // L2Relative correction
    Files_PF.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_DATA_L2L3Residual_AK4PF.txt"); // L2L3Residual correction
    JetCorrector JEC_PF(Files_PF);
    JetUncertainty JEU_PF("../../../JetEnergyCorrections/Autumn18_HI_V8_MC_Uncertainty_AK4PF.txt");

    vector<string> Files_Calo;
    Files_Calo.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_DATA_L2Relative_AK4Calo.txt");
    Files_Calo.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_DATA_L2L3Residual_AK4Calo.txt");
    JetCorrector JEC_Calo(Files_Calo);
    JetUncertainty JEU_Calo("../../../JetEnergyCorrections/Autumn18_HI_V8_MC_Uncertainty_AK4Calo.txt");
    /// print out some info
    readConfig();


    // Load random-cone UE maps (used for RC-subtracted fastJet pT)
    TFile *f_RC_maps = nullptr;

    if(CENT_SCHEME == CENT_NOMINAL){
      std::cout << "ERROR:  No background maps created for CENT_NOMINAL scheme yet (as of 2026-08-26).  Exiting...\n";
      return;
    }
    else if(CENT_SCHEME == CENT_ULTRAFINE){
      f_RC_maps = TFile::Open("/eos/cms/store/group/phys_heavyions/cbennett/maps/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-17_ultraFineCentBins.root");
    }
    else if(CENT_SCHEME == CENT_PERIPH90){
      std::cout << "ERROR:  No background maps created for CENT_PERIPH90 scheme yet (as of 2026-08-26).  Exiting...\n";
      return;
    }

    
    for(int i = 0; i < NCentralityIndices; i++){

      f_RC_maps->GetObject(Form("h_randConeEtaPhi_C%i",i), h_RC_map[i]);
      f_RC_maps->GetObject(Form("h_randConeEtaPhi_geoCorr_C%i",i), h_RC_geoCorr_map[i]);
      f_RC_maps->GetObject(Form("h_dPTEtaPhi_PF_PFCs_C%i",i), h_dPT_map[i]);
      f_RC_maps->GetObject(Form("h_dPTEtaPhi_PF_PFCs_geoCorr_C%i",i), h_dPT_geoCorr_map[i]);
      f_RC_maps->GetObject(Form("h_dPTEtaPhi_PF_PFCs_dPTAbove0_C%i",i), h_dPT_dPTAbove0_map[i]);
      f_RC_maps->GetObject(Form("h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60_C%i",i), h_dPT_PFCsPTAbove60_map[i]);
				
     
      // if(useDeltaPTMapsForBkgSub) f_RC_maps->GetObject(Form("h_dPTEtaPhi_PF_PFCs_C%i", i), h_RC_map[i]);
      // else{
      // 	if(useGeoCorrForRCMap) f_RC_maps->GetObject(Form("h_randConeEtaPhi_geoCorr_C%i", i), h_RC_map[i]);
      // 	else f_RC_maps->GetObject(Form("h_randConeEtaPhi_C%i", i), h_RC_map[i]);
	
      // }
      
      if(h_RC_map[i]) h_RC_map[i]->SetDirectory(nullptr);
      if(h_RC_geoCorr_map[i]) h_RC_geoCorr_map[i]->SetDirectory(nullptr);
      if(h_dPT_map[i]) h_dPT_map[i]->SetDirectory(nullptr);
      if(h_dPT_geoCorr_map[i]) h_dPT_geoCorr_map[i]->SetDirectory(nullptr);
      if(h_dPT_dPTAbove0_map[i]) h_dPT_dPTAbove0_map[i]->SetDirectory(nullptr);
      if(h_dPT_PFCsPTAbove60_map[i]) h_dPT_PFCsPTAbove60_map[i]->SetDirectory(nullptr);
      
    }

    for(int i = 0; i < NCentralityIndices; i++){
      if(!h_RC_map[i]) printf("WARNING: h_randConeEtaPhi_C%i missing from the RC map file\n", i);
      if(!h_RC_map[i]) printf("WARNING: h_randConeEtaPhi_C%i MISSING\n", i);
      if(!h_RC_geoCorr_map[i]) printf("WARNING: h_randConeEtaPhi_geoCorr_C%i MISSING\n", i);
      if(!h_dPT_map[i]) printf("WARNING: h_dPTEtaPhi_PF_PFCs_C%i MISSING\n", i);

    }

    // define histograms
    h_eventsBeforeSelection = new TH1D("h_eventsBeforeSelection","events before selection",2,0,1);
    h_eventsAfterSelection = new TH1D("h_eventsAfterSelection","events before selection",2,0,1);
    h_NEvents = new TH1D("h_NEvents","Number of events (L3Mu5 trigger)",100,0,10000);
    h_hiBin = new TH1D("h_hiBin","hiBin, inclusive events",200,0,200);
    h_hiBin_triggerOn = new TH1D("h_hiBin_triggerOn","hiBin, events with triggerOn",200,0,200);
    h_hiBin_jet60 = new TH1D("h_hiBin_jet60","hiBin, events with jet60",200,0,200);
    h_hiBin_jet80 = new TH1D("h_hiBin_jet80","hiBin, events with jet80",200,0,200);
    h_hiBin_jet100 = new TH1D("h_hiBin_jet100","hiBin, events with jet100",200,0,200);
    h_hiBin_jet100_withJetAbove130 = new TH1D("h_hiBin_jet100_withJetAbove130","hiBin, events with jet100 & jet with p_{T} > 130 GeV",200,0,200);
    h_hiBin_jet = new TH1D("h_hiBin_jet","hiBin, events with inclRecoJet",200,0,200);
    h_hiBin_inclRecoMuonTag = new TH1D("h_hiBin_inclRecoMuonTag","hiBin, events with inclRecoJet-inclRecoMuonTag",200,0,200);
    h_hiBin_inclRecoMuonTag_triggerOn = new TH1D("h_hiBin_inclRecoMuonTag_triggerOn","hiBin, events with inclRecoJet-inclRecoMuonTag-triggerOn",200,0,200);
    h_inclMuPt = new TH1D("h_inclMuPt","incl. muon p_{T}; muon p_{T}; Entries",NMuPtBins,muPtMin,muPtMax);
    h_nMu_triggerOn = new TH1D("h_nMu_triggerOn","nMu in triggered events",10,0,10);

    h_eventsBeforeSelection->Sumw2();
    h_eventsAfterSelection->Sumw2();
    h_NEvents->Sumw2();
    h_hiBin->Sumw2();
    h_hiBin_triggerOn->Sumw2();
    h_hiBin_jet60->Sumw2();
    h_hiBin_jet80->Sumw2();
    h_hiBin_jet100->Sumw2();
    h_hiBin_jet100_withJetAbove130->Sumw2();
    h_hiBin_jet->Sumw2();
    h_hiBin_inclRecoMuonTag->Sumw2();
    h_hiBin_inclRecoMuonTag_triggerOn->Sumw2();
    h_inclMuPt->Sumw2();
    h_nMu_triggerOn->Sumw2();

    for(int j = 0; j < NJetPtIndices; j++){
      if(j==0){
	h_muptrel_hiBin[j] = new TH2D(Form("h_muptrel_hiBin_J%i",j),Form("muptrel vs hiBin, %3.0f < pTjet < %3.0f",jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax,NhiBinBins,hiBinMin,hiBinMax) ;
      }
      else{
	h_muptrel_hiBin[j] = new TH2D(Form("h_muptrel_hiBin_J%i",j),Form("muptrel vs hiBin, %3.0f < pTjet < %3.0f",jetPtEdges[j-1],jetPtEdges[j]),NMuRelPtBins,muRelPtMin,muRelPtMax,NhiBinBins,hiBinMin,hiBinMax) ;
      }
      h_muptrel_hiBin[j]->Sumw2();
    }
  
    // for loop through the centrality indices
    for(int i = 0; i < NCentralityIndices; i++){
      // the inclusive centrality histogram	
      if(i==0){
	// ---------------------- event histograms --------------------------------
	h_vz[i] = new TH1D(Form("h_vz_C%i",i),Form("vz, inclusive events, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NVzBins,vzMin,vzMax);
	h_vz_triggerOn[i] = new TH1D(Form("h_vz_triggerOn_C%i",i),Form("vz, events with muTrigger, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NVzBins,vzMin,vzMax);
	h_vz_jet60[i] = new TH1D(Form("h_vz_jet60_C%i",i),Form("vz, events with jet60, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NVzBins,vzMin,vzMax);
	h_vz_jet80[i] = new TH1D(Form("h_vz_jet80_C%i",i),Form("vz, events with jet80, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NVzBins,vzMin,vzMax);
	h_vz_jet100[i] = new TH1D(Form("h_vz_jet100_C%i",i),Form("vz, events with jet100, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NVzBins,vzMin,vzMax);
	h_vz_jet100_withJetAbove130[i] = new TH1D(Form("h_vz_jet100_withJetAbove130_C%i",i),Form("vz, events with jet100 & jet with p_{T} > 130 GeV, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NVzBins,vzMin,vzMax);
	h_vz_jet[i] = new TH1D(Form("h_vz_jet_C%i",i),Form("vz, events with inclRecoJet, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NVzBins,vzMin,vzMax);
	h_vz_inclRecoMuonTag[i] = new TH1D(Form("h_vz_inclRecoMuonTag_C%i",i),Form("vz, events with inclRecoJet-inclRecoMuonTag, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NVzBins,vzMin,vzMax);
	h_vz_inclRecoMuonTag_triggerOn[i] = new TH1D(Form("h_vz_inclRecoMuonTag_triggerOn_C%i",i),Form("vz, events with inclRecoJet-inclRecoMuonTag-triggerOn, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NVzBins,vzMin,vzMax);
	h_NJetPerEvent[i] = new TH1D(Form("h_NJetPerEvent_C%i",i),Form("Number of incl. jets per event, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),100,0,100);
	h_NMuTaggedJetPerEvent[i] = new TH1D(Form("h_NMuTaggedJetPerEvent_C%i",i),Form("Number of #it{#mu}-tagged jets per event, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),100,0,100);
	// ----------------------------------------- incl. reco jets --------------   
	h_inclRecoJetPt[i] = new TH1D(Form("h_inclRecoJetPt_C%i",i),Form("incl. reco p_{T}^{jet}, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_inclRawJetPt[i] = new TH1D(Form("h_inclRawJetPt_C%i",i),Form("incl. raw p_{T}^{jet}, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_inclRecoJetEta[i] = new TH1D(Form("h_inclRecoJetEta_C%i",i),Form("incl. reco #eta^{jet}, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NEtaBins,etaMin,etaMax);
	h_inclRecoJetPhi[i] = new TH1D(Form("h_inclRecoJetPhi_C%i",i),Form("incl. reco #phi^{jet}, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPhiBins,phiMin,phiMax);
	h_inclRecoJetPt_inclRecoJetEta[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetEta_C%i",i),Form("incl. reco #eta^{jet} vs. incl reco p_{T}^{jet}, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	h_inclRecoJetPt_inclRecoJetPhi[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetPhi_C%i",i),Form("incl. reco #phi^{jet} vs. incl reco p_{T}^{jet}, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,NPhiBins,phiMin,phiMax);
	// ------------------------------ incl. reco jets + incl. reco muon tag ----
	h_inclRecoJetPt_inclRecoMuonTag[i] = new TH1D(Form("h_inclRecoJetPt_inclRecoMuonTag_C%i",i),Form("incl. reco p_{T}^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_inclRecoJetEta_inclRecoMuonTag[i] = new TH1D(Form("h_inclRecoJetEta_inclRecoMuonTag_C%i",i),Form("incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NEtaBins,etaMin,etaMax);
	h_inclRecoJetPhi_inclRecoMuonTag[i] = new TH1D(Form("h_inclRecoJetPhi_inclRecoMuonTag_C%i",i),Form("incl. reco #phi^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPhiBins,phiMin,phiMax);
	h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_C%i",i),Form("incl. reco #eta^{jet} vs. incl reco p_{T}^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_C%i",i),Form("incl. reco #phi^{jet} vs. incl reco p_{T}^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,NPhiBins,phiMin,phiMax);
	// ---------- events w/ incl. reco jet + incl. reco muon tag + trigger on --
	h_inclRecoJetPt_inclRecoMuonTag_triggerOn[i] = new TH1D(Form("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco p_{T}^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_inclRecoJetEta_inclRecoMuonTag_triggerOn[i] = new TH1D(Form("h_inclRecoJetEta_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NEtaBins,etaMin,etaMax);
	h_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i] = new TH1D(Form("h_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco #phi^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPhiBins,phiMin,phiMax);
	h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco #eta^{jet} vs. incl reco p_{T}^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco #phi^{jet} vs. incl reco p_{T}^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,NPhiBins,phiMin,phiMax);
	// --------------- analysis histograms ------------------------------------
	h_mupt_jetpt[i] = new TH2D(Form("h_mupt_jetpt_C%i",i),Form("jetPt vs. muPt, hiBin %i - %i", centEdges[0], centEdges[NCentralityIndices-1]),NMuPtBins,muPtMin,muPtMax,NPtBins,ptMin,ptMax);
	h_muptrel_jetpt[i] = new TH2D(Form("h_muptrel_jetpt_C%i",i),Form("jetPt vs. muRelPt, hiBin %i - %i", centEdges[0], centEdges[NCentralityIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	// muon-based 2d histograms
	h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC[i] = new TH2D(Form("h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC_C%i",i),Form("fastJet muon #it{p}_{T}^{rel} vs fastJet #it{p}_{T}, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC[i] = new TH2D(Form("h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC_C%i",i),Form("fastJet muon #Delta R vs fastJet #it{p}_{T}, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_fastJetMuonDR_inclusiveClosestFastJet[i] = new TH2D(Form("h_fastJetMuonDR_inclusiveClosestFastJet_C%i",i),Form("fastJet muon #Delta R vs fastJet #it{p}_{T}, inclusive closest fastJet, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_muonDR_inclusiveClosestJet[i] = new TH2D(Form("h_muonDR_inclusiveClosestJet_C%i",i),Form("muon #Delta R vs jet, inclusive closest jet, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_muonDR_inclusiveClosestJet_triggerOn[i] = new TH2D(Form("h_muonDR_inclusiveClosestJet_triggerOn_C%i",i),Form("muon #Delta R vs jet, inclusive closest jet, triggerOn, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_donorJetPt_noInject[i] = new TH1D(Form("h_donorJetPt_noInject_C%i",i),Form("donor fastJet #it{p}_{T}, no injection, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_donorJetPt_inject[i] = new TH1D(Form("h_donorJetPt_inject_C%i",i),Form("donor fastJet #it{p}_{T}, muon injected, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_injMuonDR_donorJetPt_noInject[i] = new TH2D(Form("h_injMuonDR_donorJetPt_noInject_C%i",i),Form("mixed #mu #Delta R vs donor fastJet #it{p}_{T}, no injection, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_injMuonPtRel_donorJetPt_noInject[i] = new TH2D(Form("h_injMuonPtRel_donorJetPt_noInject_C%i",i),Form("mixed #mu #it{p}_{T}^{rel} vs donor fastJet #it{p}_{T}, no injection, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_injMuonDR_donorJetPt_inject[i] = new TH2D(Form("h_injMuonDR_donorJetPt_inject_C%i",i),Form("mixed #mu #Delta R vs donor fastJet #it{p}_{T}, muon injected, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_injMuonPtRel_donorJetPt_inject[i] = new TH2D(Form("h_injMuonPtRel_donorJetPt_inject_C%i",i),Form("mixed #mu #it{p}_{T}^{rel} vs donor fastJet #it{p}_{T}, muon injected, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_injMuonDR_donorJetPt_injectConstit[i] = new TH2D(Form("h_injMuonDR_donorJetPt_injectConstit_C%i",i),Form("mixed #mu #Delta R vs donor fastJet #it{p}_{T}, injected, constituent-tagged, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_injMuonPtRel_donorJetPt_injectConstit[i] = new TH2D(Form("h_injMuonPtRel_donorJetPt_injectConstit_C%i",i),Form("mixed #mu #it{p}_{T}^{rel} vs donor fastJet #it{p}_{T}, injected, constituent-tagged, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_mixedMuonPtRel_recoJetPt[i] = new TH2D(Form("h_mixedMuonPtRel_recoJetPt_C%i",i),Form("mixed-event muon #it{p}_{T}^{rel} vs reco jet #it{p}_{T}, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_realMuonPtRel_mixedFastJetPt[i] = new TH2D(Form("h_realMuonPtRel_mixedFastJetPt_C%i",i),Form("real muon #it{p}_{T}^{rel} vs mixed-event fastJet #it{p}_{T}, %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("muon #it{p}_{T}^{rel} vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_mupt_recoJetPt_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_mupt_recoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("muon #it{p}_{T} vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NMuPtBins,muPtMin,muPtMax,NPtBins,ptMin,ptMax);
	h_mueta_recoJetPt_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_mueta_recoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("muon #it{#eta} vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NTrkEtaBins,trkEtaMin,trkEtaMax,NPtBins,ptMin,ptMax);
	h_muphi_recoJetPt_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_muphi_recoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("muon #it{#phi} vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NPhiBins,phiMin,phiMax,NPtBins,ptMin,ptMax);
	h_muJetDr_recoJetPt[i] = new TH2D(Form("h_muJetDr_recoJetPt_C%i",i),Form("#it{#Delta r}(muon,jet) vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_dimuonMass[i] = new TH1D(Form("h_dimuonMass_C%i",i),Form("dimuon mass, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NDimuonMassBins,dimuonMassMin,dimuonMassMax);
	h_dimuonMass_sameSign[i] = new TH1D(Form("h_dimuonMass_sameSign_C%i",i),Form("dimuon mass, %i < hiBin < %i",centEdges[0], centEdges[NCentralityIndices-1]),NDimuonMassBins,dimuonMassMin,dimuonMassMax);
	h_Jet60HLT[i] = new TH1D(Form("h_Jet60HLT_C%i",i),Form("Jet60 HLT , %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),2,0,1);
	h_Jet60HLT_Prescale[i] = new TH1D(Form("h_Jet60HLT_Prescale_C%i",i),Form("Jet60 HLT Prescale , %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),10,1,10);
	h_Jet80HLT[i] = new TH1D(Form("h_Jet80HLT_C%i",i),Form("Jet80 HLT , %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),2,0,1);
	h_Jet80HLT_Prescale[i] = new TH1D(Form("h_Jet80HLT_Prescale_C%i",i),Form("Jet80 HLT Prescale , %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),10,1,10);
	h_Jet100HLT[i] = new TH1D(Form("h_Jet100HLT_C%i",i),Form("Jet100 HLT , %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),2,0,1);
	h_Jet100HLT_Prescale[i] = new TH1D(Form("h_Jet100HLT_Prescale_C%i",i),Form("Jet100 HLT Prescale , %i < hiBin < %i",centEdges[0],centEdges[NCentralityIndices-1]),10,1,10);
	// pfCand / FastJet histograms
	h_pfPt[i] = new TH1D(Form("h_pfPt_C%i",i),Form("pfCand pT, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),100,0,100);
	h_pseudoJetPt[i] = new TH1D(Form("h_pseudoJetPt_C%i",i),Form("PseudoJet pT, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_pseudoJetPt_geoCorr[i] = new TH1D(Form("h_pseudoJetPt_geoCorr_C%i",i),Form("PseudoJet pT, geometric correction, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF[i] = new TH1D(Form("h_fastJetPt_PF_C%i",i),Form("FastJet (PF) anti-kT pT, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC[i] = new TH1D(Form("h_fastJetPt_PF_JEC_C%i",i),Form("FastJet anti-kT pT (PF, JEC), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PFCs[i] = new TH1D(Form("h_fastJetPt_PFCs_C%i",i),Form("FastJet (PFCs) anti-kT pT, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PFCs_JEC[i] = new TH1D(Form("h_fastJetPt_PFCs_JEC_C%i",i),Form("FastJet anti-kT pT (PFCs, JEC), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_RC[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_RC_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, random-cone), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_RC[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_RC_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, random-cone), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_RC_sigSel[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_RC_sigSel_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, random-cone, signal-selected same-event), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_RC_sigSel[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_RC_sigSel_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, random-cone, signal-selected same-event), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_RC_geoCorr[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_RC_geoCorr_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, random-cone, geometry correction), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, random-cone, geometry correction), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, random-cone, geometry correction, #eta-reflect), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, random-cone, geometry correction, #eta-reflect), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_dPT[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_dPT_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, dPT), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_dPT[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_dPT_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, dPT), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_dPT_geoCorr[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_dPT_geoCorr_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, dPT, geometry correction), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_dPT_geoCorr[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_dPT_geoCorr_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, dPT, geometry correction), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_dPT_dPTAbove0[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_dPT_dPTAbove0_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, dPT), dPT > 0, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_dPT_dPTAbove0[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_dPT_dPTAbove0_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, dPT), dPT > 0, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, dPT), pT(PFCs) > 60, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_dPT_PFCsPTAbove60[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_dPT_PFCsPTAbove60_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, dPT), pT(PFCs) > 60, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_nPFcand[i]   = new TH1D(Form("h_nPFcand_C%i",i),  Form("N PF cands per event, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),10000,0,10000);
	h_nPFcandCS[i] = new TH1D(Form("h_nPFcandCS_C%i",i), Form("N PFCS cands per event, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),10000,0,10000);
	h_nPFcandFastJet[i] = new TH1D(Form("h_nPFcandFastJet_C%i",i), Form("N PF cands clustered by FastJet, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),10000,0,10000);
	h_randConeEtaPhi[i] = new TProfile2D(Form("h_randConeEtaPhi_C%i",i),Form("Mean random-cone p_{T} vs (#eta,#phi), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_randConeEtaPhi_geoCorr[i] = new TProfile2D(Form("h_randConeEtaPhi_geoCorr_C%i",i),Form("Mean random-cone p_{T} (w/ geometric correction) vs (#eta,#phi), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_pseudoJetPt_sigSel[i] = new TH1D(Form("h_pseudoJetPt_sigSel_C%i",i),Form("PseudoJet pT, signal-selected same-event, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_randConeEtaPhi_sigSel[i] = new TProfile2D(Form("h_randConeEtaPhi_sigSel_C%i",i),Form("Mean random-cone p_{T} vs (#eta,#phi), signal-selected same-event, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dPTEtaPhi_PF_PFCs[i] = new TProfile2D(Form("h_dPTEtaPhi_PF_PFCs_C%i",i),Form("p_{T}^{fastJet(PF)} - p_{T}^{fastJet(PFCs)} vs (#eta,#phi), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dPTEtaPhi_PF_PFCs_geoCorr[i] = new TProfile2D(Form("h_dPTEtaPhi_PF_PFCs_geoCorr_C%i",i),Form("(p_{T}^{fastJet(PF)} - p_{T}^{fastJet(PFCs)}) #times cos(#Delta R) vs (#eta,#phi), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dPTEtaPhi_PF_PFCs_dPTAbove0[i] = new TProfile2D(Form("h_dPTEtaPhi_PF_PFCs_dPTAbove0_C%i",i),Form("p_{T}^{fastJet(PF)} - p_{T}^{fastJet(PFCs)} vs (#eta,#phi), #Delta p_{T} > 0, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60[i] = new TProfile2D(Form("h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60_C%i",i),Form("p_{T}^{fastJet(PF)} - p_{T}^{fastJet(PFCs)} (p_{T}^{PFCs} > 60 GeV)  vs (#eta,#phi), #Delta p_{T} > 0, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dRmin_PF_PFCs[i] = new TH1D(Form("h_dRmin_PF_PFCs_C%i",i), Form("min #Delta r[fastJet(PF),fastJet(PFCs)], hiBin %i-%i",centEdges[0],centEdges[NCentralityIndices-1]),NdRBins,dRBinMin,dRBinMax);
						
      }
      else{
	// ---------------------- event histograms --------------------------------
	h_vz[i] = new TH1D(Form("h_vz_C%i",i),Form("vz, inclusive events, hiBin %i - %i",centEdges[i-1],centEdges[i]),NVzBins,vzMin,vzMax);
	h_vz_triggerOn[i] = new TH1D(Form("h_vz_triggerOn_C%i",i),Form("vz, events with triggerOn, hiBin %i - %i",centEdges[i-1],centEdges[i]),NVzBins,vzMin,vzMax);
	h_vz_jet60[i] = new TH1D(Form("h_vz_jet60_C%i",i),Form("vz, events with jet60, hiBin %i - %i",centEdges[i-1],centEdges[i]),NVzBins,vzMin,vzMax);
	h_vz_jet80[i] = new TH1D(Form("h_vz_jet80_C%i",i),Form("vz, events with jet80, hiBin %i - %i",centEdges[i-1],centEdges[i]),NVzBins,vzMin,vzMax);
	h_vz_jet100[i] = new TH1D(Form("h_vz_jet100_C%i",i),Form("vz, events with jet100, hiBin %i - %i",centEdges[i-1],centEdges[i]),NVzBins,vzMin,vzMax);
	h_vz_jet100_withJetAbove130[i] = new TH1D(Form("h_vz_jet100_withJetAbove130_C%i",i),Form("vz, events with jet100 & jet with p_{T} > 130 GeV, hiBin %i - %i",centEdges[i-1],centEdges[i]),NVzBins,vzMin,vzMax);
	h_vz_jet[i] = new TH1D(Form("h_vz_jet_C%i",i),Form("vz, events with inclRecoJet, hiBin %i - %i",centEdges[i-1],centEdges[i]),NVzBins,vzMin,vzMax);
	h_vz_inclRecoMuonTag[i] = new TH1D(Form("h_vz_inclRecoMuonTag_C%i",i),Form("vz, events with inclRecoJet-inclRecoMuonTag, hiBin %i - %i",centEdges[i-1],centEdges[i]),NVzBins,vzMin,vzMax);
	h_vz_inclRecoMuonTag_triggerOn[i] = new TH1D(Form("h_vz_inclRecoMuonTag_triggerOn_C%i",i),Form("vz, events with inclRecoJet-inclRecoMuonTag-triggerOn, hiBin %i - %i",centEdges[i-1],centEdges[i]),NVzBins,vzMin,vzMax);
	h_NJetPerEvent[i] = new TH1D(Form("h_NJetPerEvent_C%i",i),Form("Number of incl. jets per event, hiBin %i - %i",centEdges[i-1],centEdges[i]),100,0,100);
	h_NMuTaggedJetPerEvent[i] = new TH1D(Form("h_NMuTaggedJetPerEvent_C%i",i),Form("Number of #it{#mu}-tagged jets per event, hiBin %i - %i",centEdges[i-1],centEdges[i]),100,0,100);
	// ----------------------------------------- incl. reco jets --------------   
	h_inclRecoJetPt[i] = new TH1D(Form("h_inclRecoJetPt_C%i",i),Form("incl. reco p_{T}^{jet}, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_inclRawJetPt[i] = new TH1D(Form("h_inclRawJetPt_C%i",i),Form("incl. raw p_{T}^{jet}, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_inclRecoJetEta[i] = new TH1D(Form("h_inclRecoJetEta_C%i",i),Form("incl. reco #eta^{jet}, hiBin %i - %i",centEdges[i-1],centEdges[i]),NEtaBins,etaMin,etaMax);
	h_inclRecoJetPhi[i] = new TH1D(Form("h_inclRecoJetPhi_C%i",i),Form("incl. reco #phi^{jet}, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPhiBins,phiMin,phiMax);
	h_inclRecoJetPt_inclRecoJetEta[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetEta_C%i",i),Form("incl. reco #eta^{jet} vs. incl reco p_{T}^{jet}, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	h_inclRecoJetPt_inclRecoJetPhi[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetPhi_C%i",i),Form("incl. reco #phi^{jet} vs. incl reco p_{T}^{jet}, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,NPhiBins,phiMin,phiMax);
	// ------------------------------ incl. reco jets + incl. reco muon tag ----
	h_inclRecoJetPt_inclRecoMuonTag[i] = new TH1D(Form("h_inclRecoJetPt_inclRecoMuonTag_C%i",i),Form("incl. reco p_{T}^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_inclRecoJetEta_inclRecoMuonTag[i] = new TH1D(Form("h_inclRecoJetEta_inclRecoMuonTag_C%i",i),Form("incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[i-1],centEdges[i]),NEtaBins,etaMin,etaMax);
	h_inclRecoJetPhi_inclRecoMuonTag[i] = new TH1D(Form("h_inclRecoJetPhi_inclRecoMuonTag_C%i",i),Form("incl. reco #phi^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPhiBins,phiMin,phiMax);
	h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_C%i",i),Form("incl. reco #eta^{jet} vs. incl reco p_{T}^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_C%i",i),Form("incl. reco #phi^{jet} vs. incl reco p_{T}^{jet}, tagged with incl. reco muon, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,NPhiBins,phiMin,phiMax);
	// ---------- events w/ incl. reco jet + incl. reco muon tag + trigger on --
	h_inclRecoJetPt_inclRecoMuonTag_triggerOn[i] = new TH1D(Form("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco p_{T}^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_inclRecoJetEta_inclRecoMuonTag_triggerOn[i] = new TH1D(Form("h_inclRecoJetEta_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[i-1],centEdges[i]),NEtaBins,etaMin,etaMax);
	h_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i] = new TH1D(Form("h_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco #phi^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPhiBins,phiMin,phiMax);
	h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco #eta^{jet} vs. incl reco p_{T}^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%i",i),Form("incl. reco #phi^{jet} vs. incl reco p_{T}^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,NPhiBins,phiMin,phiMax);
	// ----------- analysis histograms -----------------------------------------
	h_mupt_jetpt[i] = new TH2D(Form("h_mupt_jetpt_C%i",i),Form("jetPt vs. muPt, hiBin %i - %i", centEdges[i-1], centEdges[i]),NMuPtBins,muPtMin,muPtMax,NPtBins,ptMin,ptMax);
	h_muptrel_jetpt[i] = new TH2D(Form("h_muptrel_jetpt_C%i",i),Form("jetPt vs. muRelPt, hiBin %i - %i", centEdges[i-1], centEdges[i]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	// muon-based 2d histograms
	h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC[i] = new TH2D(Form("h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC_C%i",i),Form("fastJet muon #it{p}_{T}^{rel} vs fastJet #it{p}_{T}, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC[i] = new TH2D(Form("h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC_C%i",i),Form("fastJet muon #Delta R vs fastJet #it{p}_{T}, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_fastJetMuonDR_inclusiveClosestFastJet[i] = new TH2D(Form("h_fastJetMuonDR_inclusiveClosestFastJet_C%i",i),Form("fastJet muon #Delta R vs fastJet #it{p}_{T}, inclusive closest fastJet, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_muonDR_inclusiveClosestJet[i] = new TH2D(Form("h_muonDR_inclusiveClosestJet_C%i",i),Form("muon #Delta R vs jet, inclusive closest jet, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_muonDR_inclusiveClosestJet_triggerOn[i] = new TH2D(Form("h_muonDR_inclusiveClosestJet_triggerOn_C%i",i),Form("muon #Delta R vs jet, inclusive closest jet, triggerOn, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_donorJetPt_noInject[i] = new TH1D(Form("h_donorJetPt_noInject_C%i",i),Form("donor fastJet #it{p}_{T}, no injection, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_donorJetPt_inject[i] = new TH1D(Form("h_donorJetPt_inject_C%i",i),Form("donor fastJet #it{p}_{T}, muon injected, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_injMuonDR_donorJetPt_noInject[i] = new TH2D(Form("h_injMuonDR_donorJetPt_noInject_C%i",i),Form("mixed #mu #Delta R vs donor fastJet #it{p}_{T}, no injection, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_injMuonPtRel_donorJetPt_noInject[i] = new TH2D(Form("h_injMuonPtRel_donorJetPt_noInject_C%i",i),Form("mixed #mu #it{p}_{T}^{rel} vs donor fastJet #it{p}_{T}, no injection, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_injMuonDR_donorJetPt_inject[i] = new TH2D(Form("h_injMuonDR_donorJetPt_inject_C%i",i),Form("mixed #mu #Delta R vs donor fastJet #it{p}_{T}, muon injected, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_injMuonPtRel_donorJetPt_inject[i] = new TH2D(Form("h_injMuonPtRel_donorJetPt_inject_C%i",i),Form("mixed #mu #it{p}_{T}^{rel} vs donor fastJet #it{p}_{T}, muon injected, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_injMuonDR_donorJetPt_injectConstit[i] = new TH2D(Form("h_injMuonDR_donorJetPt_injectConstit_C%i",i),Form("mixed #mu #Delta R vs donor fastJet #it{p}_{T}, injected, constituent-tagged, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_injMuonPtRel_donorJetPt_injectConstit[i] = new TH2D(Form("h_injMuonPtRel_donorJetPt_injectConstit_C%i",i),Form("mixed #mu #it{p}_{T}^{rel} vs donor fastJet #it{p}_{T}, injected, constituent-tagged, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_mixedMuonPtRel_recoJetPt[i] = new TH2D(Form("h_mixedMuonPtRel_recoJetPt_C%i",i),Form("mixed-event muon #it{p}_{T}^{rel} vs reco jet #it{p}_{T}, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_realMuonPtRel_mixedFastJetPt[i] = new TH2D(Form("h_realMuonPtRel_mixedFastJetPt_C%i",i),Form("real muon #it{p}_{T}^{rel} vs mixed-event fastJet #it{p}_{T}, %i < hiBin < %i",centEdges[i-1],centEdges[i]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("muon #it{p}_{T}^{rel} vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NMuRelPtBins,muRelPtMin,muRelPtMax,NPtBins,ptMin,ptMax);
	h_mupt_recoJetPt_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_mupt_recoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("muon #it{p}_{T} vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NMuPtBins,muPtMin,muPtMax,NPtBins,ptMin,ptMax);
	h_mueta_recoJetPt_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_mueta_recoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("muon #it{#eta} vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NTrkEtaBins,trkEtaMin,trkEtaMax,NPtBins,ptMin,ptMax);
	h_muphi_recoJetPt_inclRecoMuonTag_triggerOn[i] = new TH2D(Form("h_muphi_recoJetPt_inclRecoMuonTag_triggerOn_C%i",i),Form("muon #it{#phi} vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NPhiBins,phiMin,phiMax,NPtBins,ptMin,ptMax);
	h_muJetDr_recoJetPt[i] = new TH2D(Form("h_muJetDr_recoJetPt_C%i",i),Form("#it{#Delta r}(muon,jet) vs jet #it{p}_{T}, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NdRBins,dRBinMin,dRBinMax,NPtBins,ptMin,ptMax);
	h_dimuonMass[i] = new TH1D(Form("h_dimuonMass_C%i",i),Form("dimuon mass, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NDimuonMassBins,dimuonMassMin,dimuonMassMax);
	h_dimuonMass_sameSign[i] = new TH1D(Form("h_dimuonMass_sameSign_C%i",i),Form("dimuon mass, %i < hiBin < %i",centEdges[i-1], centEdges[i]),NDimuonMassBins,dimuonMassMin,dimuonMassMax);
	h_Jet60HLT[i] = new TH1D(Form("h_Jet60HLT_C%i",i),Form("Jet60 HLT , %i < hiBin < %i",centEdges[i-1], centEdges[i]),2,0,1);
	h_Jet60HLT_Prescale[i] = new TH1D(Form("h_Jet60HLT_Prescale_C%i",i),Form("Jet60 HLT Prescale , %i < hiBin < %i",centEdges[i-1], centEdges[i]),10,1,10);
	h_Jet80HLT[i] = new TH1D(Form("h_Jet80HLT_C%i",i),Form("Jet80 HLT , %i < hiBin < %i",centEdges[i-1], centEdges[i]),2,0,1);
	h_Jet80HLT_Prescale[i] = new TH1D(Form("h_Jet80HLT_Prescale_C%i",i),Form("Jet80 HLT Prescale , %i < hiBin < %i",centEdges[i-1], centEdges[i]),10,1,10);
	h_Jet100HLT[i] = new TH1D(Form("h_Jet100HLT_C%i",i),Form("Jet100 HLT , %i < hiBin < %i",centEdges[i-1], centEdges[i]),2,0,1);
	h_Jet100HLT_Prescale[i] = new TH1D(Form("h_Jet100HLT_Prescale_C%i",i),Form("Jet100 HLT Prescale , %i < hiBin < %i",centEdges[i-1], centEdges[i]),10,1,10);
	// pfCand / FastJet histograms
	h_pfPt[i] = new TH1D(Form("h_pfPt_C%i",i),Form("pfCand pT, hiBin %i - %i",centEdges[i-1],centEdges[i]),100,0,100);
	h_pseudoJetPt[i] = new TH1D(Form("h_pseudoJetPt_C%i",i),Form("PseudoJet pT, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_pseudoJetPt_geoCorr[i] = new TH1D(Form("h_pseudoJetPt_geoCorr_C%i",i),Form("PseudoJet pT, geometry correction, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF[i] = new TH1D(Form("h_fastJetPt_PF_C%i",i),Form("FastJet (PF) anti-kT pT, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC[i] = new TH1D(Form("h_fastJetPt_PF_JEC_C%i",i),Form("FastJet (PF) anti-kT pT (JEC), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PFCs[i] = new TH1D(Form("h_fastJetPt_PFCs_C%i",i),Form("FastJet (PFCs) anti-kT pT, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PFCs_JEC[i] = new TH1D(Form("h_fastJetPt_PFCs_JEC_C%i",i),Form("FastJet (PFCs) anti-kT pT (JEC), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);


	h_fastJetPt_PF_bkgSub_RC_sigSel[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_RC_sigSel_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, random-cone, signal-selected same-event), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_RC_sigSel[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_RC_sigSel_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, random-cone, signal-selected same-event), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_RC[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_RC_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, random-cone), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_RC[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_RC_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, random-cone), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_RC_geoCorr[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_RC_geoCorr_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, random-cone, geometry correction), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, random-cone, geometry correction), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, random-cone, geometry correction, #eta-reflect), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, random-cone, geometry correction, #eta-reflect), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_dPT[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_dPT_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, dPT), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_dPT[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_dPT_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, dPT), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_dPT_geoCorr[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_dPT_geoCorr_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, dPT, geometry correction), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_dPT_geoCorr[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_dPT_geoCorr_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, dPT, geometry correction), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_dPT_dPTAbove0[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_dPT_dPTAbove0_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, dPT), dPT > 0, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_dPT_dPTAbove0[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_dPT_dPTAbove0_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, dPT), dPT > 0, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60[i] = new TH1D(Form("h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60_C%i",i),Form("FastJet anti-kT pT (PF, bkg sub, dPT), pT(PFCs) > 60, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_PF_JEC_bkgSub_dPT_PFCsPTAbove60[i] = new TH1D(Form("h_fastJetPt_PF_JEC_bkgSub_dPT_PFCsPTAbove60_C%i",i),Form("FastJet anti-kT pT (PF, JEC, bkg sub, dPT), pT(PFCs) > 60, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);




	h_nPFcand[i]   = new TH1D(Form("h_nPFcand_C%i",i),  Form("N PF cands per event, hiBin %i - %i",centEdges[i-1],centEdges[i]),10000,0,10000);
	h_nPFcandCS[i] = new TH1D(Form("h_nPFcandCS_C%i",i), Form("N PFCS cands per event, hiBin %i - %i",centEdges[i-1],centEdges[i]),10000,0,10000);
	h_nPFcandFastJet[i] = new TH1D(Form("h_nPFcandFastJet_C%i",i), Form("N PF cands clustered by FastJet, hiBin %i - %i",centEdges[i-1],centEdges[i]),10000,0,10000);
	h_randConeEtaPhi[i] = new TProfile2D(Form("h_randConeEtaPhi_C%i",i),Form("Mean random-cone p_{T} vs (#eta,#phi), hiBin %i - %i",centEdges[i-1],centEdges[i]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_randConeEtaPhi_geoCorr[i] = new TProfile2D(Form("h_randConeEtaPhi_geoCorr_C%i",i),Form("Mean random-cone p_{T} (w/ geometric correction) vs (#eta,#phi), hiBin %i - %i",centEdges[i-1],centEdges[i]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_pseudoJetPt_sigSel[i] = new TH1D(Form("h_pseudoJetPt_sigSel_C%i",i),Form("PseudoJet pT, signal-selected same-event, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_randConeEtaPhi_sigSel[i] = new TProfile2D(Form("h_randConeEtaPhi_sigSel_C%i",i),Form("Mean random-cone p_{T} vs (#eta,#phi), signal-selected same-event, hiBin %i - %i",centEdges[i-1],centEdges[i]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dPTEtaPhi_PF_PFCs[i] = new TProfile2D(Form("h_dPTEtaPhi_PF_PFCs_C%i",i),Form("p_{T}^{fastJet(PF)} - p_{T}^{fastJet(PFCs)} vs (#eta,#phi), hiBin %i - %i",centEdges[i-1],centEdges[i]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dPTEtaPhi_PF_PFCs_geoCorr[i] = new TProfile2D(Form("h_dPTEtaPhi_PF_PFCs_geoCorr_C%i",i),Form("(p_{T}^{fastJet(PF)} - p_{T}^{fastJet(PFCs)}) #times cos(#Delta R) vs (#eta,#phi), hiBin %i - %i",centEdges[i-1],centEdges[i]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dPTEtaPhi_PF_PFCs_dPTAbove0[i] = new TProfile2D(Form("h_dPTEtaPhi_PF_PFCs_dPTAbove0_C%i",i),Form("p_{T}^{fastJet(PF)} - p_{T}^{fastJet(PFCs)} vs (#eta,#phi), #Delta p_{T} > 0, hiBin %i - %i",centEdges[i-1],centEdges[i]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60[i] = new TProfile2D(Form("h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60_C%i",i),Form("p_{T}^{fastJet(PF)} - p_{T}^{fastJet(PFCs)} (p_{T}^{PFCs} > 60 GeV)  vs (#eta,#phi), #Delta p_{T} > 0, hiBin %i - %i",centEdges[i-1],centEdges[i]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
	h_dRmin_PF_PFCs[i] = new TH1D(Form("h_dRmin_PF_PFCs_C%i",i), Form("min #Delta r[fastJet(PF),fastJet(PFCs)], hiBin %i-%i",centEdges[i-1],centEdges[i]),NdRBins,dRBinMin,dRBinMax);
      }
      // sumw2 commands
      h_NJetPerEvent[i]->Sumw2();
      h_NMuTaggedJetPerEvent[i]->Sumw2();
      h_vz[i]->Sumw2();
      h_vz_triggerOn[i]->Sumw2();
      h_vz_jet60[i]->Sumw2();
      h_vz_jet80[i]->Sumw2();
      h_vz_jet100[i]->Sumw2();
      h_vz_jet100_withJetAbove130[i]->Sumw2();
      h_vz_jet[i]->Sumw2();
      h_vz_inclRecoMuonTag[i]->Sumw2();
      h_vz_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_inclRecoJetPt[i]->Sumw2();
      h_inclRawJetPt[i]->Sumw2();
      h_inclRecoJetEta[i]->Sumw2();
      h_inclRecoJetPhi[i]->Sumw2();
      h_inclRecoJetPt_inclRecoJetEta[i]->Sumw2();
      h_inclRecoJetPt_inclRecoJetPhi[i]->Sumw2();
      h_inclRecoJetPt_inclRecoMuonTag[i]->Sumw2();
      h_inclRecoJetEta_inclRecoMuonTag[i]->Sumw2();
      h_inclRecoJetPhi_inclRecoMuonTag[i]->Sumw2();
      h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag[i]->Sumw2();
      h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag[i]->Sumw2();
      h_inclRecoJetPt_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_inclRecoJetEta_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_mupt_jetpt[i]->Sumw2();
      h_muptrel_jetpt[i]->Sumw2();
      h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_mupt_recoJetPt_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_mueta_recoJetPt_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_muphi_recoJetPt_inclRecoMuonTag_triggerOn[i]->Sumw2();
      h_muJetDr_recoJetPt[i]->Sumw2();
      h_dimuonMass[i]->Sumw2();
      h_dimuonMass_sameSign[i]->Sumw2();
      h_Jet60HLT[i]->Sumw2();
      h_Jet60HLT_Prescale[i]->Sumw2();
      h_Jet80HLT[i]->Sumw2();
      h_Jet80HLT_Prescale[i]->Sumw2();
      h_Jet100HLT[i]->Sumw2();
      h_Jet100HLT_Prescale[i]->Sumw2();
      h_pfPt[i]->Sumw2();
      h_pseudoJetPt[i]->Sumw2();
      h_pseudoJetPt_geoCorr[i]->Sumw2();
      h_fastJetPt_PF[i]->Sumw2();
      h_fastJetPt_PF_JEC[i]->Sumw2();
      h_fastJetPt_PFCs[i]->Sumw2();
      h_fastJetPt_PFCs_JEC[i]->Sumw2();

      h_fastJetPt_PF_bkgSub_RC[i]->Sumw2();
      h_fastJetPt_PF_JEC_bkgSub_RC[i]->Sumw2();
      h_fastJetPt_PF_bkgSub_RC_sigSel[i]->Sumw2();
      h_fastJetPt_PF_JEC_bkgSub_RC_sigSel[i]->Sumw2();
      h_fastJetPt_PF_bkgSub_RC_geoCorr[i]->Sumw2();
      h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr[i]->Sumw2();
      h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect[i]->Sumw2();
      h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect[i]->Sumw2();
      h_fastJetPt_PF_bkgSub_dPT[i]->Sumw2();
      h_fastJetPt_PF_JEC_bkgSub_dPT[i]->Sumw2();
      h_fastJetPt_PF_bkgSub_dPT_geoCorr[i]->Sumw2();
      h_fastJetPt_PF_JEC_bkgSub_dPT_geoCorr[i]->Sumw2();
      h_fastJetPt_PF_bkgSub_dPT_dPTAbove0[i]->Sumw2();
      h_fastJetPt_PF_JEC_bkgSub_dPT_dPTAbove0[i]->Sumw2();
      h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60[i]->Sumw2();
      h_fastJetPt_PF_JEC_bkgSub_dPT_PFCsPTAbove60[i]->Sumw2();

      h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC[i]->Sumw2();
      h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC[i]->Sumw2();
      h_fastJetMuonDR_inclusiveClosestFastJet[i]->Sumw2();
      h_muonDR_inclusiveClosestJet[i]->Sumw2();
      h_muonDR_inclusiveClosestJet_triggerOn[i]->Sumw2();
      h_donorJetPt_noInject[i]->Sumw2();
      h_donorJetPt_inject[i]->Sumw2();
      h_injMuonDR_donorJetPt_noInject[i]->Sumw2();
      h_injMuonPtRel_donorJetPt_noInject[i]->Sumw2();
      h_injMuonDR_donorJetPt_inject[i]->Sumw2();
      h_injMuonPtRel_donorJetPt_inject[i]->Sumw2();
      h_injMuonDR_donorJetPt_injectConstit[i]->Sumw2();
      h_injMuonPtRel_donorJetPt_injectConstit[i]->Sumw2();
      h_mixedMuonPtRel_recoJetPt[i]->Sumw2();
      h_realMuonPtRel_mixedFastJetPt[i]->Sumw2();

      h_nPFcand[i]->Sumw2();
      h_nPFcandCS[i]->Sumw2();
      h_nPFcandFastJet[i]->Sumw2();
      h_randConeEtaPhi[i]->Sumw2();
      h_randConeEtaPhi_geoCorr[i]->Sumw2();
      h_pseudoJetPt_sigSel[i]->Sumw2();
      h_randConeEtaPhi_sigSel[i]->Sumw2();
      h_dPTEtaPhi_PF_PFCs[i]->Sumw2();
      h_dPTEtaPhi_PF_PFCs_geoCorr[i]->Sumw2();
      h_dPTEtaPhi_PF_PFCs_dPTAbove0[i]->Sumw2();
      h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60[i]->Sumw2();
      h_dRmin_PF_PFCs[i]->Sumw2();

      // loop through jet pt indices
      for(int j = 0; j < NJetPtIndices; j++){
	if(i==0 && j==0){
	  h_fastJetPtPF_dRmin[i][j] = new TH2D(Form("h_fastJetPtPF_dRmin_C%iJ%i",i,j),Form("dR(PF,PFCs) vs fastJetPt(PF), hiBin %i - %i, p_{T}^{PFCs} 20 - 500",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,NdRBins,dRBinMin,dRBinMax);
	  h_fastJetPtPF_etaPFCs[i][j] = new TH2D(Form("h_fastJetPtPF_etaPFCs_C%iJ%i",i,j),Form("#eta(PFCs) vs fastJetPt(PF), hiBin %i - %i, p_{T}^{PFCs} 20-500",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	  h_muPtOverJetPt[i][j] = new TH1D(Form("h_muPtOverJetPt_C%iJ%i",i,j),Form("p_{T}^{#mu} / p_{T}^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),100,0,1);
	  h_mupt_muptrel[i][j] = new TH2D(Form("h_mupt_muptrel_C%iJ%i",i,j),Form("muPtRel vs. muPt, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f",centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NMuPtBins,muPtMin,muPtMax,NMuRelPtBins,muRelPtMin,muRelPtMax);
	  h_inclRecoJetEta_inclRecoJetPhi[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	}
	else if(i==0){
	  h_fastJetPtPF_dRmin[i][j] = new TH2D(Form("h_fastJetPtPF_dRmin_C%iJ%i",i,j),Form("dR(PF,PFCs) vs fastJetPt(PF), hiBin %i - %i, p_{T}^{PFCs} %3.0f - %3.0f",centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NPtBins,ptMin,ptMax,NdRBins,dRBinMin,dRBinMax);
	  h_fastJetPtPF_etaPFCs[i][j] = new TH2D(Form("h_fastJetPtPF_etaPFCs_C%iJ%i",i,j),Form("#eta(PFCs) vs fastJetPt(PF), hiBin %i - %i, p_{T}^{PFCs} %3.0f - %3.0f",centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	  h_muPtOverJetPt[i][j] = new TH1D(Form("h_muPtOverJetPt_C%iJ%i",i,j),Form("p_{T}^{#mu} / p_{T}^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),100,0,1);
	  h_mupt_muptrel[i][j] = new TH2D(Form("h_mupt_muptrel_C%iJ%i",i,j),Form("muPtRel vs. muPt, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f",centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NMuPtBins,muPtMin,muPtMax,NMuRelPtBins,muRelPtMin,muRelPtMax);
	  h_inclRecoJetEta_inclRecoJetPhi[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	}
	else if(j==0){
	  h_fastJetPtPF_dRmin[i][j] = new TH2D(Form("h_fastJetPtPF_dRmin_C%iJ%i",i,j),Form("dR(PF,PFCs) vs fastJetPt(PF), hiBin %i - %i, p_{T}^{PFCs} 20 - 500",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,NdRBins,dRBinMin,dRBinMax);
	  h_fastJetPtPF_etaPFCs[i][j] = new TH2D(Form("h_fastJetPtPF_etaPFCs_C%iJ%i",i,j),Form("#eta(PFCs) vs fastJetPt(PF), hiBin %i - %i, p_{T}^{PFCs} 20 - 500",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	  h_muPtOverJetPt[i][j] = new TH1D(Form("h_muPtOverJetPt_C%iJ%i",i,j),Form("p_{T}^{#mu} / p_{T}^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),100,0,1);
	  h_mupt_muptrel[i][j] = new TH2D(Form("h_mupt_muptrel_C%iJ%i",i,j),Form("muPtRel vs. muPt, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f",centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NMuPtBins,muPtMin,muPtMax,NMuRelPtBins,muRelPtMin,muRelPtMax);
	  h_inclRecoJetEta_inclRecoJetPhi[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	}
	else{
	  h_fastJetPtPF_dRmin[i][j] = new TH2D(Form("h_fastJetPtPF_dRmin_C%iJ%i",i,j),Form("dR(PF,PFCs) vs fastJetPt(PF), hiBin %i - %i, p_{T}^{PFCs} %3.0f - %3.0f",centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NPtBins,ptMin,ptMax,NdRBins,dRBinMin,dRBinMax);
	  h_fastJetPtPF_etaPFCs[i][j] = new TH2D(Form("h_fastJetPtPF_etaPFCs_C%iJ%i",i,j),Form("#eta(PFCs) vs fastJetPt(PF), hiBin %i - %i, p_{T}^{PFCs} %3.0f - %3.0f",centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NPtBins,ptMin,ptMax,NEtaBins,etaMin,etaMax);
	  h_muPtOverJetPt[i][j] = new TH1D(Form("h_muPtOverJetPt_C%iJ%i",i,j),Form("p_{T}^{#mu} / p_{T}^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),100,0,1);
	  h_mupt_muptrel[i][j] = new TH2D(Form("h_mupt_muptrel_C%iJ%i",i,j),Form("muPtRel vs. muPt, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f",centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NMuPtBins,muPtMin,muPtMax,NMuRelPtBins,muRelPtMin,muRelPtMax);
	  h_inclRecoJetEta_inclRecoJetPhi[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	}

	h_fastJetPtPF_dRmin[i][j]->Sumw2();
	h_fastJetPtPF_etaPFCs[i][j]->Sumw2();
	h_muPtOverJetPt[i][j]->Sumw2();
	h_mupt_muptrel[i][j]->Sumw2();
	h_inclRecoJetEta_inclRecoJetPhi[i][j]->Sumw2();
	h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j]->Sumw2();
	h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j]->Sumw2();
      }

    }

  
    TFile *f = TFile::Open(input.c_str());
    cout << "	File opened!" << endl;
    auto em = new eventMap(f);
    em->isMC = isMC_status;
    em->AASetup = AASetup_status;
    cout << "	Initializing variables ... " << endl;
    em->init();
    cout << "	Loading jet..." << endl;
    if(useCaloJetsOverride) em->loadJet("akPu4CaloJetAnalyzer/t");
    else if(useFlowJetsOverride) em->loadJet("akFlowPuCs4PFJetAnalyzer/t");
    else em->loadJet("akCs4PFJetAnalyzer/t");
    cout << "	Loading muon..." << endl;
    em->loadMuon("ggHiNtuplizerGED/EventTree");
    cout << "	Loading muon triggers..." << endl;
    em->loadHLT("hltanalysis/HltTree");
    cout << "	Loading PF candidates..." << endl;
    em->loadParticleFlowAnalyzer_PFCs("pfcandAnalyzerCS");
    em->loadParticleFlowAnalyzer("pfcandAnalyzer");
    cout << "	Variables initilized!" << endl << endl ;
    int NEvents = em->evtTree->GetEntries();
    cout << "	Number of events = " << NEvents << endl;
    //int NJets = em->recoJetTree->GetEntries();
    int NJets = em->recoJetTree->GetEntries();
    cout << "     Number of jets = " << NJets << endl;

    // define event filters
    if(doSingleMuonSample){
      em->regEventFilter(NeventFilters_SingleMuon, eventFilters_SingleMuon);
    }
    else if(doMinBiasSample){
      em->regEventFilter(NeventFilters_SingleMuon, eventFilters_SingleMuon);
    }
    else if(doHardProbesSample){
      em->regEventFilter(NeventFilters_SingleMuon, eventFilters_SingleMuon);
    }
    else{};
  
    loadFitFxn_PYTHIA_JERCorrection();
    loadFitFxn_PbPb_HLT();
  
    TRandom *randomGenerator = new TRandom2();



    TFile *f_neutrino_energy_fraction_map = TFile::Open("/eos/cms/store/group/phys_heavyions/cbennett/maps/neutrino_energy_fraction_map.root");
    TH2D *neutrino_energy_fraction_map;
    TH1D *neutrino_energy_fraction_map_proj;
    f_neutrino_energy_fraction_map->GetObject("neutrino_energy_fraction_map",neutrino_energy_fraction_map);

    TFile *f_neutrino_energy_map = TFile::Open("/eos/cms/store/group/phys_heavyions/cbennett/maps/neutrino_energy_map.root");
    TH2D *neutrino_energy_map;
    f_neutrino_energy_map->GetObject("neutrino_energy_map",neutrino_energy_map);

    TFile *f_neutrino_tag_fraction = TFile::Open("/eos/cms/store/group/phys_heavyions/cbennett/maps/neutrino_tag_fraction.root");
    TH1D *neutrino_tag_fraction;
    f_neutrino_tag_fraction->GetObject("neutrino_tag_fraction",neutrino_tag_fraction);



    // jet-energy resolution fit function
    TF1 *JER_fxn = new TF1("JER_fxn","sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))",50,500);
    JER_fxn->SetParameter(0,1.26585e-01);
    JER_fxn->SetParameter(1,-9.72986e-01);
    JER_fxn->SetParameter(2,3.67352e-04);

    if(doHiBinReweightToHardProbesJet80){
      loadFitFxn_hiBin_jet80();
    }



    
    // event loop
    int eventCounter = 0;
    int evi_frac = 0;
    for(int evi = 0; evi < NEvents; evi++){


      em->getEvent(evi); // load event info from eventMap
      //cout << "em->nPFCspart = " << em->nPFCspart << endl;

      if(evi == 0) {
	std::cout << "Processing events...\n";
      }

      if((100*evi / NEvents) % 5 == 0 && (100*evi / NEvents) > evi_frac){

	cout << "evt frac: " << evi_frac << "%" << endl;

      }

      evi_frac = 100*evi / NEvents;

      h_eventsBeforeSelection->Fill(1);

      // global event cuts
      if(fabs(em->vz) > 15.0) continue;
      // event filters
      if(em->checkEventFilter()) continue;
      // hiHF cut
      if(em->hiHF > 6000) continue;

      h_eventsAfterSelection->Fill(1);

      // In data, event weight = 1
      double w = 1.0;

      double w_reweight_hiBin = 1.0;
      if(doHiBinReweightToHardProbesJet80){
	w_reweight_hiBin = fitFxn_hiBin->Eval(em->hiBin);
	w = w * w_reweight_hiBin;
      }

   
      int matchFlag[10] = {0,0,0,0,0,0,0,0,0,0};

      int matchFlagR[10] = {0,0,0,0,0,0,0,0,0,0};
	
      int CentralityIndex = getCentBin(em->hiBin);
      if(CentralityIndex < 0) continue;

      // apply min-bias trigger if activated in config
      if(applyMinBiasTrigger){
	if(em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part1_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part2_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part3_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part4_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part5_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part6_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part7_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part8_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part9_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part10_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part11_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part12_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part13_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part14_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part15_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part16_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part17_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part18_v1 == 0 &&
	   em->HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part19_v1 == 0) continue;	
      }
    
      h_vz[0]->Fill(em->vz,w);
      h_vz[CentralityIndex]->Fill(em->vz,w);
      h_hiBin->Fill(em->hiBin,w);

      h_nPFcand[0]->Fill(em->nPFpart, w);
      h_nPFcand[CentralityIndex]->Fill(em->nPFpart, w);
      h_nPFcandCS[0]->Fill(em->nPFCspart, w);
      h_nPFcandCS[CentralityIndex]->Fill(em->nPFCspart, w);

      if(em->HLT_HIL3Mu12_v1 == 1){
	h_vz_triggerOn[0]->Fill(em->vz,w);
	h_vz_triggerOn[CentralityIndex]->Fill(em->vz,w);
	h_hiBin_triggerOn->Fill(em->hiBin,w);
      }

      if(em->HLT_HICsAK4PFJet60Eta1p5_v1 == 1){
	h_vz_jet60[0]->Fill(em->vz,w);
	h_vz_jet60[CentralityIndex]->Fill(em->vz,w);
	h_hiBin_jet60->Fill(em->hiBin,w);
      }

      if(em->HLT_HICsAK4PFJet80Eta1p5_v1 == 1){
	h_vz_jet80[0]->Fill(em->vz,w);
	h_vz_jet80[CentralityIndex]->Fill(em->vz,w);
	h_hiBin_jet80->Fill(em->hiBin,w);
      }

      if(em->HLT_HICsAK4PFJet100Eta1p5_v1 == 1){
	h_vz_jet100[0]->Fill(em->vz,w);
	h_vz_jet100[CentralityIndex]->Fill(em->vz,w);
	h_hiBin_jet100->Fill(em->hiBin,w);
      }

      h_Jet60HLT[0]->Fill(em->HLT_HICsAK4PFJet60Eta1p5_v1);
      h_Jet60HLT_Prescale[0]->Fill(em->HLT_HICsAK4PFJet60Eta1p5_v1_Prescl);
      h_Jet80HLT[0]->Fill(em->HLT_HICsAK4PFJet80Eta1p5_v1);
      h_Jet80HLT_Prescale[0]->Fill(em->HLT_HICsAK4PFJet80Eta1p5_v1_Prescl);
      h_Jet100HLT[0]->Fill(em->HLT_HICsAK4PFJet100Eta1p5_v1);
      h_Jet100HLT_Prescale[0]->Fill(em->HLT_HICsAK4PFJet100Eta1p5_v1_Prescl);

      h_Jet60HLT[CentralityIndex]->Fill(em->HLT_HICsAK4PFJet60Eta1p5_v1);
      h_Jet60HLT_Prescale[CentralityIndex]->Fill(em->HLT_HICsAK4PFJet60Eta1p5_v1_Prescl);
      h_Jet80HLT[CentralityIndex]->Fill(em->HLT_HICsAK4PFJet80Eta1p5_v1);
      h_Jet80HLT_Prescale[CentralityIndex]->Fill(em->HLT_HICsAK4PFJet80Eta1p5_v1_Prescl);
      h_Jet100HLT[CentralityIndex]->Fill(em->HLT_HICsAK4PFJet100Eta1p5_v1);
      h_Jet100HLT_Prescale[CentralityIndex]->Fill(em->HLT_HICsAK4PFJet100Eta1p5_v1_Prescl);

      bool evtTriggerDecision = false;

      int triggerDecision_mu5 = em->HLT_HIL3Mu5_NHitQ10_v1;
      int triggerDecision_mu5_Prescl = em->HLT_HIL3Mu5_NHitQ10_v1_Prescl;

      int triggerDecision_mu7 = em->HLT_HIL3Mu7_NHitQ10_v1;
      int triggerDecision_mu7_Prescl = em->HLT_HIL3Mu7_NHitQ10_v1_Prescl;

      int triggerDecision_mu12 = em->HLT_HIL3Mu12_v1;
      int triggerDecision_mu12_Prescl = em->HLT_HIL3Mu12_v1_Prescl;

    
      // apply jet-trigger if activated in config
      if(applyJet60Trigger){
	if(em->HLT_HICsAK4PFJet60Eta1p5_v1 == 0) continue;
      }
      if(applyJet80Trigger){
	if(em->HLT_HICsAK4PFJet80Eta1p5_v1 == 0) continue;
      }
      if(applyJet100Trigger){
	if(em->HLT_HICsAK4PFJet100Eta1p5_v1 == 0) continue;
      }

 
      if(fillMu5){
	if(triggerIsOn(triggerDecision_mu5,triggerDecision_mu5_Prescl)){
	  evtTriggerDecision = true;
	  // h_vz_triggerOn->Fill(em->vz,w);
	  // h_hiBin_triggerOn->Fill(em->hiBin,w);
	  eventCounter++;
	  h_nMu_triggerOn->Fill(em->nMu,w);
	}
      }
      else if(fillMu7){
	if(triggerIsOn(triggerDecision_mu7,triggerDecision_mu7_Prescl)){
	  evtTriggerDecision = true;
	  // h_vz_triggerOn->Fill(em->vz,w);
	  // h_hiBin_triggerOn->Fill(em->hiBin,w);
	  eventCounter++;
	  h_nMu_triggerOn->Fill(em->nMu,w);
	}
      }
      else if(fillMu12){
	if(triggerIsOn(triggerDecision_mu12,triggerDecision_mu12_Prescl)){
	  evtTriggerDecision = true;
	  // h_vz_triggerOn->Fill(em->vz,w);
	  // h_hiBin_triggerOn->Fill(em->hiBin,w);
	  eventCounter++;
	  h_nMu_triggerOn->Fill(em->nMu,w);
	}
      }
      else{};

      // Does this event pass the signal jet-pT selection? Computed here, ahead
      // of the RECO JET LOOP further down, because the RC block right below
      // needs it. Only meaningful same-event (doEventMixing pulls candidates
      // from OTHER events, so "does THIS event have a jet" does not apply to
      // what is being mixed in). See doSignalSelectedRC in pseudoJets.h.
      bool eventHasSignalJet = false;
      if(doSignalSelectedRC && !doEventMixing){
	for(int i = 0; i < em->njet; i++){
	  double jetPt_i;
	  if(signalJetPtCutIsRaw){
	    jetPt_i = em->rawpt[i];
	  }
	  else{
	    JEC_PF.SetJetPT(em->rawpt[i]);
	    JEC_PF.SetJetEta(em->jeteta[i]);
	    JEC_PF.SetJetPhi(em->jetphi[i]);
	    jetPt_i = JEC_PF.GetCorrectedPT();
	  }
	  if(jetPt_i > signalJetPtCut){ eventHasSignalJet = true; break; }
	}
      }



      ///// PF Candidate Analyzer / Pseudo-Jet Calculator / FastJet Clustering

      int NCandidatesToSample = em->nPFpart;
      
      // pre-load mixed-event PF candidates from same-centrality events into a pool
      std::vector<double> pool_pfPt, pool_pfEta, pool_pfPhi, pool_pfId;

      // DONOR EVENT for the muon-injection study. The first centrality-matched
      // event at or after evi+1 -- which is exactly the first event this pool
      // loop accepts, so it is captured here rather than costing another
      // getEvent() pass. Its candidates are clustered on their own (baseline)
      // and again with a mixed-event muon added (injected), so the two can be
      // compared jet by jet.
      std::vector<double> donor_pfPt, donor_pfEta, donor_pfPhi;
      std::vector<int>    donor_pfId;

      if(doEventMixing){
	int eventsInPool = 0;
	int jPool = 0;
	while(eventsInPool < N_mixedEventsInPool && jPool < NEvents){
	  int mixedEventIndex = (evi + jPool + 1) % NEvents;
	  em->getEvent(mixedEventIndex);
	  if(getCentBin(em->hiBin) != CentralityIndex){ jPool++; continue; }
	  for(int l = 0; l < em->nPFpart; l++){
	    pool_pfPt.push_back(em->pfPt->at(l));
	    pool_pfEta.push_back(em->pfEta->at(l));
	    pool_pfPhi.push_back(em->pfPhi->at(l));
	    pool_pfId.push_back(em->pfId->at(l));
	    if(eventsInPool == 0){          // this event is the donor
	      donor_pfPt.push_back(em->pfPt->at(l));
	      donor_pfEta.push_back(em->pfEta->at(l));
	      donor_pfPhi.push_back(em->pfPhi->at(l));
	      donor_pfId.push_back(em->pfId->at(l));
	    }
	  }
	  eventsInPool++;
	  jPool++;
	}
	em->getEvent(evi); // restore current event
      }

      int poolSize = (int)pool_pfPt.size();
      std::mt19937 rng(std::random_device{}());
      double pi_pfcand = TMath::Pi();
      double dR_max_pfcand = 0.4;

#ifdef DO_FASTJET
      // Baseline clustering of the donor event, done ONCE per event since the
      // donor does not change with the resample index. The injected version is
      // reclustered per muon inside the mixed-muon loop below.
      //
      // NOTE: no jetTrkMax filter is applied to either collection here. That
      // filter rejects jets with trkMax/jetPt > 0.98, and injecting a 15+ GeV
      // muon into a soft combinatorial jet drives that ratio towards 1 -- it
      // would preferentially delete exactly the jets where the injection
      // matters most, and would change the jet population between the two
      // collections, which is the one thing this comparison must not do.
      std::vector<fastjet::PseudoJet> donorInputs;
      for(size_t l = 0; l < donor_pfPt.size(); l++){
	double pt = donor_pfPt[l];
	if(pt < pseudoJetCandPt_min) continue;
	double eta = donor_pfEta[l], phi = donor_pfPhi[l];
	fastjet::PseudoJet pj(pt*TMath::Cos(phi), pt*TMath::Sin(phi),
			      pt*TMath::SinH(eta), pt*TMath::CosH(eta));
	pj.set_user_info(new CandInfo((int)l, donor_pfId[l]));
	donorInputs.push_back(pj);
      }
      fastjet::JetDefinition donorJetDef(fastjet::antikt_algorithm, dR_max_pfcand);
      // unique_ptr, not a raw new/delete pair: the PseudoJets in donorJets stay
      // valid only while their ClusterSequence lives, so it has to outlast the
      // muon loop, and there is no safe single place to delete it that a future
      // "continue" in the event loop could not skip. One allocation per event
      // leaked would be fatal over a full scan.
      std::unique_ptr<fastjet::ClusterSequence> donorCS;
      std::vector<fastjet::PseudoJet> donorJets;
      if(!donorInputs.empty()){
	donorCS.reset(new fastjet::ClusterSequence(donorInputs, donorJetDef));
	donorJets = fastjet::sorted_by_pt(donorCS->inclusive_jets(0.));
      }
#endif

      for(int k = 0; k < N_mixedEventsInPool; k++){

	double randEta_k = 3.2*randomGenerator->Rndm() - 1.6;
	double randPhi_k = 2*pi_pfcand*randomGenerator->Rndm() - pi_pfcand;
	double pseudoJetPt_k = 0.;
	double pseudoJetPt_geoCorr_k = 0.;

	if(doEventMixing){
	  if(poolSize > 0){
	    std::uniform_int_distribution<int> poolDist(0, poolSize - 1);
	    for(int s = 0; s < NCandidatesToSample; s++){
	      int idx = poolDist(rng);
	      double pfPt_l  = pool_pfPt[idx];
	      double pfEta_l = pool_pfEta[idx];
	      double pfPhi_l = pool_pfPhi[idx];
	      double dR_kl = getDr(randEta_k, randPhi_k, pfEta_l, pfPhi_l);
	      if(pfPt_l > pseudoJetCandPt_min && dR_kl < dR_max_pfcand){
		h_pfPt[0]->Fill(pfPt_l, w);
		h_pfPt[CentralityIndex]->Fill(pfPt_l, w);
		pseudoJetPt_k += pfPt_l;
		pseudoJetPt_geoCorr_k += pfPt_l * TMath::Cos(dR_kl);
	      }
	    }
	  }
	}
	else{
	  for(int l = 0; l < em->nPFpart; l++){
	    double pfPt_l  = em->pfPt->at(l);
	    double pfEta_l = em->pfEta->at(l);
	    double pfPhi_l = em->pfPhi->at(l);
	    double dR_kl = getDr(randEta_k, randPhi_k, pfEta_l, pfPhi_l);
	    if(pfPt_l > pseudoJetCandPt_min && dR_kl < dR_max_pfcand){
	      h_pfPt[0]->Fill(pfPt_l, w);
	      h_pfPt[CentralityIndex]->Fill(pfPt_l, w);
	      pseudoJetPt_k += pfPt_l;
	      pseudoJetPt_geoCorr_k += pfPt_l * TMath::Cos(dR_kl);
	    }
	  }
	}

	h_pseudoJetPt[0]->Fill(pseudoJetPt_k, w);
	h_pseudoJetPt[CentralityIndex]->Fill(pseudoJetPt_k, w);

	h_pseudoJetPt_geoCorr[0]->Fill(pseudoJetPt_geoCorr_k, w);
	h_pseudoJetPt_geoCorr[CentralityIndex]->Fill(pseudoJetPt_geoCorr_k, w);

	// record cone pT at its (eta,phi) throw location for the UE map
	h_randConeEtaPhi[0]->Fill(randEta_k, randPhi_k, pseudoJetPt_k, w);
	h_randConeEtaPhi[CentralityIndex]->Fill(randEta_k, randPhi_k, pseudoJetPt_k, w);

	h_randConeEtaPhi_geoCorr[0]->Fill(randEta_k, randPhi_k, pseudoJetPt_geoCorr_k, w);
	h_randConeEtaPhi_geoCorr[CentralityIndex]->Fill(randEta_k, randPhi_k, pseudoJetPt_geoCorr_k, w);

	if(eventHasSignalJet){
	  h_pseudoJetPt_sigSel[0]->Fill(pseudoJetPt_k, w);
	  h_pseudoJetPt_sigSel[CentralityIndex]->Fill(pseudoJetPt_k, w);
	  h_randConeEtaPhi_sigSel[0]->Fill(randEta_k, randPhi_k, pseudoJetPt_k, w);
	  h_randConeEtaPhi_sigSel[CentralityIndex]->Fill(randEta_k, randPhi_k, pseudoJetPt_k, w);
	}



	
      }

#ifdef DO_FASTJET
      if(doFastJetClustering){

	// Redraw NCandidatesToSample from the SAME pool this many times, each a
	// fresh independent FastJet clustering, to multiply the statistics of
	// every histogram filled below without scanning more real events. Every
	// resample fills with weight w_resample = w/N_fastJetResamples, so the
	// TOTAL contribution of this one real event is unchanged regardless of
	// N_fastJetResamples -- only the pool-sampling noise shrinks. See
	// N_fastJetMixedEventResamples in pseudoJets.h for the full rationale,
	// including why this is a DIFFERENT convention from the pre-existing
	// random-cone resampling (h_pseudoJetPt etc.), which fills every throw
	// at the full weight w and requires dividing by N_pool downstream.
	// No effect when !doEventMixing: same-event clustering is deterministic
	// (identical candidates every time), so this runs exactly once.
	int N_fastJetResamples = doEventMixing ? N_fastJetMixedEventResamples : 1;
	double w_resample = w / N_fastJetResamples;

	for(int fjResample = 0; fjResample < N_fastJetResamples; fjResample++){

	std::vector<double> mixedEventPFCandidates_pt;
	std::vector<double> mixedEventPFCandidates_eta;
	std::vector<double> mixedEventPFCandidates_phi;
	std::vector<int> mixedEventPFCandidates_id;
	// FastJet anti-kT clustering on PF candidates
        std::vector<fastjet::PseudoJet> fjInputs;
        if(doEventMixing && poolSize > 0){
          std::uniform_int_distribution<int> poolDist(0, poolSize - 1);
          for(int s = 0; s < NCandidatesToSample; s++){
            int idx = poolDist(rng);
            double pt  = pool_pfPt[idx];
            double eta = pool_pfEta[idx];
            double phi = pool_pfPhi[idx];
	    int id = pool_pfId[idx];
            if(pt < pseudoJetCandPt_min) continue;
            double px = pt * TMath::Cos(phi);
            double py = pt * TMath::Sin(phi);
            double pz = pt * TMath::SinH(eta);
            double E  = pt * TMath::CosH(eta);
	    bool isCharged = (id == 1 || id == 2 || id == 3); // id = h, e, mu
	    fastjet::PseudoJet pseudoJet_s(px, py, pz, E);
	    //pseudoJet_s.set_user_index(isCharged ? 1 : 0);
	    pseudoJet_s.set_user_info(new CandInfo(idx, id));
            fjInputs.push_back(pseudoJet_s);
	    mixedEventPFCandidates_pt.push_back(pt);
	    mixedEventPFCandidates_eta.push_back(eta);
	    mixedEventPFCandidates_phi.push_back(phi);
	    mixedEventPFCandidates_id.push_back(id);
          }
        }
        else{
          for(int l = 0; l < em->nPFpart; l++){
            double pt  = em->pfPt->at(l);
            double eta = em->pfEta->at(l);
            double phi = em->pfPhi->at(l);
	    int id = em->pfId->at(l);
            if(pt < pseudoJetCandPt_min) continue;
            double px = pt * TMath::Cos(phi);
            double py = pt * TMath::Sin(phi);
            double pz = pt * TMath::SinH(eta);
            double E  = pt * TMath::CosH(eta);
	    bool isCharged = (id == 1 || id == 2 || id == 3); // id = h, e, mu
	    fastjet::PseudoJet pseudoJet_l(px, py, pz, E);
	    //pseudoJet_l.set_user_index(isCharged ? 1 : 0);
	    pseudoJet_l.set_user_info(new CandInfo(l, id));
            fjInputs.push_back(pseudoJet_l);
          }
        }
        // Diagnostic: how many candidates actually went into the clustering.
        // In same-event running this must equal em->nPFpart (i.e. h_nPFcand),
        // since the loop just copies the event. In mixed-event running it is the
        // number drawn from the pool, which is meant to be NCandidatesToSample =
        // em->nPFpart as well -- comparing this against h_nPFcand is the check
        // that the mixed event is being built at the right multiplicity.
        h_nPFcandFastJet[0]->Fill((int)fjInputs.size(), w_resample);
        h_nPFcandFastJet[CentralityIndex]->Fill((int)fjInputs.size(), w_resample);

        fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, dR_max_pfcand);
        fastjet::ClusterSequence cs(fjInputs, jetDef);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(0.));
        for(const auto& jet : jets){
          if(TMath::Abs(jet.eta()) > 1.6) continue;
	  std::vector<fastjet::PseudoJet> constituents = jet.constituents();
	  if(skipSingleConstituentJets){
	    if (constituents.size() < 2) continue;  // cut single-track/single-constituent jets
	  }
	  
	  double trackMaxPt = 0.0;
	  bool hasFastJetRecoMuonTag = false;
	  double fastJetMuonPt = 0.;
	  double fastJetMuonEta = 0.;
	  double fastJetMuonPhi = 0.;
	  double fastJetMuonPtRel = -999.;
	  double fastJetMuonDR = -999.;
	  for(const auto& c : constituents){
	    //if(c.user_index() == 1 && c.pt() > trackMaxPt) trackMaxPt = c.pt();
	    if(!c.has_user_info<CandInfo>()) continue;
	    const CandInfo &candinfo = c.user_info<CandInfo>();
	    if(candinfo.isCharged() && c.pt() > trackMaxPt) trackMaxPt = c.pt();
	    // c.pt() > fastJetMuonPt keeps the LEADING muon constituent. Without it
	    // this kept whichever muon happened to come last in constituent order,
	    // which is an arbitrary choice of tag when a jet contains more than one
	    // candidate above muPtCut. findRecoMuonTag, which the (real mu)
	    // templates use, resolves ties by muon index rather than pT, but at
	    // least "leading" is a defined quantity.
	    if(candinfo.getId() == 3 && c.pt() > muPtCut && evtTriggerDecision && c.pt() > fastJetMuonPt) {
	      hasFastJetRecoMuonTag = true;
	      fastJetMuonPt = c.pt();
	      fastJetMuonEta = c.eta();
	      fastJetMuonPhi = c.phi_std();
	    }
	    
	    
	  }
	           
          JEC_PF.SetJetPT(jet.pt());
          JEC_PF.SetJetEta(jet.eta());
          JEC_PF.SetJetPhi(jet.phi_std());
          double fastJetPt_JEC = JEC_PF.GetCorrectedPT();
	  
	  if(doJetTrkMaxFilter){
	    if(!passesJetTrkMaxFilter(trackMaxPt,fastJetPt_JEC)) continue;
	  }


	  
	  h_fastJetPt_PF[0]->Fill(jet.pt(), w_resample);
	  h_fastJetPt_PF[CentralityIndex]->Fill(jet.pt(), w_resample);
	  h_fastJetPt_PF_JEC[0]->Fill(fastJetPt_JEC, w_resample);
          h_fastJetPt_PF_JEC[CentralityIndex]->Fill(fastJetPt_JEC, w_resample);
          // RC-subtracted fastJet pT: subtract mean UE pT at the jet's (eta,phi) location
          if(h_RC_map[CentralityIndex]){
            double rcMeanPt = h_RC_map[CentralityIndex]->GetBinContent(
                                h_RC_map[CentralityIndex]->FindBin(jet.eta(), jet.phi_std()));

	    double rcMeanPt_geoCorr = h_RC_geoCorr_map[CentralityIndex]->GetBinContent(
                                h_RC_geoCorr_map[CentralityIndex]->FindBin(jet.eta(), jet.phi_std()));


	    double rcMeanPt_geoCorr_etaReflect = h_RC_geoCorr_map[CentralityIndex]->GetBinContent(
                                h_RC_geoCorr_map[CentralityIndex]->FindBin(-1.*jet.eta(), jet.phi_std()));


	    double dPTMeanPt = h_dPT_map[CentralityIndex]->GetBinContent(h_dPT_map[CentralityIndex]->FindBin(jet.eta(), jet.phi_std()));

	    // h_dPT_geoCorr_map is only present in maps files regenerated after the
	    // geoCorr dPT map was added to this analyzer -- guard against an
	    // older maps file that predates it (GetObject leaves the pointer
	    // null rather than throwing), instead of segfaulting on every event.
	    double dPTMeanPt_geoCorr = h_dPT_geoCorr_map[CentralityIndex] ?
	      h_dPT_geoCorr_map[CentralityIndex]->GetBinContent(h_dPT_geoCorr_map[CentralityIndex]->FindBin(jet.eta(), jet.phi_std())) : 0.;

	    double dPTMeanPt_dPTAbove0 = h_dPT_dPTAbove0_map[CentralityIndex]->GetBinContent(h_dPT_dPTAbove0_map[CentralityIndex]->FindBin(jet.eta(), jet.phi_std()));

	    double dPTMeanPt_PFCsPTAbove60 = h_dPT_PFCsPTAbove60_map[CentralityIndex]->GetBinContent(h_dPT_PFCsPTAbove60_map[CentralityIndex]->FindBin(jet.eta(), jet.phi_std()));



	    double fastJetPt_rcSub = jet.pt() - rcMeanPt;
	    double fastJetPt_rcSub_geoCorr = jet.pt() - rcMeanPt_geoCorr;
	    double fastJetPt_rcSub_geoCorr_etaReflect = jet.pt() - rcMeanPt_geoCorr_etaReflect;
	    double fastJetPt_dPTSub = jet.pt() - dPTMeanPt;
	    double fastJetPt_dPTSub_geoCorr = jet.pt() - dPTMeanPt_geoCorr;
	    double fastJetPt_dPTSub_dPTAbove0 = jet.pt() - dPTMeanPt_dPTAbove0;
	    double fastJetPt_dPTSub_PFCsPTAbove60 = jet.pt() - dPTMeanPt_PFCsPTAbove60;

	    JEC_PF.SetJetEta(jet.eta());
	    JEC_PF.SetJetPhi(jet.phi_std());
	    
	    JEC_PF.SetJetPT(fastJetPt_rcSub);
	    double fastJetPt_JEC_rcSub = JEC_PF.GetCorrectedPT();

	    JEC_PF.SetJetPT(fastJetPt_rcSub_geoCorr);
	    double fastJetPt_JEC_rcSub_geoCorr = JEC_PF.GetCorrectedPT();

	    JEC_PF.SetJetPT(fastJetPt_rcSub_geoCorr_etaReflect);
	    double fastJetPt_JEC_rcSub_geoCorr_etaReflect = JEC_PF.GetCorrectedPT();

	    JEC_PF.SetJetPT(fastJetPt_dPTSub);
	    double fastJetPt_JEC_dPTSub = JEC_PF.GetCorrectedPT();

	    JEC_PF.SetJetPT(fastJetPt_dPTSub_geoCorr);
	    double fastJetPt_JEC_dPTSub_geoCorr = JEC_PF.GetCorrectedPT();

	    JEC_PF.SetJetPT(fastJetPt_dPTSub_dPTAbove0);
	    double fastJetPt_JEC_dPTSub_dPTAbove0 = JEC_PF.GetCorrectedPT();

	    JEC_PF.SetJetPT(fastJetPt_dPTSub_PFCsPTAbove60);
	    double fastJetPt_JEC_dPTSub_PFCsPTAbove60 = JEC_PF.GetCorrectedPT();


	    

	    if(fastJetPt_JEC_rcSub > 20.){
	      h_fastJetPt_PF_bkgSub_RC[0]->Fill(fastJetPt_rcSub, w_resample);
              h_fastJetPt_PF_bkgSub_RC[CentralityIndex]->Fill(fastJetPt_rcSub, w_resample);

	      if(hasFastJetRecoMuonTag){
		// JEC-corrected, background-subtracted jet pT throughout: for the
		// ptRel itself (which is built on the jet's momentum, so the scale
		// matters) and for the Y axis of both histograms. These are the
		// (fake mu, fake jet) term of the ptRel decomposition and are
		// subtracted from h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn,
		// whose jet pT is JEC-corrected -- on the raw rcSub scale a given
		// pT window would select systematically harder jets here than in the
		// data it corrects.
		//
		// The neighbouring 1D h_fastJetPt_PF_* spectra deliberately keep BOTH
		// scales as a matched pair, distinguished by _JEC_ in the name. These
		// two carry no such marker and had no JEC twin, so raw was simply the
		// wrong scale rather than a documented choice.
		fastJetMuonPtRel = getPtRel(fastJetMuonPt,fastJetMuonEta,fastJetMuonPhi,fastJetPt_JEC_rcSub,jet.eta(),jet.phi_std());
		fastJetMuonDR = getDr(fastJetMuonEta,fastJetMuonPhi,jet.eta(),jet.phi_std());
		h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC[0]->Fill(fastJetMuonPtRel,fastJetPt_JEC_rcSub,w_resample);
		h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC[CentralityIndex]->Fill(fastJetMuonPtRel,fastJetPt_JEC_rcSub,w_resample);
		h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC[0]->Fill(fastJetMuonDR,fastJetPt_JEC_rcSub,w_resample);
		h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC[CentralityIndex]->Fill(fastJetMuonDR,fastJetPt_JEC_rcSub,w_resample);
	      }

	      h_fastJetPt_PF_JEC_bkgSub_RC[0]->Fill(fastJetPt_JEC_rcSub, w_resample);
              h_fastJetPt_PF_JEC_bkgSub_RC[CentralityIndex]->Fill(fastJetPt_JEC_rcSub, w_resample);

	      if(eventHasSignalJet){
	        h_fastJetPt_PF_bkgSub_RC_sigSel[0]->Fill(fastJetPt_rcSub, w_resample);
                h_fastJetPt_PF_bkgSub_RC_sigSel[CentralityIndex]->Fill(fastJetPt_rcSub, w_resample);

	        h_fastJetPt_PF_JEC_bkgSub_RC_sigSel[0]->Fill(fastJetPt_JEC_rcSub, w_resample);
                h_fastJetPt_PF_JEC_bkgSub_RC_sigSel[CentralityIndex]->Fill(fastJetPt_JEC_rcSub, w_resample);
	      }
	    }

	    if(fastJetPt_JEC_rcSub_geoCorr > 20.){
	      h_fastJetPt_PF_bkgSub_RC_geoCorr[0]->Fill(fastJetPt_rcSub_geoCorr, w_resample);
              h_fastJetPt_PF_bkgSub_RC_geoCorr[CentralityIndex]->Fill(fastJetPt_rcSub_geoCorr, w_resample);

	      h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr[0]->Fill(fastJetPt_JEC_rcSub_geoCorr, w_resample);
              h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr[CentralityIndex]->Fill(fastJetPt_JEC_rcSub_geoCorr, w_resample);
	    }

	    if(fastJetPt_JEC_rcSub_geoCorr_etaReflect > 20.){
	      h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect[0]->Fill(fastJetPt_rcSub_geoCorr_etaReflect, w_resample);
              h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect[CentralityIndex]->Fill(fastJetPt_rcSub_geoCorr_etaReflect, w_resample);

	      h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect[0]->Fill(fastJetPt_JEC_rcSub_geoCorr_etaReflect, w_resample);
              h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect[CentralityIndex]->Fill(fastJetPt_JEC_rcSub_geoCorr_etaReflect, w_resample);
	    }

	    if(fastJetPt_JEC_dPTSub > 20.){
	      h_fastJetPt_PF_bkgSub_dPT[0]->Fill(fastJetPt_dPTSub,w_resample);
	      h_fastJetPt_PF_bkgSub_dPT[CentralityIndex]->Fill(fastJetPt_dPTSub,w_resample);

	      h_fastJetPt_PF_JEC_bkgSub_dPT[0]->Fill(fastJetPt_JEC_dPTSub,w_resample);
	      h_fastJetPt_PF_JEC_bkgSub_dPT[CentralityIndex]->Fill(fastJetPt_JEC_dPTSub,w_resample);
	    }

	    if(h_dPT_geoCorr_map[CentralityIndex] && fastJetPt_JEC_dPTSub_geoCorr > 20.){
	      h_fastJetPt_PF_bkgSub_dPT_geoCorr[0]->Fill(fastJetPt_dPTSub_geoCorr,w_resample);
	      h_fastJetPt_PF_bkgSub_dPT_geoCorr[CentralityIndex]->Fill(fastJetPt_dPTSub_geoCorr,w_resample);

	      h_fastJetPt_PF_JEC_bkgSub_dPT_geoCorr[0]->Fill(fastJetPt_JEC_dPTSub_geoCorr,w_resample);
	      h_fastJetPt_PF_JEC_bkgSub_dPT_geoCorr[CentralityIndex]->Fill(fastJetPt_JEC_dPTSub_geoCorr,w_resample);
	    }

	    if(fastJetPt_JEC_dPTSub_dPTAbove0 > 20.){
	      h_fastJetPt_PF_bkgSub_dPT_dPTAbove0[0]->Fill(fastJetPt_dPTSub_dPTAbove0,w_resample);
	      h_fastJetPt_PF_bkgSub_dPT_dPTAbove0[CentralityIndex]->Fill(fastJetPt_dPTSub_dPTAbove0,w_resample);

	      h_fastJetPt_PF_JEC_bkgSub_dPT_dPTAbove0[0]->Fill(fastJetPt_JEC_dPTSub_dPTAbove0,w_resample);
	      h_fastJetPt_PF_JEC_bkgSub_dPT_dPTAbove0[CentralityIndex]->Fill(fastJetPt_JEC_dPTSub_dPTAbove0,w_resample);
	    }

	    if(fastJetPt_JEC_dPTSub_PFCsPTAbove60 > 20.){
	      h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60[0]->Fill(fastJetPt_dPTSub_PFCsPTAbove60,w_resample);
	      h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60[CentralityIndex]->Fill(fastJetPt_dPTSub_PFCsPTAbove60,w_resample);

	      h_fastJetPt_PF_JEC_bkgSub_dPT_PFCsPTAbove60[0]->Fill(fastJetPt_JEC_dPTSub_PFCsPTAbove60,w_resample);
	      h_fastJetPt_PF_JEC_bkgSub_dPT_PFCsPTAbove60[CentralityIndex]->Fill(fastJetPt_JEC_dPTSub_PFCsPTAbove60,w_resample);
	    }
	    
          }
        }

	if(doEventMixing){
	  for(int i = 0; i < mixedEventPFCandidates_id.size(); i++){

	    if(mixedEventPFCandidates_id.at(i) != 3) continue; // skip if candidate is not a muon
	    if(mixedEventPFCandidates_pt.at(i) < muPtCut) continue;
	    if(fabs(mixedEventPFCandidates_eta.at(i)) > 2.) continue;
	    
	    double fastJetMuonDR_i = 999.;
	    double jetPt_JEC_rcSub_match_i = -999.;

	    for(const auto& jet : jets){

	      double muonJetDR_ij = 999.;

	      if(fabs(jet.eta()) > 1.6) continue;
	      std::vector<fastjet::PseudoJet> constituents = jet.constituents();
	      double trackMaxPt = 0.0;
	      for(const auto& c : constituents){
		if(!c.has_user_info<CandInfo>()) continue;
		const CandInfo &candinfo = c.user_info<CandInfo>();
		if(candinfo.isCharged() && c.pt() > trackMaxPt) trackMaxPt = c.pt();
	      }
	      JEC_PF.SetJetPT(jet.pt());
	      JEC_PF.SetJetEta(jet.eta());
	      JEC_PF.SetJetPhi(jet.phi_std());
	      double jetPt_JEC = JEC_PF.GetCorrectedPT();
	      if(!h_RC_map[CentralityIndex]) continue;
	      double rcMeanPt = h_RC_map[CentralityIndex]->GetBinContent(
                                h_RC_map[CentralityIndex]->FindBin(jet.eta(), jet.phi_std()));
	      double jetPt_rcSub = jet.pt() - rcMeanPt;
	      JEC_PF.SetJetPT(jetPt_rcSub);
	      double jetPt_JEC_rcSub = JEC_PF.GetCorrectedPT();

	      if(jetPt_JEC_rcSub < 20.) continue;
	      
	      if(doJetTrkMaxFilter){
		if(!passesJetTrkMaxFilter(trackMaxPt,jetPt_JEC)) continue;
	      }

	      muonJetDR_ij = getDr(mixedEventPFCandidates_eta.at(i),mixedEventPFCandidates_phi.at(i),jet.eta(),jet.phi_std());
	      if(muonJetDR_ij < fastJetMuonDR_i) {
		fastJetMuonDR_i = muonJetDR_ij;
		jetPt_JEC_rcSub_match_i = jetPt_JEC_rcSub;
	      }

	    }
	    h_fastJetMuonDR_inclusiveClosestFastJet[0]->Fill(fastJetMuonDR_i,jetPt_JEC_rcSub_match_i,w_resample);
	    h_fastJetMuonDR_inclusiveClosestFastJet[CentralityIndex]->Fill(fastJetMuonDR_i,jetPt_JEC_rcSub_match_i,w_resample);

	    double fastJetMuonDR_recoJet_i = 999.;
	    double recoJet_match_i = 0.;
	    // jet axis of the closest-reco-jet match, kept so ptRel can be formed
	    // below -- jet pT alone is not enough
	    double recoJetEta_match_i = -999.;
	    double recoJetPhi_match_i = -999.;

	    // loop through jets
	    for(int j = 0; j < em->njet; j++){

	      double jetPt_j = 0.;
	      if(useCaloJetsOverride){
		JEC_Calo.SetJetPT(em->rawpt[j]);
		JEC_Calo.SetJetEta(em->jeteta[j]);
		JEC_Calo.SetJetPhi(em->jetphi[j]);
		jetPt_j = JEC_Calo.GetCorrectedPT();
	      }
	      else{
		JEC_PF.SetJetPT(em->rawpt[j]);
		JEC_PF.SetJetEta(em->jeteta[j]);
		JEC_PF.SetJetPhi(em->jetphi[j]);
		jetPt_j = JEC_PF.GetCorrectedPT();
	      }

	      double jetEta_j = em->jeteta[j];
	      double jetPhi_j = em->jetphi[j];
	      double jetTrkMax_j = em->jetTrkMax[j];

	      if(jetPt_j < 20.) continue;
	      if(fabs(jetEta_j) > etaMax) continue;
	      if(doJetTrkMaxFilter){
		if(!passesJetTrkMaxFilter(jetTrkMax_j,jetPt_j)) continue;
	      }

	      double dR_ij = getDr(mixedEventPFCandidates_eta.at(i),mixedEventPFCandidates_phi.at(i),jetEta_j,jetPhi_j);
	      if(dR_ij < fastJetMuonDR_recoJet_i) {
		fastJetMuonDR_recoJet_i = dR_ij;
		recoJet_match_i = jetPt_j;
		recoJetEta_match_i = jetEta_j;
		recoJetPhi_match_i = jetPhi_j;
	      }

	      
	    
	    }

	    h_muonDR_inclusiveClosestJet[0]->Fill(fastJetMuonDR_recoJet_i,recoJet_match_i,w_resample);
	    h_muonDR_inclusiveClosestJet[CentralityIndex]->Fill(fastJetMuonDR_recoJet_i,recoJet_match_i,w_resample);
	    // same quantities, same weight, gated so this can be compared against
	    // the trigger-gated ptRel template below
	    if(evtTriggerDecision){
	      h_muonDR_inclusiveClosestJet_triggerOn[0]->Fill(fastJetMuonDR_recoJet_i,recoJet_match_i,w_resample);
	      h_muonDR_inclusiveClosestJet_triggerOn[CentralityIndex]->Fill(fastJetMuonDR_recoJet_i,recoJet_match_i,w_resample);
	    }

	    // (fake mu, real jet) template: a mixed-event PF muon that lands close
	    // enough to a REAL reco jet to have been tagged. The dR < epsilon_mm
	    // requirement is what makes this a background estimate for the data
	    // tagging rate rather than just a nearest-neighbour distance -- the
	    // h_muonDR histogram above deliberately has no such cut.
	    //
	    // evtTriggerDecision matches the guard on the data fill: without it
	    // this template accumulates in events the data histogram never sees,
	    // and the per-event rates being subtracted are normalised over
	    // different event samples. Note it gates only this fill, NOT the
	    // enclosing muon loop -- the two h_*muonDR_inclusiveClosest*
	    // histograms above are pre-existing diagnostics with their own
	    // (untriggered) convention, and changing them is out of scope here.
	    if(evtTriggerDecision && fastJetMuonDR_recoJet_i < epsilon_mm && recoJet_match_i > 0.){
	      double mixedMuonPtRel_i = getPtRel(mixedEventPFCandidates_pt.at(i),
						 mixedEventPFCandidates_eta.at(i),
						 mixedEventPFCandidates_phi.at(i),
						 recoJet_match_i,
						 recoJetEta_match_i,
						 recoJetPhi_match_i);
	      h_mixedMuonPtRel_recoJetPt[0]->Fill(mixedMuonPtRel_i,recoJet_match_i,w_resample);
	      h_mixedMuonPtRel_recoJetPt[CentralityIndex]->Fill(mixedMuonPtRel_i,recoJet_match_i,w_resample);
	    }

	    // ---- muon-injection study ------------------------------------
	    // Same mixed-event muon, matched to donor-event fastJets twice:
	    // once against jets clustered without it, once after adding it to
	    // the donor candidates and reclustering. See the histogram
	    // declarations for what the difference measures.
	    //
	    // No jetTrkMax filter here, deliberately: it cuts trkMax/jetPt >
	    // 0.98, and a 15+ GeV muon injected into a soft combinatorial jet
	    // approaches that, so it would preferentially remove the jets the
	    // study is about and change the population between the two
	    // collections.
	    if(!donorJets.empty() && h_RC_map[CentralityIndex]){

	      double muPt_inj  = mixedEventPFCandidates_pt.at(i);
	      double muEta_inj = mixedEventPFCandidates_eta.at(i);
	      double muPhi_inj = mixedEventPFCandidates_phi.at(i);

	      // (a) BASELINE: donor jets as clustered WITHOUT the muon
	      double drB = 999., ptB = -999., etaB = -999., phiB = -999.;
	      for(const auto& dj : donorJets){
		if(fabs(dj.eta()) > 1.6) continue;
		double rc = h_RC_map[CentralityIndex]->GetBinContent(
			      h_RC_map[CentralityIndex]->FindBin(dj.eta(), dj.phi_std()));
		JEC_PF.SetJetPT(dj.pt() - rc);
		JEC_PF.SetJetEta(dj.eta());
		JEC_PF.SetJetPhi(dj.phi_std());
		double djPt = JEC_PF.GetCorrectedPT();
		if(djPt < 20.) continue;
		h_donorJetPt_noInject[0]->Fill(djPt, w_resample);
		h_donorJetPt_noInject[CentralityIndex]->Fill(djPt, w_resample);
		double dr = getDr(muEta_inj, muPhi_inj, dj.eta(), dj.phi_std());
		if(dr < drB){ drB = dr; ptB = djPt; etaB = dj.eta(); phiB = dj.phi_std(); }
	      }
	      if(ptB > 0.){
		h_injMuonDR_donorJetPt_noInject[0]->Fill(drB, ptB, w_resample);
		h_injMuonDR_donorJetPt_noInject[CentralityIndex]->Fill(drB, ptB, w_resample);
		if(drB < epsilon_mm){
		  double pr = getPtRel(muPt_inj, muEta_inj, muPhi_inj, ptB, etaB, phiB);
		  h_injMuonPtRel_donorJetPt_noInject[0]->Fill(pr, ptB, w_resample);
		  h_injMuonPtRel_donorJetPt_noInject[CentralityIndex]->Fill(pr, ptB, w_resample);
		}
	      }

	      // (b) INJECTED: add the muon to the donor candidates, recluster
	      std::vector<fastjet::PseudoJet> injInputs = donorInputs;
	      fastjet::PseudoJet muPJ(muPt_inj*TMath::Cos(muPhi_inj),
				      muPt_inj*TMath::Sin(muPhi_inj),
				      muPt_inj*TMath::SinH(muEta_inj),
				      muPt_inj*TMath::CosH(muEta_inj));
	      muPJ.set_user_info(new CandInfo(-1, 3));  // index -1 marks the injected muon
	      injInputs.push_back(muPJ);
	      fastjet::ClusterSequence injCS(injInputs, donorJetDef);
	      std::vector<fastjet::PseudoJet> injJets = fastjet::sorted_by_pt(injCS.inclusive_jets(0.));

	      double drI = 999., ptI = -999., etaI = -999., phiI = -999.;
	      double drO = -999., ptO = -999., etaO = -999., phiO = -999.;   // owning jet
	      for(const auto& ij : injJets){
		if(fabs(ij.eta()) > 1.6) continue;
		double rc = h_RC_map[CentralityIndex]->GetBinContent(
			      h_RC_map[CentralityIndex]->FindBin(ij.eta(), ij.phi_std()));
		JEC_PF.SetJetPT(ij.pt() - rc);
		JEC_PF.SetJetEta(ij.eta());
		JEC_PF.SetJetPhi(ij.phi_std());
		double ijPt = JEC_PF.GetCorrectedPT();
		if(ijPt < 20.) continue;
		h_donorJetPt_inject[0]->Fill(ijPt, w_resample);
		h_donorJetPt_inject[CentralityIndex]->Fill(ijPt, w_resample);
		double dr = getDr(muEta_inj, muPhi_inj, ij.eta(), ij.phi_std());
		if(dr < drI){ drI = dr; ptI = ijPt; etaI = ij.eta(); phiI = ij.phi_std(); }
		// is this the jet anti-kT gave the injected muon to?
		for(const auto& c : ij.constituents()){
		  if(!c.has_user_info<CandInfo>()) continue;
		  if(c.user_info<CandInfo>().getIndex() == -1){
		    drO = dr; ptO = ijPt; etaO = ij.eta(); phiO = ij.phi_std();
		    break;
		  }
		}
	      }
	      if(ptI > 0.){
		h_injMuonDR_donorJetPt_inject[0]->Fill(drI, ptI, w_resample);
		h_injMuonDR_donorJetPt_inject[CentralityIndex]->Fill(drI, ptI, w_resample);
		if(drI < epsilon_mm){
		  double pr = getPtRel(muPt_inj, muEta_inj, muPhi_inj, ptI, etaI, phiI);
		  h_injMuonPtRel_donorJetPt_inject[0]->Fill(pr, ptI, w_resample);
		  h_injMuonPtRel_donorJetPt_inject[CentralityIndex]->Fill(pr, ptI, w_resample);
		}
	      }
	      // No dR cut on the constituent-tagged case: the owning jet is the
	      // owning jet however far its axis ends up from the muon. It must
	      // still clear the same 20 GeV and |eta| cuts as the others, so a
	      // muon whose jet falls below threshold contributes to neither.
	      if(ptO > 0.){
		h_injMuonDR_donorJetPt_injectConstit[0]->Fill(drO, ptO, w_resample);
		h_injMuonDR_donorJetPt_injectConstit[CentralityIndex]->Fill(drO, ptO, w_resample);
		double pr = getPtRel(muPt_inj, muEta_inj, muPhi_inj, ptO, etaO, phiO);
		h_injMuonPtRel_donorJetPt_injectConstit[0]->Fill(pr, ptO, w_resample);
		h_injMuonPtRel_donorJetPt_injectConstit[CentralityIndex]->Fill(pr, ptO, w_resample);
	      }
	    }

	  }

	  // (real mu, fake jet) template: a genuine reco muon tagged to a
	  // MIXED-EVENT fastJet, i.e. a combinatorial jet that has nothing to do
	  // with this muon. findRecoMuonTag is reused verbatim so the muon
	  // selection (tight quality ID, pT window, |eta|<2, W-decay filter,
	  // dR < epsilon_mm) is byte-identical to the data tagging -- writing the
	  // cuts out by hand here would risk the template and the data drifting
	  // apart silently.
	  //
	  // matchFlagR_mixedJet is a SEPARATE array from the real-jet matching:
	  // sharing it would let a muon consumed here suppress a later real-jet
	  // tag in the same event (or vice versa) and bias both.
	  //
	  // Sized to em->nMu rather than the fixed [10] used by the data path.
	  // A fixed 10 would need a guard skipping events with more muons, and
	  // such an event would then contribute nothing here while still
	  // contributing to the data histogram this is subtracted from -- biasing
	  // the template low in exactly the busiest events. (The pre-existing
	  // matchFlagR[10] in the reco-jet loop has the same latent overflow and
	  // is deliberately left alone here.)
	  //
	  // evtTriggerDecision matches the guard on the data fill, so both are
	  // normalised over the same event sample.
	  if(evtTriggerDecision){
	    std::vector<int> matchFlagR_mixedJet(em->nMu > 0 ? em->nMu : 1, 0);
	    for(const auto& jet : jets){
		if(fabs(jet.eta()) > 1.6) continue;
		if(!h_RC_map[CentralityIndex]) continue;

		// trackMaxPt over the jet's own constituents, then the same
		// jetTrkMax filter the data, T3 and T4 jet loops all apply --
		// without it this template accepts jets the other three reject.
		// Applied against the JEC'd unsubtracted pT, matching how the
		// neighbouring fastJet loops call it.
		std::vector<fastjet::PseudoJet> constituents_m = jet.constituents();
		double trackMaxPt_m = 0.0;
		for(const auto& c : constituents_m){
		  if(!c.has_user_info<CandInfo>()) continue;
		  const CandInfo &candinfo_m = c.user_info<CandInfo>();
		  if(candinfo_m.isCharged() && c.pt() > trackMaxPt_m) trackMaxPt_m = c.pt();
		}
		JEC_PF.SetJetPT(jet.pt());
		JEC_PF.SetJetEta(jet.eta());
		JEC_PF.SetJetPhi(jet.phi_std());
		double jetPt_JEC_m = JEC_PF.GetCorrectedPT();
		if(doJetTrkMaxFilter){
		  if(!passesJetTrkMaxFilter(trackMaxPt_m,jetPt_JEC_m)) continue;
		}

		double rcMeanPt_m = h_RC_map[CentralityIndex]->GetBinContent(
				      h_RC_map[CentralityIndex]->FindBin(jet.eta(), jet.phi_std()));
		double jetPt_rcSub_m = jet.pt() - rcMeanPt_m;
		JEC_PF.SetJetPT(jetPt_rcSub_m);
		JEC_PF.SetJetEta(jet.eta());
		JEC_PF.SetJetPhi(jet.phi_std());
		double jetPt_JEC_rcSub_m = JEC_PF.GetCorrectedPT();
		if(jetPt_JEC_rcSub_m < 20.) continue;

		// JEC-corrected pT, NOT the raw rcSub one, for both the ptRel and
		// the Y axis. This template is subtracted from
		// h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn, whose jet pT is
		// JEC-corrected; filling on the raw scale would mean a "80-100 GeV"
		// window selected systematically harder jets here than in the data
		// it corrects, and ptRel itself would be built on a different jet
		// momentum. Same scale as every other term of the decomposition.
		double muPtRel_m = -999., muPt_m = -999., muEta_m = -999., muPhi_m = -999., muJetDr_m = -999.;
		if(findRecoMuonTag(em, jetPt_JEC_rcSub_m, jet.eta(), jet.phi_std(), matchFlagR_mixedJet.data(),
				   muPtRel_m, muPt_m, muEta_m, muPhi_m, muJetDr_m)){
		  h_realMuonPtRel_mixedFastJetPt[0]->Fill(muPtRel_m,jetPt_JEC_rcSub_m,w_resample);
		  h_realMuonPtRel_mixedFastJetPt[CentralityIndex]->Fill(muPtRel_m,jetPt_JEC_rcSub_m,w_resample);
		}
	    }
	  }

	}

	// FastJet anti-kT clustering on PFCs candidates
        std::vector<fastjet::PseudoJet> fjInputs_PFCs;
        
        
	for(int l = 0; l < em->nPFCspart; l++){
	  double pt  = em->pfCsPt->at(l);
	  double eta = em->pfCsEta->at(l);
	  double phi = em->pfCsPhi->at(l);
	  double id = em->pfCsId->at(l);
	  if(pt < pseudoJetCandPt_min) continue;
	  double px = pt * TMath::Cos(phi);
	  double py = pt * TMath::Sin(phi);
	  double pz = pt * TMath::SinH(eta);
	  double E  = pt * TMath::CosH(eta);
	  bool isCharged = (id == 1 || id == 2 || id == 3); // id = h, e, mu,
	  fastjet::PseudoJet pseudoJet_l(px, py, pz, E);
	  pseudoJet_l.set_user_index(isCharged ? 1 : 0);
	  fjInputs_PFCs.push_back(pseudoJet_l);
	}
        
	fastjet::ClusterSequence cs_PFCs(fjInputs_PFCs, jetDef);
        std::vector<fastjet::PseudoJet> jets_PFCs = fastjet::sorted_by_pt(cs_PFCs.inclusive_jets(0.));
        for(const auto& jet_PFCs : jets_PFCs){
          if(TMath::Abs(jet_PFCs.eta()) > 1.6) continue;
	  std::vector<fastjet::PseudoJet> constituents_PFCs = jet_PFCs.constituents();
	  if(skipSingleConstituentJets){
	    if (constituents_PFCs.size() < 2) continue;  // cut single-track/single-constituent jets
	  }
	  double trackMaxPt = 0.0;
	  for(const auto& c : constituents_PFCs){
	    if(c.user_index() == 1 && c.pt() > trackMaxPt) trackMaxPt = c.pt();
	  }
	  
          JEC_PF.SetJetPT(jet_PFCs.pt());
          JEC_PF.SetJetEta(jet_PFCs.eta());
          JEC_PF.SetJetPhi(jet_PFCs.phi_std());
          double fastJetPt_PFCs_JEC = JEC_PF.GetCorrectedPT();
	  int jetPtIndex_PFCs = -1;
	  if(fastJetPt_PFCs_JEC < 20.) continue;
	  if(doJetTrkMaxFilter){
	    // do track-max pT cut
	    if(!passesJetTrkMaxFilter(trackMaxPt,fastJetPt_PFCs_JEC)) continue;
	  }
	  jetPtIndex_PFCs = getJetPtBin(fastJetPt_PFCs_JEC);
	  
	  h_fastJetPt_PFCs[0]->Fill(jet_PFCs.pt(), w_resample);
          h_fastJetPt_PFCs[CentralityIndex]->Fill(jet_PFCs.pt(), w_resample);
          h_fastJetPt_PFCs_JEC[0]->Fill(fastJetPt_PFCs_JEC, w_resample);
          h_fastJetPt_PFCs_JEC[CentralityIndex]->Fill(fastJetPt_PFCs_JEC, w_resample);

	  double dR_min = 999.0;
	  double dPT = -999.0;
	  double matchedPtPF = -999.0;
	  // match h_fastJetPt_PFCs to h_fastJetPt_PF
	  for(const auto& jet_PF : jets){
	    std::vector<fastjet::PseudoJet> constituents_PF = jet_PF.constituents();
	    if(skipSingleConstituentJets){
	      if (constituents_PF.size() < 2) continue;  // cut single-track/single-constituent jets
	    }
	    double dR = getDr(jet_PF.eta(),jet_PF.phi_std(),jet_PFCs.eta(),jet_PFCs.phi_std());
	    if(dR < dR_min){
	      dR_min = dR;
	      dPT = jet_PF.pt() - jet_PFCs.pt();
	      matchedPtPF = jet_PF.pt();
	    }
	  }
	  h_dRmin_PF_PFCs[0]->Fill(dR_min,w_resample);
	  h_dRmin_PF_PFCs[CentralityIndex]->Fill(dR_min,w_resample);
	  h_dPTEtaPhi_PF_PFCs[0]->Fill(jet_PFCs.eta(),jet_PFCs.phi_std(),dPT);
	  h_dPTEtaPhi_PF_PFCs[CentralityIndex]->Fill(jet_PFCs.eta(),jet_PFCs.phi_std(),dPT);
	  double dPT_geoCorr = dPT * TMath::Cos(dR_min);
	  h_dPTEtaPhi_PF_PFCs_geoCorr[0]->Fill(jet_PFCs.eta(),jet_PFCs.phi_std(),dPT_geoCorr);
	  h_dPTEtaPhi_PF_PFCs_geoCorr[CentralityIndex]->Fill(jet_PFCs.eta(),jet_PFCs.phi_std(),dPT_geoCorr);
	  if(dPT > 0){
	    h_dPTEtaPhi_PF_PFCs_dPTAbove0[0]->Fill(jet_PFCs.eta(),jet_PFCs.phi_std(),dPT);
	    h_dPTEtaPhi_PF_PFCs_dPTAbove0[CentralityIndex]->Fill(jet_PFCs.eta(),jet_PFCs.phi_std(),dPT);
	  }
	  if(jet_PFCs.pt() > 60.){
	    h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60[0]->Fill(jet_PFCs.eta(),jet_PFCs.phi_std(),dPT);
	    h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60[CentralityIndex]->Fill(jet_PFCs.eta(),jet_PFCs.phi_std(),dPT);
	  }

	  h_fastJetPtPF_dRmin[0][0]->Fill(matchedPtPF,dR_min,w_resample);
	  h_fastJetPtPF_dRmin[CentralityIndex][0]->Fill(matchedPtPF,dR_min,w_resample);
	  if(jetPtIndex_PFCs > 0){
	    h_fastJetPtPF_dRmin[0][jetPtIndex_PFCs]->Fill(matchedPtPF,dR_min,w_resample);
	    h_fastJetPtPF_dRmin[CentralityIndex][jetPtIndex_PFCs]->Fill(matchedPtPF,dR_min,w_resample);
	  }

	  h_fastJetPtPF_etaPFCs[0][0]->Fill(matchedPtPF,jet_PFCs.eta(),w_resample);
	  h_fastJetPtPF_etaPFCs[CentralityIndex][0]->Fill(matchedPtPF,jet_PFCs.eta(),w_resample);
	  if(jetPtIndex_PFCs > 0){
	    h_fastJetPtPF_etaPFCs[0][jetPtIndex_PFCs]->Fill(matchedPtPF,jet_PFCs.eta(),w_resample);
	    h_fastJetPtPF_etaPFCs[CentralityIndex][jetPtIndex_PFCs]->Fill(matchedPtPF,jet_PFCs.eta(),w_resample);
	  }

        }

	} // end for(fjResample ...)









	
      }
#endif

      bool eventHasGoodJet = false;
      bool eventHasInclRecoMuonTag = false;
      bool eventHasInclRecoMuonTagPlusTrigger = false;
      bool eventHasMatchedRecoMuonTag = false;
      bool eventHasMatchedRecoMuonTagPlusTrigger = false;
      int inclJetCounter = 0;
      int muTaggedJetCounter = 0;

      //cout << "nMu = " << em->nMu << endl;
    

      double w_trig = w;

      // if(applyMu12TriggerEfficiencyCorrection){
      //   if(CentralityIndex == 4) w_trig = w / fitFxn_PbPb_HLT_C4->Eval(leadingMuonPt);
      //   else if(CentralityIndex == 3) w_trig = w / fitFxn_PbPb_HLT_C3->Eval(leadingMuonPt);
      //   else if(CentralityIndex == 2) w_trig = w / fitFxn_PbPb_HLT_C2->Eval(leadingMuonPt);
      //   else if(CentralityIndex == 1) w_trig = w / fitFxn_PbPb_HLT_C1->Eval(leadingMuonPt);
      //   else{};
      // }


      double leadingRecoJetPt = 0.0;
      // RECO JET LOOP
      for(int i = 0; i < em->njet ; i++){

	if(onlyOneMuonTaggedJetPerEvent){
	  if(eventHasInclRecoMuonTag) continue;
	}
      
	// JET VARIABLES

	double x = 0.;
	if(useCaloJetsOverride){
	  JEC_Calo.SetJetPT(em->rawpt[i]);
	  JEC_Calo.SetJetEta(em->jeteta[i]);
	  JEC_Calo.SetJetPhi(em->jetphi[i]);
	  x = JEC_Calo.GetCorrectedPT();
	}
	else{
	  JEC_PF.SetJetPT(em->rawpt[i]);
	  JEC_PF.SetJetEta(em->jeteta[i]);
	  JEC_PF.SetJetPhi(em->jetphi[i]);
	  x = JEC_PF.GetCorrectedPT();
	}
	
	double rawJetPt_i = em->rawpt[i];
	//double x = em->jetpt[i]; // use built-in JEC
	double y = em->jeteta[i]; // recoJetEta
	double z = em->jetphi[i]; // recoJetPhi
	double jetTrkMax_i = em->jetTrkMax[i]; 


	// double sigma = 0.663*JER_fxn->Eval(x);
	// double mu = 1.0;

	// double smear = randomGenerator->Gaus(mu,sigma);
	//cout << "smear = " << smear << endl;
      
	//x = x*smear;
	if(doJetTrkMaxFilter){
	  if(!passesJetTrkMaxFilter(jetTrkMax_i,x)) continue;
	}      


	//cout << "rawPt = " << em->rawpt[i] << "  |  jetPt = " << em->jetpt[i] << "  |  corrPt = " << x << endl;
	if(doEtaPhiMask){
	  if(etaPhiMask(y,z)) continue;
	}

	// JEU must follow the jet collection, exactly as JEC does above. Both
	// branches previously passed JEU_Calo, so a PF jet was assigned the
	// AK4Calo uncertainty and JEU_PF was declared but never used anywhere.
	// Harmless while apply_JEU_shift_up/down are both false -- applyJEU_JER
	// only reads the JEU inside those branches -- but wrong the moment a JEU
	// systematic is switched on.
	//
	// The else branch covers akCs4PF and, via useFlowJetsOverride,
	// akFlowPuCs4PF; both are PF collections, so JEU_PF is right for each.
	if(useCaloJetsOverride){
	  x = applyJEU_JER(x, JEU_Calo, JER_fxn, randomGenerator,
			   fitFxn_PYTHIA_JERCorrection,
			   neutrino_tag_fraction, neutrino_energy_map);
	}
	else{
	  x = applyJEU_JER(x, JEU_PF, JER_fxn, randomGenerator,
			   fitFxn_PYTHIA_JERCorrection,
			   neutrino_tag_fraction, neutrino_energy_map);
	}
	if(x < 0) continue;

	//cout << "Event " << evi << ", jet " << i << endl;
	//cout << "~~~~  jetPt = " << em->jetpt[i] << ", corrJetPt = " << x << endl;

	double muPtRel = -1.0;
	double muPt = -1.0;
	double muEta = -1.0;
	double muPhi = -1.0;
	double muJetDr = -1.0;

	// jet kinematic cuts
	if(TMath::Abs(y) > etaMax || x < jetPtCut) continue;

	if(x > leadingRecoJetPt) leadingRecoJetPt = x;

	eventHasGoodJet = true;

	int jetPtIndex = getJetPtBin(x);

	//cout << "nMu = " << em->nMu << endl;

	// look for recoMuon match to recoJet
	bool hasInclRecoMuonTag = findRecoMuonTag(em, x, y, z, matchFlagR,
	                                          muPtRel, muPt, muEta, muPhi, muJetDr);

	if(applyMu12TriggerEfficiencyCorrection){
	  if(CentralityIndex == 4) w_trig = w / fitFxn_PbPb_HLT_C4->Eval(muPt);
	  else if(CentralityIndex == 3) w_trig = w / fitFxn_PbPb_HLT_C3->Eval(muPt);
	  else if(CentralityIndex == 2) w_trig = w / fitFxn_PbPb_HLT_C2->Eval(muPt);
	  else if(CentralityIndex == 1) w_trig = w / fitFxn_PbPb_HLT_C1->Eval(muPt);
	  else{};
	}

	// Fill the jet/event histograms

	inclJetCounter++;
	h_inclRecoJetPt[0]->Fill(x,w);
	h_inclRecoJetPt[CentralityIndex]->Fill(x,w);

	h_inclRawJetPt[0]->Fill(rawJetPt_i,w);
	h_inclRawJetPt[CentralityIndex]->Fill(rawJetPt_i,w);

	h_inclRecoJetEta[0]->Fill(y,w);
	h_inclRecoJetEta[CentralityIndex]->Fill(y,w);

	h_inclRecoJetPhi[0]->Fill(z,w);
	h_inclRecoJetPhi[CentralityIndex]->Fill(z,w);

	h_inclRecoJetPt_inclRecoJetEta[0]->Fill(x,y,w);
	h_inclRecoJetPt_inclRecoJetEta[CentralityIndex]->Fill(x,y,w);

	h_inclRecoJetPt_inclRecoJetPhi[0]->Fill(x,z,w);
	h_inclRecoJetPt_inclRecoJetPhi[CentralityIndex]->Fill(x,z,w);

	h_inclRecoJetEta_inclRecoJetPhi[0][0]->Fill(y,z,w);
	h_inclRecoJetEta_inclRecoJetPhi[CentralityIndex][0]->Fill(y,z,w);
	if(jetPtIndex > 0){
	  h_inclRecoJetEta_inclRecoJetPhi[0][jetPtIndex]->Fill(y,z,w);
	  h_inclRecoJetEta_inclRecoJetPhi[CentralityIndex][jetPtIndex]->Fill(y,z,w);
	}
     
	if(hasInclRecoMuonTag){

	  eventHasInclRecoMuonTag = true;
	  muTaggedJetCounter++;
	
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[0][0]->Fill(y,z,w);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[CentralityIndex][0]->Fill(y,z,w);
	  if(jetPtIndex > 0){
	    h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[0][jetPtIndex]->Fill(y,z,w);
	    h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[CentralityIndex][jetPtIndex]->Fill(y,z,w);
	  }

	  h_inclRecoJetPt_inclRecoMuonTag[0]->Fill(x,w);
	  h_inclRecoJetPt_inclRecoMuonTag[CentralityIndex]->Fill(x,w);

	  h_inclRecoJetEta_inclRecoMuonTag[0]->Fill(y,w);
	  h_inclRecoJetEta_inclRecoMuonTag[CentralityIndex]->Fill(y,w);

	  h_inclRecoJetPhi_inclRecoMuonTag[0]->Fill(z,w);
	  h_inclRecoJetPhi_inclRecoMuonTag[CentralityIndex]->Fill(z,w);

	  h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag[0]->Fill(x,y,w);
	  h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag[CentralityIndex]->Fill(x,y,w);

	  h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag[0]->Fill(x,z,w);
	  h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag[CentralityIndex]->Fill(x,z,w);
       
	  if(evtTriggerDecision){
	 
	    eventHasInclRecoMuonTagPlusTrigger = true;

	    h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[0][0]->Fill(y,z,w_trig);
	    h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[CentralityIndex][0]->Fill(y,z,w_trig);
	    if(jetPtIndex > 0){
	      h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[0][jetPtIndex]->Fill(y,z,w_trig);
	      h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[CentralityIndex][jetPtIndex]->Fill(y,z,w_trig);
	    }

	    h_inclRecoJetPt_inclRecoMuonTag_triggerOn[0]->Fill(x,w_trig);
	    h_inclRecoJetPt_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(x,w_trig);

	    h_inclRecoJetEta_inclRecoMuonTag_triggerOn[0]->Fill(y,w_trig);
	    h_inclRecoJetEta_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(y,w_trig);

	    h_inclRecoJetPhi_inclRecoMuonTag_triggerOn[0]->Fill(z,w_trig);
	    h_inclRecoJetPhi_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(z,w_trig);

	    h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn[0]->Fill(x,y,w_trig);
	    h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(x,y,w_trig);

	    h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn[0]->Fill(x,z,w_trig);
	    h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(x,z,w_trig);

	    h_mupt_jetpt[0]->Fill(muPt,x,w_trig);
	    h_mupt_jetpt[CentralityIndex]->Fill(muPt,x,w_trig);

	    h_muptrel_jetpt[0]->Fill(muPtRel,x,w_trig);
	    h_muptrel_jetpt[CentralityIndex]->Fill(muPtRel,x,w_trig);

	  
	    h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn[0]->Fill(muPtRel,x,w_trig);
	    h_mupt_recoJetPt_inclRecoMuonTag_triggerOn[0]->Fill(muPt,x,w_trig);
	    h_mueta_recoJetPt_inclRecoMuonTag_triggerOn[0]->Fill(muEta,x,w_trig);
	    h_muphi_recoJetPt_inclRecoMuonTag_triggerOn[0]->Fill(muPhi,x,w_trig);
	    h_muJetDr_recoJetPt[0]->Fill(muJetDr,x,w_trig);

	    h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(muPtRel,x,w_trig);
	    h_mupt_recoJetPt_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(muPt,x,w_trig);
	    h_mueta_recoJetPt_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(muEta,x,w_trig);
	    h_muphi_recoJetPt_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(muPhi,x,w_trig);
	    h_muJetDr_recoJetPt[CentralityIndex]->Fill(muJetDr,x,w_trig);

	    h_muPtOverJetPt[0][0]->Fill(muPt/x,w_trig);
	    h_mupt_muptrel[0][0]->Fill(muPt,muPtRel,w_trig);
	 
	    if(jetPtIndex > 0){
	    
	      h_muPtOverJetPt[0][jetPtIndex]->Fill(muPt/x,w_trig);
	      h_muPtOverJetPt[CentralityIndex][0]->Fill(muPt/x,w_trig);
	      h_muPtOverJetPt[CentralityIndex][jetPtIndex]->Fill(muPt/x,w_trig);

	      h_mupt_muptrel[0][jetPtIndex]->Fill(muPt,muPtRel,w_trig);
	      h_mupt_muptrel[CentralityIndex][0]->Fill(muPt,muPtRel,w_trig);
	      h_mupt_muptrel[CentralityIndex][jetPtIndex]->Fill(muPt,muPtRel,w_trig);

	      h_muptrel_hiBin[0]->Fill(muPtRel,em->hiBin,w_trig);
	      h_muptrel_hiBin[jetPtIndex]->Fill(muPtRel,em->hiBin,w_trig);
	   
	    }
	  } 

	}

      }
      // END recoJet LOOP

      h_NJetPerEvent[0]->Fill(inclJetCounter);
      h_NJetPerEvent[CentralityIndex]->Fill(inclJetCounter);

      h_NMuTaggedJetPerEvent[0]->Fill(muTaggedJetCounter);
      h_NMuTaggedJetPerEvent[CentralityIndex]->Fill(muTaggedJetCounter);

    
    
      if(eventHasGoodJet && leadingRecoJetPt > 80){

	if(em->HLT_HICsAK4PFJet100Eta1p5_v1 == 1 && leadingRecoJetPt > 130){
	  h_vz_jet100_withJetAbove130[0]->Fill(em->vz,w);
	  h_vz_jet100_withJetAbove130[CentralityIndex]->Fill(em->vz,w);
	  h_hiBin_jet100_withJetAbove130->Fill(em->hiBin,w);
	}

	h_vz_jet[0]->Fill(em->vz,w);
	h_vz_jet[CentralityIndex]->Fill(em->vz,w);
	
	h_hiBin_jet->Fill(em->hiBin,w);

	if(eventHasInclRecoMuonTag){

	  h_vz_inclRecoMuonTag[0]->Fill(em->vz,w);
	  h_vz_inclRecoMuonTag[CentralityIndex]->Fill(em->vz,w);

	  h_hiBin_inclRecoMuonTag->Fill(em->hiBin,w);

	  if(eventHasInclRecoMuonTagPlusTrigger){

	    h_vz_inclRecoMuonTag_triggerOn[0]->Fill(em->vz,w_trig);
	    h_vz_inclRecoMuonTag_triggerOn[CentralityIndex]->Fill(em->vz,w_trig);

	    h_hiBin_inclRecoMuonTag_triggerOn->Fill(em->hiBin,w_trig);
	 
	  }
       

	}

      }


      // RECO MUON LOOP

      int loopMuonTrigger = 0;
      int loopJetTrigger = 0;
    
      double loopMuonPtCut = 0.0;
      double loopJetPtCut = 0.0;

      if(applyJet60Trigger){

	loopJetTrigger = em->HLT_HICsAK4PFJet60Eta1p5_v1;
	loopJetPtCut = 100.0;

      }
      else if(applyJet80Trigger){

	loopJetTrigger = em->HLT_HICsAK4PFJet80Eta1p5_v1;
	loopJetPtCut = 130.0;

      }
      else if(applyJet100Trigger){

	loopJetTrigger = em->HLT_HICsAK4PFJet100Eta1p5_v1;
	loopJetPtCut = 150.0;

      }
      else{
	loopJetTrigger = 1;
	loopJetPtCut = 0.;
      }
    
      if(fillMu5) loopMuonTrigger = triggerDecision_mu5;
      else if(fillMu7) loopMuonTrigger = triggerDecision_mu7;
      else if(fillMu12) loopMuonTrigger = triggerDecision_mu12;
      else{};

      double leadingMuonPt = 0.0;
      double etaCut_Zloop = 2.4;
    
      //if(triggerIsOn(loopMuonTrigger,1) && triggerIsOn(loopJetTrigger,1) && leadingRecoJetPt > loopJetPtCut){
      if(triggerIsOn(loopMuonTrigger,1)){
	for(int m = 0; m < em->nMu; m++){

	  double muPt_m = em->muPt->at(m);
	  double muEta_m = em->muEta->at(m);
	  double muPhi_m = em->muPhi->at(m);

	  //cout << "(muPt, muEta, muPhi) = (" << muPt_m << ", " << muEta_m << ", " << muPhi_m << ")" << endl;

	  // skip if muon has already been matched to a jet in this event
	  // muon kinematic cuts
	  if(muPt_m < muPtCut || muPt_m > muPtMaxCut || fabs(muEta_m) > etaCut_Zloop) continue;
	  // muon quality cuts
	  if(fillMu12){
	    if(!isQualityMuon_tight(em->muChi2NDF->at(m),
				    em->muInnerD0->at(m),
				    em->muInnerDz->at(m),
				    em->muMuonHits->at(m),
				    em->muPixelHits->at(m),
				    em->muIsGlobal->at(m),
				    em->muIsPF->at(m),
				    em->muStations->at(m),
				    em->muTrkLayers->at(m))) continue; // skip if muon doesnt pass quality cuts
	  }

	  else if(fillMu5 || fillMu7){
	    if(!isQualityMuon_hybridSoft(em->muChi2NDF->at(m),
					 em->muInnerD0->at(m),
					 em->muInnerDz->at(m),
					 em->muPixelHits->at(m),
					 em->muIsTracker->at(m),
					 em->muIsGlobal->at(m),
					 em->muTrkLayers->at(m))) continue; // skip if muon doesnt pass quality cuts
	  }
	  else{};

	  if(muPt_m > leadingMuonPt) leadingMuonPt = muPt_m;

	  h_inclMuPt->Fill(muPt_m,w);

	  for(int k = m+1; k < em->nMu; k++){
	    //for(int k = 0; k < em->nMu; k++){ // double count error (for debugging)

	    double muPt_k = em->muPt->at(k);
	    double muEta_k = em->muEta->at(k);
	    double muPhi_k = em->muPhi->at(k);

	    if(muPt_k < muPtCut || muPt_k > muPtMaxCut || fabs(muEta_k) > etaCut_Zloop) continue;

	    if(fillMu12){
	      if(!isQualityMuon_tight(em->muChi2NDF->at(k),
				      em->muInnerD0->at(k),
				      em->muInnerDz->at(k),
				      em->muMuonHits->at(k),
				      em->muPixelHits->at(k),
				      em->muIsGlobal->at(k),
				      em->muIsPF->at(k),
				      em->muStations->at(k),
				      em->muTrkLayers->at(k))) continue; // skip if muon doesnt pass quality cuts
	    }

	    else if(fillMu5 || fillMu7){
	      if(!isQualityMuon_hybridSoft(em->muChi2NDF->at(k),
					   em->muInnerD0->at(k),
					   em->muInnerDz->at(k),
					   em->muPixelHits->at(k),
					   em->muIsTracker->at(k),
					   em->muIsGlobal->at(k),
					   em->muTrkLayers->at(k))) continue; // skip if muon doesnt pass quality cuts
	    }
	    else{};

	    double w_mk = w;

	    if(CentralityIndex == 4) w_mk = w / (1. - (1. - fitFxn_PbPb_HLT_C4->Eval(muPt_m))*(1. - fitFxn_PbPb_HLT_C4->Eval(muPt_k)));
	    else if(CentralityIndex == 3) w_mk = w / (1. - (1. - fitFxn_PbPb_HLT_C3->Eval(muPt_m))*(1. - fitFxn_PbPb_HLT_C3->Eval(muPt_k)));
	    else if(CentralityIndex == 2) w_mk = w / (1. - (1. - fitFxn_PbPb_HLT_C2->Eval(muPt_m))*(1. - fitFxn_PbPb_HLT_C2->Eval(muPt_k)));
	    else if(CentralityIndex == 1) w_mk = w / (1. - (1. - fitFxn_PbPb_HLT_C1->Eval(muPt_m))*(1. - fitFxn_PbPb_HLT_C1->Eval(muPt_k)));
	    else{};

	    if(em->muCharge->at(m)*em->muCharge->at(k) == -1){

	      //h_dimuonMass[0]->Fill(calculateDimuonMass(muPt_m,muEta_m,muPhi_m,muPt_k,muEta_k,muPhi_k,em->runNumber,em->lumiSection,em->evtNumber),w_mk);
	      h_dimuonMass[0]->Fill(calculateDimuonMass(muPt_m,muEta_m,muPhi_m,muPt_k,muEta_k,muPhi_k),w_mk);
	      h_dimuonMass[CentralityIndex]->Fill(calculateDimuonMass(muPt_m,muEta_m,muPhi_m,muPt_k,muEta_k,muPhi_k),w_mk);
	  
	    }

	    else if(em->muCharge->at(m)*em->muCharge->at(k) == 1){

	      h_dimuonMass_sameSign[0]->Fill(calculateDimuonMass(muPt_m,muEta_m,muPhi_m,muPt_k,muEta_k,muPhi_k),w_mk);
	      h_dimuonMass_sameSign[CentralityIndex]->Fill(calculateDimuonMass(muPt_m,muEta_m,muPhi_m,muPt_k,muEta_k,muPhi_k),w_mk);
	  
	    }
	
	  }
	}
      }


  

    

    } // end event loop

    h_NEvents->Fill(eventCounter);

 
    delete f;
    // WRITE
    auto wf = TFile::Open(output,"recreate");

    h_eventsBeforeSelection->Write();
    h_eventsAfterSelection->Write();
    h_NEvents->Write();
    h_hiBin->Write();
    h_hiBin_triggerOn->Write();
    h_hiBin_jet60->Write();
    h_hiBin_jet80->Write();
    h_hiBin_jet100->Write();
    h_hiBin_jet->Write();
    h_hiBin_inclRecoMuonTag->Write();
    h_hiBin_inclRecoMuonTag_triggerOn->Write();
    h_inclMuPt->Write();
    h_nMu_triggerOn->Write();

    for(int i = 0; i < NCentralityIndices; i++){

      h_NJetPerEvent[i]->Write();
      h_NMuTaggedJetPerEvent[i]->Write();
    
      h_vz[i]->Write();
      h_vz_triggerOn[i]->Write();
      h_vz_jet60[i]->Write();
      h_vz_jet80[i]->Write();
      h_vz_jet100[i]->Write();
      h_vz_jet100_withJetAbove130[i]->Write();
      h_vz_jet[i]->Write();
      h_vz_inclRecoMuonTag[i]->Write();
      h_vz_inclRecoMuonTag_triggerOn[i]->Write();
      h_dimuonMass[i]->Write();
      h_dimuonMass_sameSign[i]->Write();
   
      h_inclRecoJetPt[i]->Write();
      h_inclRawJetPt[i]->Write();
      h_inclRecoJetEta[i]->Write();
      h_inclRecoJetPhi[i]->Write();
      h_inclRecoJetPt_inclRecoJetEta[i]->Write();
      h_inclRecoJetPt_inclRecoJetPhi[i]->Write();
   
      h_inclRecoJetPt_inclRecoMuonTag[i]->Write();
      h_inclRecoJetEta_inclRecoMuonTag[i]->Write();
      h_inclRecoJetPhi_inclRecoMuonTag[i]->Write();
      h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag[i]->Write();
      h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag[i]->Write();
   
      h_inclRecoJetPt_inclRecoMuonTag_triggerOn[i]->Write();
      h_inclRecoJetEta_inclRecoMuonTag_triggerOn[i]->Write();
      h_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i]->Write();
      h_inclRecoJetPt_inclRecoJetEta_inclRecoMuonTag_triggerOn[i]->Write();
      h_inclRecoJetPt_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i]->Write();

      h_mupt_jetpt[i]->Write();
      h_muptrel_jetpt[i]->Write();

      h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn[i]->Write();
      h_mupt_recoJetPt_inclRecoMuonTag_triggerOn[i]->Write();
      h_mueta_recoJetPt_inclRecoMuonTag_triggerOn[i]->Write();
      h_muphi_recoJetPt_inclRecoMuonTag_triggerOn[i]->Write();
      h_muJetDr_recoJetPt[i]->Write();
      h_Jet60HLT[i]->Write();
      h_Jet60HLT_Prescale[i]->Write();
      h_Jet80HLT[i]->Write();
      h_Jet80HLT_Prescale[i]->Write();
      h_Jet100HLT[i]->Write();
      h_Jet100HLT_Prescale[i]->Write();

      h_pfPt[i]->Write();
      h_pseudoJetPt[i]->Write();
      h_pseudoJetPt_geoCorr[i]->Write();
      h_fastJetPt_PF[i]->Write();
      h_fastJetPt_PF_JEC[i]->Write();
      h_fastJetPt_PFCs[i]->Write();
      h_fastJetPt_PFCs_JEC[i]->Write();
      h_fastJetPt_PF_bkgSub_RC[i]->Write();
      h_fastJetPt_PF_JEC_bkgSub_RC[i]->Write();
      h_fastJetPt_PF_bkgSub_RC_sigSel[i]->Write();
      h_fastJetPt_PF_JEC_bkgSub_RC_sigSel[i]->Write();
      h_fastJetPt_PF_bkgSub_RC_geoCorr[i]->Write();
      h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr[i]->Write();
      h_fastJetPt_PF_bkgSub_RC_geoCorr_etaReflect[i]->Write();
      h_fastJetPt_PF_JEC_bkgSub_RC_geoCorr_etaReflect[i]->Write();
      h_fastJetPt_PF_bkgSub_dPT[i]->Write();
      h_fastJetPt_PF_JEC_bkgSub_dPT[i]->Write();
      h_fastJetPt_PF_bkgSub_dPT_geoCorr[i]->Write();
      h_fastJetPt_PF_JEC_bkgSub_dPT_geoCorr[i]->Write();
      h_fastJetPt_PF_bkgSub_dPT_dPTAbove0[i]->Write();
      h_fastJetPt_PF_JEC_bkgSub_dPT_dPTAbove0[i]->Write();
      h_fastJetPt_PF_bkgSub_dPT_PFCsPTAbove60[i]->Write();
      h_fastJetPt_PF_JEC_bkgSub_dPT_PFCsPTAbove60[i]->Write();

      h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC[i]->Write();
      h_fastJetMuonDR_fastJetPt_PF_bkgSub_RC[i]->Write();
      h_fastJetMuonDR_inclusiveClosestFastJet[i]->Write();
      h_muonDR_inclusiveClosestJet[i]->Write();
      h_muonDR_inclusiveClosestJet_triggerOn[i]->Write();
      h_donorJetPt_noInject[i]->Write();
      h_donorJetPt_inject[i]->Write();
      h_injMuonDR_donorJetPt_noInject[i]->Write();
      h_injMuonPtRel_donorJetPt_noInject[i]->Write();
      h_injMuonDR_donorJetPt_inject[i]->Write();
      h_injMuonPtRel_donorJetPt_inject[i]->Write();
      h_injMuonDR_donorJetPt_injectConstit[i]->Write();
      h_injMuonPtRel_donorJetPt_injectConstit[i]->Write();
      h_mixedMuonPtRel_recoJetPt[i]->Write();
      h_realMuonPtRel_mixedFastJetPt[i]->Write();
      
      h_nPFcand[i]->Write();
      h_nPFcandCS[i]->Write();
      h_nPFcandFastJet[i]->Write();
      h_randConeEtaPhi[i]->Write();
      h_randConeEtaPhi_geoCorr[i]->Write();
      h_pseudoJetPt_sigSel[i]->Write();
      h_randConeEtaPhi_sigSel[i]->Write();
      h_dPTEtaPhi_PF_PFCs[i]->Write();
      h_dPTEtaPhi_PF_PFCs_geoCorr[i]->Write();
      h_dPTEtaPhi_PF_PFCs_dPTAbove0[i]->Write();
      h_dPTEtaPhi_PF_PFCs_PFCsPTAbove60[i]->Write();
      h_dRmin_PF_PFCs[i]->Write();

      for(int j = 0; j < NJetPtIndices; j++){

	h_inclRecoJetEta_inclRecoJetPhi[i][j]->Write();
	h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j]->Write();
	h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j]->Write();
	h_muPtOverJetPt[i][j]->Write();
	h_mupt_muptrel[i][j]->Write();
	h_fastJetPtPF_dRmin[i][j]->Write();
	h_fastJetPtPF_etaPFCs[i][j]->Write();
     
      }
    }

    for(int j = 0; j < NJetPtIndices; j++){
      h_muptrel_hiBin[j]->Write();
    }

    wf->Close();
    return;
    // END WRITE

  }

}
