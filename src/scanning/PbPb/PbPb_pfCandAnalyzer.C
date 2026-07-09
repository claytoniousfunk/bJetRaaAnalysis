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
#include "TSystem.h"
#include "TRandom2.h"

// FastJet
#ifdef DO_FASTJET
#include "fastjet/ClusterSequence.hh"
#endif

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
TH1D     *h_fastJetPt[NCentralityIndices];
TH1D     *h_fastJetPt_JEC[NCentralityIndices];
TH1D     *h_nPFcand[NCentralityIndices];
TH1D     *h_nPFcandCS[NCentralityIndices];
// Random-cone eta/phi map: mean cone pT at each (eta,phi) cell.
// Cell size ~0.1x0.1 — sub-cone resolution, enough stats per cell.
// Used to build a local UE background map for FastJet jet subtraction.
const int    NRC_EtaBins = 32;   // -1.6 to 1.6, width 0.1
const int    NRC_PhiBins = 64;   // -pi to pi,   width ~0.098
TProfile2D  *h_randConeEtaPhi[NCentralityIndices];
TH1D        *h_fastJetPt_bkgSub_RC[NCentralityIndices];
// RC eta/phi maps loaded from external file at run time
TProfile2D  *h_RC_map[NCentralityIndices];

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
  if(doSingleMuonSample) inputFileList = "../../../fileNames/fileNames_HISingleMuon_HIRun2018A-04Apr2019-v1.txt";
  else if(doMinBiasSample) inputFileList = "../../../fileNames/fileNames_HIMinimumBias0_Part1_withTracksAndPFCandidates.txt";
  else if(doHardProbesSample) inputFileList = "../../../fileNames/fileNames_HIHardProbes_HIRun2018A-04Apr2019-v1_103X_dataRun2_Prompt_fixEcalADCToGeV_v1.txt";
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
						   pseudoJetCandPt_min,
						   doEventMixing,
						   doConstituentSubtraction);


    TString suffixEdit = CENT_SCHEME_SUFFIX;

    TString output = Form("%s%s%s/PbPb_pfCandAnalyzer_output_%i.root",outputBaseDir.Data(),outputDatasetName.Data(),suffixEdit.Data(),group);

    std::cout << "output dataset = " << output << std::endl;

    if(gSystem->AccessPathName(Form("%s%s%s",outputBaseDir.Data(),outputDatasetName.Data(),suffixEdit.Data()))){
      std::cout << "\033[1;31m Output directory not found: \033[0m " << Form("%s%s%s",outputBaseDir.Data(),outputDatasetName.Data(),suffixEdit.Data()) << std::endl;
      return;
    }
  

  
    // JET ENERGY CORRECTIONS
    vector<string> Files;
    Files.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_DATA_L2Relative_AK4PF.txt"); // L2Relative correction
    Files.push_back("../../../JetEnergyCorrections/Autumn18_HI_V8_DATA_L2L3Residual_AK4PF.txt"); // L2L3Residual correction
    JetCorrector JEC(Files);
    JetUncertainty JEU("../../../JetEnergyCorrections/Autumn18_HI_V8_MC_Uncertainty_AK4PF.txt");
    /// print out some info
    readConfig();


    // Load random-cone UE maps (used for RC-subtracted fastJet pT)
    TFile *f_RC_maps = TFile::Open("/eos/cms/store/group/phys_heavyions/cbennett/maps/PbPb_MinBias_Part1_randConeEtaPhi_2026-7-8_ultraFineCentBins.root");
    for(int i = 0; i < NCentralityIndices; i++){
      f_RC_maps->GetObject(Form("h_randConeEtaPhi_C%d", i), h_RC_map[i]);
      if(h_RC_map[i]) h_RC_map[i]->SetDirectory(nullptr);
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
	h_fastJetPt[i] = new TH1D(Form("h_fastJetPt_C%i",i),Form("FastJet anti-kT pT, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_JEC[i] = new TH1D(Form("h_fastJetPt_JEC_C%i",i),Form("FastJet anti-kT pT (JEC), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_fastJetPt_bkgSub_RC[i] = new TH1D(Form("h_fastJetPt_bkgSub_RC_C%i",i),Form("FastJet anti-kT pT (bkg sub, random-cone), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NPtBins,ptMin,ptMax);
	h_nPFcand[i]   = new TH1D(Form("h_nPFcand_C%i",i),  Form("N PF cands per event, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),10000,0,10000);
	h_nPFcandCS[i] = new TH1D(Form("h_nPFcandCS_C%i",i), Form("N PFCS cands per event, hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),10000,0,10000);
	h_randConeEtaPhi[i] = new TProfile2D(Form("h_randConeEtaPhi_C%i",i),Form("Mean random-cone p_{T} vs (#eta,#phi), hiBin %i - %i",centEdges[0],centEdges[NCentralityIndices-1]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
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
	h_fastJetPt[i] = new TH1D(Form("h_fastJetPt_C%i",i),Form("FastJet anti-kT pT, hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_JEC[i] = new TH1D(Form("h_fastJetPt_JEC_C%i",i),Form("FastJet anti-kT pT (JEC), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_fastJetPt_bkgSub_RC[i] = new TH1D(Form("h_fastJetPt_bkgSub_RC_C%i",i),Form("FastJet anti-kT pT (bkg sub, random-cone), hiBin %i - %i",centEdges[i-1],centEdges[i]),NPtBins,ptMin,ptMax);
	h_nPFcand[i]   = new TH1D(Form("h_nPFcand_C%i",i),  Form("N PF cands per event, hiBin %i - %i",centEdges[i-1],centEdges[i]),10000,0,10000);
	h_nPFcandCS[i] = new TH1D(Form("h_nPFcandCS_C%i",i), Form("N PFCS cands per event, hiBin %i - %i",centEdges[i-1],centEdges[i]),10000,0,10000);
	h_randConeEtaPhi[i] = new TProfile2D(Form("h_randConeEtaPhi_C%i",i),Form("Mean random-cone p_{T} vs (#eta,#phi), hiBin %i - %i",centEdges[i-1],centEdges[i]),NRC_EtaBins,etaMin,etaMax,NRC_PhiBins,phiMin,phiMax);
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
      h_fastJetPt[i]->Sumw2();
      h_fastJetPt_JEC[i]->Sumw2();
      h_fastJetPt_bkgSub_RC[i]->Sumw2();
      h_nPFcand[i]->Sumw2();
      h_nPFcandCS[i]->Sumw2();
      h_randConeEtaPhi[i]->Sumw2();

      // loop through jet pt indices
      for(int j = 0; j < NJetPtIndices; j++){
	if(i==0 && j==0){
	  h_muPtOverJetPt[i][j] = new TH1D(Form("h_muPtOverJetPt_C%iJ%i",i,j),Form("p_{T}^{#mu} / p_{T}^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),100,0,1);
	  h_mupt_muptrel[i][j] = new TH2D(Form("h_mupt_muptrel_C%iJ%i",i,j),Form("muPtRel vs. muPt, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f",centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NMuPtBins,muPtMin,muPtMax,NMuRelPtBins,muRelPtMin,muRelPtMax);
	  h_inclRecoJetEta_inclRecoJetPhi[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	}
	else if(i==0){	
	  h_muPtOverJetPt[i][j] = new TH1D(Form("h_muPtOverJetPt_C%iJ%i",i,j),Form("p_{T}^{#mu} / p_{T}^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),100,0,1);
	  h_mupt_muptrel[i][j] = new TH2D(Form("h_mupt_muptrel_C%iJ%i",i,j),Form("muPtRel vs. muPt, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f",centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NMuPtBins,muPtMin,muPtMax,NMuRelPtBins,muRelPtMin,muRelPtMax);
	  h_inclRecoJetEta_inclRecoJetPhi[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[0],centEdges[NCentralityIndices-1],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	}
	else if(j==0){
	  h_muPtOverJetPt[i][j] = new TH1D(Form("h_muPtOverJetPt_C%iJ%i",i,j),Form("p_{T}^{#mu} / p_{T}^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),100,0,1);
	  h_mupt_muptrel[i][j] = new TH2D(Form("h_mupt_muptrel_C%iJ%i",i,j),Form("muPtRel vs. muPt, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f",centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NMuPtBins,muPtMin,muPtMax,NMuRelPtBins,muRelPtMin,muRelPtMax);
	  h_inclRecoJetEta_inclRecoJetPhi[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	}
	else{
	  h_muPtOverJetPt[i][j] = new TH1D(Form("h_muPtOverJetPt_C%iJ%i",i,j),Form("p_{T}^{#mu} / p_{T}^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),100,0,1);
	  h_mupt_muptrel[i][j] = new TH2D(Form("h_mupt_muptrel_C%iJ%i",i,j),Form("muPtRel vs. muPt, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f",centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NMuPtBins,muPtMin,muPtMax,NMuRelPtBins,muRelPtMin,muRelPtMax);
	  h_inclRecoJetEta_inclRecoJetPhi[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	  h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j] = new TH2D(Form("h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn_C%iJ%i",i,j),Form("incl. reco #phi^{jet} vs. incl. reco #eta^{jet}, tagged with incl. reco muon, trigger ON, hiBin %i - %i, p_{T}^{jet} %3.0f - %3.0f", centEdges[i-1],centEdges[i],jetPtEdges[j-1],jetPtEdges[j]),NEtaBins,etaMin,etaMax,NPhiBins,phiMin,phiMax);
	}

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
    em->loadJet("akCs4PFJetAnalyzer/t");
    cout << "	Loading muon..." << endl;
    em->loadMuon("ggHiNtuplizerGED/EventTree");
    cout << "	Loading muon triggers..." << endl;
    em->loadHLT("hltanalysis/HltTree");
    cout << "	Loading PF candidates..." << endl;
    if(doConstituentSubtraction) em->loadParticleFlowAnalyzer("pfcandAnalyzerCS");
    else                         em->loadParticleFlowAnalyzer("pfcandAnalyzer");
    cout << "	Variables initilized!" << endl << endl ;
    // Open a separate file handle so these branch addresses don't overwrite
    // the ones em->loadParticleFlowAnalyzer() set for em->nPFpart/pfPt/etc.
    TFile* f_pfRaw = TFile::Open(input.c_str());
    TTree* pfTreeStd_raw = (TTree*) f_pfRaw->Get("pfcandAnalyzer/pfTree");
    TTree* pfTreeCS_raw  = (TTree*) f_pfRaw->Get("pfcandAnalyzerCS/pfTree");
    int nPFcand_std = 0, nPFcand_cs = 0;
    if(pfTreeStd_raw) pfTreeStd_raw->SetBranchAddress("nPFpart", &nPFcand_std);
    if(pfTreeCS_raw)  pfTreeCS_raw ->SetBranchAddress("nPFpart", &nPFcand_cs);
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


    // event loop
    int eventCounter = 0;
    int evi_frac = 0;
    for(int evi = 0; evi < NEvents; evi++){


      em->getEvent(evi); // load event info from eventMap
      if(pfTreeStd_raw) pfTreeStd_raw->GetEntry(evi);
      if(pfTreeCS_raw)  pfTreeCS_raw ->GetEntry(evi);

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

      h_nPFcand[0]->Fill(nPFcand_std, w);
      h_nPFcand[CentralityIndex]->Fill(nPFcand_std, w);
      h_nPFcandCS[0]->Fill(nPFcand_cs, w);
      h_nPFcandCS[CentralityIndex]->Fill(nPFcand_cs, w);

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

    
   
      ///// PF Candidate Analyzer / Pseudo-Jet Calculator / FastJet Clustering

      int NCandidatesToSample = em->nPFpart;

      // pre-load mixed-event PF candidates from same-centrality events into a pool
      std::vector<double> pool_pfPt, pool_pfEta, pool_pfPhi;
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

      for(int k = 0; k < N_mixedEventsInPool; k++){

	double randEta_k = 3.2*randomGenerator->Rndm() - 1.6;
	double randPhi_k = 2*pi_pfcand*randomGenerator->Rndm() - pi_pfcand;
	double pseudoJetPt_k = 0.;

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
	    }
	  }
	}

	h_pseudoJetPt[0]->Fill(pseudoJetPt_k, w);
	h_pseudoJetPt[CentralityIndex]->Fill(pseudoJetPt_k, w);

	// record cone pT at its (eta,phi) throw location for the UE map
	h_randConeEtaPhi[0]->Fill(randEta_k, randPhi_k, pseudoJetPt_k, w);
	h_randConeEtaPhi[CentralityIndex]->Fill(randEta_k, randPhi_k, pseudoJetPt_k, w);

      }

      // FastJet anti-kT clustering on PF candidates
#ifdef DO_FASTJET
      if(doFastJetClustering){
        std::vector<fastjet::PseudoJet> fjInputs;
        if(doEventMixing && poolSize > 0){
          std::uniform_int_distribution<int> poolDist(0, poolSize - 1);
          for(int s = 0; s < NCandidatesToSample; s++){
            int idx = poolDist(rng);
            double pt  = pool_pfPt[idx];
            double eta = pool_pfEta[idx];
            double phi = pool_pfPhi[idx];
            if(pt < pseudoJetCandPt_min) continue;
            double px = pt * TMath::Cos(phi);
            double py = pt * TMath::Sin(phi);
            double pz = pt * TMath::SinH(eta);
            double E  = pt * TMath::CosH(eta);
            fjInputs.push_back(fastjet::PseudoJet(px, py, pz, E));
          }
        }
        else{
          for(int l = 0; l < em->nPFpart; l++){
            double pt  = em->pfPt->at(l);
            double eta = em->pfEta->at(l);
            double phi = em->pfPhi->at(l);
            if(pt < pseudoJetCandPt_min) continue;
            double px = pt * TMath::Cos(phi);
            double py = pt * TMath::Sin(phi);
            double pz = pt * TMath::SinH(eta);
            double E  = pt * TMath::CosH(eta);
            fjInputs.push_back(fastjet::PseudoJet(px, py, pz, E));
          }
        }
        fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, dR_max_pfcand);
        fastjet::ClusterSequence cs(fjInputs, jetDef);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(0.));
        for(const auto& jet : jets){
          if(TMath::Abs(jet.eta()) > 1.6) continue;
          h_fastJetPt[0]->Fill(jet.pt(), w);
          h_fastJetPt[CentralityIndex]->Fill(jet.pt(), w);
          JEC.SetJetPT(jet.pt());
          JEC.SetJetEta(jet.eta());
          JEC.SetJetPhi(jet.phi());
          double fastJetPt_JEC = JEC.GetCorrectedPT();
          h_fastJetPt_JEC[0]->Fill(fastJetPt_JEC, w);
          h_fastJetPt_JEC[CentralityIndex]->Fill(fastJetPt_JEC, w);
          // RC-subtracted fastJet pT: subtract mean UE pT at the jet's (eta,phi) location
          if(h_RC_map[CentralityIndex]){
            double rcMeanPt = h_RC_map[CentralityIndex]->GetBinContent(
                                h_RC_map[CentralityIndex]->FindBin(jet.eta(), jet.phi()));
            double fastJetPt_rcSub = jet.pt() - rcMeanPt;
	        
            if(fastJetPt_rcSub > 0){
	      std::cout << "event : " << evi << ", cent = " << CentralityIndex << ", jet.pt() = " << jet.pt() << ", rcMeanPt = " << rcMeanPt << ", jet.eta() = " << jet.eta() << ", jet.phi() = " << jet.phi() <<   "\n";
              h_fastJetPt_bkgSub_RC[0]->Fill(fastJetPt_rcSub, w);
              h_fastJetPt_bkgSub_RC[CentralityIndex]->Fill(fastJetPt_rcSub, w);
            }
          }
        }
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
		
	JEC.SetJetPT(em->rawpt[i]);
	JEC.SetJetEta(em->jeteta[i]);
	JEC.SetJetPhi(em->jetphi[i]);

	double x = JEC.GetCorrectedPT();  // use manual JEC
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

	x = applyJEU_JER(x, JEU, JER_fxn, randomGenerator,
	                 fitFxn_PYTHIA_JERCorrection,
	                 neutrino_tag_fraction, neutrino_energy_map);
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
      h_fastJetPt[i]->Write();
      h_fastJetPt_JEC[i]->Write();
      h_fastJetPt_bkgSub_RC[i]->Write();
      h_nPFcand[i]->Write();
      h_nPFcandCS[i]->Write();
      h_randConeEtaPhi[i]->Write();

      for(int j = 0; j < NJetPtIndices; j++){

	h_inclRecoJetEta_inclRecoJetPhi[i][j]->Write();
	h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag[i][j]->Write();
	h_inclRecoJetEta_inclRecoJetPhi_inclRecoMuonTag_triggerOn[i][j]->Write();
	h_muPtOverJetPt[i][j]->Write();
	h_mupt_muptrel[i][j]->Write();
     
      }
    }

    for(int j = 0; j < NJetPtIndices; j++){
      h_muptrel_hiBin[j]->Write();
    }

    wf->Close();
    f_pfRaw->Close();
    return;
    // END WRITE

  }

}
