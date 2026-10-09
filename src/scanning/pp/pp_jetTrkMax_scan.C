// Leading-track (and all-track) ptRel in pp jet60 data, from the 2025 skims of the
// HighEGJet forest that stores the jet trackMax direction (the current pp forests do
// not). 2026-10-09 (user): + track pT axis, for an ATLAS-style data-driven light
// template (arXiv:2204.13530: track-jet pairs in inclusive-jet data, track pT
// reweighted to the light-jet muon pT); see src/plots/bPurity/studyDataDrivenLight_pp_PF.C.
//   h_jetTrkMaxPtRel_trkPt_recoJetPt   leading track: ptRel x track pT x jet pT
// (The skims hold filterTree, evtTree, hltTree and jetTree only -- no track or muon
// tree, checked 2026-10-09 -- so the ATLAS all-track version is not possible here.)
// On the ptRel x muon pT x jet pT binning of the muon templates
// (headers/AnalysisSetup/ptRelMuPt3D.h: track pT from 15 GeV, the muon threshold).
// JEC: the current pp PF set (Spring18 ppRef V6, as pp_scan.C), was Fall17 (2025-11-24).
// Run per skim group:  root -l -b -q 'pp_jetTrkMax_scan.C(<group>)'
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

// event map
#include "../../../eventMap/eventMap.h"
#include "../../../headers/AnalysisSetup/ptRelMuPt3D.h"   // ptRel x (track) pT x jet pT binning
// jet corrector
#include "../../../JetEnergyCorrections/JetCorrector.h"
// general analysis variables
//#include "../../../headers/AnalysisSetupV2p2.h"
#include "../../../headers/AnalysisSetupV2p3.h"
// eta-phi mask function
#include "../../../headers/functions/etaPhiMask.h"
// getDr function
#include "../../../headers/functions/getDr.h"
// getJetPtBin function
#include "../../../headers/functions/getJetPtBin.h"
// getCentBin function
#include "../../../headers/functions/getCentBin_v2.h"
// getPtRel function
#include "../../../headers/functions/getPtRel.h"
// isQualityMuon_hybridSoft function
#include "../../../headers/functions/isQualityMuon_hybridSoft.h"
// isQualityMuon_tight function
#include "../../../headers/functions/isQualityMuon_tight.h"
// isWDecayMuon function
#include "../../../headers/functions/isWDecayMuon.h"
// triggerIsOn function
#include "../../../headers/functions/triggerIsOn.h"
// pthat filter function
#include "../../../headers/functions/passesLeadingGenJetPthatFilter.h"
// JetTrkMax filter function
#include "../../../headers/functions/jet_filter/passesJetTrkMaxFilter.h"
// print introduction
#include "../../../headers/introductions/printIntroduction_pp_scan_V3p7.h"
// analysis config
#include "../../../headers/config/config_pp.h"
//#include "../../../headers/config/config_pp_MB.h"
// read config
#include "../../../headers/config/readConfig.h"
// initialize histograms
TH1D *h_jetPt;
TH1D *h_jetTrkMaxPt[NJetPtIndices];
TH1D *h_jetTrkMaxPtOverJetPt[NJetPtIndices];
TH1D *h_jetTrkMaxEta[NJetPtIndices];
TH1D *h_jetTrkMaxPhi[NJetPtIndices];
TH1D *h_jetTrkMaxDR[NJetPtIndices];
TH1D *h_jetTrkMaxPtRel[NJetPtIndices];
TH3D *h_jetTrkMaxPtRel_trkPt_recoJetPt;   // leading track


///////////////////////  start the program
void pp_jetTrkMax_scan(int group = 1){

  TString input = Form("/eos/cms/store/group/phys_heavyions/cbennett/skims/output_skims_pp_HighEGJet_withJetTrackMaxInfo/pp_skim_output_%i.root",group);
  TString output = Form("/eos/cms/store/group/phys_heavyions/cbennett/scanningOutput/output_pp_jetTrkMax_trkPtAxis_jet60_ppRefJEC_2026-10-9/pp_scan_output_%i.root",group);


  // JET ENERGY CORRECTIONS
  vector<string> Files;
  Files.push_back("../../../JetEnergyCorrections/Spring18_ppRef5TeV_V6_DATA_L2Relative_AK4PF.txt"); // L2Relative correction (as pp_scan.C, 2026-10-09)
  Files.push_back("../../../JetEnergyCorrections/Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt"); // L2L3Residual correction
  // Files.push_back("../../../JetEnergyCorrections/Fall17_17Nov2017F_V6_DATA_L2Relative_AK4PF.txt"); // the 2025-11-24 scan
  // Files.push_back("../../../JetEnergyCorrections/Fall17_17Nov2017F_V6_DATA_L2L3Residual_AK4PF.txt");
  JetCorrector JEC(Files);
  /// >>>>>>>>>>>>>>> print out some info
  printIntroduction_pp_scan_V3p7();
  readConfig();
  /////////////  Define histograms
  h_jetPt = new TH1D("h_jetPt","jetPt",NPtBins,ptMin,ptMax);

  h_jetPt->Sumw2();
  h_jetTrkMaxPtRel_trkPt_recoJetPt = bookPtRelMuPt3D("h_jetTrkMaxPtRel_trkPt_recoJetPt", "leading track p_{T}^{rel} x track p_{T} x jet p_{T}");

  // loop through jet pT indices
  for(int j = 0; j < NJetPtIndices; j++){

    if(j==0){
      
      h_jetTrkMaxPt[j] = new TH1D(Form("h_jetTrkMaxPt_J%i",j),Form("jetTrkMaxPt, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),500,0,500);
      h_jetTrkMaxPtOverJetPt[j] = new TH1D(Form("h_jetTrkMaxPtOverJetPt_J%i",j),Form("jetTrkMaxPt / jetPt, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),100,0,1);
      h_jetTrkMaxEta[j] = new TH1D(Form("h_jetTrkMaxEta_J%i",j),Form("jetTrkMaxEta, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NTrkEtaBins,trkEtaMin,trkEtaMax);
      h_jetTrkMaxPhi[j] = new TH1D(Form("h_jetTrkMaxPhi_J%i",j),Form("jetTrkMaxPhi, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NPhiBins,phiMin,phiMax);
      h_jetTrkMaxDR[j] = new TH1D(Form("h_jetTrkMaxDR_J%i",j),Form("jetTrkMaxDR, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NdRBins,dRBinMin,dRBinMax);
      h_jetTrkMaxPtRel[j] = new TH1D(Form("h_jetTrkMaxPtRel_J%i",j),Form("jetTrkMaxPtRel, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[0],jetPtEdges[NJetPtIndices-1]),NMuRelPtBins,muRelPtMin,muRelPtMax);

    }

    else{

      h_jetTrkMaxPt[j] = new TH1D(Form("h_jetTrkMaxPt_J%i",j),Form("jetTrkMaxPt, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[j-1],jetPtEdges[j]),500,0,500);
      h_jetTrkMaxPtOverJetPt[j] = new TH1D(Form("h_jetTrkMaxPtOverJetPt_J%i",j),Form("jetTrkMaxPt / jetPt, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[j-1],jetPtEdges[j]),100,0,1);
      h_jetTrkMaxEta[j] = new TH1D(Form("h_jetTrkMaxEta_J%i",j),Form("jetTrkMaxEta, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[j-1],jetPtEdges[j]),NTrkEtaBins,trkEtaMin,trkEtaMax);
      h_jetTrkMaxPhi[j] = new TH1D(Form("h_jetTrkMaxPhi_J%i",j),Form("jetTrkMaxPhi, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[j-1],jetPtEdges[j]),NPhiBins,phiMin,phiMax);
      h_jetTrkMaxDR[j] = new TH1D(Form("h_jetTrkMaxDR_J%i",j),Form("jetTrkMaxDR, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[j-1],jetPtEdges[j]),NdRBins,dRBinMin,dRBinMax);
      h_jetTrkMaxPtRel[j] = new TH1D(Form("h_jetTrkMaxPtRel_J%i",j),Form("jetTrkMaxPtRel, p_{T}^{jet} %3.0f - %3.0f",jetPtEdges[j-1],jetPtEdges[j]),NMuRelPtBins,muRelPtMin,muRelPtMax);

    }

    h_jetTrkMaxPt[j]->Sumw2();
    h_jetTrkMaxPtOverJetPt[j]->Sumw2();
    h_jetTrkMaxEta[j]->Sumw2();
    h_jetTrkMaxPhi[j]->Sumw2();
    h_jetTrkMaxDR[j]->Sumw2();
    h_jetTrkMaxPtRel[j]->Sumw2();
    
  }


  TFile *f = TFile::Open(input);
  cout << "	File opened!" << endl;
  auto em = new eventMap(f);
  em->isMC = isMC_status;
  em->AASetup = AASetup_status;
  cout << "	Initializing variables ... " << endl;
  em->init("evtTree");   // the skims' event tree (eventMap::init() opens the forest's hiEvtAnalyzer/HiTree)
  cout << "	Loading jet..." << endl;
  em->loadJet(jetTreeString);
  // no muon tree in the skims (and no muons are used here)
  cout << "	Loading muon triggers..." << endl;
  em->loadHLT(hltString);   // was loadMuonTrigger (renamed in eventMap.h)
  // no track tree in the skims: only the jet trackMax branches
  // data: no gen particles (the old em->loadGenParticle() no longer compiles)
  cout << "	Variables initilized!" << endl << endl ;
  int NEvents = em->evtTree->GetEntries();
  cout << "	Number of events = " << NEvents << endl;


  // event filters: already applied in the skims (their filterTree is not the forest's
  // skimanalysis/HltTree that regEventFilter opens, and unregistered filters would
  // reject every event)

  
  // event loop
  int evi_frac = 0;
  for(int evi = 0; evi < NEvents; evi++){

    if(evi == 0) cout << "Processing events..." << endl;

    em->getEvent(evi); // load event info from eventMap
    // em->muonTriggerTree->GetEntry(evi);
    // em->jetEvtTree->GetEntry(evi);
    // if(em->njet > 0){
    //   em->recoJetTree->GetEntry(evi);
    // }
    // em->muonEvtTree->GetEntry(evi);
    // if(em->nMu > 0){
    //   em->muonTree->GetEntry(evi);
    // }

    if((100*evi / NEvents) % 5 == 0 && (100*evi / NEvents) > evi_frac){

      cout << "evt frac: " << evi_frac << "%" << endl;

    }

    evi_frac = 100*evi / NEvents;

    // global event cuts
    if(fabs(em->vz) > 15.0) continue;

    // event filters

    if(em->HLT_HIAK4PFJet60_v1 == 0) continue;
    //if(em->HLT_HIAK4PFJet80_v1 == 0) continue;
    //if(em->HLT_HIAK4PFJet100_v1 == 0) continue;

    // In data, event weight = 1
    double w = 1.0;

    // RECO JET LOOP
    for(int i = 0; i < em->njet ; i++){

      // JET VARIABLES
	
      JEC.SetJetPT(em->rawpt[i]);
      JEC.SetJetEta(em->jeteta[i]);
      JEC.SetJetPhi(em->jetphi[i]);

      double recoJetPt_i = JEC.GetCorrectedPT();  // apply manual JEC
      double recoJetEta_i = em->jeteta[i]; // recoJetEta
      double recoJetPhi_i = em->jetphi[i]; // recoJetPhi
      double jetTrkMax_i = em->jetTrkMax[i];
      double jetTrkMaxEta_i = em->jetTrkMaxEta[i];
      double jetTrkMaxPhi_i = em->jetTrkMaxPhi[i];
      double jetTrkMaxDR_i = em->jetTrkMaxDR[i];
      double jetTrkMaxPtRel_i = getPtRel(jetTrkMax_i,jetTrkMaxEta_i,jetTrkMaxPhi_i,recoJetPt_i,recoJetEta_i,recoJetPhi_i);

      if(doJetTrkMaxFilter){
	if(!passesJetTrkMaxFilter(jetTrkMax_i,recoJetPt_i)) continue;
      }
     
      if(doEtaPhiMask){
	if(etaPhiMask(recoJetEta_i,recoJetPhi_i)) continue;
      }

      //cout << "x" << endl;
     
     		
      if(TMath::Abs(recoJetEta_i) > 1.6 || recoJetPt_i < 80.) continue;

      if(jetTrkMax_i < 14.) continue;
      h_jetTrkMaxPtRel_trkPt_recoJetPt->Fill(jetTrkMaxPtRel_i, jetTrkMax_i, recoJetPt_i, w);

      int jetPtIndex = getJetPtBin(recoJetPt_i);
      
      if(jetPtIndex < 0) continue;

      h_jetPt->Fill(recoJetPt_i,w);
      
      h_jetTrkMaxPt[0]->Fill(jetTrkMax_i,w);
      h_jetTrkMaxPtOverJetPt[0]->Fill(jetTrkMax_i/recoJetPt_i,w);
      h_jetTrkMaxEta[0]->Fill(jetTrkMaxEta_i,w);
      h_jetTrkMaxPhi[0]->Fill(jetTrkMaxPhi_i,w);
      h_jetTrkMaxDR[0]->Fill(jetTrkMaxDR_i,w);
      h_jetTrkMaxPtRel[0]->Fill(jetTrkMaxPtRel_i,w);

      h_jetTrkMaxPt[jetPtIndex]->Fill(jetTrkMax_i,w);
      h_jetTrkMaxPtOverJetPt[jetPtIndex]->Fill(jetTrkMax_i/recoJetPt_i,w);
      h_jetTrkMaxEta[jetPtIndex]->Fill(jetTrkMaxEta_i,w);
      h_jetTrkMaxPhi[jetPtIndex]->Fill(jetTrkMaxPhi_i,w);
      h_jetTrkMaxDR[jetPtIndex]->Fill(jetTrkMaxDR_i,w);
      h_jetTrkMaxPtRel[jetPtIndex]->Fill(jetTrkMaxPtRel_i,w);

  
    }
    // END recoJet LOOP

 

  } // end event loop


 
  delete f;
  // WRITE
  auto wf = TFile::Open(output,"recreate");
  // >>>>>>>>>> write histograms

  h_jetPt->Write();
  h_jetTrkMaxPtRel_trkPt_recoJetPt->Write();

  for(int j = 0; j < NJetPtIndices; j++){

    h_jetTrkMaxPt[j]->Write();
    h_jetTrkMaxPtOverJetPt[j]->Write();
    h_jetTrkMaxEta[j]->Write();
    h_jetTrkMaxPhi[j]->Write();
    h_jetTrkMaxDR[j]->Write();
    h_jetTrkMaxPtRel[j]->Write();
      
  }


  wf->Close();
  return;
  // END WRITE

}
