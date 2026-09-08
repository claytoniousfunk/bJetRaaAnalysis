#include "../../headers/functions/divideByBinwidth.h"
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

bool doFakeJetSubtraction = false;
bool doUnfolding = false;
bool useMinBiasOnly = false;

TString fakeJetsPath  = "../../rootFiles/fakeJets/fakeJets.root";
TString outputRootPath = "./rootFiles/JetsPerZ/histograms_JetsPerZ_lightJets_rebinned.root";

int N_iter_pp = 1;
int N_iter_C4 = 1;
int N_iter_C3 = 1;
int N_iter_C2 = 1;
int N_iter_C1 = 1;

// scanning files
TFile *file_pp_HighEGJet_jet60, *file_pp_HighEGJet_jet80, *file_pp_HighEGJet_jet100, *file_pp_SingleMuon, *file_pp_MinBias;
TFile *file_PbPb_HardProbes_jet60, *file_PbPb_HardProbes_jet80, *file_PbPb_HardProbes_jet100, *file_PbPb_SingleMuon, *file_PbPb_MinBias;
TFile *file_PYTHIA_response, *file_PH_response;

// pub data files
TFile *f_pubATLAS_C1, *f_pubATLAS_C2, *f_pubATLAS_C3, *f_pubATLAS_C4;
TH1D *h_pubATLAS_C1, *h_pubATLAS_C1_err, *h_pubATLAS_C1_orig;
TH1D *h_pubATLAS_C2, *h_pubATLAS_C2_err, *h_pubATLAS_C2_orig;
TH1D *h_pubATLAS_C3, *h_pubATLAS_C3_err, *h_pubATLAS_C3_orig;
TH1D *h_pubATLAS_C4, *h_pubATLAS_C4_err, *h_pubATLAS_C4_orig;
TH1D *h_pubATLAS_spectra_pp, *h_pubATLAS_spectra_pp_err, *h_pubATLAS_spectra_pp_orig;
TH1D *h_pubATLAS_spectra_C1, *h_pubATLAS_spectra_C1_err, *h_pubATLAS_spectra_C1_orig;
TDirectoryFile *dirATLAS_C1, *dirATLAS_C2, *dirATLAS_C3, *dirATLAS_C4;
TDirectoryFile *dirATLAS_spectra_pp, *dirATLAS_spectra_C1;
TFile *f_pubCMS_C1, *f_pubCMS_C2, *f_pubCMS_C3, *f_pubCMS_C4;
TH1D *h_pubCMS_C1, *h_pubCMS_C1_err, *h_pubCMS_C1_orig;
TH1D *h_pubCMS_C2, *h_pubCMS_C2_err, *h_pubCMS_C2_orig;
TH1D *h_pubCMS_C3, *h_pubCMS_C3_err, *h_pubCMS_C3_orig;
TH1D *h_pubCMS_C4, *h_pubCMS_C4_err, *h_pubCMS_C4_orig;
TDirectoryFile *dirCMS_C1, *dirCMS_C2, *dirCMS_C3, *dirCMS_C4;
TFile *f_pubATLAS_spectra_pp, *f_pubATLAS_spectra_C1;

// loaded histograms
TH1D *h_pp, *h_C4, *h_C3, *h_C2, *h_C1;
TH1D *h_pp_jetMB, *h_C4_jetMB, *h_C3_jetMB, *h_C2_jetMB, *h_C1_jetMB;
TH1D *h_pp_jet60, *h_C4_jet60, *h_C3_jet60, *h_C2_jet60, *h_C1_jet60;
TH1D *h_pp_jet80, *h_C4_jet80, *h_C3_jet80, *h_C2_jet80, *h_C1_jet80;
TH1D *h_pp_jet100, *h_C4_jet100, *h_C3_jet100, *h_C2_jet100, *h_C1_jet100;
TH1D *hiBin_pp, *hiBin_PbPb;
TH1D *h_vz_pp, *h_vz_C4, *h_vz_C3, *h_vz_C2, *h_vz_C1;
TH1D *r_C4, *r_C3, *r_C2, *r_C1;
double N_coll_avg_C4, N_coll_avg_C3, N_coll_avg_C2, N_coll_avg_C1;
TH1D *h_vz_MinBias_C4, *h_vz_MinBias_C3, *h_vz_MinBias_C2, *h_vz_MinBias_C1;
TH1D *h_vz_MinBias_Jet100_C4, *h_vz_MinBias_Jet100_C3, *h_vz_MinBias_Jet100_C2, *h_vz_MinBias_Jet100_C1;
TH1D *h_vz_MinBias_Jet100_withJetAbove130_C4, *h_vz_MinBias_Jet100_withJetAbove130_C3, *h_vz_MinBias_Jet100_withJetAbove130_C2, *h_vz_MinBias_Jet100_withJetAbove130_C1;
double N_evt_MinBias_C4, N_evt_MinBias_C3, N_evt_MinBias_C2, N_evt_MinBias_C1;
double N_evt_MinBias_Jet100_C4, N_evt_MinBias_Jet100_C3, N_evt_MinBias_Jet100_C2, N_evt_MinBias_Jet100_C1;
double N_evt_MinBias_Jet100_withJetAbove130_C4, N_evt_MinBias_Jet100_withJetAbove130_C3, N_evt_MinBias_Jet100_withJetAbove130_C2, N_evt_MinBias_Jet100_withJetAbove130_C1;
double minBiasEquivalentEventConversion_C4, minBiasEquivalentEventConversion_C3, minBiasEquivalentEventConversion_C2, minBiasEquivalentEventConversion_C1;
double triggerLossCorrection_C4, triggerLossCorrection_C3, triggerLossCorrection_C2, triggerLossCorrection_C1;

TH1D *h_dimuonMass_pp, *h_dimuonMass_C4, *h_dimuonMass_C3, *h_dimuonMass_C2, *h_dimuonMass_C1;

TFile *file_fakeJets;
TH1D *h_fakeJets_C4, *h_fakeJets_C3, *h_fakeJets_C2, *h_fakeJets_C1;

TH1D *h_C1_clone, *h_fakeJets_C1_clone, *h_C1_sub, *r_C1_sub;
TH1D *h_C2_clone, *h_fakeJets_C2_clone, *h_C2_sub, *r_C2_sub;
TH1D *h_C3_clone, *h_fakeJets_C3_clone, *h_C3_sub, *r_C3_sub;
TH1D *h_C4_clone, *h_fakeJets_C4_clone, *h_C4_sub, *r_C4_sub;

TCanvas *canv_fakeJets_C1, *canv_fakeJets_C2, *canv_fakeJets_C3, *canv_fakeJets_C4;
TPad *pad_fakeJets_C1_upper, *pad_fakeJets_C1_lower;
TPad *pad_fakeJets_C2_upper, *pad_fakeJets_C2_lower;
TPad *pad_fakeJets_C3_upper, *pad_fakeJets_C3_lower;
TPad *pad_fakeJets_C4_upper, *pad_fakeJets_C4_lower;
TLegend *leg_fakeJets_C1, *leg_fakeJets_C2, *leg_fakeJets_C3, *leg_fakeJets_C4;

TH2D *h_response_pp, *h_response_C4, *h_response_C3, *h_response_C2, *h_response_C1;
TH1D *h_meas_pp, *h_meas_C4, *h_meas_C3, *h_meas_C2, *h_meas_C1;
TH1D *h_truth_pp, *h_truth_C4, *h_truth_C3, *h_truth_C2, *h_truth_C1;

RooUnfoldResponse* response_pp = nullptr;
RooUnfoldResponse* response_C4 = nullptr;
RooUnfoldResponse* response_C3 = nullptr;
RooUnfoldResponse* response_C2 = nullptr;
RooUnfoldResponse* response_C1 = nullptr;

RooUnfoldBayes* unfold_pp = nullptr;
RooUnfoldBayes* unfold_C4 = nullptr;
RooUnfoldBayes* unfold_C3 = nullptr;
RooUnfoldBayes* unfold_C2 = nullptr;
RooUnfoldBayes* unfold_C1 = nullptr;

TH1D *h_pp_unfold, *h_C4_unfold, *h_C3_unfold, *h_C2_unfold, *h_C1_unfold;
TLegend *leg;

const int N_edge = 21;
double newAxis[N_edge] = {60,65,70,75,80,85,90,95,100,110,120,130,140,150,160,180,200,240,280,350,500};

double NZ_pp, NZ_C4, NZ_C3, NZ_C2, NZ_C1;

TH1D *R_C4, *R_C3, *R_C2, *R_C1; // systematics from NZ counts

TCanvas *canv_prime;
TPad *pad_prime;


void openScanFiles(){
  
  file_pp_HighEGJet_jet60 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/pp_HighEGJet_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_2026-3-10.root");
  file_pp_HighEGJet_jet80 = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_HighEGJet_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-27.root");
  file_pp_HighEGJet_jet100 = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_HighEGJet_Jet100HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-27.root");
  //file_pp_SingleMuon = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/pp_SingleMuon_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-3-16.root");
  file_pp_SingleMuon = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/pp_SingleMuon_mu12_pTmu-15to999_tight_deltaR-40_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root");
  file_pp_MinBias = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root");

  file_PbPb_HardProbes_jet60 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet60HLT_mu12_pTmu-15to999_tight_WDecayFilter_2026-3-11.root");
  file_PbPb_HardProbes_jet80 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet80HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-5.root");
  file_PbPb_HardProbes_jet100 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet100HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-5.root");
  file_PbPb_SingleMuon = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root");
  file_PbPb_MinBias = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_MinBias_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root");
  //file_PbPb_MinBias = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_noRhoModification_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-25_ultraFineCentBins.root");

  file_PYTHIA_response = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIA/PYTHIA_DiJet_response_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-8-3.root");
  file_PH_response = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/response/PYTHIAHYDJET_response_DiJet_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-4-6_evenEvents.root");
  
  return;
}

void openPubDataFiles(){

  f_pubATLAS_C1 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1673184-v1-Table_19.root");
  f_pubATLAS_C2 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1673184-v1-Table_20.root");
  f_pubATLAS_C3 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1673184-v1-Table_22.root");
  f_pubATLAS_C4 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1673184-v1-Table_24.root");
  
  f_pubCMS_C1 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1848440-v1-Figure_7-3.root");
  f_pubCMS_C2 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1848440-v1-Figure_7-3.root");
  f_pubCMS_C3 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1848440-v1-Figure_7-3.root");
  f_pubCMS_C4 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1848440-v1-Figure_7-3.root");

  f_pubATLAS_spectra_pp = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1673184-v1-Table_4.root");
  f_pubATLAS_spectra_C1 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/HEPData/HEPData-ins1673184-v1-Table_11.root");

  dirATLAS_C1 = (TDirectoryFile*) f_pubATLAS_C1->Get("Table 19");
  dirATLAS_C2 = (TDirectoryFile*) f_pubATLAS_C2->Get("Table 20");
  dirATLAS_C3 = (TDirectoryFile*) f_pubATLAS_C3->Get("Table 22");
  dirATLAS_C4 = (TDirectoryFile*) f_pubATLAS_C4->Get("Table 24");

  dirCMS_C1 = (TDirectoryFile*) f_pubCMS_C1->Get("Figure 7-3");
  dirCMS_C2 = (TDirectoryFile*) f_pubCMS_C2->Get("Figure 7-3");
  dirCMS_C3 = (TDirectoryFile*) f_pubCMS_C3->Get("Figure 7-3");
  dirCMS_C4 = (TDirectoryFile*) f_pubCMS_C4->Get("Figure 7-3");

  dirATLAS_spectra_pp = (TDirectoryFile*) f_pubATLAS_spectra_pp->Get("Table 4");
  dirATLAS_spectra_C1 = (TDirectoryFile*) f_pubATLAS_spectra_C1->Get("Table 11");

  return;
  
}

void stylizeHistograms(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4){

  double lw  = 2;
  double ms = 2;
  
  h0->SetLineColor(kBlack);
  h0->SetMarkerColor(kBlack);
  h0->SetLineWidth(lw);
  h0->SetMarkerSize(ms);
  h0->SetMarkerStyle(24);
  h0->SetStats(0);
  h0->SetTitle("");

  h4->SetLineColor(kBlack);
  h4->SetMarkerColor(kBlack);
  h4->SetLineWidth(lw);
  h4->SetMarkerSize(ms);
  h4->SetMarkerStyle(33);
  h4->SetStats(0);
  h4->SetTitle("");

  h3->SetLineColor(kBlue-4);
  h3->SetMarkerColor(kBlue-4);
  h3->SetLineWidth(lw);
  h3->SetMarkerSize(ms);
  h3->SetMarkerStyle(21);
  h3->SetStats(0);
  h3->SetTitle("");
  
  h2->SetLineColor(kGreen+2);
  h2->SetMarkerColor(kGreen+2);
  h2->SetLineWidth(lw);
  h2->SetMarkerSize(ms);
  h2->SetMarkerStyle(34);
  h2->SetStats(0);
  h2->SetTitle("");

  h1->SetLineColor(kRed-4);
  h1->SetMarkerColor(kRed-4);
  h1->SetLineWidth(lw);
  h1->SetMarkerSize(ms);
  h1->SetMarkerStyle(47);
  h1->SetStats(0);
  h1->SetTitle("");

  return;

}

void stylizeSystematics(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4){

  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0);
  h3->SetMarkerSize(0);
  h4->SetMarkerSize(0);

  double alpha = 0.5;
  h1->SetFillColorAlpha(kRed-4,alpha);
  h2->SetFillColorAlpha(kGreen+2,alpha);
  h3->SetFillColorAlpha(kBlue-4,alpha);
  h4->SetFillColorAlpha(kGray,alpha);

}


void stylizePubHistograms_ATLAS(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4){

  double lw  = 2;
  double ms = 2;

  h1->SetLineColor(kRed-4);
  h1->SetMarkerColor(kRed-4);
  h1->SetLineWidth(lw);
  h1->SetMarkerSize(ms);
  h1->SetMarkerStyle(24);
  h1->SetStats(0);
  h1->SetTitle("");

  h2->SetLineColor(kGreen+2);
  h2->SetMarkerColor(kGreen+2);
  h2->SetLineWidth(lw);
  h2->SetMarkerSize(ms);
  h2->SetMarkerStyle(24);
  h2->SetStats(0);
  h2->SetTitle("");

  h3->SetLineColor(kBlue-4);
  h3->SetMarkerColor(kBlue-4);
  h3->SetLineWidth(lw);
  h3->SetMarkerSize(ms);
  h3->SetMarkerStyle(24);
  h3->SetStats(0);
  h3->SetTitle("");

  h4->SetLineColor(kBlack);
  h4->SetMarkerColor(kBlack);
  h4->SetLineWidth(lw);
  h4->SetMarkerSize(ms);
  h4->SetMarkerStyle(24);
  h4->SetStats(0);
  h4->SetTitle("");

}

void stylizePubHistograms_CMS(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4){

  double lw  = 2;
  double ms = 2;

  h1->SetLineColor(kRed-4);
  h1->SetMarkerColor(kRed-4);
  h1->SetLineWidth(lw);
  h1->SetMarkerSize(ms);
  h1->SetMarkerStyle(46);
  h1->SetStats(0);
  h1->SetTitle("");

  h2->SetLineColor(kGreen+2);
  h2->SetMarkerColor(kGreen+2);
  h2->SetLineWidth(lw);
  h2->SetMarkerSize(ms);
  h2->SetMarkerStyle(28);
  h2->SetStats(0);
  h2->SetTitle("");

  h3->SetLineColor(kBlue-4);
  h3->SetMarkerColor(kBlue-4);
  h3->SetLineWidth(lw);
  h3->SetMarkerSize(ms);
  h3->SetMarkerStyle(25);
  h3->SetStats(0);
  h3->SetTitle("");

  h4->SetLineColor(kBlack);
  h4->SetMarkerColor(kBlack);
  h4->SetLineWidth(lw);
  h4->SetMarkerSize(ms);
  h4->SetMarkerStyle(27);
  h4->SetStats(0);
  h4->SetTitle("");

}


void buildPubDataHistograms(){

  // ATLAS
  
  h_pubATLAS_C1 = (TH1D*) dirATLAS_C1->Get("Hist1D_y1");
  h_pubATLAS_C1_orig = (TH1D*) dirATLAS_C1->Get("Hist1D_y1");
  h_pubATLAS_C1_err = (TH1D*) dirATLAS_C1->Get("Hist1D_y1_e1plus");
  for(int i = 0; i < h_pubATLAS_C1->GetSize(); i++){
    h_pubATLAS_C1->SetBinContent(i,h_pubATLAS_C1_orig->GetBinContent(i));
    h_pubATLAS_C1->SetBinError(i,h_pubATLAS_C1_err->GetBinContent(i));
  }

  h_pubATLAS_C2 = (TH1D*) dirATLAS_C2->Get("Hist1D_y1");
  h_pubATLAS_C2_orig = (TH1D*) dirATLAS_C2->Get("Hist1D_y1");
  h_pubATLAS_C2_err = (TH1D*) dirATLAS_C2->Get("Hist1D_y1_e1plus");
  for(int i = 0; i < h_pubATLAS_C2->GetSize(); i++){
    h_pubATLAS_C2->SetBinContent(i,h_pubATLAS_C2_orig->GetBinContent(i));
    h_pubATLAS_C2->SetBinError(i,h_pubATLAS_C2_err->GetBinContent(i));
  }

  h_pubATLAS_C3 = (TH1D*) dirATLAS_C3->Get("Hist1D_y1");
  h_pubATLAS_C3_orig = (TH1D*) dirATLAS_C3->Get("Hist1D_y1");
  h_pubATLAS_C3_err = (TH1D*) dirATLAS_C3->Get("Hist1D_y1_e1plus");
  for(int i = 0; i < h_pubATLAS_C3->GetSize(); i++){
    h_pubATLAS_C3->SetBinContent(i,h_pubATLAS_C3_orig->GetBinContent(i));
    h_pubATLAS_C3->SetBinError(i,h_pubATLAS_C3_err->GetBinContent(i));
  }

  h_pubATLAS_C4 = (TH1D*) dirATLAS_C4->Get("Hist1D_y1");
  h_pubATLAS_C4_orig = (TH1D*) dirATLAS_C4->Get("Hist1D_y1");
  h_pubATLAS_C4_err = (TH1D*) dirATLAS_C4->Get("Hist1D_y1_e1plus");
  for(int i = 0; i < h_pubATLAS_C4->GetSize(); i++){
    h_pubATLAS_C4->SetBinContent(i,h_pubATLAS_C4_orig->GetBinContent(i));
    h_pubATLAS_C4->SetBinError(i,h_pubATLAS_C4_err->GetBinContent(i));
  }

  h_pubATLAS_spectra_pp = (TH1D*) dirATLAS_spectra_pp->Get("Hist1D_y1");
  h_pubATLAS_spectra_pp_orig = (TH1D*) dirATLAS_spectra_pp->Get("Hist1D_y1");
  h_pubATLAS_spectra_pp_err = (TH1D*) dirATLAS_spectra_pp->Get("Hist1D_y1_e1plus");
  for(int i = 0; i < h_pubATLAS_spectra_pp->GetSize(); i++){
    h_pubATLAS_spectra_pp->SetBinContent(i,h_pubATLAS_spectra_pp_orig->GetBinContent(i));
    h_pubATLAS_spectra_pp->SetBinError(i,h_pubATLAS_spectra_pp_err->GetBinContent(i));
  }

  h_pubATLAS_spectra_C1 = (TH1D*) dirATLAS_spectra_C1->Get("Hist1D_y1");
  h_pubATLAS_spectra_C1_orig = (TH1D*) dirATLAS_spectra_C1->Get("Hist1D_y1");
  h_pubATLAS_spectra_C1_err = (TH1D*) dirATLAS_spectra_C1->Get("Hist1D_y1_e1plus");
  for(int i = 0; i < h_pubATLAS_spectra_C1->GetSize(); i++){
    h_pubATLAS_spectra_C1->SetBinContent(i,h_pubATLAS_spectra_C1_orig->GetBinContent(i));
    h_pubATLAS_spectra_C1->SetBinError(i,h_pubATLAS_spectra_C1_err->GetBinContent(i));
  }

  // CMS

  h_pubCMS_C1 = (TH1D*) dirCMS_C1->Get("Hist1D_y1");
  h_pubCMS_C1_orig = (TH1D*) dirCMS_C1->Get("Hist1D_y1");
  h_pubCMS_C1_err = (TH1D*) dirCMS_C1->Get("Hist1D_y1_e1plus");
  for(int i = 0; i < h_pubCMS_C1->GetSize(); i++){
    h_pubCMS_C1->SetBinContent(i,h_pubCMS_C1_orig->GetBinContent(i));
    h_pubCMS_C1->SetBinError(i,h_pubCMS_C1_err->GetBinContent(i));
  }

  h_pubCMS_C2 = (TH1D*) dirCMS_C2->Get("Hist1D_y2");
  h_pubCMS_C2_orig = (TH1D*) dirCMS_C2->Get("Hist1D_y2");
  h_pubCMS_C2_err = (TH1D*) dirCMS_C2->Get("Hist1D_y2_e1plus");
  for(int i = 0; i < h_pubCMS_C2->GetSize(); i++){
    h_pubCMS_C2->SetBinContent(i,h_pubCMS_C2_orig->GetBinContent(i));
    h_pubCMS_C2->SetBinError(i,h_pubCMS_C2_err->GetBinContent(i));
  }

  h_pubCMS_C3 = (TH1D*) dirCMS_C3->Get("Hist1D_y3");
  h_pubCMS_C3_orig = (TH1D*) dirCMS_C3->Get("Hist1D_y3");
  h_pubCMS_C3_err = (TH1D*) dirCMS_C3->Get("Hist1D_y3_e1plus");
  for(int i = 0; i < h_pubCMS_C3->GetSize(); i++){
    h_pubCMS_C3->SetBinContent(i,h_pubCMS_C3_orig->GetBinContent(i));
    h_pubCMS_C3->SetBinError(i,h_pubCMS_C3_err->GetBinContent(i));
  }

  h_pubCMS_C4 = (TH1D*) dirCMS_C4->Get("Hist1D_y4");
  h_pubCMS_C4_orig = (TH1D*) dirCMS_C4->Get("Hist1D_y4");
  h_pubCMS_C4_err = (TH1D*) dirCMS_C4->Get("Hist1D_y4_e1plus");
  for(int i = 0; i < h_pubCMS_C4->GetSize(); i++){
    h_pubCMS_C4->SetBinContent(i,h_pubCMS_C4_orig->GetBinContent(i));
    h_pubCMS_C4->SetBinError(i,h_pubCMS_C4_err->GetBinContent(i));
  }


  stylizePubHistograms_ATLAS(h_pubATLAS_C1, h_pubATLAS_C2, h_pubATLAS_C3, h_pubATLAS_C4);
  stylizePubHistograms_CMS(h_pubCMS_C1, h_pubCMS_C2, h_pubCMS_C3, h_pubCMS_C4);

  return;

}

void fetchHistograms(){

  file_pp_SingleMuon->GetObject("h_dimuonMass",h_dimuonMass_pp);
  file_PbPb_SingleMuon->GetObject("h_dimuonMass_C4",h_dimuonMass_C4);
  file_PbPb_SingleMuon->GetObject("h_dimuonMass_C3",h_dimuonMass_C3);
  file_PbPb_SingleMuon->GetObject("h_dimuonMass_C2",h_dimuonMass_C2);
  file_PbPb_SingleMuon->GetObject("h_dimuonMass_C1",h_dimuonMass_C1);

  file_PbPb_MinBias->GetObject("h_vz_C4",h_vz_MinBias_C4);
  file_PbPb_MinBias->GetObject("h_vz_C3",h_vz_MinBias_C3);
  file_PbPb_MinBias->GetObject("h_vz_C2",h_vz_MinBias_C2);
  file_PbPb_MinBias->GetObject("h_vz_C1",h_vz_MinBias_C1);

  file_PbPb_MinBias->GetObject("h_vz_jet100_withJetAbove130_C4",h_vz_MinBias_Jet100_withJetAbove130_C4);
  file_PbPb_MinBias->GetObject("h_vz_jet100_withJetAbove130_C3",h_vz_MinBias_Jet100_withJetAbove130_C3);
  file_PbPb_MinBias->GetObject("h_vz_jet100_withJetAbove130_C2",h_vz_MinBias_Jet100_withJetAbove130_C2);
  file_PbPb_MinBias->GetObject("h_vz_jet100_withJetAbove130_C1",h_vz_MinBias_Jet100_withJetAbove130_C1);

  file_PbPb_MinBias->GetObject("h_vz_jet100_C4",h_vz_MinBias_Jet100_C4);
  file_PbPb_MinBias->GetObject("h_vz_jet100_C3",h_vz_MinBias_Jet100_C3);
  file_PbPb_MinBias->GetObject("h_vz_jet100_C2",h_vz_MinBias_Jet100_C2);
  file_PbPb_MinBias->GetObject("h_vz_jet100_C1",h_vz_MinBias_Jet100_C1);

  triggerLossCorrection_C4 = h_vz_MinBias_Jet100_C4->Integral() / h_vz_MinBias_Jet100_withJetAbove130_C4->Integral();
  triggerLossCorrection_C3 = h_vz_MinBias_Jet100_C3->Integral() / h_vz_MinBias_Jet100_withJetAbove130_C3->Integral();
  triggerLossCorrection_C2 = h_vz_MinBias_Jet100_C2->Integral() / h_vz_MinBias_Jet100_withJetAbove130_C2->Integral();
  triggerLossCorrection_C1 = h_vz_MinBias_Jet100_C1->Integral() / h_vz_MinBias_Jet100_withJetAbove130_C1->Integral();

  minBiasEquivalentEventConversion_C4 = h_vz_MinBias_C4->Integral() / h_vz_MinBias_Jet100_C4->Integral();
  minBiasEquivalentEventConversion_C3 = h_vz_MinBias_C3->Integral() / h_vz_MinBias_Jet100_C3->Integral();
  minBiasEquivalentEventConversion_C2 = h_vz_MinBias_C2->Integral() / h_vz_MinBias_Jet100_C2->Integral();
  minBiasEquivalentEventConversion_C1 = h_vz_MinBias_C1->Integral() / h_vz_MinBias_Jet100_C1->Integral();

  file_pp_HighEGJet_jet60->GetObject("h_inclRecoJetPt",h_pp_jet60);
  file_pp_HighEGJet_jet80->GetObject("h_inclRecoJetPt",h_pp_jet80);
  file_pp_HighEGJet_jet100->GetObject("h_inclRecoJetPt",h_pp_jet100);
  file_pp_MinBias->GetObject("h_inclRecoJetPt",h_pp_jetMB);

  file_pp_HighEGJet_jet100->GetObject("h_vz",h_vz_pp);

  file_PbPb_MinBias->GetObject("h_inclRecoJetPt_C4",h_C4_jetMB);
  file_PbPb_MinBias->GetObject("h_inclRecoJetPt_C3",h_C3_jetMB);
  file_PbPb_MinBias->GetObject("h_inclRecoJetPt_C2",h_C2_jetMB);
  file_PbPb_MinBias->GetObject("h_inclRecoJetPt_C1",h_C1_jetMB);

  file_PbPb_HardProbes_jet60->GetObject("h_inclRecoJetPt_C4",h_C4_jet60);
  file_PbPb_HardProbes_jet60->GetObject("h_inclRecoJetPt_C3",h_C3_jet60);
  file_PbPb_HardProbes_jet60->GetObject("h_inclRecoJetPt_C2",h_C2_jet60);
  file_PbPb_HardProbes_jet60->GetObject("h_inclRecoJetPt_C1",h_C1_jet60);

  file_PbPb_HardProbes_jet80->GetObject("h_inclRecoJetPt_C4",h_C4_jet80);
  file_PbPb_HardProbes_jet80->GetObject("h_inclRecoJetPt_C3",h_C3_jet80);
  file_PbPb_HardProbes_jet80->GetObject("h_inclRecoJetPt_C2",h_C2_jet80);
  file_PbPb_HardProbes_jet80->GetObject("h_inclRecoJetPt_C1",h_C1_jet80);

  file_PbPb_HardProbes_jet100->GetObject("h_inclRecoJetPt_C4",h_C4_jet100);
  file_PbPb_HardProbes_jet100->GetObject("h_inclRecoJetPt_C3",h_C3_jet100);
  file_PbPb_HardProbes_jet100->GetObject("h_inclRecoJetPt_C2",h_C2_jet100);
  file_PbPb_HardProbes_jet100->GetObject("h_inclRecoJetPt_C1",h_C1_jet100);

  file_fakeJets = TFile::Open(fakeJetsPath);

  file_fakeJets->GetObject("h_fakeJets_C4", h_fakeJets_C4);
  file_fakeJets->GetObject("h_fakeJets_C3", h_fakeJets_C3);
  file_fakeJets->GetObject("h_fakeJets_C2", h_fakeJets_C2);
  file_fakeJets->GetObject("h_fakeJets_C1", h_fakeJets_C1);

  file_PbPb_HardProbes_jet100->GetObject("h_vz_jet100_C4",h_vz_C4);
  file_PbPb_HardProbes_jet100->GetObject("h_vz_jet100_C3",h_vz_C3);
  file_PbPb_HardProbes_jet100->GetObject("h_vz_jet100_C2",h_vz_C2);
  file_PbPb_HardProbes_jet100->GetObject("h_vz_jet100_C1",h_vz_C1);

  return;
  
}


void normalizeMinBiasDataByNEvent(){

  h_C4_jetMB->Scale(1./h_vz_MinBias_C4->Integral());
  h_C3_jetMB->Scale(1./h_vz_MinBias_C3->Integral());
  h_C2_jetMB->Scale(1./h_vz_MinBias_C2->Integral());
  h_C1_jetMB->Scale(1./h_vz_MinBias_C1->Integral());

  return;

}


void cloneMinBiasDataHistogramsForDrawing(){

  h_C1_clone = (TH1D*) h_C1_jetMB->Clone("h_C1_clone");
  h_fakeJets_C1_clone = (TH1D*) h_fakeJets_C1->Clone("h_fakeJets_C1_clone");
  h_C1_sub = (TH1D*) h_C1_jetMB->Clone("h_C1_sub");

  h_C2_clone = (TH1D*) h_C2_jetMB->Clone("h_C2_clone");
  h_fakeJets_C2_clone = (TH1D*) h_fakeJets_C2->Clone("h_fakeJets_C2_clone");
  h_C2_sub = (TH1D*) h_C2_jetMB->Clone("h_C2_sub");

  h_C3_clone = (TH1D*) h_C3_jetMB->Clone("h_C3_clone");
  h_fakeJets_C3_clone = (TH1D*) h_fakeJets_C3->Clone("h_fakeJets_C3_clone");
  h_C3_sub = (TH1D*) h_C3_jetMB->Clone("h_C3_sub");

  h_C4_clone = (TH1D*) h_C4_jetMB->Clone("h_C4_clone");
  h_fakeJets_C4_clone = (TH1D*) h_fakeJets_C4->Clone("h_fakeJets_C4_clone");
  h_C4_sub = (TH1D*) h_C4_jetMB->Clone("h_C4_sub");
  
}

void subtractFakeJetsFromMinBiasDataForDrawing(){

  h_C1_sub->Add(h_fakeJets_C1,-1);
  h_C2_sub->Add(h_fakeJets_C2,-1);
  h_C3_sub->Add(h_fakeJets_C3,-1);
  h_C4_sub->Add(h_fakeJets_C4,-1);

}


void drawFakeJetsVsMinBiasData(){
  
  h_C1_clone->SetLineColor(kRed-4);
  h_fakeJets_C1_clone->SetLineColor(kGray+1);
  h_C1_sub->SetLineColor(kRed-4);
  h_fakeJets_C1_clone->SetLineStyle(1);
  h_C1_sub->SetLineStyle(7);
  h_C1_clone->SetLineWidth(2);
  h_fakeJets_C1_clone->SetLineWidth(2);
  h_C1_sub->SetLineWidth(2);
  r_C1_sub = (TH1D*) h_C1_sub->Clone("r_C1_sub");
  r_C1_sub->Divide(h_C1_sub,h_C1_clone,1,1,"B");
  canv_fakeJets_C1 = new TCanvas("canv_fakeJets_C1","canv_fakeJets_C1",700,700);
  canv_fakeJets_C1->cd();
  pad_fakeJets_C1_upper = new TPad("pad_fakeJets_C1_upper","pad_fakeJets_C1_upper",0,0.4,1,1);
  pad_fakeJets_C1_lower = new TPad("pad_fakeJets_C1_lower","pad_fakeJets_C1_lower",0,0,1,0.4);
  pad_fakeJets_C1_upper->SetLeftMargin(0.15);
  pad_fakeJets_C1_lower->SetLeftMargin(0.15);
  pad_fakeJets_C1_upper->SetBottomMargin(0.);
  pad_fakeJets_C1_lower->SetBottomMargin(0.18);
  pad_fakeJets_C1_upper->SetTopMargin(0.15);
  pad_fakeJets_C1_lower->SetTopMargin(0.);
  pad_fakeJets_C1_upper->SetLogy();
  pad_fakeJets_C1_upper->Draw();
  pad_fakeJets_C1_lower->Draw();
  pad_fakeJets_C1_upper->cd();
  h_C1_clone->SetStats(0);
  h_C1_clone->SetTitle("");
  h_C1_clone->GetYaxis()->SetLabelSize(0.045);
  h_C1_clone->GetYaxis()->SetTitleSize(0.055);
  h_C1_clone->GetYaxis()->SetTitle("#frac{1}{#it{N}^{evt}} #frac{d#it{N}^{jet}}{d#it{p}_{T}} [GeV^{-1}]");
  h_C1_clone->Draw("hist");
  h_fakeJets_C1_clone->Draw("hist same");
  h_C1_sub->Draw("hist same");
  leg_fakeJets_C1 = new TLegend(0.3,0.7,0.7,0.84);
  leg_fakeJets_C1->SetBorderSize(0);
  leg_fakeJets_C1->SetTextSize(0.045);
  leg_fakeJets_C1->AddEntry(h_C1_clone,"PbPb MinBias 0-10%");
  leg_fakeJets_C1->AddEntry(h_fakeJets_C1_clone,"fake-jets, data (mixed-event) 0-10%");
  leg_fakeJets_C1->AddEntry(h_C1_sub,"Corrected PbPb MinBias 0-10%");
  leg_fakeJets_C1->Draw();
  pad_fakeJets_C1_lower->cd();
  r_C1_sub->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r_C1_sub->GetYaxis()->SetTitle("Corrected / Nominal");
  r_C1_sub->GetXaxis()->SetLabelSize(0.06);
  r_C1_sub->GetXaxis()->SetTitleSize(0.07);
  r_C1_sub->GetYaxis()->SetLabelSize(0.06);
  r_C1_sub->GetYaxis()->SetTitleSize(0.07);
  r_C1_sub->SetTitle("");
  r_C1_sub->SetStats(0);
  r_C1_sub->GetYaxis()->SetRangeUser(0,1);
  r_C1_sub->Draw("hist");

  h_C2_clone->SetLineColor(kRed-4);
  h_fakeJets_C2_clone->SetLineColor(kGray+1);
  h_C2_sub->SetLineColor(kRed-4);
  h_fakeJets_C2_clone->SetLineStyle(1);
  h_C2_sub->SetLineStyle(7);
  h_C2_clone->SetLineWidth(2);
  h_fakeJets_C2_clone->SetLineWidth(2);
  h_C2_sub->SetLineWidth(2);
  r_C2_sub = (TH1D*) h_C2_sub->Clone("r_C2_sub");
  r_C2_sub->Divide(h_C2_sub,h_C2_clone,1,1,"B");
  canv_fakeJets_C2 = new TCanvas("canv_fakeJets_C2","canv_fakeJets_C2",700,700);
  canv_fakeJets_C2->cd();
  pad_fakeJets_C2_upper = new TPad("pad_fakeJets_C2_upper","pad_fakeJets_C2_upper",0,0.4,1,1);
  pad_fakeJets_C2_lower = new TPad("pad_fakeJets_C2_lower","pad_fakeJets_C2_lower",0,0,1,0.4);
  pad_fakeJets_C2_upper->SetLeftMargin(0.15);
  pad_fakeJets_C2_lower->SetLeftMargin(0.15);
  pad_fakeJets_C2_upper->SetBottomMargin(0.);
  pad_fakeJets_C2_lower->SetBottomMargin(0.18);
  pad_fakeJets_C2_upper->SetTopMargin(0.15);
  pad_fakeJets_C2_lower->SetTopMargin(0.);
  pad_fakeJets_C2_upper->SetLogy();
  pad_fakeJets_C2_upper->Draw();
  pad_fakeJets_C2_lower->Draw();
  pad_fakeJets_C2_upper->cd();
  h_C2_clone->SetStats(0);
  h_C2_clone->SetTitle("");
  h_C2_clone->GetYaxis()->SetLabelSize(0.045);
  h_C2_clone->GetYaxis()->SetTitleSize(0.055);
  h_C2_clone->GetYaxis()->SetTitle("#frac{1}{#it{N}^{evt}} #frac{d#it{N}^{jet}}{d#it{p}_{T}} [GeV^{-1}]");
  h_C2_clone->Draw("hist");
  h_fakeJets_C2_clone->Draw("hist same");
  h_C2_sub->Draw("hist same");
  leg_fakeJets_C2 = new TLegend(0.3,0.7,0.7,0.84);
  leg_fakeJets_C2->SetBorderSize(0);
  leg_fakeJets_C2->SetTextSize(0.045);
  leg_fakeJets_C2->AddEntry(h_C2_clone,"PbPb MinBias 10-30%");
  leg_fakeJets_C2->AddEntry(h_fakeJets_C2_clone,"fake-jets, data (mixed-event) 10-30%");
  leg_fakeJets_C2->AddEntry(h_C2_sub,"Corrected PbPb MinBias 10-30%");
  leg_fakeJets_C2->Draw();
  pad_fakeJets_C2_lower->cd();
  r_C2_sub->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r_C2_sub->GetYaxis()->SetTitle("Corrected / Nominal");
  r_C2_sub->GetXaxis()->SetLabelSize(0.06);
  r_C2_sub->GetXaxis()->SetTitleSize(0.07);
  r_C2_sub->GetYaxis()->SetLabelSize(0.06);
  r_C2_sub->GetYaxis()->SetTitleSize(0.07);
  r_C2_sub->SetTitle("");
  r_C2_sub->SetStats(0);
  r_C2_sub->Draw("hist");

  h_C3_clone->SetLineColor(kRed-4);
  h_fakeJets_C3_clone->SetLineColor(kGray+1);
  h_C3_sub->SetLineColor(kRed-4);
  h_fakeJets_C3_clone->SetLineStyle(1);
  h_C3_sub->SetLineStyle(7);
  h_C3_clone->SetLineWidth(2);
  h_fakeJets_C3_clone->SetLineWidth(2);
  h_C3_sub->SetLineWidth(2);
  r_C3_sub = (TH1D*) h_C3_sub->Clone("r_C3_sub");
  r_C3_sub->Divide(h_C3_sub,h_C3_clone,1,1,"B");
  canv_fakeJets_C3 = new TCanvas("canv_fakeJets_C3","canv_fakeJets_C3",700,700);
  canv_fakeJets_C3->cd();
  pad_fakeJets_C3_upper = new TPad("pad_fakeJets_C3_upper","pad_fakeJets_C3_upper",0,0.4,1,1);
  pad_fakeJets_C3_lower = new TPad("pad_fakeJets_C3_lower","pad_fakeJets_C3_lower",0,0,1,0.4);
  pad_fakeJets_C3_upper->SetLeftMargin(0.15);
  pad_fakeJets_C3_lower->SetLeftMargin(0.15);
  pad_fakeJets_C3_upper->SetBottomMargin(0.);
  pad_fakeJets_C3_lower->SetBottomMargin(0.18);
  pad_fakeJets_C3_upper->SetTopMargin(0.15);
  pad_fakeJets_C3_lower->SetTopMargin(0.);
  pad_fakeJets_C3_upper->SetLogy();
  pad_fakeJets_C3_upper->Draw();
  pad_fakeJets_C3_lower->Draw();
  pad_fakeJets_C3_upper->cd();
  h_C3_clone->SetStats(0);
  h_C3_clone->SetTitle("");
  h_C3_clone->GetYaxis()->SetLabelSize(0.045);
  h_C3_clone->GetYaxis()->SetTitleSize(0.055);
  h_C3_clone->GetYaxis()->SetTitle("#frac{1}{#it{N}^{evt}} #frac{d#it{N}^{jet}}{d#it{p}_{T}} [GeV^{-1}]");
  h_C3_clone->Draw("hist");
  h_fakeJets_C3_clone->Draw("hist same");
  h_C3_sub->Draw("hist same");
  leg_fakeJets_C3 = new TLegend(0.3,0.7,0.7,0.84);
  leg_fakeJets_C3->SetBorderSize(0);
  leg_fakeJets_C3->SetTextSize(0.045);
  leg_fakeJets_C3->AddEntry(h_C3_clone,"PbPb MinBias 30-50%");
  leg_fakeJets_C3->AddEntry(h_fakeJets_C3_clone,"fake-jets, data (mixed-event) 30-50%");
  leg_fakeJets_C3->AddEntry(h_C3_sub,"Corrected PbPb MinBias 30-50%");
  leg_fakeJets_C3->Draw();
  pad_fakeJets_C3_lower->cd();
  r_C3_sub->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r_C3_sub->GetYaxis()->SetTitle("Corrected / Nominal");
  r_C3_sub->GetXaxis()->SetLabelSize(0.06);
  r_C3_sub->GetXaxis()->SetTitleSize(0.07);
  r_C3_sub->GetYaxis()->SetLabelSize(0.06);
  r_C3_sub->GetYaxis()->SetTitleSize(0.07);
  r_C3_sub->SetTitle("");
  r_C3_sub->SetStats(0);
  r_C3_sub->Draw("hist");

  h_C4_clone->SetLineColor(kRed-4);
  h_fakeJets_C4_clone->SetLineColor(kGray+1);
  h_C4_sub->SetLineColor(kRed-4);
  h_fakeJets_C4_clone->SetLineStyle(1);
  h_C4_sub->SetLineStyle(7);
  h_C4_clone->SetLineWidth(2);
  h_fakeJets_C4_clone->SetLineWidth(2);
  h_C4_sub->SetLineWidth(2);
  r_C4_sub = (TH1D*) h_C4_sub->Clone("r_C4_sub");
  r_C4_sub->Divide(h_C4_sub,h_C4_clone,1,1,"B");
  canv_fakeJets_C4 = new TCanvas("canv_fakeJets_C4","canv_fakeJets_C4",700,700);
  canv_fakeJets_C4->cd();
  pad_fakeJets_C4_upper = new TPad("pad_fakeJets_C4_upper","pad_fakeJets_C4_upper",0,0.4,1,1);
  pad_fakeJets_C4_lower = new TPad("pad_fakeJets_C4_lower","pad_fakeJets_C4_lower",0,0,1,0.4);
  pad_fakeJets_C4_upper->SetLeftMargin(0.15);
  pad_fakeJets_C4_lower->SetLeftMargin(0.15);
  pad_fakeJets_C4_upper->SetBottomMargin(0.);
  pad_fakeJets_C4_lower->SetBottomMargin(0.18);
  pad_fakeJets_C4_upper->SetTopMargin(0.15);
  pad_fakeJets_C4_lower->SetTopMargin(0.);
  pad_fakeJets_C4_upper->SetLogy();
  pad_fakeJets_C4_upper->Draw();
  pad_fakeJets_C4_lower->Draw();
  pad_fakeJets_C4_upper->cd();
  h_C4_clone->SetStats(0);
  h_C4_clone->SetTitle("");
  h_C4_clone->GetYaxis()->SetLabelSize(0.045);
  h_C4_clone->GetYaxis()->SetTitleSize(0.055);
  h_C4_clone->GetYaxis()->SetTitle("#frac{1}{#it{N}^{evt}} #frac{d#it{N}^{jet}}{d#it{p}_{T}} [GeV^{-1}]");
  h_C4_clone->Draw("hist");
  h_fakeJets_C4_clone->Draw("hist same");
  h_C4_sub->Draw("hist same");
  leg_fakeJets_C4 = new TLegend(0.3,0.7,0.7,0.84);
  leg_fakeJets_C4->SetBorderSize(0);
  leg_fakeJets_C4->SetTextSize(0.045);
  leg_fakeJets_C4->AddEntry(h_C4_clone,"PbPb MinBias 50-80%");
  leg_fakeJets_C4->AddEntry(h_fakeJets_C4_clone,"fake-jets, data (mixed-event) 50-80%");
  leg_fakeJets_C4->AddEntry(h_C4_sub,"Corrected PbPb MinBias 50-80%");
  leg_fakeJets_C4->Draw();
  pad_fakeJets_C4_lower->cd();
  r_C4_sub->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r_C4_sub->GetYaxis()->SetTitle("Corrected / Nominal");
  r_C4_sub->GetXaxis()->SetLabelSize(0.06);
  r_C4_sub->GetXaxis()->SetTitleSize(0.07);
  r_C4_sub->GetYaxis()->SetLabelSize(0.06);
  r_C4_sub->GetYaxis()->SetTitleSize(0.07);
  r_C4_sub->SetTitle("");
  r_C4_sub->SetStats(0);
  r_C4_sub->Draw("hist");

  canv_fakeJets_C1->SaveAs("../../figures/JetsPerZ/fakeJets_C1.pdf");
  canv_fakeJets_C2->SaveAs("../../figures/JetsPerZ/fakeJets_C2.pdf");
  canv_fakeJets_C3->SaveAs("../../figures/JetsPerZ/fakeJets_C3.pdf");
  canv_fakeJets_C4->SaveAs("../../figures/JetsPerZ/fakeJets_C4.pdf");


}

TH1D* stitchSamples(TH1D *h_jetMB, TH1D *h_fakeJets, TH1D *h_jet60, TH1D *h_jet80, TH1D *h_jet100,
		    bool isPP, bool isPbPb){

  TH1D *h_return = (TH1D*) h_jet100->Clone("h_return");

  // double jetMB_pTmin = 60.;
  // double jet60_pTmin = 200;
  // double jet80_pTmin = 300;
  // double jet100_pTmin = 500;
  // double smallShift = 0.01;

  double jetMB_pTmin = 60.;
  double jet60_pTmin = 150;
  double jet80_pTmin = 150;
  double jet100_pTmin = 200;
  double smallShift = 0.01;

  if(!isPP){
    jet60_pTmin = jet80_pTmin;
  }




  double N_jet100 = h_jet100->Integral(h_jet100->FindBin(jet100_pTmin + smallShift),h_jet100->FindBin(500. - smallShift));
  double N_jet80 = h_jet80->Integral(h_jet80->FindBin(jet100_pTmin + smallShift),h_jet80->FindBin(500. - smallShift));
  double N_jet60 = h_jet60->Integral(h_jet60->FindBin(jet80_pTmin + smallShift),h_jet60->FindBin(jet100_pTmin - smallShift));

  // For PbPb, jet60_pTmin == jet80_pTmin (Jet60 sample not used), so the
  // [jet60_pTmin, jet80_pTmin] overlap range is empty and gives N_jetMB = 0.
  // Normalize MinBias directly to Jet80.
  //
  // The window used to sit at [100,150], on the claim that Jet80 is ~100%
  // efficient above ~90 GeV.  It is not: the MinBias/Jet80 ratio varies by 2.5x
  // across [100,150] in C1 (1.7x in C2, 1.4x in C3/C4), i.e. the window sat
  // inside the trigger turn-on, and the turn-on is slowest in central events
  // where the background is largest.  That made k too small for C1 in
  // particular, scaled the whole sub-150 GeV region down by a
  // centrality-dependent factor, and inverted the PbPb/pp ordering: at 108 GeV
  // C1 came out BELOW C2 (0.71 vs 0.78) even though the raw per-event MinBias
  // spectra, which involve no stitching at all, give C1 > C2 > C3 > C4.
  //
  // MinBias/Jet80 reaches its plateau around 130 GeV.  Normalizing over
  // [130,170] restores the monotonic ordering and reproduces the shape seen in
  // the stitching-free MinBias comparison: at 88 GeV that gives C1/C2 = 1.41,
  // [130,170] gives 1.33, the old window gave 1.13.  The cost is the MinBias
  // normalization statistics: 0.44% -> 0.96% in C1, 1.60% -> 3.20% in C4.
  double mbNormLo = 130.;
  double mbNormHi = 170.;
  double N_jetMB     = h_jetMB->Integral(h_jetMB->FindBin(mbNormLo + smallShift), h_jetMB->FindBin(mbNormHi - smallShift));
  double N_fakeJets  = h_fakeJets->Integral(h_fakeJets->FindBin(mbNormLo + smallShift), h_fakeJets->FindBin(mbNormHi - smallShift));

  TH1D *h_jet80_scaled = (TH1D*) h_jet80->Clone("h_jet80_scaled");
  h_jet80_scaled->Scale(N_jet100 / N_jet80);

  double N_jet80_scaled = h_jet80_scaled->Integral(h_jet80_scaled->FindBin(jet80_pTmin + smallShift),h_jet80_scaled->FindBin(jet100_pTmin - smallShift));

  TH1D *h_jet60_scaled = (TH1D*) h_jet60->Clone("h_jet60_scaled");
  h_jet60_scaled->Scale(N_jet80_scaled / N_jet60);

  // For PbPb: normalize MinBias directly to Jet80 over [jetMB_pTmin, jet80_pTmin].
  // For pp: normalize MinBias to (scaled) Jet60 over [jet60_pTmin, jet80_pTmin] as before.
  double N_jet80_scaled_mbRange = isPbPb
    ? h_jet80_scaled->Integral(h_jet80_scaled->FindBin(mbNormLo + smallShift), h_jet80_scaled->FindBin(mbNormHi - smallShift))
    : h_jet60_scaled->Integral(h_jet60_scaled->FindBin(mbNormLo + smallShift), h_jet60_scaled->FindBin(mbNormHi - smallShift));

  TH1D *h_jetMB_scaled = (TH1D*) h_jetMB->Clone("h_jetMB_scaled");
  h_jetMB_scaled->Scale(N_jet80_scaled_mbRange / N_jetMB);

  TH1D *h_fakeJets_scaled = (TH1D*) h_fakeJets->Clone("h_fakeJets_scaled");
  h_fakeJets_scaled->Scale(N_jet80_scaled_mbRange / N_jetMB);
    
  // subtract mixed-event FastJet fake jets from MinBias spectrum
  if(isPbPb && doFakeJetSubtraction) h_jetMB_scaled->Add(h_fakeJets_scaled,-1);

  if(useMinBiasOnly){
    // Return the MinBias spectrum over the whole pT range instead of handing the
    // high-pT region to the triggered samples -- but keep the normalisation the
    // stitch just computed for it. That factor is what puts MinBias into
    // jet100-trigger units, so the per-Z scaling downstream stays valid; without
    // it the ratio carries a stray 1/N_evt and is meaningless in absolute terms.
    for(int i = 0; i < h_jetMB_scaled->GetSize(); i++){
      h_return->SetBinContent(i, h_jetMB_scaled->GetBinContent(i));
      h_return->SetBinError  (i, h_jetMB_scaled->GetBinError(i));
    }
    return h_return;
  }
  
  for(int i = 0; i < h_jet60->GetSize(); i++){
    double pT = h_jet60->GetBinCenter(i);
    if(pT < jet60_pTmin){
      h_return->SetBinContent(i,h_jetMB_scaled->GetBinContent(i));
      h_return->SetBinError(i,h_jetMB_scaled->GetBinError(i));
    }
    else if(pT > jet60_pTmin && pT < jet80_pTmin){
      h_return->SetBinContent(i,h_jet60_scaled->GetBinContent(i));
      h_return->SetBinError(i,h_jet60_scaled->GetBinError(i));
    }
    else if(pT > jet80_pTmin && pT < jet100_pTmin){
      h_return->SetBinContent(i,h_jet80_scaled->GetBinContent(i));
      h_return->SetBinError(i,h_jet80_scaled->GetBinError(i));
    }
    else if(pT > jet100_pTmin){
      h_return->SetBinContent(i,h_jet100->GetBinContent(i));
      h_return->SetBinError(i,h_jet100->GetBinError(i));
    }
    else{};
  }

  return h_return;
  
}

void stitchHistograms(){

  h_pp = (TH1D*) stitchSamples(h_pp_jetMB, h_fakeJets_C4, h_pp_jet60, h_pp_jet80, h_pp_jet100, 1, 0);
  h_C4 = (TH1D*) stitchSamples(h_C4_jetMB, h_fakeJets_C4, h_C4_jet60, h_C4_jet80, h_C4_jet100, 0, 1);
  h_C3 = (TH1D*) stitchSamples(h_C3_jetMB, h_fakeJets_C3, h_C3_jet60, h_C3_jet80, h_C3_jet100, 0, 1);
  h_C2 = (TH1D*) stitchSamples(h_C2_jetMB, h_fakeJets_C2, h_C2_jet60, h_C2_jet80, h_C2_jet100, 0, 1);
  h_C1 = (TH1D*) stitchSamples(h_C1_jetMB, h_fakeJets_C1, h_C1_jet60, h_C1_jet80, h_C1_jet100, 0, 1);

  return;

}

void fetchHistogramsMC(){

  file_PYTHIA_response->GetObject("h_matchedRecoJetPt_genJetPt_allJets",h_response_pp);
  file_PH_response->GetObject("h_matchedRecoJetPt_genJetPt_allJets_C4",h_response_C4);
  file_PH_response->GetObject("h_matchedRecoJetPt_genJetPt_allJets_C3",h_response_C3);
  file_PH_response->GetObject("h_matchedRecoJetPt_genJetPt_allJets_C2",h_response_C2);
  file_PH_response->GetObject("h_matchedRecoJetPt_genJetPt_allJets_C1",h_response_C1);

  return;

}

void getMeasuredAndTruthHistograms(){

  h_meas_pp = (TH1D*) h_response_pp->ProjectionX("h_meas_pp");
  h_meas_C4 = (TH1D*) h_response_C4->ProjectionX("h_meas_C4");
  h_meas_C3 = (TH1D*) h_response_C3->ProjectionX("h_meas_C3");
  h_meas_C2 = (TH1D*) h_response_C2->ProjectionX("h_meas_C2");
  h_meas_C1 = (TH1D*) h_response_C1->ProjectionX("h_meas_C1");

  h_truth_pp = (TH1D*) h_response_pp->ProjectionY("h_truth_pp");
  h_truth_C4 = (TH1D*) h_response_C4->ProjectionY("h_truth_C4");
  h_truth_C3 = (TH1D*) h_response_C3->ProjectionY("h_truth_C3");
  h_truth_C2 = (TH1D*) h_response_C2->ProjectionY("h_truth_C2");
  h_truth_C1 = (TH1D*) h_response_C1->ProjectionY("h_truth_C1");

  return;

}

void defineRooUnfoldObjects(){

  response_pp = new RooUnfoldResponse(h_meas_pp, h_truth_pp, h_response_pp, "response_pp", "pp response");
  response_C4 = new RooUnfoldResponse(h_meas_C4, h_truth_C4, h_response_C4, "response_C4", "C4 response");
  response_C3 = new RooUnfoldResponse(h_meas_C3, h_truth_C3, h_response_C3, "response_C3", "C3 response");
  response_C2 = new RooUnfoldResponse(h_meas_C2, h_truth_C2, h_response_C2, "response_C2", "C2 response");
  response_C1 = new RooUnfoldResponse(h_meas_C1, h_truth_C1, h_response_C1, "response_C1", "C1 response");

  unfold_pp = new RooUnfoldBayes(response_pp, h_pp, N_iter_pp);
  unfold_C4 = new RooUnfoldBayes(response_C4, h_C4, N_iter_C4);
  unfold_C3 = new RooUnfoldBayes(response_C3, h_C3, N_iter_C3);
  unfold_C2 = new RooUnfoldBayes(response_C2, h_C2, N_iter_C2);
  unfold_C1 = new RooUnfoldBayes(response_C1, h_C1, N_iter_C1);

  return;
  
}

void unfoldSpectra(){

  if(doUnfolding){
    h_pp_unfold = (TH1D*) unfold_pp->Hunfold();
    h_C4_unfold = (TH1D*) unfold_C4->Hunfold();
    h_C3_unfold = (TH1D*) unfold_C3->Hunfold();
    h_C2_unfold = (TH1D*) unfold_C2->Hunfold();
    h_C1_unfold = (TH1D*) unfold_C1->Hunfold();
  }
  else{
    h_pp_unfold = (TH1D*) h_pp->Clone("h_pp_noUnfold");
    h_C4_unfold = (TH1D*) h_C4->Clone("h_C4_noUnfold");
    h_C3_unfold = (TH1D*) h_C3->Clone("h_C3_noUnfold");
    h_C2_unfold = (TH1D*) h_C2->Clone("h_C2_noUnfold");
    h_C1_unfold = (TH1D*) h_C1->Clone("h_C1_noUnfold");
  }

  stylizeHistograms(h_pp_unfold,h_C1_unfold,h_C2_unfold,h_C3_unfold,h_C4_unfold);

  return;

}

void rebinUnfoldedHistograms(){

  h_pp_unfold = (TH1D*) h_pp_unfold->Rebin(N_edge-1,"h_pp_unfold",newAxis);
  h_C4_unfold = (TH1D*) h_C4_unfold->Rebin(N_edge-1,"h_C4_unfold",newAxis);
  h_C3_unfold = (TH1D*) h_C3_unfold->Rebin(N_edge-1,"h_C3_unfold",newAxis);
  h_C2_unfold = (TH1D*) h_C2_unfold->Rebin(N_edge-1,"h_C2_unfold",newAxis);
  h_C1_unfold = (TH1D*) h_C1_unfold->Rebin(N_edge-1,"h_C1_unfold",newAxis);

  return;
  
}

void normalizeUnfoldedHistogramsByPTAndEtaBinWidth(){
  
   divideByBinwidth(h_pp_unfold); // pT normalization
   h_pp_unfold->Scale(1./3.2); // eta normalization

   divideByBinwidth(h_C4_unfold);
   h_C4_unfold->Scale(1./3.2);

   divideByBinwidth(h_C3_unfold);
   h_C3_unfold->Scale(1./3.2);

   divideByBinwidth(h_C2_unfold);
   h_C2_unfold->Scale(1./3.2);

   divideByBinwidth(h_C1_unfold);
   h_C1_unfold->Scale(1./3.2);

   return;

}

void calculateZBosonYield(){

  NZ_pp = h_dimuonMass_pp->Integral(h_dimuonMass_pp->GetXaxis()->FindBin(75),h_dimuonMass_pp->GetXaxis()->FindBin(105));
  NZ_C4 = h_dimuonMass_C4->Integral(h_dimuonMass_C4->GetXaxis()->FindBin(75),h_dimuonMass_C4->GetXaxis()->FindBin(105));
  NZ_C3 = h_dimuonMass_C3->Integral(h_dimuonMass_C3->GetXaxis()->FindBin(75),h_dimuonMass_C3->GetXaxis()->FindBin(105));
  NZ_C2 = h_dimuonMass_C2->Integral(h_dimuonMass_C2->GetXaxis()->FindBin(75),h_dimuonMass_C2->GetXaxis()->FindBin(105));
  NZ_C1 = h_dimuonMass_C1->Integral(h_dimuonMass_C1->GetXaxis()->FindBin(75),h_dimuonMass_C1->GetXaxis()->FindBin(105));

  // scale by muon reconstruction efficiency
  NZ_pp = NZ_pp/(0.9708*0.9708);
  NZ_C4 = NZ_C4/(0.9778*0.9778);
  NZ_C3 = NZ_C3/(0.9856*0.9856);
  NZ_C2 = NZ_C2/(0.9069*0.9069);
  NZ_C1 = NZ_C1/(0.8627*0.8627);

  return;

}

void normalizeUnfoldedHistogramsByZBosonYield(){

  h_pp_unfold->Scale(1./NZ_pp);
  h_C4_unfold->Scale(1./NZ_C4);
  h_C3_unfold->Scale(1./NZ_C3);
  h_C2_unfold->Scale(1./NZ_C2);
  h_C1_unfold->Scale(1./NZ_C1);

  return;

}

void calculateUnfoldedPbPbToPPRatio(){

  r_C4 = (TH1D*) h_C4_unfold->Clone("r_C4");
  r_C4->Divide(h_C4_unfold,h_pp_unfold,1,1,"");

  r_C3 = (TH1D*) h_C3_unfold->Clone("r_C3");
  r_C3->Divide(h_C3_unfold,h_pp_unfold,1,1,"");

  r_C2 = (TH1D*) h_C2_unfold->Clone("r_C2");
  r_C2->Divide(h_C2_unfold,h_pp_unfold,1,1,"");

  r_C1 = (TH1D*) h_C1_unfold->Clone("r_C1");
  r_C1->Divide(h_C1_unfold,h_pp_unfold,1,1,"");

  return;

}

void calculateSystematicFromZBosonYield(){

  R_C4 = (TH1D*) r_C4->Clone("R_C4");
  R_C3 = (TH1D*) r_C3->Clone("R_C3");
  R_C2 = (TH1D*) r_C2->Clone("R_C2");
  R_C1 = (TH1D*) r_C1->Clone("R_C1");

  for(int i = 0; i < R_C4->GetSize(); i++){
    R_C4->SetBinError(i,1./TMath::Sqrt(NZ_C4));
    R_C3->SetBinError(i,1./TMath::Sqrt(NZ_C3));
    R_C2->SetBinError(i,1./TMath::Sqrt(NZ_C2));
    R_C1->SetBinError(i,1./TMath::Sqrt(NZ_C1));
  }

  stylizeSystematics(R_C1,R_C2,R_C3,R_C4);

  return;

}

void drawJetsPerZ(){
  
  leg = new TLegend(0.3,0.7,0.68,0.88);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_C4_unfold,"PbPb 50-80%","p");  
  leg->AddEntry(h_C3_unfold,"PbPb 30-50%","p");
  leg->AddEntry(h_C2_unfold,"PbPb 10-30%","p");
  leg->AddEntry(h_C1_unfold,"PbPb 0-10%","p");


  canv_prime = new TCanvas("canv_prime","canv_prime",700,700);
  canv_prime->cd();
  pad_prime = new TPad("pad_prime","pad_prime",0,0,1,1);
  pad_prime->SetLeftMargin(0.2);
  pad_prime->SetBottomMargin(0.15);
  //  pad->SetLogy();
  //pad->SetLogx();
  pad_prime->Draw();
  pad_prime->cd();
  R_C1->SetTitle("");
  R_C1->SetStats(0);
  R_C1->GetYaxis()->SetTitleOffset(2.0);
  //R_C1->GetYaxis()->SetTitle("#it{R}_{AA}^{jet}");
  R_C1->GetYaxis()->SetTitle("#frac{1}{#it{N}_{Z}^{PbPb}} #frac{d#it{N}_{jet}^{PbPb}}{d#it{p}_{T}} #scale[3.0]{/} #frac{1}{#it{N}_{Z}^{pp}} #frac{d#it{N}_{jet}^{pp}}{d#it{p}_{T}}");
  R_C1->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  //R_C1->GetYaxis()->SetRangeUser(0,1.);
  R_C1->GetYaxis()->SetRangeUser(0,2.);
  R_C1->GetXaxis()->SetRangeUser(newAxis[0],newAxis[N_edge-1]);
  R_C1->Draw("e2");
  r_C1->Draw("same");
  R_C2->Draw("e2 same");
  r_C2->Draw("same");
  R_C3->Draw("e2 same");
  r_C3->Draw("same");
  R_C4->Draw("e2 same");
  r_C4->Draw("same");
  TLine *li = new TLine();
  li->SetLineStyle(7);
  li->DrawLine(newAxis[0],1,newAxis[N_edge-1],1);
  leg->Draw();

  canv_prime->SaveAs("../../figures/JetsPerZ/JetsPerZ_lightJets.pdf");

  return;
  
}

int main(){

  openScanFiles();
  openPubDataFiles();
  buildPubDataHistograms();
  fetchHistograms();
  normalizeMinBiasDataByNEvent();
  cloneMinBiasDataHistogramsForDrawing();
  subtractFakeJetsFromMinBiasDataForDrawing();
  drawFakeJetsVsMinBiasData();
  stitchHistograms();
  fetchHistogramsMC();
  getMeasuredAndTruthHistograms();
  defineRooUnfoldObjects();
  unfoldSpectra();
  rebinUnfoldedHistograms();
  normalizeUnfoldedHistogramsByPTAndEtaBinWidth();
  calculateZBosonYield();
  normalizeUnfoldedHistogramsByZBosonYield();
  calculateUnfoldedPbPbToPPRatio();
  calculateSystematicFromZBosonYield();
  drawJetsPerZ();

  return -1;

}
