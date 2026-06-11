
#include "../../../headers/functions/divideByBinwidth.h"

void stylizeHistograms(TH1D *h1, TH1D *h2, TH1D *h3){

  h1->SetLineColor(kBlack);
  h2->SetLineColor(kRed-4);
  h3->SetLineColor(kBlue-4);

  h1->SetLineWidth(2.0);
  h2->SetLineWidth(2.0);
  h3->SetLineWidth(2.0);


}



void subleadingPFCand(){

  TFile *f_PFPT_10 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/HYDJET/canonical/pfCandAnalyzer/HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_2026-6-4_ultraFineCentBins_pfCandAnalyzer_PFPT-10.root");

  TFile *f_PFPT_15 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/HYDJET/canonical/pfCandAnalyzer/HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_2026-6-4_ultraFineCentBins_pfCandAnalyzer_PFPT-15.root");

  TFile *f_PFPT_20 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/HYDJET/canonical/pfCandAnalyzer/HYDJET_pThat-unweighted_mu12_pTmu-15to999_tight_hiBinShift-0_jetTrkMaxFilter_WDecayFilter_2026-6-4_ultraFineCentBins_pfCandAnalyzer_PFPT-20.root");

  TH2D *H_recoJetPt_PFPT_10_C16, *H_recoJetPt_PFPT_10_C15, *H_recoJetPt_PFPT_10_C14, *H_recoJetPt_PFPT_10_C13, *H_recoJetPt_PFPT_10_C12, *H_recoJetPt_PFPT_10_C11, *H_recoJetPt_PFPT_10_C10, *H_recoJetPt_PFPT_10_C9, *H_recoJetPt_PFPT_10_C8, *H_recoJetPt_PFPT_10_C7, *H_recoJetPt_PFPT_10_C6, *H_recoJetPt_PFPT_10_C5, *H_recoJetPt_PFPT_10_C4, *H_recoJetPt_PFPT_10_C3, *H_recoJetPt_PFPT_10_C2, *H_recoJetPt_PFPT_10_C1;
  TH2D *H_recoJetPt_PFPT_15_C16, *H_recoJetPt_PFPT_15_C15, *H_recoJetPt_PFPT_15_C14, *H_recoJetPt_PFPT_15_C13, *H_recoJetPt_PFPT_15_C12, *H_recoJetPt_PFPT_15_C11, *H_recoJetPt_PFPT_15_C10, *H_recoJetPt_PFPT_15_C9, *H_recoJetPt_PFPT_15_C8, *H_recoJetPt_PFPT_15_C7, *H_recoJetPt_PFPT_15_C6, *H_recoJetPt_PFPT_15_C5, *H_recoJetPt_PFPT_15_C4, *H_recoJetPt_PFPT_15_C3, *H_recoJetPt_PFPT_15_C2, *H_recoJetPt_PFPT_15_C1;
  TH2D *H_recoJetPt_PFPT_20_C16, *H_recoJetPt_PFPT_20_C15, *H_recoJetPt_PFPT_20_C14, *H_recoJetPt_PFPT_20_C13, *H_recoJetPt_PFPT_20_C12, *H_recoJetPt_PFPT_20_C11, *H_recoJetPt_PFPT_20_C10, *H_recoJetPt_PFPT_20_C9, *H_recoJetPt_PFPT_20_C8, *H_recoJetPt_PFPT_20_C7, *H_recoJetPt_PFPT_20_C6, *H_recoJetPt_PFPT_20_C5, *H_recoJetPt_PFPT_20_C4, *H_recoJetPt_PFPT_20_C3, *H_recoJetPt_PFPT_20_C2, *H_recoJetPt_PFPT_20_C1;

  TH2D *H_recoJetPtNoSubleading_PFPT_10_C16, *H_recoJetPtNoSubleading_PFPT_10_C15, *H_recoJetPtNoSubleading_PFPT_10_C14, *H_recoJetPtNoSubleading_PFPT_10_C13, *H_recoJetPtNoSubleading_PFPT_10_C12, *H_recoJetPtNoSubleading_PFPT_10_C11, *H_recoJetPtNoSubleading_PFPT_10_C10, *H_recoJetPtNoSubleading_PFPT_10_C9, *H_recoJetPtNoSubleading_PFPT_10_C8, *H_recoJetPtNoSubleading_PFPT_10_C7, *H_recoJetPtNoSubleading_PFPT_10_C6, *H_recoJetPtNoSubleading_PFPT_10_C5, *H_recoJetPtNoSubleading_PFPT_10_C4, *H_recoJetPtNoSubleading_PFPT_10_C3, *H_recoJetPtNoSubleading_PFPT_10_C2, *H_recoJetPtNoSubleading_PFPT_10_C1;
  TH2D *H_recoJetPtNoSubleading_PFPT_15_C16, *H_recoJetPtNoSubleading_PFPT_15_C15, *H_recoJetPtNoSubleading_PFPT_15_C14, *H_recoJetPtNoSubleading_PFPT_15_C13, *H_recoJetPtNoSubleading_PFPT_15_C12, *H_recoJetPtNoSubleading_PFPT_15_C11, *H_recoJetPtNoSubleading_PFPT_15_C10, *H_recoJetPtNoSubleading_PFPT_15_C9, *H_recoJetPtNoSubleading_PFPT_15_C8, *H_recoJetPtNoSubleading_PFPT_15_C7, *H_recoJetPtNoSubleading_PFPT_15_C6, *H_recoJetPtNoSubleading_PFPT_15_C5, *H_recoJetPtNoSubleading_PFPT_15_C4, *H_recoJetPtNoSubleading_PFPT_15_C3, *H_recoJetPtNoSubleading_PFPT_15_C2, *H_recoJetPtNoSubleading_PFPT_15_C1;
  TH2D *H_recoJetPtNoSubleading_PFPT_20_C16, *H_recoJetPtNoSubleading_PFPT_20_C15, *H_recoJetPtNoSubleading_PFPT_20_C14, *H_recoJetPtNoSubleading_PFPT_20_C13, *H_recoJetPtNoSubleading_PFPT_20_C12, *H_recoJetPtNoSubleading_PFPT_20_C11, *H_recoJetPtNoSubleading_PFPT_20_C10, *H_recoJetPtNoSubleading_PFPT_20_C9, *H_recoJetPtNoSubleading_PFPT_20_C8, *H_recoJetPtNoSubleading_PFPT_20_C7, *H_recoJetPtNoSubleading_PFPT_20_C6, *H_recoJetPtNoSubleading_PFPT_20_C5, *H_recoJetPtNoSubleading_PFPT_20_C4, *H_recoJetPtNoSubleading_PFPT_20_C3, *H_recoJetPtNoSubleading_PFPT_20_C2, *H_recoJetPtNoSubleading_PFPT_20_C1;

  TH1D *h_recoJetPt_PFPT_10_C16, *h_recoJetPt_PFPT_10_C15, *h_recoJetPt_PFPT_10_C14, *h_recoJetPt_PFPT_10_C13, *h_recoJetPt_PFPT_10_C12, *h_recoJetPt_PFPT_10_C11, *h_recoJetPt_PFPT_10_C10, *h_recoJetPt_PFPT_10_C9, *h_recoJetPt_PFPT_10_C8, *h_recoJetPt_PFPT_10_C7, *h_recoJetPt_PFPT_10_C6, *h_recoJetPt_PFPT_10_C5, *h_recoJetPt_PFPT_10_C4, *h_recoJetPt_PFPT_10_C3, *h_recoJetPt_PFPT_10_C2, *h_recoJetPt_PFPT_10_C1;
  TH1D *h_recoJetPt_PFPT_15_C16, *h_recoJetPt_PFPT_15_C15, *h_recoJetPt_PFPT_15_C14, *h_recoJetPt_PFPT_15_C13, *h_recoJetPt_PFPT_15_C12, *h_recoJetPt_PFPT_15_C11, *h_recoJetPt_PFPT_15_C10, *h_recoJetPt_PFPT_15_C9, *h_recoJetPt_PFPT_15_C8, *h_recoJetPt_PFPT_15_C7, *h_recoJetPt_PFPT_15_C6, *h_recoJetPt_PFPT_15_C5, *h_recoJetPt_PFPT_15_C4, *h_recoJetPt_PFPT_15_C3, *h_recoJetPt_PFPT_15_C2, *h_recoJetPt_PFPT_15_C1;
  TH1D *h_recoJetPt_PFPT_20_C16, *h_recoJetPt_PFPT_20_C15, *h_recoJetPt_PFPT_20_C14, *h_recoJetPt_PFPT_20_C13, *h_recoJetPt_PFPT_20_C12, *h_recoJetPt_PFPT_20_C11, *h_recoJetPt_PFPT_20_C10, *h_recoJetPt_PFPT_20_C9, *h_recoJetPt_PFPT_20_C8, *h_recoJetPt_PFPT_20_C7, *h_recoJetPt_PFPT_20_C6, *h_recoJetPt_PFPT_20_C5, *h_recoJetPt_PFPT_20_C4, *h_recoJetPt_PFPT_20_C3, *h_recoJetPt_PFPT_20_C2, *h_recoJetPt_PFPT_20_C1;

  TH1D *h_recoJetPtNoSubleading_PFPT_10_C16, *h_recoJetPtNoSubleading_PFPT_10_C15, *h_recoJetPtNoSubleading_PFPT_10_C14, *h_recoJetPtNoSubleading_PFPT_10_C13, *h_recoJetPtNoSubleading_PFPT_10_C12, *h_recoJetPtNoSubleading_PFPT_10_C11, *h_recoJetPtNoSubleading_PFPT_10_C10, *h_recoJetPtNoSubleading_PFPT_10_C9, *h_recoJetPtNoSubleading_PFPT_10_C8, *h_recoJetPtNoSubleading_PFPT_10_C7, *h_recoJetPtNoSubleading_PFPT_10_C6, *h_recoJetPtNoSubleading_PFPT_10_C5, *h_recoJetPtNoSubleading_PFPT_10_C4, *h_recoJetPtNoSubleading_PFPT_10_C3, *h_recoJetPtNoSubleading_PFPT_10_C2, *h_recoJetPtNoSubleading_PFPT_10_C1;
  TH1D *h_recoJetPtNoSubleading_PFPT_15_C16, *h_recoJetPtNoSubleading_PFPT_15_C15, *h_recoJetPtNoSubleading_PFPT_15_C14, *h_recoJetPtNoSubleading_PFPT_15_C13, *h_recoJetPtNoSubleading_PFPT_15_C12, *h_recoJetPtNoSubleading_PFPT_15_C11, *h_recoJetPtNoSubleading_PFPT_15_C10, *h_recoJetPtNoSubleading_PFPT_15_C9, *h_recoJetPtNoSubleading_PFPT_15_C8, *h_recoJetPtNoSubleading_PFPT_15_C7, *h_recoJetPtNoSubleading_PFPT_15_C6, *h_recoJetPtNoSubleading_PFPT_15_C5, *h_recoJetPtNoSubleading_PFPT_15_C4, *h_recoJetPtNoSubleading_PFPT_15_C3, *h_recoJetPtNoSubleading_PFPT_15_C2, *h_recoJetPtNoSubleading_PFPT_15_C1;
  TH1D *h_recoJetPtNoSubleading_PFPT_20_C16, *h_recoJetPtNoSubleading_PFPT_20_C15, *h_recoJetPtNoSubleading_PFPT_20_C14, *h_recoJetPtNoSubleading_PFPT_20_C13, *h_recoJetPtNoSubleading_PFPT_20_C12, *h_recoJetPtNoSubleading_PFPT_20_C11, *h_recoJetPtNoSubleading_PFPT_20_C10, *h_recoJetPtNoSubleading_PFPT_20_C9, *h_recoJetPtNoSubleading_PFPT_20_C8, *h_recoJetPtNoSubleading_PFPT_20_C7, *h_recoJetPtNoSubleading_PFPT_20_C6, *h_recoJetPtNoSubleading_PFPT_20_C5, *h_recoJetPtNoSubleading_PFPT_20_C4, *h_recoJetPtNoSubleading_PFPT_20_C3, *h_recoJetPtNoSubleading_PFPT_20_C2, *h_recoJetPtNoSubleading_PFPT_20_C1;

  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C16",H_recoJetPt_PFPT_10_C16);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C15",H_recoJetPt_PFPT_10_C15);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C14",H_recoJetPt_PFPT_10_C14);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C13",H_recoJetPt_PFPT_10_C13);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C12",H_recoJetPt_PFPT_10_C12);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C11",H_recoJetPt_PFPT_10_C11);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C10",H_recoJetPt_PFPT_10_C10);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C9",H_recoJetPt_PFPT_10_C9);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C8",H_recoJetPt_PFPT_10_C8);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C7",H_recoJetPt_PFPT_10_C7);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C6",H_recoJetPt_PFPT_10_C6);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C5",H_recoJetPt_PFPT_10_C5);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C4",H_recoJetPt_PFPT_10_C4);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C3",H_recoJetPt_PFPT_10_C3);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C2",H_recoJetPt_PFPT_10_C2);
  f_PFPT_10->GetObject("h_inclRecoJetPt_flavor_C1",H_recoJetPt_PFPT_10_C1);

  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C16",H_recoJetPtNoSubleading_PFPT_10_C16);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C15",H_recoJetPtNoSubleading_PFPT_10_C15);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C14",H_recoJetPtNoSubleading_PFPT_10_C14);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C13",H_recoJetPtNoSubleading_PFPT_10_C13);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C12",H_recoJetPtNoSubleading_PFPT_10_C12);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C11",H_recoJetPtNoSubleading_PFPT_10_C11);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C10",H_recoJetPtNoSubleading_PFPT_10_C10);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C9",H_recoJetPtNoSubleading_PFPT_10_C9);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C8",H_recoJetPtNoSubleading_PFPT_10_C8);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C7",H_recoJetPtNoSubleading_PFPT_10_C7);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C6",H_recoJetPtNoSubleading_PFPT_10_C6);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C5",H_recoJetPtNoSubleading_PFPT_10_C5);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C4",H_recoJetPtNoSubleading_PFPT_10_C4);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C3",H_recoJetPtNoSubleading_PFPT_10_C3);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C2",H_recoJetPtNoSubleading_PFPT_10_C2);
  f_PFPT_10->GetObject("h_inclRecoJetPtNoSubleading_flavor_C1",H_recoJetPtNoSubleading_PFPT_10_C1);

  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C16",H_recoJetPt_PFPT_15_C16);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C15",H_recoJetPt_PFPT_15_C15);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C14",H_recoJetPt_PFPT_15_C14);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C13",H_recoJetPt_PFPT_15_C13);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C12",H_recoJetPt_PFPT_15_C12);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C11",H_recoJetPt_PFPT_15_C11);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C10",H_recoJetPt_PFPT_15_C10);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C9",H_recoJetPt_PFPT_15_C9);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C8",H_recoJetPt_PFPT_15_C8);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C7",H_recoJetPt_PFPT_15_C7);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C6",H_recoJetPt_PFPT_15_C6);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C5",H_recoJetPt_PFPT_15_C5);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C4",H_recoJetPt_PFPT_15_C4);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C3",H_recoJetPt_PFPT_15_C3);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C2",H_recoJetPt_PFPT_15_C2);
  f_PFPT_15->GetObject("h_inclRecoJetPt_flavor_C1",H_recoJetPt_PFPT_15_C1);

  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C16",H_recoJetPtNoSubleading_PFPT_15_C16);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C15",H_recoJetPtNoSubleading_PFPT_15_C15);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C14",H_recoJetPtNoSubleading_PFPT_15_C14);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C13",H_recoJetPtNoSubleading_PFPT_15_C13);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C12",H_recoJetPtNoSubleading_PFPT_15_C12);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C11",H_recoJetPtNoSubleading_PFPT_15_C11);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C10",H_recoJetPtNoSubleading_PFPT_15_C10);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C9",H_recoJetPtNoSubleading_PFPT_15_C9);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C8",H_recoJetPtNoSubleading_PFPT_15_C8);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C7",H_recoJetPtNoSubleading_PFPT_15_C7);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C6",H_recoJetPtNoSubleading_PFPT_15_C6);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C5",H_recoJetPtNoSubleading_PFPT_15_C5);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C4",H_recoJetPtNoSubleading_PFPT_15_C4);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C3",H_recoJetPtNoSubleading_PFPT_15_C3);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C2",H_recoJetPtNoSubleading_PFPT_15_C2);
  f_PFPT_15->GetObject("h_inclRecoJetPtNoSubleading_flavor_C1",H_recoJetPtNoSubleading_PFPT_15_C1);

  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C16",H_recoJetPt_PFPT_20_C16);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C15",H_recoJetPt_PFPT_20_C15);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C14",H_recoJetPt_PFPT_20_C14);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C13",H_recoJetPt_PFPT_20_C13);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C12",H_recoJetPt_PFPT_20_C12);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C11",H_recoJetPt_PFPT_20_C11);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C10",H_recoJetPt_PFPT_20_C10);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C9",H_recoJetPt_PFPT_20_C9);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C8",H_recoJetPt_PFPT_20_C8);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C7",H_recoJetPt_PFPT_20_C7);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C6",H_recoJetPt_PFPT_20_C6);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C5",H_recoJetPt_PFPT_20_C5);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C4",H_recoJetPt_PFPT_20_C4);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C3",H_recoJetPt_PFPT_20_C3);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C2",H_recoJetPt_PFPT_20_C2);
  f_PFPT_20->GetObject("h_inclRecoJetPt_flavor_C1",H_recoJetPt_PFPT_20_C1);

  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C16",H_recoJetPtNoSubleading_PFPT_20_C16);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C15",H_recoJetPtNoSubleading_PFPT_20_C15);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C14",H_recoJetPtNoSubleading_PFPT_20_C14);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C13",H_recoJetPtNoSubleading_PFPT_20_C13);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C12",H_recoJetPtNoSubleading_PFPT_20_C12);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C11",H_recoJetPtNoSubleading_PFPT_20_C11);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C10",H_recoJetPtNoSubleading_PFPT_20_C10);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C9",H_recoJetPtNoSubleading_PFPT_20_C9);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C8",H_recoJetPtNoSubleading_PFPT_20_C8);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C7",H_recoJetPtNoSubleading_PFPT_20_C7);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C6",H_recoJetPtNoSubleading_PFPT_20_C6);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C5",H_recoJetPtNoSubleading_PFPT_20_C5);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C4",H_recoJetPtNoSubleading_PFPT_20_C4);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C3",H_recoJetPtNoSubleading_PFPT_20_C3);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C2",H_recoJetPtNoSubleading_PFPT_20_C2);
  f_PFPT_20->GetObject("h_inclRecoJetPtNoSubleading_flavor_C1",H_recoJetPtNoSubleading_PFPT_20_C1);
  

  h_recoJetPt_PFPT_10_C16 = (TH1D*) H_recoJetPt_PFPT_10_C16->ProjectionX("h_recoJetPt_PFPT_10_C16");
  h_recoJetPt_PFPT_10_C15 = (TH1D*) H_recoJetPt_PFPT_10_C15->ProjectionX("h_recoJetPt_PFPT_10_C15");
  h_recoJetPt_PFPT_10_C14 = (TH1D*) H_recoJetPt_PFPT_10_C14->ProjectionX("h_recoJetPt_PFPT_10_C14");
  h_recoJetPt_PFPT_10_C13 = (TH1D*) H_recoJetPt_PFPT_10_C13->ProjectionX("h_recoJetPt_PFPT_10_C13");
  h_recoJetPt_PFPT_10_C12 = (TH1D*) H_recoJetPt_PFPT_10_C12->ProjectionX("h_recoJetPt_PFPT_10_C12");
  h_recoJetPt_PFPT_10_C11 = (TH1D*) H_recoJetPt_PFPT_10_C11->ProjectionX("h_recoJetPt_PFPT_10_C11");
  h_recoJetPt_PFPT_10_C10 = (TH1D*) H_recoJetPt_PFPT_10_C10->ProjectionX("h_recoJetPt_PFPT_10_C10");
  h_recoJetPt_PFPT_10_C9 = (TH1D*) H_recoJetPt_PFPT_10_C9->ProjectionX("h_recoJetPt_PFPT_10_C9");
  h_recoJetPt_PFPT_10_C8 = (TH1D*) H_recoJetPt_PFPT_10_C8->ProjectionX("h_recoJetPt_PFPT_10_C8");
  h_recoJetPt_PFPT_10_C7 = (TH1D*) H_recoJetPt_PFPT_10_C7->ProjectionX("h_recoJetPt_PFPT_10_C7");
  h_recoJetPt_PFPT_10_C6 = (TH1D*) H_recoJetPt_PFPT_10_C6->ProjectionX("h_recoJetPt_PFPT_10_C6");
  h_recoJetPt_PFPT_10_C5 = (TH1D*) H_recoJetPt_PFPT_10_C5->ProjectionX("h_recoJetPt_PFPT_10_C5");
  h_recoJetPt_PFPT_10_C4 = (TH1D*) H_recoJetPt_PFPT_10_C4->ProjectionX("h_recoJetPt_PFPT_10_C4");
  h_recoJetPt_PFPT_10_C3 = (TH1D*) H_recoJetPt_PFPT_10_C3->ProjectionX("h_recoJetPt_PFPT_10_C3");
  h_recoJetPt_PFPT_10_C2 = (TH1D*) H_recoJetPt_PFPT_10_C2->ProjectionX("h_recoJetPt_PFPT_10_C2");
  h_recoJetPt_PFPT_10_C1 = (TH1D*) H_recoJetPt_PFPT_10_C1->ProjectionX("h_recoJetPt_PFPT_10_C1");

  h_recoJetPtNoSubleading_PFPT_10_C16 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C16->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C16");
  h_recoJetPtNoSubleading_PFPT_10_C15 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C15->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C15");
  h_recoJetPtNoSubleading_PFPT_10_C14 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C14->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C14");
  h_recoJetPtNoSubleading_PFPT_10_C13 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C13->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C13");
  h_recoJetPtNoSubleading_PFPT_10_C12 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C12->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C12");
  h_recoJetPtNoSubleading_PFPT_10_C11 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C11->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C11");
  h_recoJetPtNoSubleading_PFPT_10_C10 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C10->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C10");
  h_recoJetPtNoSubleading_PFPT_10_C9 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C9->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C9");
  h_recoJetPtNoSubleading_PFPT_10_C8 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C8->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C8");
  h_recoJetPtNoSubleading_PFPT_10_C7 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C7->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C7");
  h_recoJetPtNoSubleading_PFPT_10_C6 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C6->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C6");
  h_recoJetPtNoSubleading_PFPT_10_C5 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C5->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C5");
  h_recoJetPtNoSubleading_PFPT_10_C4 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C4->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C4");
  h_recoJetPtNoSubleading_PFPT_10_C3 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C3->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C3");
  h_recoJetPtNoSubleading_PFPT_10_C2 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C2->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C2");
  h_recoJetPtNoSubleading_PFPT_10_C1 = (TH1D*) H_recoJetPtNoSubleading_PFPT_10_C1->ProjectionX("h_recoJetPtNoSubleading_PFPT_10_C1");


  h_recoJetPt_PFPT_15_C16 = (TH1D*) H_recoJetPt_PFPT_15_C16->ProjectionX("h_recoJetPt_PFPT_15_C16");
  h_recoJetPt_PFPT_15_C15 = (TH1D*) H_recoJetPt_PFPT_15_C15->ProjectionX("h_recoJetPt_PFPT_15_C15");
  h_recoJetPt_PFPT_15_C14 = (TH1D*) H_recoJetPt_PFPT_15_C14->ProjectionX("h_recoJetPt_PFPT_15_C14");
  h_recoJetPt_PFPT_15_C13 = (TH1D*) H_recoJetPt_PFPT_15_C13->ProjectionX("h_recoJetPt_PFPT_15_C13");
  h_recoJetPt_PFPT_15_C12 = (TH1D*) H_recoJetPt_PFPT_15_C12->ProjectionX("h_recoJetPt_PFPT_15_C12");
  h_recoJetPt_PFPT_15_C11 = (TH1D*) H_recoJetPt_PFPT_15_C11->ProjectionX("h_recoJetPt_PFPT_15_C11");
  h_recoJetPt_PFPT_15_C10 = (TH1D*) H_recoJetPt_PFPT_15_C10->ProjectionX("h_recoJetPt_PFPT_15_C10");
  h_recoJetPt_PFPT_15_C9 = (TH1D*) H_recoJetPt_PFPT_15_C9->ProjectionX("h_recoJetPt_PFPT_15_C9");
  h_recoJetPt_PFPT_15_C8 = (TH1D*) H_recoJetPt_PFPT_15_C8->ProjectionX("h_recoJetPt_PFPT_15_C8");
  h_recoJetPt_PFPT_15_C7 = (TH1D*) H_recoJetPt_PFPT_15_C7->ProjectionX("h_recoJetPt_PFPT_15_C7");
  h_recoJetPt_PFPT_15_C6 = (TH1D*) H_recoJetPt_PFPT_15_C6->ProjectionX("h_recoJetPt_PFPT_15_C6");
  h_recoJetPt_PFPT_15_C5 = (TH1D*) H_recoJetPt_PFPT_15_C5->ProjectionX("h_recoJetPt_PFPT_15_C5");
  h_recoJetPt_PFPT_15_C4 = (TH1D*) H_recoJetPt_PFPT_15_C4->ProjectionX("h_recoJetPt_PFPT_15_C4");
  h_recoJetPt_PFPT_15_C3 = (TH1D*) H_recoJetPt_PFPT_15_C3->ProjectionX("h_recoJetPt_PFPT_15_C3");
  h_recoJetPt_PFPT_15_C2 = (TH1D*) H_recoJetPt_PFPT_15_C2->ProjectionX("h_recoJetPt_PFPT_15_C2");
  h_recoJetPt_PFPT_15_C1 = (TH1D*) H_recoJetPt_PFPT_15_C1->ProjectionX("h_recoJetPt_PFPT_15_C1");

  h_recoJetPtNoSubleading_PFPT_15_C16 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C16->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C16");
  h_recoJetPtNoSubleading_PFPT_15_C15 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C15->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C15");
  h_recoJetPtNoSubleading_PFPT_15_C14 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C14->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C14");
  h_recoJetPtNoSubleading_PFPT_15_C13 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C13->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C13");
  h_recoJetPtNoSubleading_PFPT_15_C12 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C12->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C12");
  h_recoJetPtNoSubleading_PFPT_15_C11 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C11->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C11");
  h_recoJetPtNoSubleading_PFPT_15_C10 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C10->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C10");
  h_recoJetPtNoSubleading_PFPT_15_C9 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C9->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C9");
  h_recoJetPtNoSubleading_PFPT_15_C8 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C8->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C8");
  h_recoJetPtNoSubleading_PFPT_15_C7 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C7->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C7");
  h_recoJetPtNoSubleading_PFPT_15_C6 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C6->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C6");
  h_recoJetPtNoSubleading_PFPT_15_C5 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C5->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C5");
  h_recoJetPtNoSubleading_PFPT_15_C4 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C4->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C4");
  h_recoJetPtNoSubleading_PFPT_15_C3 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C3->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C3");
  h_recoJetPtNoSubleading_PFPT_15_C2 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C2->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C2");
  h_recoJetPtNoSubleading_PFPT_15_C1 = (TH1D*) H_recoJetPtNoSubleading_PFPT_15_C1->ProjectionX("h_recoJetPtNoSubleading_PFPT_15_C1");


  h_recoJetPt_PFPT_20_C16 = (TH1D*) H_recoJetPt_PFPT_20_C16->ProjectionX("h_recoJetPt_PFPT_20_C16");
  h_recoJetPt_PFPT_20_C15 = (TH1D*) H_recoJetPt_PFPT_20_C15->ProjectionX("h_recoJetPt_PFPT_20_C15");
  h_recoJetPt_PFPT_20_C14 = (TH1D*) H_recoJetPt_PFPT_20_C14->ProjectionX("h_recoJetPt_PFPT_20_C14");
  h_recoJetPt_PFPT_20_C13 = (TH1D*) H_recoJetPt_PFPT_20_C13->ProjectionX("h_recoJetPt_PFPT_20_C13");
  h_recoJetPt_PFPT_20_C12 = (TH1D*) H_recoJetPt_PFPT_20_C12->ProjectionX("h_recoJetPt_PFPT_20_C12");
  h_recoJetPt_PFPT_20_C11 = (TH1D*) H_recoJetPt_PFPT_20_C11->ProjectionX("h_recoJetPt_PFPT_20_C11");
  h_recoJetPt_PFPT_20_C10 = (TH1D*) H_recoJetPt_PFPT_20_C10->ProjectionX("h_recoJetPt_PFPT_20_C10");
  h_recoJetPt_PFPT_20_C9 = (TH1D*) H_recoJetPt_PFPT_20_C9->ProjectionX("h_recoJetPt_PFPT_20_C9");
  h_recoJetPt_PFPT_20_C8 = (TH1D*) H_recoJetPt_PFPT_20_C8->ProjectionX("h_recoJetPt_PFPT_20_C8");
  h_recoJetPt_PFPT_20_C7 = (TH1D*) H_recoJetPt_PFPT_20_C7->ProjectionX("h_recoJetPt_PFPT_20_C7");
  h_recoJetPt_PFPT_20_C6 = (TH1D*) H_recoJetPt_PFPT_20_C6->ProjectionX("h_recoJetPt_PFPT_20_C6");
  h_recoJetPt_PFPT_20_C5 = (TH1D*) H_recoJetPt_PFPT_20_C5->ProjectionX("h_recoJetPt_PFPT_20_C5");
  h_recoJetPt_PFPT_20_C4 = (TH1D*) H_recoJetPt_PFPT_20_C4->ProjectionX("h_recoJetPt_PFPT_20_C4");
  h_recoJetPt_PFPT_20_C3 = (TH1D*) H_recoJetPt_PFPT_20_C3->ProjectionX("h_recoJetPt_PFPT_20_C3");
  h_recoJetPt_PFPT_20_C2 = (TH1D*) H_recoJetPt_PFPT_20_C2->ProjectionX("h_recoJetPt_PFPT_20_C2");
  h_recoJetPt_PFPT_20_C1 = (TH1D*) H_recoJetPt_PFPT_20_C1->ProjectionX("h_recoJetPt_PFPT_20_C1");

  h_recoJetPtNoSubleading_PFPT_20_C16 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C16->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C16");
  h_recoJetPtNoSubleading_PFPT_20_C15 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C15->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C15");
  h_recoJetPtNoSubleading_PFPT_20_C14 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C14->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C14");
  h_recoJetPtNoSubleading_PFPT_20_C13 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C13->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C13");
  h_recoJetPtNoSubleading_PFPT_20_C12 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C12->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C12");
  h_recoJetPtNoSubleading_PFPT_20_C11 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C11->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C11");
  h_recoJetPtNoSubleading_PFPT_20_C10 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C10->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C10");
  h_recoJetPtNoSubleading_PFPT_20_C9 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C9->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C9");
  h_recoJetPtNoSubleading_PFPT_20_C8 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C8->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C8");
  h_recoJetPtNoSubleading_PFPT_20_C7 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C7->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C7");
  h_recoJetPtNoSubleading_PFPT_20_C6 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C6->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C6");
  h_recoJetPtNoSubleading_PFPT_20_C5 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C5->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C5");
  h_recoJetPtNoSubleading_PFPT_20_C4 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C4->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C4");
  h_recoJetPtNoSubleading_PFPT_20_C3 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C3->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C3");
  h_recoJetPtNoSubleading_PFPT_20_C2 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C2->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C2");
  h_recoJetPtNoSubleading_PFPT_20_C1 = (TH1D*) H_recoJetPtNoSubleading_PFPT_20_C1->ProjectionX("h_recoJetPtNoSubleading_PFPT_20_C1");


  TH1D *h_recoJetPt_PFPT_10_Cm4, *h_recoJetPt_PFPT_10_Cm3, *h_recoJetPt_PFPT_10_Cm2, *h_recoJetPt_PFPT_10_Cm1;
  TH1D *h_recoJetPtNoSubleading_PFPT_10_Cm4, *h_recoJetPtNoSubleading_PFPT_10_Cm3, *h_recoJetPtNoSubleading_PFPT_10_Cm2, *h_recoJetPtNoSubleading_PFPT_10_Cm1;

  TH1D *h_recoJetPt_PFPT_15_Cm4, *h_recoJetPt_PFPT_15_Cm3, *h_recoJetPt_PFPT_15_Cm2, *h_recoJetPt_PFPT_15_Cm1;
  TH1D *h_recoJetPtNoSubleading_PFPT_15_Cm4, *h_recoJetPtNoSubleading_PFPT_15_Cm3, *h_recoJetPtNoSubleading_PFPT_15_Cm2, *h_recoJetPtNoSubleading_PFPT_15_Cm1;

  TH1D *h_recoJetPt_PFPT_20_Cm4, *h_recoJetPt_PFPT_20_Cm3, *h_recoJetPt_PFPT_20_Cm2, *h_recoJetPt_PFPT_20_Cm1;
  TH1D *h_recoJetPtNoSubleading_PFPT_20_Cm4, *h_recoJetPtNoSubleading_PFPT_20_Cm3, *h_recoJetPtNoSubleading_PFPT_20_Cm2, *h_recoJetPtNoSubleading_PFPT_20_Cm1;

  h_recoJetPt_PFPT_10_Cm4 = (TH1D*) h_recoJetPt_PFPT_10_C16->Clone("h_recoJetPt_PFPT_10_Cm4"); // 50-80
  h_recoJetPt_PFPT_10_Cm4->Add(h_recoJetPt_PFPT_10_C15);
  h_recoJetPt_PFPT_10_Cm4->Add(h_recoJetPt_PFPT_10_C14);
  h_recoJetPt_PFPT_10_Cm4->Add(h_recoJetPt_PFPT_10_C13);
  h_recoJetPt_PFPT_10_Cm4->Add(h_recoJetPt_PFPT_10_C12);
  h_recoJetPt_PFPT_10_Cm4->Add(h_recoJetPt_PFPT_10_C11);

  h_recoJetPt_PFPT_10_Cm3 = (TH1D*) h_recoJetPt_PFPT_10_C10->Clone("h_recoJetPt_PFPT_10_Cm3");  // 30-50
  h_recoJetPt_PFPT_10_Cm3->Add(h_recoJetPt_PFPT_10_C9);
  h_recoJetPt_PFPT_10_Cm3->Add(h_recoJetPt_PFPT_10_C8);
  h_recoJetPt_PFPT_10_Cm3->Add(h_recoJetPt_PFPT_10_C7);

  h_recoJetPt_PFPT_10_Cm2 = (TH1D*) h_recoJetPt_PFPT_10_C6->Clone("h_recoJetPt_PFPT_10_Cm2");  // 10-30
  h_recoJetPt_PFPT_10_Cm2->Add(h_recoJetPt_PFPT_10_C5);
  h_recoJetPt_PFPT_10_Cm2->Add(h_recoJetPt_PFPT_10_C4);
  h_recoJetPt_PFPT_10_Cm2->Add(h_recoJetPt_PFPT_10_C3);

  h_recoJetPt_PFPT_10_Cm1 = (TH1D*) h_recoJetPt_PFPT_10_C2->Clone("h_recoJetPt_PFPT_10_Cm1");  // 0-10
  h_recoJetPt_PFPT_10_Cm1->Add(h_recoJetPt_PFPT_10_C1);

  h_recoJetPtNoSubleading_PFPT_10_Cm4 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_C16->Clone("h_recoJetPtNoSubleading_PFPT_10_Cm4"); // 50-80
  h_recoJetPtNoSubleading_PFPT_10_Cm4->Add(h_recoJetPtNoSubleading_PFPT_10_C15);
  h_recoJetPtNoSubleading_PFPT_10_Cm4->Add(h_recoJetPtNoSubleading_PFPT_10_C14);
  h_recoJetPtNoSubleading_PFPT_10_Cm4->Add(h_recoJetPtNoSubleading_PFPT_10_C13);
  h_recoJetPtNoSubleading_PFPT_10_Cm4->Add(h_recoJetPtNoSubleading_PFPT_10_C12);
  h_recoJetPtNoSubleading_PFPT_10_Cm4->Add(h_recoJetPtNoSubleading_PFPT_10_C11);

  h_recoJetPtNoSubleading_PFPT_10_Cm3 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_C10->Clone("h_recoJetPtNoSubleading_PFPT_10_Cm3");  // 30-50
  h_recoJetPtNoSubleading_PFPT_10_Cm3->Add(h_recoJetPtNoSubleading_PFPT_10_C9);
  h_recoJetPtNoSubleading_PFPT_10_Cm3->Add(h_recoJetPtNoSubleading_PFPT_10_C8);
  h_recoJetPtNoSubleading_PFPT_10_Cm3->Add(h_recoJetPtNoSubleading_PFPT_10_C7);

  h_recoJetPtNoSubleading_PFPT_10_Cm2 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_C6->Clone("h_recoJetPtNoSubleading_PFPT_10_Cm2");  // 10-30
  h_recoJetPtNoSubleading_PFPT_10_Cm2->Add(h_recoJetPtNoSubleading_PFPT_10_C5);
  h_recoJetPtNoSubleading_PFPT_10_Cm2->Add(h_recoJetPtNoSubleading_PFPT_10_C4);
  h_recoJetPtNoSubleading_PFPT_10_Cm2->Add(h_recoJetPtNoSubleading_PFPT_10_C3);

  h_recoJetPtNoSubleading_PFPT_10_Cm1 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_C2->Clone("h_recoJetPtNoSubleading_PFPT_10_Cm1");  // 0-10
  h_recoJetPtNoSubleading_PFPT_10_Cm1->Add(h_recoJetPtNoSubleading_PFPT_10_C1);


  h_recoJetPt_PFPT_15_Cm4 = (TH1D*) h_recoJetPt_PFPT_15_C16->Clone("h_recoJetPt_PFPT_15_Cm4"); // 50-80
  h_recoJetPt_PFPT_15_Cm4->Add(h_recoJetPt_PFPT_15_C15);
  h_recoJetPt_PFPT_15_Cm4->Add(h_recoJetPt_PFPT_15_C14);
  h_recoJetPt_PFPT_15_Cm4->Add(h_recoJetPt_PFPT_15_C13);
  h_recoJetPt_PFPT_15_Cm4->Add(h_recoJetPt_PFPT_15_C12);
  h_recoJetPt_PFPT_15_Cm4->Add(h_recoJetPt_PFPT_15_C11);

  h_recoJetPt_PFPT_15_Cm3 = (TH1D*) h_recoJetPt_PFPT_15_C10->Clone("h_recoJetPt_PFPT_15_Cm3");  // 30-50
  h_recoJetPt_PFPT_15_Cm3->Add(h_recoJetPt_PFPT_15_C9);
  h_recoJetPt_PFPT_15_Cm3->Add(h_recoJetPt_PFPT_15_C8);
  h_recoJetPt_PFPT_15_Cm3->Add(h_recoJetPt_PFPT_15_C7);

  h_recoJetPt_PFPT_15_Cm2 = (TH1D*) h_recoJetPt_PFPT_15_C6->Clone("h_recoJetPt_PFPT_15_Cm2");  // 10-30
  h_recoJetPt_PFPT_15_Cm2->Add(h_recoJetPt_PFPT_15_C5);
  h_recoJetPt_PFPT_15_Cm2->Add(h_recoJetPt_PFPT_15_C4);
  h_recoJetPt_PFPT_15_Cm2->Add(h_recoJetPt_PFPT_15_C3);

  h_recoJetPt_PFPT_15_Cm1 = (TH1D*) h_recoJetPt_PFPT_15_C2->Clone("h_recoJetPt_PFPT_15_Cm1");  // 0-10
  h_recoJetPt_PFPT_15_Cm1->Add(h_recoJetPt_PFPT_15_C1);

  h_recoJetPtNoSubleading_PFPT_15_Cm4 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_C16->Clone("h_recoJetPtNoSubleading_PFPT_15_Cm4"); // 50-80
  h_recoJetPtNoSubleading_PFPT_15_Cm4->Add(h_recoJetPtNoSubleading_PFPT_15_C15);
  h_recoJetPtNoSubleading_PFPT_15_Cm4->Add(h_recoJetPtNoSubleading_PFPT_15_C14);
  h_recoJetPtNoSubleading_PFPT_15_Cm4->Add(h_recoJetPtNoSubleading_PFPT_15_C13);
  h_recoJetPtNoSubleading_PFPT_15_Cm4->Add(h_recoJetPtNoSubleading_PFPT_15_C12);
  h_recoJetPtNoSubleading_PFPT_15_Cm4->Add(h_recoJetPtNoSubleading_PFPT_15_C11);

  h_recoJetPtNoSubleading_PFPT_15_Cm3 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_C10->Clone("h_recoJetPtNoSubleading_PFPT_15_Cm3");  // 30-50
  h_recoJetPtNoSubleading_PFPT_15_Cm3->Add(h_recoJetPtNoSubleading_PFPT_15_C9);
  h_recoJetPtNoSubleading_PFPT_15_Cm3->Add(h_recoJetPtNoSubleading_PFPT_15_C8);
  h_recoJetPtNoSubleading_PFPT_15_Cm3->Add(h_recoJetPtNoSubleading_PFPT_15_C7);

  h_recoJetPtNoSubleading_PFPT_15_Cm2 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_C6->Clone("h_recoJetPtNoSubleading_PFPT_15_Cm2");  // 10-30
  h_recoJetPtNoSubleading_PFPT_15_Cm2->Add(h_recoJetPtNoSubleading_PFPT_15_C5);
  h_recoJetPtNoSubleading_PFPT_15_Cm2->Add(h_recoJetPtNoSubleading_PFPT_15_C4);
  h_recoJetPtNoSubleading_PFPT_15_Cm2->Add(h_recoJetPtNoSubleading_PFPT_15_C3);

  h_recoJetPtNoSubleading_PFPT_15_Cm1 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_C2->Clone("h_recoJetPtNoSubleading_PFPT_15_Cm1");  // 0-10
  h_recoJetPtNoSubleading_PFPT_15_Cm1->Add(h_recoJetPtNoSubleading_PFPT_15_C1);



  h_recoJetPt_PFPT_20_Cm4 = (TH1D*) h_recoJetPt_PFPT_20_C16->Clone("h_recoJetPt_PFPT_20_Cm4"); // 50-80
  h_recoJetPt_PFPT_20_Cm4->Add(h_recoJetPt_PFPT_20_C15);
  h_recoJetPt_PFPT_20_Cm4->Add(h_recoJetPt_PFPT_20_C14);
  h_recoJetPt_PFPT_20_Cm4->Add(h_recoJetPt_PFPT_20_C13);
  h_recoJetPt_PFPT_20_Cm4->Add(h_recoJetPt_PFPT_20_C12);
  h_recoJetPt_PFPT_20_Cm4->Add(h_recoJetPt_PFPT_20_C11);

  h_recoJetPt_PFPT_20_Cm3 = (TH1D*) h_recoJetPt_PFPT_20_C10->Clone("h_recoJetPt_PFPT_20_Cm3");  // 30-50
  h_recoJetPt_PFPT_20_Cm3->Add(h_recoJetPt_PFPT_20_C9);
  h_recoJetPt_PFPT_20_Cm3->Add(h_recoJetPt_PFPT_20_C8);
  h_recoJetPt_PFPT_20_Cm3->Add(h_recoJetPt_PFPT_20_C7);

  h_recoJetPt_PFPT_20_Cm2 = (TH1D*) h_recoJetPt_PFPT_20_C6->Clone("h_recoJetPt_PFPT_20_Cm2");  // 10-30
  h_recoJetPt_PFPT_20_Cm2->Add(h_recoJetPt_PFPT_20_C5);
  h_recoJetPt_PFPT_20_Cm2->Add(h_recoJetPt_PFPT_20_C4);
  h_recoJetPt_PFPT_20_Cm2->Add(h_recoJetPt_PFPT_20_C3);

  h_recoJetPt_PFPT_20_Cm1 = (TH1D*) h_recoJetPt_PFPT_20_C2->Clone("h_recoJetPt_PFPT_20_Cm1");  // 0-10
  h_recoJetPt_PFPT_20_Cm1->Add(h_recoJetPt_PFPT_20_C1);

  h_recoJetPtNoSubleading_PFPT_20_Cm4 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_C16->Clone("h_recoJetPtNoSubleading_PFPT_20_Cm4"); // 50-80
  h_recoJetPtNoSubleading_PFPT_20_Cm4->Add(h_recoJetPtNoSubleading_PFPT_20_C15);
  h_recoJetPtNoSubleading_PFPT_20_Cm4->Add(h_recoJetPtNoSubleading_PFPT_20_C14);
  h_recoJetPtNoSubleading_PFPT_20_Cm4->Add(h_recoJetPtNoSubleading_PFPT_20_C13);
  h_recoJetPtNoSubleading_PFPT_20_Cm4->Add(h_recoJetPtNoSubleading_PFPT_20_C12);
  h_recoJetPtNoSubleading_PFPT_20_Cm4->Add(h_recoJetPtNoSubleading_PFPT_20_C11);

  h_recoJetPtNoSubleading_PFPT_20_Cm3 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_C10->Clone("h_recoJetPtNoSubleading_PFPT_20_Cm3");  // 30-50
  h_recoJetPtNoSubleading_PFPT_20_Cm3->Add(h_recoJetPtNoSubleading_PFPT_20_C9);
  h_recoJetPtNoSubleading_PFPT_20_Cm3->Add(h_recoJetPtNoSubleading_PFPT_20_C8);
  h_recoJetPtNoSubleading_PFPT_20_Cm3->Add(h_recoJetPtNoSubleading_PFPT_20_C7);

  h_recoJetPtNoSubleading_PFPT_20_Cm2 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_C6->Clone("h_recoJetPtNoSubleading_PFPT_20_Cm2");  // 10-30
  h_recoJetPtNoSubleading_PFPT_20_Cm2->Add(h_recoJetPtNoSubleading_PFPT_20_C5);
  h_recoJetPtNoSubleading_PFPT_20_Cm2->Add(h_recoJetPtNoSubleading_PFPT_20_C4);
  h_recoJetPtNoSubleading_PFPT_20_Cm2->Add(h_recoJetPtNoSubleading_PFPT_20_C3);

  h_recoJetPtNoSubleading_PFPT_20_Cm1 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_C2->Clone("h_recoJetPtNoSubleading_PFPT_20_Cm1");  // 0-10
  h_recoJetPtNoSubleading_PFPT_20_Cm1->Add(h_recoJetPtNoSubleading_PFPT_20_C1);



  // rebin
  const int N_edge = 5;
  double newAxis[N_edge] = {20,30,50,100,500};
  
  h_recoJetPt_PFPT_10_Cm4 = (TH1D*) h_recoJetPt_PFPT_10_Cm4->Rebin(N_edge-1,"h_recoJetPt_PFPT_10_Cm4",newAxis);
  h_recoJetPt_PFPT_10_Cm3 = (TH1D*) h_recoJetPt_PFPT_10_Cm3->Rebin(N_edge-1,"h_recoJetPt_PFPT_10_Cm3",newAxis);
  h_recoJetPt_PFPT_10_Cm2 = (TH1D*) h_recoJetPt_PFPT_10_Cm2->Rebin(N_edge-1,"h_recoJetPt_PFPT_10_Cm2",newAxis);
  h_recoJetPt_PFPT_10_Cm1 = (TH1D*) h_recoJetPt_PFPT_10_Cm1->Rebin(N_edge-1,"h_recoJetPt_PFPT_10_Cm1",newAxis);
  divideByBinwidth(h_recoJetPt_PFPT_10_Cm4);
  divideByBinwidth(h_recoJetPt_PFPT_10_Cm3);
  divideByBinwidth(h_recoJetPt_PFPT_10_Cm2);
  divideByBinwidth(h_recoJetPt_PFPT_10_Cm1);

  h_recoJetPtNoSubleading_PFPT_10_Cm4 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_Cm4->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_10_Cm4",newAxis);
  h_recoJetPtNoSubleading_PFPT_10_Cm3 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_Cm3->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_10_Cm3",newAxis);
  h_recoJetPtNoSubleading_PFPT_10_Cm2 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_Cm2->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_10_Cm2",newAxis);
  h_recoJetPtNoSubleading_PFPT_10_Cm1 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_Cm1->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_10_Cm1",newAxis);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_10_Cm4);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_10_Cm3);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_10_Cm2);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_10_Cm1);


  h_recoJetPt_PFPT_15_Cm4 = (TH1D*) h_recoJetPt_PFPT_15_Cm4->Rebin(N_edge-1,"h_recoJetPt_PFPT_15_Cm4",newAxis);
  h_recoJetPt_PFPT_15_Cm3 = (TH1D*) h_recoJetPt_PFPT_15_Cm3->Rebin(N_edge-1,"h_recoJetPt_PFPT_15_Cm3",newAxis);
  h_recoJetPt_PFPT_15_Cm2 = (TH1D*) h_recoJetPt_PFPT_15_Cm2->Rebin(N_edge-1,"h_recoJetPt_PFPT_15_Cm2",newAxis);
  h_recoJetPt_PFPT_15_Cm1 = (TH1D*) h_recoJetPt_PFPT_15_Cm1->Rebin(N_edge-1,"h_recoJetPt_PFPT_15_Cm1",newAxis);
  divideByBinwidth(h_recoJetPt_PFPT_15_Cm4);
  divideByBinwidth(h_recoJetPt_PFPT_15_Cm3);
  divideByBinwidth(h_recoJetPt_PFPT_15_Cm2);
  divideByBinwidth(h_recoJetPt_PFPT_15_Cm1);

  h_recoJetPtNoSubleading_PFPT_15_Cm4 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_Cm4->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_15_Cm4",newAxis);
  h_recoJetPtNoSubleading_PFPT_15_Cm3 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_Cm3->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_15_Cm3",newAxis);
  h_recoJetPtNoSubleading_PFPT_15_Cm2 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_Cm2->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_15_Cm2",newAxis);
  h_recoJetPtNoSubleading_PFPT_15_Cm1 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_Cm1->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_15_Cm1",newAxis);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_15_Cm4);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_15_Cm3);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_15_Cm2);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_15_Cm1);


  h_recoJetPt_PFPT_20_Cm4 = (TH1D*) h_recoJetPt_PFPT_20_Cm4->Rebin(N_edge-1,"h_recoJetPt_PFPT_20_Cm4",newAxis);
  h_recoJetPt_PFPT_20_Cm3 = (TH1D*) h_recoJetPt_PFPT_20_Cm3->Rebin(N_edge-1,"h_recoJetPt_PFPT_20_Cm3",newAxis);
  h_recoJetPt_PFPT_20_Cm2 = (TH1D*) h_recoJetPt_PFPT_20_Cm2->Rebin(N_edge-1,"h_recoJetPt_PFPT_20_Cm2",newAxis);
  h_recoJetPt_PFPT_20_Cm1 = (TH1D*) h_recoJetPt_PFPT_20_Cm1->Rebin(N_edge-1,"h_recoJetPt_PFPT_20_Cm1",newAxis);
  divideByBinwidth(h_recoJetPt_PFPT_20_Cm4);
  divideByBinwidth(h_recoJetPt_PFPT_20_Cm3);
  divideByBinwidth(h_recoJetPt_PFPT_20_Cm2);
  divideByBinwidth(h_recoJetPt_PFPT_20_Cm1);

  h_recoJetPtNoSubleading_PFPT_20_Cm4 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_Cm4->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_20_Cm4",newAxis);
  h_recoJetPtNoSubleading_PFPT_20_Cm3 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_Cm3->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_20_Cm3",newAxis);
  h_recoJetPtNoSubleading_PFPT_20_Cm2 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_Cm2->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_20_Cm2",newAxis);
  h_recoJetPtNoSubleading_PFPT_20_Cm1 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_Cm1->Rebin(N_edge-1,"h_recoJetPtNoSubleading_PFPT_20_Cm1",newAxis);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_20_Cm4);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_20_Cm3);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_20_Cm2);
  divideByBinwidth(h_recoJetPtNoSubleading_PFPT_20_Cm1);
  

  

  TH1D *r_PFPT_10_Cm4, *r_PFPT_10_Cm3, *r_PFPT_10_Cm2, *r_PFPT_10_Cm1;
  TH1D *r_PFPT_15_Cm4, *r_PFPT_15_Cm3, *r_PFPT_15_Cm2, *r_PFPT_15_Cm1;
  TH1D *r_PFPT_20_Cm4, *r_PFPT_20_Cm3, *r_PFPT_20_Cm2, *r_PFPT_20_Cm1;

  
  r_PFPT_10_Cm4 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_Cm4->Clone("r_PFPT_10_Cm4");
  r_PFPT_10_Cm4->Divide(h_recoJetPtNoSubleading_PFPT_10_Cm4,h_recoJetPt_PFPT_10_Cm4, 1,1,"B");

  r_PFPT_10_Cm3 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_Cm3->Clone("r_PFPT_10_Cm3");
  r_PFPT_10_Cm3->Divide(h_recoJetPtNoSubleading_PFPT_10_Cm3,h_recoJetPt_PFPT_10_Cm3, 1,1,"B");

  r_PFPT_10_Cm2 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_Cm2->Clone("r_PFPT_10_Cm2");
  r_PFPT_10_Cm2->Divide(h_recoJetPtNoSubleading_PFPT_10_Cm2,h_recoJetPt_PFPT_10_Cm2, 1,1,"B");

  r_PFPT_10_Cm1 = (TH1D*) h_recoJetPtNoSubleading_PFPT_10_Cm1->Clone("r_PFPT_10_Cm1");
  r_PFPT_10_Cm1->Divide(h_recoJetPtNoSubleading_PFPT_10_Cm1,h_recoJetPt_PFPT_10_Cm1, 1,1,"B");


  r_PFPT_15_Cm4 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_Cm4->Clone("r_PFPT_15_Cm4");
  r_PFPT_15_Cm4->Divide(h_recoJetPtNoSubleading_PFPT_15_Cm4,h_recoJetPt_PFPT_15_Cm4, 1,1,"B");

  r_PFPT_15_Cm3 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_Cm3->Clone("r_PFPT_15_Cm3");
  r_PFPT_15_Cm3->Divide(h_recoJetPtNoSubleading_PFPT_15_Cm3,h_recoJetPt_PFPT_15_Cm3, 1,1,"B");

  r_PFPT_15_Cm2 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_Cm2->Clone("r_PFPT_15_Cm2");
  r_PFPT_15_Cm2->Divide(h_recoJetPtNoSubleading_PFPT_15_Cm2,h_recoJetPt_PFPT_15_Cm2, 1,1,"B");

  r_PFPT_15_Cm1 = (TH1D*) h_recoJetPtNoSubleading_PFPT_15_Cm1->Clone("r_PFPT_15_Cm1");
  r_PFPT_15_Cm1->Divide(h_recoJetPtNoSubleading_PFPT_15_Cm1,h_recoJetPt_PFPT_15_Cm1, 1,1,"B");


  r_PFPT_20_Cm4 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_Cm4->Clone("r_PFPT_20_Cm4");
  r_PFPT_20_Cm4->Divide(h_recoJetPtNoSubleading_PFPT_20_Cm4,h_recoJetPt_PFPT_20_Cm4, 1,1,"B");

  r_PFPT_20_Cm3 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_Cm3->Clone("r_PFPT_20_Cm3");
  r_PFPT_20_Cm3->Divide(h_recoJetPtNoSubleading_PFPT_20_Cm3,h_recoJetPt_PFPT_20_Cm3, 1,1,"B");

  r_PFPT_20_Cm2 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_Cm2->Clone("r_PFPT_20_Cm2");
  r_PFPT_20_Cm2->Divide(h_recoJetPtNoSubleading_PFPT_20_Cm2,h_recoJetPt_PFPT_20_Cm2, 1,1,"B");

  r_PFPT_20_Cm1 = (TH1D*) h_recoJetPtNoSubleading_PFPT_20_Cm1->Clone("r_PFPT_20_Cm1");
  r_PFPT_20_Cm1->Divide(h_recoJetPtNoSubleading_PFPT_20_Cm1,h_recoJetPt_PFPT_20_Cm1, 1,1,"B");


  stylizeHistograms(r_PFPT_10_Cm4, r_PFPT_15_Cm4, r_PFPT_20_Cm4);
  stylizeHistograms(r_PFPT_10_Cm3, r_PFPT_15_Cm3, r_PFPT_20_Cm3);
  stylizeHistograms(r_PFPT_10_Cm2, r_PFPT_15_Cm2, r_PFPT_20_Cm2);
  stylizeHistograms(r_PFPT_10_Cm1, r_PFPT_15_Cm1, r_PFPT_20_Cm1);

  TLegend *leg = new TLegend(0.65,0.65,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(r_PFPT_10_Cm4,"#it{p}_{T}^{PF} > 10 GeV");
  leg->AddEntry(r_PFPT_15_Cm4,"#it{p}_{T}^{PF} > 15 GeV");
  leg->AddEntry(r_PFPT_20_Cm4,"#it{p}_{T}^{PF} > 20 GeV");
  
  TCanvas *canv_Cm4 = new TCanvas("canv_Cm4","canv_Cm4",700,700);
  canv_Cm4->cd();
  TPad *pad_Cm4 = new TPad("pad_Cm4","pad_Cm4",0,0,1,1);
  pad_Cm4->SetLeftMargin(0.15);
  pad_Cm4->SetBottomMargin(0.15);
  pad_Cm4->Draw();
  pad_Cm4->cd();
  r_PFPT_10_Cm4->SetTitleSize(0.06);
  r_PFPT_10_Cm4->SetTitle("HYDJET 50-80%");
  r_PFPT_10_Cm4->SetStats(0);
  r_PFPT_10_Cm4->GetYaxis()->SetTitleSize(0.055);
  r_PFPT_10_Cm4->GetYaxis()->SetLabelSize(0.045);
  r_PFPT_10_Cm4->GetYaxis()->SetTitle("No-subleading / All");
  r_PFPT_10_Cm4->GetYaxis()->SetRangeUser(0,1.05);
  r_PFPT_10_Cm4->GetXaxis()->SetTitleSize(0.055);
  r_PFPT_10_Cm4->GetXaxis()->SetLabelSize(0.045);
  r_PFPT_10_Cm4->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r_PFPT_10_Cm4->Draw();
  r_PFPT_15_Cm4->Draw("same");
  r_PFPT_20_Cm4->Draw("same");
  leg->Draw();
  canv_Cm4->SaveAs("../../../figures/subleadingPFCand/canv_Cm4.pdf");

  TCanvas *canv_Cm3 = new TCanvas("canv_Cm3","canv_Cm3",700,700);
  canv_Cm3->cd();
  TPad *pad_Cm3 = new TPad("pad_Cm3","pad_Cm3",0,0,1,1);
  pad_Cm3->SetLeftMargin(0.15);
  pad_Cm3->SetBottomMargin(0.15);
  pad_Cm3->Draw();
  pad_Cm3->cd();
  r_PFPT_10_Cm3->SetTitleSize(0.06);
  r_PFPT_10_Cm3->SetTitle("HYDJET 30-50%");
  r_PFPT_10_Cm3->SetStats(0);
  r_PFPT_10_Cm3->GetYaxis()->SetTitleSize(0.055);
  r_PFPT_10_Cm3->GetYaxis()->SetLabelSize(0.045);
  r_PFPT_10_Cm3->GetYaxis()->SetTitle("No-subleading / All");
  r_PFPT_10_Cm3->GetYaxis()->SetRangeUser(0,1.05);
  r_PFPT_10_Cm3->GetXaxis()->SetTitleSize(0.055);
  r_PFPT_10_Cm3->GetXaxis()->SetLabelSize(0.045);
  r_PFPT_10_Cm3->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r_PFPT_10_Cm3->Draw();
  r_PFPT_15_Cm3->Draw("same");
  r_PFPT_20_Cm3->Draw("same");
  leg->Draw();
  canv_Cm3->SaveAs("../../../figures/subleadingPFCand/canv_Cm3.pdf");

  TCanvas *canv_Cm2 = new TCanvas("canv_Cm2","canv_Cm2",700,700);
  canv_Cm2->cd();
  TPad *pad_Cm2 = new TPad("pad_Cm2","pad_Cm2",0,0,1,1);
  pad_Cm2->SetLeftMargin(0.15);
  pad_Cm2->SetBottomMargin(0.15);
  pad_Cm2->Draw();
  pad_Cm2->cd();
  r_PFPT_10_Cm2->SetTitleSize(0.06);
  r_PFPT_10_Cm2->SetTitle("HYDJET 10-30%");
  r_PFPT_10_Cm2->SetStats(0);
  r_PFPT_10_Cm2->GetYaxis()->SetTitleSize(0.055);
  r_PFPT_10_Cm2->GetYaxis()->SetLabelSize(0.045);
  r_PFPT_10_Cm2->GetYaxis()->SetTitle("No-subleading / All");
  r_PFPT_10_Cm2->GetYaxis()->SetRangeUser(0,1.05);
  r_PFPT_10_Cm2->GetXaxis()->SetTitleSize(0.055);
  r_PFPT_10_Cm2->GetXaxis()->SetLabelSize(0.045);
  r_PFPT_10_Cm2->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r_PFPT_10_Cm2->Draw();
  r_PFPT_15_Cm2->Draw("same");
  r_PFPT_20_Cm2->Draw("same");
  leg->Draw();
  canv_Cm2->SaveAs("../../../figures/subleadingPFCand/canv_Cm2.pdf");

  TCanvas *canv_Cm1 = new TCanvas("canv_Cm1","canv_Cm1",700,700);
  canv_Cm1->cd();
  TPad *pad_Cm1 = new TPad("pad_Cm1","pad_Cm1",0,0,1,1);
  pad_Cm1->SetLeftMargin(0.15);
  pad_Cm1->SetBottomMargin(0.15);
  pad_Cm1->Draw();
  pad_Cm1->cd();
  r_PFPT_10_Cm1->SetTitleSize(0.06);
  r_PFPT_10_Cm1->SetTitle("HYDJET 0-10%");
  r_PFPT_10_Cm1->SetStats(0);
  r_PFPT_10_Cm1->GetYaxis()->SetTitleSize(0.055);
  r_PFPT_10_Cm1->GetYaxis()->SetLabelSize(0.045);
  r_PFPT_10_Cm1->GetYaxis()->SetTitle("No-subleading / All");
  r_PFPT_10_Cm1->GetYaxis()->SetRangeUser(0,1.05);
  r_PFPT_10_Cm1->GetXaxis()->SetTitleSize(0.055);
  r_PFPT_10_Cm1->GetXaxis()->SetLabelSize(0.045);
  r_PFPT_10_Cm1->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  r_PFPT_10_Cm1->Draw();
  r_PFPT_15_Cm1->Draw("same");
  r_PFPT_20_Cm1->Draw("same");
  leg->Draw();
  canv_Cm1->SaveAs("../../../figures/subleadingPFCand/canv_Cm1.pdf");  
  


}
