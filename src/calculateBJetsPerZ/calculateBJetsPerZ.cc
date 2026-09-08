#include "../../headers/functions/divideByBinwidth.h"
// Roofit
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"

TFile *file_PbPb_MinBias, *file_pp_MinBias;
TFile *file_PbPb_SingleMuon, *file_pp_SingleMuon;
TFile *file_PbPb_HardProbes_jet100, *file_PbPb_HardProbes_jet80, *file_PbPb_HardProbes_jet60;
TFile *file_pp_HighEGJet_jet100, *file_pp_HighEGJet_jet80, *file_pp_HighEGJet_jet60;

TFile *file_PYTHIA, *file_PH_DiJet, *file_PH_MuJet, *file_PH_BJet;

TH2D *H_PYTHIA, *H_PH_DiJet, *H_PH_MuJet, *H_PH_BJet;
TH2D *H_pp, *H_PbPb;


TH1D *h_data;
TH1D *h_b, *h_bGS, *h_c, *h_d, *h_u, *h_s, *h_g, *h_x, *h_l;
double cFixVal = 0.0;
TH1D *h_b_draw, *h_c_draw, *h_l_draw, *h_x_draw, *h_draw, *h_roo;

const int M = 51;
double muRelPtAxis[M] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0}; // M = 51



void openTemplateFitterFiles(){

  file_pp_HighEGJet_jet60 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/pp_HighEGJet_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_2026-3-10.root");
  file_pp_HighEGJet_jet80 = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_HighEGJet_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-27.root");
  file_pp_HighEGJet_jet100 = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_HighEGJet_Jet100HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-27.root");
  file_pp_MinBias = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/pp_MinBias_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-7-6.root");
  file_pp_SingleMuon = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/pp_SingleMuon_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_2026-2-12.root");


  file_PbPb_HardProbes_jet60 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet60HLT_mu12_pTmu-15to999_tight_WDecayFilter_2026-3-11.root");
  file_PbPb_HardProbes_jet80 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet80HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-5.root");
  file_PbPb_HardProbes_jet100 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet100HLT_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-5.root");
  file_PbPb_SingleMuon = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15to999_tight_mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-5-4.root");
  file_PbPb_MinBias = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_MinBias_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root");
  file_PbPb_SingleMuon = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_SingleMuon_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_WDecayFilter_2026-2-12.root");

  file_PYTHIA = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIA/latest/PYTHIA_DiJet_pThat-30_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_2026-2-11.root");
  file_PH_DiJet = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/PYTHIAHYDJET_DiJet_pThat-30_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_vzReweight_hiBinReweight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_weightCut_2026-2-12.root");
  file_PH_MuJet = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/PYTHIAHYDJET_MuJet_pThat-30_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_vzReweight_hiBinReweight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_weightCut_2026-2-12.root");
  file_PH_BJet = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/PYTHIAHYDJET_BJet_pThat-30_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_vzReweight_hiBinReweight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_weightCut_2026-2-12.root");

  return;

}

double func_temp_1(double *x, double *par){
  double xx = x[0];
  int bin = h_b->FindBin(xx);
  double b = (par[0])*h_b->GetBinContent(bin);
  double e = (1-par[0])*h_l->GetBinContent(bin);
  return b+e;
}

double func_temp_2(double *x, double *par){
  double xx = x[0];
  int bin = h_b->FindBin(xx);
  double b = (par[0])*h_b->GetBinContent(bin);
  double c = (par[1])*h_c->GetBinContent(bin);
  double e = (1-par[0]-par[1])*h_l->GetBinContent(bin);
  return b+e+c;
}


void rebinHistograms(){
  
  h_roo = (TH1D*) h_roo->Rebin(M-1,"h_roo",muRelPtAxis);
  h_data = (TH1D*) h_data->Rebin(M-1,"h_data",muRelPtAxis);
  h_l = (TH1D*) h_l->Rebin(M-1,"h_l",muRelPtAxis);
  h_b = (TH1D*) h_b->Rebin(M-1,"h_b",muRelPtAxis);
  h_c = (TH1D*) h_c->Rebin(M-1,"h_c",muRelPtAxis);
  h_x = (TH1D*) h_x->Rebin(M-1,"h_x",muRelPtAxis);

  divideByBinwidth(h_roo);
  divideByBinwidth(h_data);
  divideByBinwidth(h_l);
  divideByBinwidth(h_b);
  divideByBinwidth(h_c);
  divideByBinwidth(h_x);

  return;
}


void stylizeHistograms(bool do2templateFit, bool do3templateFit){

  h_draw->SetTitle("");
  h_draw->SetStats(0);
  h_draw->SetMarkerStyle(8);
  h_draw->SetMarkerColor(kBlack);
  h_draw->SetLineColor(kBlack);

  h_b_draw->SetStats(0);
  h_c_draw->SetStats(0);
  h_l_draw->SetStats(0);
  h_x_draw->SetStats(0);

  double alphaVal = 0.6;
  double markerVal = 0.0;
  double markerVal2 = 1.0;

  h_b_draw->SetMarkerStyle(20);
  h_b_draw->SetMarkerSize(markerVal2);
  h_b_draw->SetMarkerColor(kRed);
  h_b_draw->SetLineColor(kRed);
  h_b_draw->SetFillColorAlpha(kRed,alphaVal);

  h_c_draw->SetMarkerStyle(20);
  h_c_draw->SetMarkerSize(markerVal2);
  h_c_draw->SetMarkerColor(kGreen+1);
  h_c_draw->SetLineColor(kGreen+1);
  h_c_draw->SetFillColorAlpha(kGreen+1,alphaVal);

  h_l_draw->SetMarkerStyle(20);
  h_l_draw->SetMarkerSize(markerVal2);
  if(do2templateFit){
    h_l_draw->SetMarkerColor(kBlue+2);
    h_l_draw->SetLineColor(kBlue+2);
    h_l_draw->SetFillColorAlpha(kBlue+2,alphaVal);
  }
  else{
    h_l_draw->SetMarkerColor(kBlue-4);
    h_l_draw->SetLineColor(kBlue-4);
    h_l_draw->SetFillColorAlpha(kBlue-4,alphaVal);
  }

  h_x_draw->SetMarkerStyle(20);
  h_x_draw->SetMarkerSize(markerVal2);
  h_x_draw->SetMarkerColor(kPink+1);
  h_x_draw->SetLineColor(kPink+1);
  h_x_draw->SetFillColorAlpha(kPink+1,alphaVal);

  return;

}

void stylizeBJetHistograms(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4){

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



double templateFitter(bool isData = 1,
		    bool ispp   = 0,
		    bool isC1   = 0,
		    bool isC2   = 1,
		    bool isC3   = 0,
		    bool isC4   = 0,		      
		    double low_jetPt = 100.,
		    double high_jetPt = 120.,
		    bool mergeCtemplates = true,
		    bool mergeBtemplates = true,
		    bool do2templateFit  = true,
		    bool do3templateFit  = false,
		    double low_x  = 0.0,
		    double high_x = 5.0,
		    double c_multiplier = 1.0,
		    double bGS_multiplier = 1.0,
		    int returnValueIndex = 1){


  TFile *f_data, *f_mc, *f_mc_bJet, *f_mc_muJet, *f_fake;
  TH2D *H_data, *H_mc, *H_mc_bJet, *H_mc_muJet, *H_fake;

  TH2D *H_mc_b, *H_mc_bGS, *H_mc_c, *H_mc_d, *H_mc_u, *H_mc_s, *H_mc_g, *H_mc_x;
  TH2D *X_mc_b, *X_mc_bGS, *X_mc_c, *X_mc_d, *X_mc_u, *X_mc_s, *X_mc_g, *X_mc_x;
  TH2D *H_mc_bJet_b, *H_mc_bJet_bGS, *H_mc_muJet_b, *H_mc_muJet_bGS, *H_mc_muJet_c;

  
  TH1D *x_b, *x_bGS, *x_c, *x_d, *x_u, *x_s, *x_g, *x_x, *x_l;
  TH1D *h_bJet_b, *h_bJet_bGS, *h_muJet_b, *h_muJet_bGS, *h_muJet_c;
  TH1D *h_fake;
  TH1D *h_vz_data, *h_vz_fake;

  int centBin = 0;

  
  
  if(ispp){

    f_mc = file_PYTHIA;
    
    if(!isData){
      f_data = f_mc;
      f_data->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_allJets_T0",H_data);
    }
    else{
      f_data = file_pp_SingleMuon;
      f_data->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn",H_data);
    }

    f_mc->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bJets_T0",H_mc_b);
    f_mc->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bGSJets_T0",H_mc_bGS);
    f_mc->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_cJets_T0",H_mc_c);
    f_mc->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_dJets_T0",H_mc_d);
    f_mc->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_uJets_T0",H_mc_u);
    f_mc->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_sJets_T0",H_mc_s);
    f_mc->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_gJets_T0",H_mc_g);
    f_mc->GetObject("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_xJets_T0",H_mc_x);

  }
  else{
    if(isC4) centBin = 4;
    else if(isC3) centBin = 3;
    else if(isC2) centBin = 2;
    else if(isC1) centBin = 1;

    f_mc = file_PH_DiJet;
    f_mc_bJet = file_PH_BJet;
    f_mc_muJet = file_PH_MuJet;

    f_fake = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_SingleMuon_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-27_coarseBins_partial.root");

    if(!isData){
      f_data = f_mc;
      f_data->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_allJets_C%iT0",centBin),H_data);
    }
    else{
      f_data = file_PbPb_SingleMuon;
      f_data->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_C%i",centBin),H_data);
    }

    f_data->GetObject(Form("h_vz_C%i",centBin),h_vz_data);
    f_fake->GetObject(Form("h_vz_C%i",centBin),h_vz_fake);
    
    f_mc->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bJets_C%iT0",centBin),H_mc_b);
    f_mc->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bGSJets_C%iT0",centBin),H_mc_bGS);
    f_mc->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_cJets_C%iT0",centBin),H_mc_c);
    f_mc->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_uJets_C%iT0",centBin),H_mc_u);
    f_mc->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_dJets_C%iT0",centBin),H_mc_d);
    f_mc->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_sJets_C%iT0",centBin),H_mc_s);
    f_mc->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_gJets_C%iT0",centBin),H_mc_g);
    f_mc->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_xJets_C%iT0",centBin),H_mc_x);

    f_mc_bJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bJets_C%iT0",centBin),H_mc_bJet_b);
    f_mc_bJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bGSJets_C%iT0",centBin),H_mc_bJet_bGS);
    f_mc_muJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bJets_C%iT0",centBin),H_mc_muJet_b);
    f_mc_muJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_bGSJets_C%iT0",centBin),H_mc_muJet_bGS);
    f_mc_muJet->GetObject(Form("h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_cJets_C%iT0",centBin),H_mc_muJet_c);

    f_fake->GetObject(Form("h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC_C%i",centBin),H_fake);

  }

  X_mc_u = (TH2D*) H_mc_u->Clone("X_mc_u");
  X_mc_d = (TH2D*) H_mc_d->Clone("X_mc_d");
  X_mc_s = (TH2D*) H_mc_s->Clone("X_mc_s");
  X_mc_g = (TH2D*) H_mc_g->Clone("X_mc_g");
  X_mc_x = (TH2D*) H_mc_x->Clone("X_mc_x");
  
  
  TH1D *binFinder = H_data->ProjectionY();
  TAxis *yaxis = binFinder->GetXaxis();
  double smallShift = 0.01;
  

  h_data = (TH1D*) H_data->ProjectionX("h_data",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
  

  if(!ispp){
    h_fake = (TH1D*) H_fake->ProjectionX("h_fake",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
    h_data->Scale(1./h_vz_data->Integral());
    h_fake->Scale(1./h_vz_fake->Integral());
    h_data->Add(h_fake,-1);
    h_data->Scale(h_vz_data->Integral());
  }

  
  h_b = (TH1D*) H_mc_b->ProjectionX("h_b",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
  h_bGS = (TH1D*) H_mc_bGS->ProjectionX("h_bGS",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
  h_c = (TH1D*) H_mc_c->ProjectionX("h_c",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
  h_d = (TH1D*) H_mc_d->ProjectionX("h_d",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
  h_u = (TH1D*) H_mc_u->ProjectionX("h_u",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
  h_s = (TH1D*) H_mc_s->ProjectionX("h_s",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
  h_g = (TH1D*) H_mc_g->ProjectionX("h_g",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
  h_x = (TH1D*) H_mc_x->ProjectionX("h_x",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));

  h_l = (TH1D*) h_u->Clone("h_l");
  h_l->Add(h_d);
  h_l->Add(h_s);
  h_l->Add(h_g);

  x_l = (TH1D*) h_l->Clone("x_l");
  x_d = (TH1D*) h_d->Clone("x_d");
  x_u = (TH1D*) h_u->Clone("x_u");
  x_s = (TH1D*) h_s->Clone("x_s");
  x_g = (TH1D*) h_g->Clone("x_g");
  x_x = (TH1D*) h_x->Clone("x_x");

  // scaling before merging the c-template.
  double c_truth = h_c->Integral() / (x_l->Integral() + h_b->Integral() + h_c->Integral() + h_bGS->Integral());
  double l_truth = x_l->Integral() / (x_l->Integral() + h_b->Integral() + h_c->Integral() + h_bGS->Integral());
  double b_truth = (h_b->Integral()+h_bGS->Integral()) / (x_l->Integral() + h_b->Integral() + h_c->Integral() + h_bGS->Integral());
  cout << "\n   c_truth = " << c_truth << "\n   l_truth = " << l_truth << "\n   b_truth = " << b_truth << "\n";

  cFixVal = c_multiplier*c_truth;

  if(!ispp){

    h_bJet_b = (TH1D*) H_mc_bJet_b->ProjectionX("h_bJet_b",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
    h_bJet_bGS = (TH1D*) H_mc_bJet_bGS->ProjectionX("h_bJet_bGS",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
    h_muJet_b = (TH1D*) H_mc_muJet_b->ProjectionX("h_muJet_b",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
    h_muJet_bGS = (TH1D*) H_mc_muJet_bGS->ProjectionX("h_muJet_bGS",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));
    h_muJet_c = (TH1D*) H_mc_muJet_c->ProjectionX("h_muJet_c",binFinder->FindBin(low_jetPt + smallShift),binFinder->FindBin(high_jetPt - smallShift));

  }

  ////////// Calculate exact count ///////////////
		
  double b_frac = 0.0;
  double eb_frac = 0.0;
  double eps = 0.001;
  if(!isData){
		
    TH1D *h_sum = (TH1D*) h_data->Clone("h_sum");

    double en_tot;
    double n_tot = h_sum->IntegralAndError(h_data->FindBin(low_x+eps),h_data->FindBin(high_x-eps),en_tot,"");
    cout << "n_tot = " << n_tot << " pm " << en_tot << endl;
    double en_b, en_bGS;
    double n_b = h_b->IntegralAndError(h_b->FindBin(low_x+eps),h_b->FindBin(high_x-eps),en_b,"");
    double n_bGS = h_bGS->IntegralAndError(h_bGS->FindBin(low_x+eps),h_bGS->FindBin(high_x-eps),en_bGS,"");
    en_b = (n_b + n_bGS) * TMath::Sqrt((en_b/n_b)*(en_b/n_b) + (en_bGS/n_bGS)*(en_bGS/n_bGS)  );
    n_b += n_bGS;
    

    
    cout << "n_b = " << n_b << " pm " << en_b << endl;
    b_frac = n_b / n_tot;
    eb_frac = b_frac*(TMath::Sqrt(TMath::Power(en_b/n_b,2)+TMath::Power(en_tot/n_tot,2)));


    cout << "exact count b-fraction = " << b_frac << "+/- " << eb_frac << endl;
    
  }

  if(!ispp){
    
    if(mergeCtemplates) h_c->Add(h_muJet_c);
    if(mergeBtemplates){
      h_b->Add(h_muJet_b);
      h_b->Add(h_bJet_b);
      h_bGS->Add(h_muJet_bGS);
      h_bGS->Add(h_bJet_bGS);
    }

  }

  double N_h_b = h_b->Integral();
  double N_h_bGS = h_bGS->Integral();
  double N_h_bTot = N_h_b + N_h_bGS;

  double f_b_truth = N_h_b / N_h_bTot;
  double f_bGS_truth = N_h_bGS / N_h_bTot; // by construction, f_b_truth + f_bGS_truth = 1
  double GS_enhancement_shift = bGS_multiplier - 1.; // GS_enhancement_shift = percent to increase the GS jets by
  double GS_enhancement_factor = 1. + (GS_enhancement_shift / f_bGS_truth); // multiply this by the current GS fraction and it will result in an increase given by GS_enhancement_shift;

   h_bGS->Scale(N_h_b*(f_bGS_truth + GS_enhancement_shift) / (N_h_bGS*(1-(f_bGS_truth + GS_enhancement_shift))));
   h_b->Add(h_bGS);

   if(do2templateFit){
    h_c->Scale(c_multiplier * c_truth  / h_c->Integral());
    h_l->Scale(l_truth / h_l->Integral());

    h_l->Add(h_c);
   }

   h_l->Scale(1./h_l->Integral());
   h_b->Scale(1./h_b->Integral());
   h_c->Scale(1./h_c->Integral());

   h_roo = (TH1D*) h_data->Clone("h_roo");
   h_data->Scale(1./h_data->Integral());
   

   rebinHistograms();

   TH1D *fitRatio = (TH1D*) h_data->Clone("fitRatio");

   h_b_draw = (TH1D*) h_b->Clone("h_b_draw");
   h_c_draw = (TH1D*) h_c->Clone("h_c_draw");
   h_l_draw = (TH1D*) h_l->Clone("h_l_draw");
   h_x_draw = (TH1D*) h_x->Clone("h_x_draw");
   h_draw = (TH1D*) h_data->Clone("h_draw");
   

   TF1 *func;

   if(do2templateFit){
     func = new TF1("func",func_temp_1,low_x,high_x,1);
     func->SetParName(0,"b");
   }
   else if(do3templateFit){
     func = new TF1("func",func_temp_2,low_x,high_x,2);
     //func->FixParameter(1,cFixVal);
     func->SetParName(0,"b");
     func->SetParName(1,"c");
   }
   else{
     std::cout << "\n\n ERROR! PLEASE SELECT 2- OR 3-TEMPLATE FIT! \n\n";
     return -1.;
   }

   func->SetParameter(0,0.5);
   func->SetParameter(1,0.4);

   // configure RooFit
   RooRealVar x("x","Observable",muRelPtAxis[0],muRelPtAxis[M-1]);
   x.setRange("fit_region",low_x,high_x);
   RooDataHist dataHist("dataHist","data histogram",x,RooFit::Import(*h_roo));
   RooDataHist bHist("bHist","MC b-jet histogram",x,RooFit::Import(*h_b));
   RooDataHist cHist("cHist","MC c-jet histogram",x,RooFit::Import(*h_c));
   RooDataHist lHist("lHist","MC light+c-jet histogram",x,RooFit::Import(*h_l));
   RooHistPdf bPdf("bPdf","MC b-jet PDF",x,bHist,0);
   RooHistPdf cPdf("cPdf","MC c-jet PDF",x,cHist,0);
   RooHistPdf lPdf("lPdf","MC light+c-jet PDF",x,lHist,0);
   RooRealVar fb("fb","b fraction",0.5,0.0,1.0);
   RooRealVar fl("fl","light fraction",0.5,0.0,1.0);
   RooFormulaVar fc("fc", "1 - fb - fl", "1 - fb - fl", RooArgList(fb, fl));
   RooAddPdf model_2param("model_2param","b+c+light templates",RooArgList(bPdf,lPdf),RooArgList(fb));
   RooAddPdf model_3param("model_3param","b+c+light templates",RooArgList(bPdf,cPdf,lPdf),RooArgList(fb,fc));
   if(do2templateFit){
     model_2param.fitTo(dataHist,
			RooFit::Minimizer("Minuit","migrad"),
			RooFit::SumW2Error(true),
			RooFit::Range("fit_region"));
   }
   else if(do3templateFit){
     model_3param.fitTo(dataHist,
			RooFit::Minimizer("Minuit","migrad"),
			RooFit::SumW2Error(true),
			RooFit::Range("fit_region"));
   }
   else{
     cout << "Please choose either a 2- or 3-parameter fit" << endl;
     return -1.;
   }

   double bFractionFromRooFit = fb.getVal();
   double errBFractionFromRooFit = fb.getError();
   double lFractionFromRooFit = fl.getVal();
   double errLFractionFromRooFit = fl.getError();
   double cFractionFromRooFit = 0.;
   double errCFractionFromRooFit = 0.;
   if(do3templateFit){
     cFractionFromRooFit = fc.getVal();
     errCFractionFromRooFit = cFractionFromRooFit*sqrt((errBFractionFromRooFit*errBFractionFromRooFit)/(bFractionFromRooFit*bFractionFromRooFit) + (errLFractionFromRooFit*errLFractionFromRooFit)/(lFractionFromRooFit*lFractionFromRooFit));
   }

   double p0 = bFractionFromRooFit;
   double ep0 = errBFractionFromRooFit;
   double p1 = cFractionFromRooFit;
   double ep1 = errCFractionFromRooFit;
   double p2 = lFractionFromRooFit;
   double ep2 = errLFractionFromRooFit;

   std::cout << "b-fraction from RooFit = " << bFractionFromRooFit << "+/- " << errBFractionFromRooFit << "\n";
   std::cout << "light-fraction from RooFit = " << fl.getVal() << "+/- " << fl.getError() << "\n";
   std::cout << "c-fraction from RooFit = " << fc.getVal() << "\n";

   stylizeHistograms(do2templateFit,do3templateFit);

   TH1D *h_l_sc = (TH1D*) h_l_draw->Clone("h_l_sc");
   TH1D *h_c_sc = (TH1D*) h_c_draw->Clone("h_c_sc");
   TH1D *h_b_sc = (TH1D*) h_b_draw->Clone("h_b_sc");
   TH1D *h_x_sc = (TH1D*) h_x_draw->Clone("h_x_sc");

   if(do2templateFit){
     h_b_sc->Scale(p0);
     h_l_sc->Scale(1-p0);    
   }
   else if(do3templateFit){
     h_b_sc->Scale(p0);
     h_c_sc->Scale(p1);
     h_l_sc->Scale(1-p0-p1);    
   }

   TH1D *h_l_sc2 = (TH1D*) h_l_sc->Clone("h_l_sc2");
   TH1D *h_c_sc2 = (TH1D*) h_c_sc->Clone("h_c_sc2");
   TH1D *h_b_sc2 = (TH1D*) h_b_sc->Clone("h_b_sc2");
   TH1D *h_x_sc2 = (TH1D*) h_x_sc->Clone("h_x_sc2");

   THStack *h_stack = new THStack("h_stack","");
   THStack *h_stack2 = new THStack("h_stack2","");

   h_stack->Add(h_b_sc);
   if(do3templateFit) h_stack->Add(h_c_sc);
   h_stack->Add(h_l_sc);

   h_stack2->Add(h_b_sc2);
   if(do3templateFit) h_stack2->Add(h_c_sc2);
   h_stack2->Add(h_l_sc2);

   TH1D *stack_histo_rep = (TH1D*) h_b_sc->Clone("stack_histo_rep");
   if(do3templateFit) stack_histo_rep->Add(h_c_sc);
   stack_histo_rep->Add(h_l_sc);


   TCanvas *c1 = new TCanvas("c1","c1",500,500);
   c1->cd();
   TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
   TPad *pad2 = new TPad("pad2","pad2",0.,0.0,1.0,0.3);
   pad1->Draw();
   pad1->cd();
   pad1->SetLeftMargin(0.2);
   pad1->SetBottomMargin(0.04);	
   h_stack->Draw();
   h_stack->GetXaxis()->SetLimits(0.0,3.0);
   h_stack->GetXaxis()->SetLabelSize(0);
   h_stack->GetYaxis()->SetTitleSize(0.08);
   h_stack->GetYaxis()->SetLabelSize(0.05);
   h_stack->GetYaxis()->SetTitle("Entries");
   h_stack->Draw("hist");
   h_stack2->Draw("same");
   h_draw->Draw("ep same");
   h_stack->SetMinimum(0.0);
   h_stack->SetMaximum(1.6);
   auto legend = new TLegend(0.68,0.65,0.88,0.87);
   // configure fit legend
   if(isData){
     legend->AddEntry(h_draw,"Data","p");
   }
   else{
     legend->AddEntry(h_draw,"Simulation","p");
   }
   legend->AddEntry(h_b_sc,"#font[52]{b} jets","f");
   if(do3templateFit){
     legend->AddEntry(h_c_sc,"#font[52]{c} jets","f");
     legend->AddEntry(h_l_sc,"light jets","f");
   }
   else legend->AddEntry(h_l_sc,"light+#font[52]{c} jets","f");
   legend->SetBorderSize(0);
   legend->Draw();

   TLatex *la = new TLatex();
   
   double x_t0 = 0.5;
   double x_t1 = 0.22;
   double x_t = 0.45;
   double y_t1 = 0.92;
   double y_t2 = 0.92;
   double y_t3 = 0.7;
   double y_t4 = 0.6;
   double y_t5 = 0.6;
   double y_t6 = 0.5;
   double y_t7 = 0.50;
   double y_t8 = 0.4;
   double y_t9 = 0.8;

   la->SetTextFont(62);
   la->SetTextSize(0.055);
   la->DrawLatexNDC(x_t1,y_t1,"CMS #scale[0.8]{#font[52]{Preliminary}}");
   la->SetTextFont(42);
   la->SetTextSize(0.047);
   if(isData){
     if(ispp){
       la->DrawLatexNDC(x_t0,y_t2,"#font[42]{pp 5.02 TeV (301 pb^{-1})}");
     }	  
     else if(isC1){
       la->DrawLatexNDC(x_t0,y_t2,"PbPb 0-10% 5.02 TeV (1689 #mub^{-1})");
     }
     else if(isC2){
       la->DrawLatexNDC(x_t0,y_t2,"PbPb 10-30% 5.02 TeV (1689 #mub^{-1})");
     }
     else if(isC3){
       la->DrawLatexNDC(x_t0,y_t2,"PbPb 30-50% 5.02 TeV (1689 #mub^{-1})");
     }
     else if(isC4){
       la->DrawLatexNDC(x_t0,y_t2,"PbPb 50-80% 5.02 TeV (1689 #mub^{-1})");
     }
     else{};
   }
   if(!isData) {
    if(ispp) la->DrawLatexNDC(x_t0,y_t2,"PYTHIA (5.02 TeV)");
    else if(isC1) la->DrawLatexNDC(x_t0,y_t2,"PYTHIA+HYDJET 0-10% (5.02 TeV)");
    else if(isC2) la->DrawLatexNDC(x_t0,y_t2,"PYTHIA+HYDJET 10-30% (5.02 TeV)");
    else if(isC3) la->DrawLatexNDC(x_t0,y_t2,"PYTHIA+HYDJET 30-50% (5.02 TeV)");
    else if(isC4) la->DrawLatexNDC(x_t0,y_t2,"PYTHIA+HYDJET 50-80% (5.02 TeV)");
    else{};
  }

   la->SetTextSize(0.045);
   la->DrawLatexNDC(0.23,0.71,Form("%3.0f < #font[52]{p}_{T}^{jet} < %3.0f GeV, |#it{#eta}^{jet}| < 1.6",low_jetPt,high_jetPt));
  la->DrawLatexNDC(0.23,0.65,"#font[52]{p}_{T}^{#mu} > 15 GeV, |#it{#eta}^{#it{#mu}}| < 2.0");	
  la->DrawLatexNDC(0.48,0.52,Form("#it{f}_{#it{b}} = %.4f #pm %.4f", p0, ep0));
  if(!isData) la->DrawLatexNDC(0.48,0.45,Form("#it{f}_{#it{b}}^{true} = %.4f #pm %.4f", b_frac, eb_frac));
  if(ispp){
    la->DrawLatexNDC(0.23,0.83,"Anti-#font[52]{k}_{T} PF jets, tight muon ID");
    la->DrawLatexNDC(0.23,0.77,"#font[82]{HLT_HIL3Mu12_v1}");
  }
  else{
    la->DrawLatexNDC(0.23,0.83,"Anti-#font[52]{k}_{T} CsPF jets, tight muon ID");
    la->DrawLatexNDC(0.23,0.77,"#font[82]{HLT_HIL3Mu12_v1}");
  }

  c1->cd();
  pad2->SetBottomMargin(0.4);
  pad2->SetLeftMargin(0.2);
  pad2->SetTopMargin(0.04);
  pad2->Draw();
  pad2->cd();

  fitRatio->Divide(fitRatio,stack_histo_rep,1,1,"");
  fitRatio->SetStats(0);
  fitRatio->SetTitle("");
  fitRatio->SetMarkerStyle(8);
  fitRatio->SetMarkerColor(kBlack);
  fitRatio->SetLineColor(kBlack);
  fitRatio->SetFillColorAlpha(kBlack,0.5);
  fitRatio->Draw("e1");
  fitRatio->GetYaxis()->SetTitle("Points / Fit");
  fitRatio->GetYaxis()->SetTitleSize(0.14);
  fitRatio->GetYaxis()->SetLabelSize(0.12);
  fitRatio->GetYaxis()->SetTitleOffset(0.5);
  fitRatio->GetXaxis()->SetTitle("muon #font[52]{p}_{T}^{rel} [GeV]");
  fitRatio->GetXaxis()->SetTitleOffset(1.1);
  fitRatio->GetXaxis()->SetTitleSize(.14);
  fitRatio->GetXaxis()->SetLabelSize(0.11);
  //fitRatio->GetXaxis()->SetRangeUser(0.0,4.0);
  fitRatio->SetMinimum(0.6);
  fitRatio->SetMaximum(1.4);
  fitRatio->GetYaxis()->SetNdivisions(404);
  TLine *line2 = new TLine(muRelPtAxis[0],1.0,muRelPtAxis[M-1],1.0);
  line2->SetLineStyle(7);
  line2->Draw();
  TLine *line1 = new TLine(muRelPtAxis[0],1.2,muRelPtAxis[M-1],1.2);
  line1->SetLineStyle(2);
  line1->Draw();
  TLine *line3 = new TLine(muRelPtAxis[0],0.8,muRelPtAxis[M-1],0.8);
  line3->SetLineStyle(2);
  line3->Draw();

  //// calculate chi2 by hand
  double chi2prime = 0.0;
  double chi2primeOverNDOF = 0.0;
  double observed_i = 0.0;
  double expected_i = 1.0;
  double weight_i = 0.0;

  for(int i = 1; i < M; i++){
   
    //observed_i = h0_inclR->GetBinContent(i);
    //expected_i = stack_histo_rep->GetBinContent(i);
    observed_i = fitRatio->GetBinContent(i);
    weight_i = fitRatio->GetBinError(i);
    //cout << "observed(i) = " << observed_i << endl;
    //chi2prime += (observed_i - expected_i)*(observed_i - expected_i) / (expected_i * expected_i);
    chi2prime += (observed_i - expected_i)*(observed_i - expected_i) / (weight_i * weight_i);
    
  }

  chi2primeOverNDOF = chi2prime / (1. * (M-1));

  pad1->cd();

  int Njet = h_data->GetEntries();

  if(isData) {
    la->DrawLatexNDC(0.48,0.45,Form("#chi^{2}/ndof = %3.1f",chi2primeOverNDOF));
    //t_b->DrawLatexNDC(0.48,0.45,Form("#chi^{2}/ndof = %3.1f",chi2/ndof));
    la->DrawLatexNDC(0.48,0.38,Form("#it{N}^{jet} = %i",Njet));
    //t_b->DrawLatexNDC(0.53,0.31,"#it{g}#rightarrow#it{b}#bar{#it{b}} jets increased by 20%");
    //t_b->DrawLatexNDC(0.53,0.24,"JER in MC smeared by 20%");
  }
  
  else{
   
    la->DrawLatexNDC(0.48,0.38,Form("#chi^{2}/ndof = %3.1f",chi2primeOverNDOF));
    //t_b->DrawLatexNDC(0.53,0.31,"#it{g}#rightarrow#it{b}#bar{#it{b}} jets increased by 20%");
  }

  if(ispp) c1->SaveAs(Form("../../figures/templateFits/pp/templateFit_bJetEnhanced_%3.0f-%3.0f.pdf",low_jetPt,high_jetPt));
  else c1->SaveAs(Form("../../figures/templateFits/PbPb/C%i/templateFit_bJetEnhanced_%3.0f-%3.0f.pdf",centBin,low_jetPt,high_jetPt));


  if(returnValueIndex == 1) return p0;
  else if(returnValueIndex == 2) return ep0;
  else if(returnValueIndex == 3) return chi2primeOverNDOF;
  else return -1.;

}


const int N_jetPtAxisEdges = 8;
double jetPtAxisEdges[N_jetPtAxisEdges] = {80,90,100,120,150,200,300,500};
string output_file_string = "";


void calculateBPurity(){

  double jetPtCenters[N_jetPtAxisEdges-1];
  double jetPtWidths[N_jetPtAxisEdges-1];

  for(int i = 0; i < N_jetPtAxisEdges-1; i++){
    jetPtCenters[i] = (jetPtAxisEdges[i+1] + jetPtAxisEdges[i])/2.;
    jetPtWidths[i] = (jetPtAxisEdges[i+1] - jetPtAxisEdges[i])/2;
  }

  double results_pp[N_jetPtAxisEdges-1], err_results_pp[N_jetPtAxisEdges-1];
  double results_C4[N_jetPtAxisEdges-1], err_results_C4[N_jetPtAxisEdges-1];
  double results_C3[N_jetPtAxisEdges-1], err_results_C3[N_jetPtAxisEdges-1];
  double results_C2[N_jetPtAxisEdges-1], err_results_C2[N_jetPtAxisEdges-1];
  double results_C1[N_jetPtAxisEdges-1], err_results_C1[N_jetPtAxisEdges-1];

  int do_data = 1;
  
  int do_mergeB = 1;
  int do_mergeC = 1;

  // one has to 1, one has to be 0 for the fit integers
  int do_2tempFit = 1;
  int do_3tempFit = 0; 

  
  double pTrel_fit_low  = 0.0;
  double pTrel_fit_high = 5.0;

  double c_multiplier = 1.0;
  double bGS_multiplier = 1.175;

  int do_JER_smear = 0;
  int do_JEU_up = 0;
  int do_JEU_down = 0;

  int do_centShiftUp = 0;
  int do_centShiftDown = 0;

  output_file_string = Form("../../rootFiles/bPurityResults/bPurityResults_DATA-%i_mergeB-%i_mergeC-%i_pTrel-%1.1fto%1.1f_cMult-%1.1f_bMult-%1.3f_JERsmear-%i_JEUShiftUp-%i_JEUShiftDown-%i_centShiftUp-%i_centShiftDown-%i.root",do_data,do_mergeB,do_mergeC,pTrel_fit_low,pTrel_fit_high,c_multiplier,bGS_multiplier,do_JER_smear,do_JEU_up,do_JEU_down,do_centShiftUp,do_centShiftDown);

  for(int i = 0; i < N_jetPtAxisEdges-1; i++){
    results_pp[i] = templateFitter(do_data, 1,0,0,0,0, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 1);
    err_results_pp[i] = templateFitter(do_data, 1,0,0,0,0, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 2);

    results_C4[i] = templateFitter(do_data, 0,0,0,0,1, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 1);
    err_results_C4[i] = templateFitter(do_data, 0,0,0,0,1, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 2);
    
    results_C3[i] = templateFitter(do_data, 0,0,0,1,0, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 1);
    err_results_C3[i] = templateFitter(do_data, 0,0,0,1,0, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 2);

    results_C2[i] = templateFitter(do_data, 0,0,1,0,0, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 1);
    err_results_C2[i] = templateFitter(do_data, 0,0,1,0,0, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 2);

    results_C1[i] = templateFitter(do_data, 0,1,0,0,0, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 1);
    err_results_C1[i] = templateFitter(do_data, 0,1,0,0,0, jetPtAxisEdges[i],jetPtAxisEdges[i+1], do_mergeC,do_mergeB, do_2tempFit,do_3tempFit, pTrel_fit_low,pTrel_fit_high, c_multiplier, bGS_multiplier, 2);

    
  }

  /////////////////////////////////////////////////////////////   plot the b-jet purity results
  
  TCanvas *canv = new TCanvas("canv","canv",600,600);
  canv->cd();
  TPad *pad = new TPad("pad","pad",0,0,1,1);
  pad->SetLeftMargin(0.15);
  pad->SetBottomMargin(0.15);
  pad->Draw();
  pad->cd();
	
  TMultiGraph *m = new TMultiGraph();
  TGraphErrors *g0 = new TGraphErrors(N_jetPtAxisEdges-1,jetPtCenters,results_pp,jetPtWidths,err_results_pp);
  TGraphErrors *g1 = new TGraphErrors(N_jetPtAxisEdges-1,jetPtCenters,results_C1,jetPtWidths,err_results_C1);
  TGraphErrors *g2 = new TGraphErrors(N_jetPtAxisEdges-1,jetPtCenters,results_C2,jetPtWidths,err_results_C2);
  TGraphErrors *g3 = new TGraphErrors(N_jetPtAxisEdges-1,jetPtCenters,results_C3,jetPtWidths,err_results_C3);
  TGraphErrors *g4 = new TGraphErrors(N_jetPtAxisEdges-1,jetPtCenters,results_C4,jetPtWidths,err_results_C4);

  double lw  = 2;
  double ms = 2;
  
  g0->SetLineColor(kBlack);
  g0->SetLineWidth(lw);
  g0->SetMarkerColor(kBlack);
  g0->SetMarkerStyle(24);
  g0->SetMarkerSize(ms);

  g4->SetLineColor(kBlack);
  g4->SetLineWidth(lw);
  g4->SetMarkerColor(kBlack);
  g4->SetMarkerStyle(33);
  g4->SetMarkerSize(ms);

  g3->SetLineColor(kBlue-4);
  g3->SetLineWidth(lw);
  g3->SetMarkerColor(kBlue-4);
  g3->SetMarkerStyle(21);
  g3->SetMarkerSize(ms);

  g2->SetLineColor(kGreen+2);
  g2->SetLineWidth(lw);
  g2->SetMarkerColor(kGreen+2);
  g2->SetMarkerStyle(34);
  g2->SetMarkerSize(ms);
  
  g1->SetLineColor(kRed-4);
  g1->SetLineWidth(lw);
  g1->SetMarkerColor(kRed-4);
  g1->SetMarkerStyle(47);
  g1->SetMarkerSize(ms);

  /////////////////////////////  generate histogram representations for the b-jet fraction
  
  TH1D *bFracResults_pp = new TH1D("bFracResults_pp","bFracResults_pp",N_jetPtAxisEdges-1,jetPtAxisEdges);
  TH1D *bFracResults_C4 = new TH1D("bFracResults_C4","bFracResults_C4",N_jetPtAxisEdges-1,jetPtAxisEdges);
  TH1D *bFracResults_C3 = new TH1D("bFracResults_C3","bFracResults_C3",N_jetPtAxisEdges-1,jetPtAxisEdges);
  TH1D *bFracResults_C2 = new TH1D("bFracResults_C2","bFracResults_C2",N_jetPtAxisEdges-1,jetPtAxisEdges);
  TH1D *bFracResults_C1 = new TH1D("bFracResults_C1","bFracResults_C1",N_jetPtAxisEdges-1,jetPtAxisEdges);
	
  for(int i = 0; i < N_jetPtAxisEdges + 1; i++){

    if(i==0){
      
      bFracResults_pp->SetBinContent(i,0);
      bFracResults_pp->SetBinError(i,0);
      bFracResults_C4->SetBinContent(i,0);
      bFracResults_C4->SetBinError(i,0);
      bFracResults_C3->SetBinContent(i,0);
      bFracResults_C3->SetBinError(i,0);
      bFracResults_C2->SetBinContent(i,0);
      bFracResults_C2->SetBinError(i,0);
      bFracResults_C1->SetBinContent(i,0);
      bFracResults_C1->SetBinError(i,0);

    }
    else{

      bFracResults_pp->SetBinContent(i,results_pp[i-1]);
      bFracResults_pp->SetBinError(i,err_results_pp[i-1]);
      bFracResults_C4->SetBinContent(i,results_C4[i-1]);
      bFracResults_C4->SetBinError(i,err_results_C4[i-1]);
      bFracResults_C3->SetBinContent(i,results_C3[i-1]);
      bFracResults_C3->SetBinError(i,err_results_C3[i-1]);
      bFracResults_C2->SetBinContent(i,results_C2[i-1]);
      bFracResults_C2->SetBinError(i,err_results_C2[i-1]);
      bFracResults_C1->SetBinContent(i,results_C1[i-1]);
      bFracResults_C1->SetBinError(i,err_results_C1[i-1]);
			
    }

  }

  m->Add(g0);
  m->Add(g1);
  m->Add(g2);
  m->Add(g3);
  m->Add(g4);
 

  m->GetYaxis()->SetTitleSize(0.06);
  m->GetXaxis()->SetTitleSize(0.06);
  m->GetYaxis()->SetLabelSize(0.045);
  m->GetXaxis()->SetLabelSize(0.045);
  m->GetYaxis()->SetTitle("#font[52]{b}-jet purity");
  //m->GetYaxis()->SetTitle("#it{f}_{#it{b}}");
  m->GetXaxis()->SetTitle("#font[52]{p}_{T}^{recoJet} [GeV]");

  m->GetYaxis()->SetRangeUser(0,1);

  m->Draw("AP");
  m->SetTitle("Nominal templates");

  TLegend *leg = new TLegend(0.46,0.75,0.75,0.88);
  
  leg->AddEntry(g0,"pp","p");
  leg->AddEntry(g4,"PbPb 50-80%","p");
  leg->AddEntry(g3,"PbPb 30-50%","p");
  leg->AddEntry(g2,"PbPb 10-30%","p");
  leg->AddEntry(g1,"PbPb 0-10%","p");
  
  leg->SetBorderSize(0);
  leg->SetTextSize(0.032);
  leg->Draw();

  canv->SaveAs("../../figures/bPurity/bPurity.pdf");

  TFile *bFracResults = TFile::Open(output_file_string.c_str(),"recreate");
	
  bFracResults_pp->Write();
  bFracResults_C4->Write();
  bFracResults_C3->Write();
  bFracResults_C2->Write();
  bFracResults_C1->Write();

  bFracResults->Close();

  return;

}

void correctionFactorsCalculator(bool ispp = 1, bool isC4 = 0, bool isC3 = 0, bool isC2 = 0, bool isC1 = 0){

  TFile *f1, *f2, *f0;

  if(ispp){
    f0 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIA/latest/response/PYTHIA_DiJet_response_pThat-15_mu12_pTmu-15_tight_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_2026-1-26.root");
    f1 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIA/latest/PYTHIA_DiJet_pThat-15_mu12_pTmu-15_tight_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_2026-2-2.root");
    f2 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/pp/latest/pp_HighEGJet_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_2026-3-5.root");
  }

  else{
    f0 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/response/PYTHIAHYDJET_response_DiJet_pThat-15_mu12_pTmu-14_tight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_2026-1-19.root");
    //f1 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/PYTHIAHYDJET_DiJet_pThat-30_mu12_pTmu-15_tight_mu12TriggerEfficiencyCorrection_vzReweight_hiBinReweight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_weightCut_2026-2-12.root");
    f1 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/PYTHIAHYDJET_DiJet_pThat-15_mu12_pTmu-15_tight_hiBinShift-10_leadingXjetDumpFilter_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_2026-1-14.root");
    //f1 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/final/allTemplates/PYTHIAHYDJET_DiJet_withGS_scan_mu12_tight_pTmu-14_pThat-25_hiHFcut_jetTrkMaxFilter_vzReweight_hiBinReweight_noWeightCut_removeHYDJETjet0p45_fineCentBins_projectableTemplates_allTemplates.root");
    f2 = TFile::Open("/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_HardProbes_Jet60HLT_mu12_pTmu-15to999_tight_2026-3-5.root");
  }

  TH2D *h0, *h1, *h2, *h3, *h4; // 2d flavor-pt maps
  TH2D *hh0;
  // CENT BIN 0
  TH1D *h0_incl; // inclusive projection
  TH1D *h0_b, *h0_bGS, *h0_c, *h0_d, *h0_g, *h0_s, *h0_u; // projections of each quark flavor
  TH1D *h0_bbar, *h0_cbar, *h0_dbar, *h0_sbar, *h0_ubar; // projections of each antiquark flavor
  TH1D *h0_ghost;
  TH1D *h0_noFlavor;


  TH2D *j0, *j1, *j2, *j3, *j4; // 2d flavor-eta maps
  TH2D *jj0, *jj1, *jj2, *jj3;
  // CENT BIN 0
  TH1D *j0_incl;
  TH1D *j0_b, *j0_bGS, *j0_c, *j0_d, *j0_g, *j0_s, *j0_u; // projections of each quark flavor
  TH1D *j0_bbar, *j0_cbar, *j0_dbar, *j0_sbar, *j0_ubar; // projections of each antiquark flavor
  TH1D *j0_ghost;
  TH1D *j0_noFlavor;

  TH1D *j1_incl;
  TH1D *j1_b, *j1_bGS, *j1_c, *j1_d, *j1_g, *j1_s, *j1_u; // projections of each quark flavor
  TH1D *j1_bbar, *j1_cbar, *j1_dbar, *j1_sbar, *j1_ubar; // projections of each antiquark flavor
  TH1D *j1_ghost;
  TH1D *j1_noFlavor;

  TH1D *j2_incl;
  TH1D *j2_b, *j2_bGS, *j2_c, *j2_d, *j2_g, *j2_s, *j2_u; // projections of each quark flavor
  TH1D *j2_bbar, *j2_cbar, *j2_dbar, *j2_sbar, *j2_ubar; // projections of each antiquark flavor
  TH1D *j2_ghost;
  TH1D *j2_noFlavor;
	
	
  TH1D *j3_incl;
  TH1D *j3_b, *j3_bGS, *j3_c, *j3_d, *j3_g, *j3_s, *j3_u; // projections of each quark flavor
  TH1D *j3_bbar, *j3_cbar, *j3_dbar, *j3_sbar, *j3_ubar; // projections of each antiquark flavor
  TH1D *j3_ghost;
  TH1D *j3_noFlavor;

  TH1D *J1_incl, *J2_incl;  // histograms in data
  TH1D *JJ1_incl, *JJ2_incl;

  if(ispp){
    //f0->GetObject("h_inclGenJetPt_flavor", h0);
    //f0->GetObject("h_inclGenJetPt_inclGenMuonTag_flavor",j0);
    f1->GetObject("h_inclRecoJetPt_flavor", h0);
    f1->GetObject("h_inclRecoJetPt_inclGenMuonTag_flavor",j0);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_flavor",j1);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor",j2);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor",j3);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag",J1_incl);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn",J2_incl);
  }
  else if(isC4){
    //f0->GetObject("h_inclGenJetPt_flavor_C4", h0);
    //f0->GetObject("h_inclGenJetPt_inclGenMuonTag_flavor_C4",j0);
    f1->GetObject("h_inclRecoJetPt_flavor_C4", h0);
    f1->GetObject("h_inclRecoJetPt_inclGenMuonTag_flavor_C4",j0);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_flavor_C4",j1);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor_C4",j2);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor_C4",j3);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag_C4",J1_incl);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C4",J2_incl);
  }
  else if(isC3){
    //f0->GetObject("h_inclGenJetPt_flavor_C3", h0);
    //f0->GetObject("h_inclGenJetPt_inclGenMuonTag_flavor_C3",j0);
    f1->GetObject("h_inclRecoJetPt_flavor_C3", h0);
    f1->GetObject("h_inclRecoJetPt_inclGenMuonTag_flavor_C3",j0);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_flavor_C3",j1);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor_C3",j2);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor_C3",j3);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag_C3",J1_incl);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C3",J2_incl);
  }
  else if(isC2){
    //f0->GetObject("h_inclGenJetPt_flavor_C2", h0);
    //f0->GetObject("h_inclGenJetPt_inclGenMuonTag_flavor_C2",j0);
    f1->GetObject("h_inclRecoJetPt_flavor_C2", h0);
    f1->GetObject("h_inclRecoJetPt_inclGenMuonTag_flavor_C2",j0);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_flavor_C2",j1);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor_C2",j2);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor_C2",j3);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag_C2",J1_incl);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C2",J2_incl);
  }
  else if(isC1){
    //f0->GetObject("h_inclGenJetPt_flavor_C1", h0);
    //f0->GetObject("h_inclGenJetPt_inclGenMuonTag_flavor_C1",j0);
    f1->GetObject("h_inclRecoJetPt_flavor_C1", h0);
    f1->GetObject("h_inclRecoJetPt_inclGenMuonTag_flavor_C1",j0);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_flavor_C1",j1);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor_C1",j2);
    f1->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_flavor_C1",j3);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag_C1",J1_incl);
    f2->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C1",J2_incl);
  }

  TH1D *binFinder = h0->ProjectionY();
  TAxis *yaxis = binFinder->GetXaxis();
  double smallShift = 0.01;

  h0_incl = (TH1D*) h0->ProjectionX("h0_incl",binFinder->FindBin(-5+smallShift),binFinder->FindBin(22-smallShift));
  h0_b = (TH1D*) h0->ProjectionX("h0_b",binFinder->FindBin(5+smallShift),binFinder->FindBin(6-smallShift));
  h0_bGS = (TH1D*) h0->ProjectionX("h0_bGS",binFinder->FindBin(17+smallShift),binFinder->FindBin(18-smallShift));
  h0_c = (TH1D*) h0->ProjectionX("h0_c",binFinder->FindBin(4+smallShift),binFinder->FindBin(5-smallShift));
  h0_d = (TH1D*) h0->ProjectionX("h0_d",binFinder->FindBin(1+smallShift),binFinder->FindBin(2-smallShift));
  h0_g = (TH1D*) h0->ProjectionX("h0_g",binFinder->FindBin(21+smallShift),binFinder->FindBin(22-smallShift));
  h0_s = (TH1D*) h0->ProjectionX("h0_s",binFinder->FindBin(3+smallShift),binFinder->FindBin(4-smallShift));
  h0_u = (TH1D*) h0->ProjectionX("h0_u",binFinder->FindBin(2+smallShift),binFinder->FindBin(3-smallShift));
  h0_bbar = (TH1D*) h0->ProjectionX("h0_bbar",binFinder->FindBin(-5+smallShift),binFinder->FindBin(-4-smallShift));
  h0_cbar = (TH1D*) h0->ProjectionX("h0_cbar",binFinder->FindBin(-4+smallShift),binFinder->FindBin(-3-smallShift));
  h0_dbar = (TH1D*) h0->ProjectionX("h0_dbar",binFinder->FindBin(-1+smallShift),binFinder->FindBin(0-smallShift));
  h0_sbar = (TH1D*) h0->ProjectionX("h0_sbar",binFinder->FindBin(-3+smallShift),binFinder->FindBin(-2-smallShift));
  h0_ubar = (TH1D*) h0->ProjectionX("h0_ubar",binFinder->FindBin(-2+smallShift),binFinder->FindBin(-1-smallShift));
  h0_ghost = (TH1D*) h0->ProjectionX("h0_ghost",binFinder->FindBin(0+smallShift),binFinder->FindBin(1-smallShift));
  h0_noFlavor = (TH1D*) h0->ProjectionX("h0_noFlavor",binFinder->FindBin(19+smallShift),binFinder->FindBin(20-smallShift));

  j0_incl = (TH1D*) j0->ProjectionX("j0_incl",binFinder->FindBin(-5+smallShift),binFinder->FindBin(22-smallShift));
  j0_b = (TH1D*) j0->ProjectionX("j0_b",binFinder->FindBin(5+smallShift),binFinder->FindBin(6-smallShift));
  j0_bGS = (TH1D*) j0->ProjectionX("j0_bGS",binFinder->FindBin(17+smallShift),binFinder->FindBin(18-smallShift));
  j0_c = (TH1D*) j0->ProjectionX("j0_c",binFinder->FindBin(4+smallShift),binFinder->FindBin(5-smallShift));
  j0_d = (TH1D*) j0->ProjectionX("j0_d",binFinder->FindBin(1+smallShift),binFinder->FindBin(2-smallShift));
  j0_g = (TH1D*) j0->ProjectionX("j0_g",binFinder->FindBin(21+smallShift),binFinder->FindBin(22-smallShift));
  j0_s = (TH1D*) j0->ProjectionX("j0_s",binFinder->FindBin(3+smallShift),binFinder->FindBin(4-smallShift));
  j0_u = (TH1D*) j0->ProjectionX("j0_u",binFinder->FindBin(2+smallShift),binFinder->FindBin(3-smallShift));
  j0_bbar = (TH1D*) j0->ProjectionX("j0_bbar",binFinder->FindBin(-5+smallShift),binFinder->FindBin(-4-smallShift));
  j0_cbar = (TH1D*) j0->ProjectionX("j0_cbar",binFinder->FindBin(-4+smallShift),binFinder->FindBin(-3-smallShift));
  j0_dbar = (TH1D*) j0->ProjectionX("j0_dbar",binFinder->FindBin(-1+smallShift),binFinder->FindBin(0-smallShift));
  j0_sbar = (TH1D*) j0->ProjectionX("j0_sbar",binFinder->FindBin(-3+smallShift),binFinder->FindBin(-2-smallShift));
  j0_ubar = (TH1D*) j0->ProjectionX("j0_ubar",binFinder->FindBin(-2+smallShift),binFinder->FindBin(-1-smallShift));
  j0_ghost = (TH1D*) j0->ProjectionX("j0_ghost",binFinder->FindBin(0+smallShift),binFinder->FindBin(1-smallShift));
  j0_noFlavor = (TH1D*) j0->ProjectionX("j0_noFlavor",binFinder->FindBin(19+smallShift),binFinder->FindBin(20-smallShift));
  

  j1_incl = (TH1D*) j1->ProjectionX("j1_incl",binFinder->FindBin(-5+smallShift),binFinder->FindBin(22-smallShift));
  j1_b = (TH1D*) j1->ProjectionX("j1_b",binFinder->FindBin(5+smallShift),binFinder->FindBin(6-smallShift));
  j1_bGS = (TH1D*) j1->ProjectionX("j1_bGS",binFinder->FindBin(17+smallShift),binFinder->FindBin(18-smallShift));
  j1_c = (TH1D*) j1->ProjectionX("j1_c",binFinder->FindBin(4+smallShift),binFinder->FindBin(5-smallShift));
  j1_d = (TH1D*) j1->ProjectionX("j1_d",binFinder->FindBin(1+smallShift),binFinder->FindBin(2-smallShift));
  j1_g = (TH1D*) j1->ProjectionX("j1_g",binFinder->FindBin(21+smallShift),binFinder->FindBin(22-smallShift));
  j1_s = (TH1D*) j1->ProjectionX("j1_s",binFinder->FindBin(3+smallShift),binFinder->FindBin(4-smallShift));
  j1_u = (TH1D*) j1->ProjectionX("j1_u",binFinder->FindBin(2+smallShift),binFinder->FindBin(3-smallShift));
  j1_bbar = (TH1D*) j1->ProjectionX("j1_bbar",binFinder->FindBin(-5+smallShift),binFinder->FindBin(-4-smallShift));
  j1_cbar = (TH1D*) j1->ProjectionX("j1_cbar",binFinder->FindBin(-4+smallShift),binFinder->FindBin(-3-smallShift));
  j1_dbar = (TH1D*) j1->ProjectionX("j1_dbar",binFinder->FindBin(-1+smallShift),binFinder->FindBin(0-smallShift));
  j1_sbar = (TH1D*) j1->ProjectionX("j1_sbar",binFinder->FindBin(-3+smallShift),binFinder->FindBin(-2-smallShift));
  j1_ubar = (TH1D*) j1->ProjectionX("j1_ubar",binFinder->FindBin(-2+smallShift),binFinder->FindBin(-1-smallShift));
  j1_ghost = (TH1D*) j1->ProjectionX("j1_ghost",binFinder->FindBin(0+smallShift),binFinder->FindBin(1-smallShift));
  j1_noFlavor = (TH1D*) j1->ProjectionX("j1_noFlavor",binFinder->FindBin(19+smallShift),binFinder->FindBin(20-smallShift));

  j2_incl = (TH1D*) j2->ProjectionX("j2_incl",binFinder->FindBin(-5+smallShift),binFinder->FindBin(22-smallShift));
  j2_b = (TH1D*) j2->ProjectionX("j2_b",binFinder->FindBin(5+smallShift),binFinder->FindBin(6-smallShift));
  j2_bGS = (TH1D*) j2->ProjectionX("j2_bGS",binFinder->FindBin(17+smallShift),binFinder->FindBin(18-smallShift));
  j2_c = (TH1D*) j2->ProjectionX("j2_c",binFinder->FindBin(4+smallShift),binFinder->FindBin(5-smallShift));
  j2_d = (TH1D*) j2->ProjectionX("j2_d",binFinder->FindBin(1+smallShift),binFinder->FindBin(2-smallShift));
  j2_g = (TH1D*) j2->ProjectionX("j2_g",binFinder->FindBin(21+smallShift),binFinder->FindBin(22-smallShift));
  j2_s = (TH1D*) j2->ProjectionX("j2_s",binFinder->FindBin(3+smallShift),binFinder->FindBin(4-smallShift));
  j2_u = (TH1D*) j2->ProjectionX("j2_u",binFinder->FindBin(2+smallShift),binFinder->FindBin(3-smallShift));
  j2_bbar = (TH1D*) j2->ProjectionX("j2_bbar",binFinder->FindBin(-5+smallShift),binFinder->FindBin(-4-smallShift));
  j2_cbar = (TH1D*) j2->ProjectionX("j2_cbar",binFinder->FindBin(-4+smallShift),binFinder->FindBin(-3-smallShift));
  j2_dbar = (TH1D*) j2->ProjectionX("j2_dbar",binFinder->FindBin(-1+smallShift),binFinder->FindBin(0-smallShift));
  j2_sbar = (TH1D*) j2->ProjectionX("j2_sbar",binFinder->FindBin(-3+smallShift),binFinder->FindBin(-2-smallShift));
  j2_ubar = (TH1D*) j2->ProjectionX("j2_ubar",binFinder->FindBin(-2+smallShift),binFinder->FindBin(-1-smallShift));
  j2_ghost = (TH1D*) j2->ProjectionX("j2_ghost",binFinder->FindBin(0+smallShift),binFinder->FindBin(1-smallShift));
  j2_noFlavor = (TH1D*) j2->ProjectionX("j2_noFlavor",binFinder->FindBin(19+smallShift),binFinder->FindBin(20-smallShift));

  j3_incl = (TH1D*) j3->ProjectionX("j3_incl",binFinder->FindBin(-5+smallShift),binFinder->FindBin(22-smallShift));
  j3_b = (TH1D*) j3->ProjectionX("j3_b",binFinder->FindBin(5+smallShift),binFinder->FindBin(6-smallShift));
  j3_bGS = (TH1D*) j3->ProjectionX("j3_bGS",binFinder->FindBin(17+smallShift),binFinder->FindBin(18-smallShift));
  j3_c = (TH1D*) j3->ProjectionX("j3_c",binFinder->FindBin(4+smallShift),binFinder->FindBin(5-smallShift));
  j3_d = (TH1D*) j3->ProjectionX("j3_d",binFinder->FindBin(1+smallShift),binFinder->FindBin(2-smallShift));
  j3_g = (TH1D*) j3->ProjectionX("j3_g",binFinder->FindBin(21+smallShift),binFinder->FindBin(22-smallShift));
  j3_s = (TH1D*) j3->ProjectionX("j3_s",binFinder->FindBin(3+smallShift),binFinder->FindBin(4-smallShift));
  j3_u = (TH1D*) j3->ProjectionX("j3_u",binFinder->FindBin(2+smallShift),binFinder->FindBin(3-smallShift));
  j3_bbar = (TH1D*) j3->ProjectionX("j3_bbar",binFinder->FindBin(-5+smallShift),binFinder->FindBin(-4-smallShift));
  j3_cbar = (TH1D*) j3->ProjectionX("j3_cbar",binFinder->FindBin(-4+smallShift),binFinder->FindBin(-3-smallShift));
  j3_dbar = (TH1D*) j3->ProjectionX("j3_dbar",binFinder->FindBin(-1+smallShift),binFinder->FindBin(0-smallShift));
  j3_sbar = (TH1D*) j3->ProjectionX("j3_sbar",binFinder->FindBin(-3+smallShift),binFinder->FindBin(-2-smallShift));
  j3_ubar = (TH1D*) j3->ProjectionX("j3_ubar",binFinder->FindBin(-2+smallShift),binFinder->FindBin(-1-smallShift));
  j3_ghost = (TH1D*) j3->ProjectionX("j3_ghost",binFinder->FindBin(0+smallShift),binFinder->FindBin(1-smallShift));
  j3_noFlavor = (TH1D*) j3->ProjectionX("j3_noFlavor",binFinder->FindBin(19+smallShift),binFinder->FindBin(20-smallShift));

  h0_b->Add(h0_bbar);
  h0_b->Add(h0_bGS);
  h0_c->Add(h0_cbar);
  h0_d->Add(h0_dbar);
  h0_s->Add(h0_sbar);
  h0_u->Add(h0_ubar);

  j0_b->Add(j0_bbar);
  j0_b->Add(j0_bGS);
  j0_c->Add(j0_cbar);
  j0_d->Add(j0_dbar);
  j0_s->Add(j0_sbar);
  j0_u->Add(j0_ubar);

	
  j1_b->Add(j1_bbar);
  j1_b->Add(j1_bGS);
  j1_c->Add(j1_cbar);
  j1_d->Add(j1_dbar);
  j1_s->Add(j1_sbar);
  j1_u->Add(j1_ubar);
	
  j2_b->Add(j2_bbar);
  j2_b->Add(j2_bGS);
  j2_c->Add(j2_cbar);
  j2_d->Add(j2_dbar);
  j2_s->Add(j2_sbar);
  j2_u->Add(j2_ubar);
	
  j3_b->Add(j3_bbar);
  j3_b->Add(j3_bGS);
  j3_c->Add(j3_cbar);
  j3_d->Add(j3_dbar);
  j3_s->Add(j3_sbar);
  j3_u->Add(j3_ubar);

  TH1D *h0_l;
  h0_l = (TH1D*) h0_u->Clone("h0_l");
  h0_l->Add(h0_d);
  h0_l->Add(h0_s);
  h0_l->Add(h0_g);
  h0_l->Add(h0_ghost);

  TH1D *j0_l;
  j0_l = (TH1D*) j0_u->Clone("j0_l");
  j0_l->Add(j0_d);
  j0_l->Add(j0_s);
  j0_l->Add(j0_g);
  j0_l->Add(j0_ghost);

  TH1D *j1_l;
  j1_l = (TH1D*) j1_u->Clone("j1_l");
  j1_l->Add(j1_d);
  j1_l->Add(j1_s);
  j1_l->Add(j1_g);
  j1_l->Add(j1_ghost);

  TH1D *j2_l;
  j2_l = (TH1D*) j2_u->Clone("j2_l");
  j2_l->Add(j2_d);
  j2_l->Add(j2_s);
  j2_l->Add(j2_g);
  j2_l->Add(j2_ghost);

  TH1D *j3_l;
  j3_l = (TH1D*) j3_u->Clone("j3_l");
  j3_l->Add(j3_d);
  j3_l->Add(j3_s);
  j3_l->Add(j3_g);
  j3_l->Add(j3_ghost);

  h0_incl = (TH1D*) h0_incl->Rebin(N_jetPtAxisEdges-1,"h0_incl",jetPtAxisEdges);
  h0_b = (TH1D*) h0_b->Rebin(N_jetPtAxisEdges-1,"h0_b",jetPtAxisEdges);
  h0_c = (TH1D*) h0_c->Rebin(N_jetPtAxisEdges-1,"h0_c",jetPtAxisEdges);
  h0_l = (TH1D*) h0_l->Rebin(N_jetPtAxisEdges-1,"h0_l",jetPtAxisEdges);
  h0_ghost = (TH1D*) h0_ghost->Rebin(N_jetPtAxisEdges-1,"h0_ghost",jetPtAxisEdges);
  h0_noFlavor = (TH1D*) h0_noFlavor->Rebin(N_jetPtAxisEdges-1,"h0_noFlavor",jetPtAxisEdges);
  j0_incl = (TH1D*) j0_incl->Rebin(N_jetPtAxisEdges-1,"j0_incl",jetPtAxisEdges);
  j0_b = (TH1D*) j0_b->Rebin(N_jetPtAxisEdges-1,"j0_b",jetPtAxisEdges);
  j0_c = (TH1D*) j0_c->Rebin(N_jetPtAxisEdges-1,"j0_c",jetPtAxisEdges);
  j0_l = (TH1D*) j0_l->Rebin(N_jetPtAxisEdges-1,"j0_l",jetPtAxisEdges);
  j0_ghost = (TH1D*) j0_ghost->Rebin(N_jetPtAxisEdges-1,"j0_ghost",jetPtAxisEdges);
  j0_noFlavor = (TH1D*) j0_noFlavor->Rebin(N_jetPtAxisEdges-1,"j0_noFlavor",jetPtAxisEdges);
  j1_incl = (TH1D*) j1_incl->Rebin(N_jetPtAxisEdges-1,"j1_incl",jetPtAxisEdges);
  j1_b = (TH1D*) j1_b->Rebin(N_jetPtAxisEdges-1,"j1_b",jetPtAxisEdges);
  j1_c = (TH1D*) j1_c->Rebin(N_jetPtAxisEdges-1,"j1_c",jetPtAxisEdges);
  j1_l = (TH1D*) j1_l->Rebin(N_jetPtAxisEdges-1,"j1_l",jetPtAxisEdges);
  j1_ghost = (TH1D*) j1_ghost->Rebin(N_jetPtAxisEdges-1,"j1_ghost",jetPtAxisEdges);
  j1_noFlavor = (TH1D*) j1_noFlavor->Rebin(N_jetPtAxisEdges-1,"j1_noFlavor",jetPtAxisEdges);
  j2_incl = (TH1D*) j2_incl->Rebin(N_jetPtAxisEdges-1,"j2_incl",jetPtAxisEdges);
  j2_b = (TH1D*) j2_b->Rebin(N_jetPtAxisEdges-1,"j2_b",jetPtAxisEdges);
  j2_c = (TH1D*) j2_c->Rebin(N_jetPtAxisEdges-1,"j2_c",jetPtAxisEdges);
  j2_l = (TH1D*) j2_l->Rebin(N_jetPtAxisEdges-1,"j2_l",jetPtAxisEdges);
  j2_ghost = (TH1D*) j2_ghost->Rebin(N_jetPtAxisEdges-1,"j2_ghost",jetPtAxisEdges);
  j2_noFlavor = (TH1D*) j2_noFlavor->Rebin(N_jetPtAxisEdges-1,"j2_noFlavor",jetPtAxisEdges);
  j3_incl = (TH1D*) j3_incl->Rebin(N_jetPtAxisEdges-1,"j3_incl",jetPtAxisEdges);
  j3_b = (TH1D*) j3_b->Rebin(N_jetPtAxisEdges-1,"j3_b",jetPtAxisEdges);
  j3_c = (TH1D*) j3_c->Rebin(N_jetPtAxisEdges-1,"j3_c",jetPtAxisEdges);
  j3_l = (TH1D*) j3_l->Rebin(N_jetPtAxisEdges-1,"j3_l",jetPtAxisEdges);
  j3_ghost = (TH1D*) j3_ghost->Rebin(N_jetPtAxisEdges-1,"j3_ghost",jetPtAxisEdges);
  j3_noFlavor = (TH1D*) j3_noFlavor->Rebin(N_jetPtAxisEdges-1,"j3_noFlavor",jetPtAxisEdges);
  J1_incl = (TH1D*) J1_incl->Rebin(N_jetPtAxisEdges-1,"J1_incl",jetPtAxisEdges);
  J2_incl = (TH1D*) J2_incl->Rebin(N_jetPtAxisEdges-1,"J2_incl",jetPtAxisEdges);

  divideByBinwidth(h0_incl);
  divideByBinwidth(h0_b);
  divideByBinwidth(h0_c);
  divideByBinwidth(h0_l);
  divideByBinwidth(h0_ghost);
  divideByBinwidth(h0_noFlavor);

  divideByBinwidth(j0_incl);
  divideByBinwidth(j0_b);
  divideByBinwidth(j0_c);
  divideByBinwidth(j0_l);
  divideByBinwidth(j0_ghost);
  divideByBinwidth(j0_noFlavor);

  divideByBinwidth(j1_incl);
  divideByBinwidth(j1_b);
  divideByBinwidth(j1_c);
  divideByBinwidth(j1_l);
  divideByBinwidth(j1_ghost);
  divideByBinwidth(j1_noFlavor);

  divideByBinwidth(j2_incl);
  divideByBinwidth(j2_b);
  divideByBinwidth(j2_c);
  divideByBinwidth(j2_l);
  divideByBinwidth(j2_ghost);
  divideByBinwidth(j2_noFlavor);

  divideByBinwidth(j3_incl);
  divideByBinwidth(j3_b);
  divideByBinwidth(j3_c);
  divideByBinwidth(j3_l);
  divideByBinwidth(j3_ghost);
  divideByBinwidth(j3_noFlavor);

  divideByBinwidth(J1_incl);
  divideByBinwidth(J2_incl);

  h0_b->SetLineColor(kRed);
  j0_b->SetLineColor(kGreen+1);
  j1_b->SetLineColor(kBlue);
  j2_b->SetLineColor(kCyan+2);
  j3_b->SetLineColor(kMagenta);

  h0_b->SetMarkerColor(kRed);
  j0_b->SetMarkerColor(kGreen+1);
  j1_b->SetMarkerColor(kBlue);
  j2_b->SetMarkerColor(kCyan+2);
  j3_b->SetMarkerColor(kMagenta);

	
  h0_c->SetLineColor(kRed);
  j0_c->SetLineColor(kGreen+1);
  j1_c->SetLineColor(kBlue);
  j2_c->SetLineColor(kCyan+2);
  j3_c->SetLineColor(kMagenta);

  h0_c->SetMarkerColor(kRed);
  j0_c->SetMarkerColor(kGreen+1);
  j1_c->SetMarkerColor(kBlue);
  j2_c->SetMarkerColor(kCyan+2);
  j3_c->SetMarkerColor(kMagenta);
	
  h0_l->SetLineColor(kRed);
  j0_l->SetLineColor(kGreen+1);
  j1_l->SetLineColor(kBlue);
  j2_l->SetLineColor(kCyan+2);
  j3_l->SetLineColor(kMagenta);

  double markSize = 1.2;
  int markStyle1 = 24;
  int markStyle2 = 28;
  int markStyle3 = 32;
  int markStyle4 = 46;
  int markStyle5 = 35;

  h0_b->SetMarkerSize(markSize);
  j0_b->SetMarkerSize(markSize);
  j1_b->SetMarkerSize(markSize);	
  j2_b->SetMarkerSize(markSize);
  j3_b->SetMarkerSize(markSize);

  h0_c->SetMarkerSize(markSize);
  j0_c->SetMarkerSize(markSize);
  j1_c->SetMarkerSize(markSize);	
  j2_c->SetMarkerSize(markSize);
  j3_c->SetMarkerSize(markSize);

  h0_l->SetMarkerSize(markSize);
  j0_l->SetMarkerSize(markSize);
  j1_l->SetMarkerSize(markSize);	
  j2_l->SetMarkerSize(markSize);
  j3_l->SetMarkerSize(markSize);

  h0_b->SetMarkerStyle(markStyle1);
  j0_b->SetMarkerStyle(markStyle2);
  j1_b->SetMarkerStyle(markStyle3);
  j2_b->SetMarkerStyle(markStyle4);
  j3_b->SetMarkerStyle(markStyle5);

  h0_c->SetMarkerStyle(markStyle1);
  j0_c->SetMarkerStyle(markStyle2);
  j1_c->SetMarkerStyle(markStyle3);
  j2_c->SetMarkerStyle(markStyle4);
  j3_c->SetMarkerStyle(markStyle5);

  h0_l->SetMarkerStyle(markStyle1);
  j0_l->SetMarkerStyle(markStyle2);
  j1_l->SetMarkerStyle(markStyle3);
  j2_l->SetMarkerStyle(markStyle4);
  j3_l->SetMarkerStyle(markStyle5);

  TH1D *r0_b = (TH1D*) j0_b->Clone("r0_b");
  r0_b->Divide(j0_b,h0_b,1,1,"B");
  TH1D *r1_b = (TH1D*) j1_b->Clone("r1_b");
  r1_b->Divide(j1_b,h0_b,1,1,"B");
  TH1D *r2_b = (TH1D*) j2_b->Clone("r2_b");
  r2_b->Divide(j2_b,h0_b,1,1,"B");
  TH1D *r3_b = (TH1D*) j3_b->Clone("r3_b");
  r3_b->Divide(j3_b,h0_b,1,1,"B");
	
  TH1D *r0_c = (TH1D*) j0_c->Clone("r0_c");
  r0_c->Divide(j0_c,h0_c,1,1,"B");
  TH1D *r1_c = (TH1D*) j1_c->Clone("r1_c");
  r1_c->Divide(j1_c,h0_c,1,1,"B");
  TH1D *r2_c = (TH1D*) j2_c->Clone("r2_c");
  r2_c->Divide(j2_c,h0_c,1,1,"B");
  TH1D *r3_c = (TH1D*) j3_c->Clone("r3_c");
  r3_c->Divide(j3_c,h0_c,1,1,"B");

  TH1D *r0_l = (TH1D*) j0_l->Clone("r0_l");
  r0_l->Divide(j0_l,h0_l,1,1,"B");
  TH1D *r1_l = (TH1D*) j1_l->Clone("r1_l");
  r1_l->Divide(j1_l,h0_l,1,1,"B");
  TH1D *r2_l = (TH1D*) j2_l->Clone("r2_l");
  r2_l->Divide(j2_l,h0_l,1,1,"B");
  TH1D *r3_l = (TH1D*) j3_l->Clone("r3_l");
  r3_l->Divide(j3_l,h0_l,1,1,"B");
  TH1D *r_J2 = (TH1D*) J2_incl->Clone("r_J2");
  r_J2->Divide(J2_incl,J1_incl,1,1,"B");

  TLatex *la = new TLatex();
  double laTextSize = 0.1;
  double laTitleSize = 0.05;
  int laTextFont = 42;
  la->SetTextFont(laTextFont);

  TH1D *corrFactor_1, *corrFactor_2, *corrFactor_3, *corrFactor_3_data;
  corrFactor_1 = (TH1D*) r0_b->Clone("corrFactor_1");
  corrFactor_2 = (TH1D*) j1_b->Clone("corrFactor_2");
  corrFactor_2->Divide(j1_b,j0_b,1,1,"B");
  corrFactor_3 = (TH1D*) r2_b->Clone("corrFactor_3");
  corrFactor_3->Divide(r2_b,r1_b,1,1,"B");
  corrFactor_3_data = (TH1D*) r_J2->Clone("corrFactor_3_data");
  corrFactor_1->SetTitle("bJet Correction Factor 1");

  double tSize = 0.05;
  double lSize = 0.038;

  TCanvas *canv_bJetCorr_1 = new TCanvas("canv_bJetCorr_1","bJet correction factor 1",700,700);
  canv_bJetCorr_1->cd();
  TPad *pad1_bJetCorr_1 = new TPad("pad1_bJetCorr_1","pad1_bJetCorr_1",0,0,1,1);
  pad1_bJetCorr_1->SetLeftMargin(0.2);
  pad1_bJetCorr_1->SetBottomMargin(0.2);
  pad1_bJetCorr_1->Draw();
  pad1_bJetCorr_1->cd();
  corrFactor_1->GetXaxis()->SetTitleSize(tSize);
  corrFactor_1->GetYaxis()->SetTitleSize(tSize);
  corrFactor_1->GetXaxis()->SetLabelSize(lSize);
  corrFactor_1->GetYaxis()->SetLabelSize(lSize);
  corrFactor_1->GetYaxis()->SetRangeUser(0,0.2);
  corrFactor_1->GetYaxis()->SetTitle("Correction");
  corrFactor_1->Draw();

  TCanvas *canv_bJetCorr_2 = new TCanvas("canv_bJetCorr_2","bJet correction factor 2",700,700);
  canv_bJetCorr_2->cd();
  TPad *pad1_bJetCorr_2 = new TPad("pad1_bJetCorr_2","pad1_bJetCorr_2",0,0,1,1);
  pad1_bJetCorr_2->SetLeftMargin(0.2);
  pad1_bJetCorr_2->SetBottomMargin(0.2);
  pad1_bJetCorr_2->Draw();
  pad1_bJetCorr_2->cd();
  corrFactor_2->GetXaxis()->SetTitleSize(tSize);
  corrFactor_2->GetYaxis()->SetTitleSize(tSize);
  corrFactor_2->GetXaxis()->SetLabelSize(lSize);
  corrFactor_2->GetYaxis()->SetLabelSize(lSize);
  corrFactor_2->SetStats(0);
  corrFactor_2->GetYaxis()->SetTitle("Correction");
  corrFactor_2->GetXaxis()->SetTitle("#font[52]{p}_{T}^{recoJet} [GeV]");
  corrFactor_2->SetTitle("bJet Correction Factor 2");
  corrFactor_2->GetYaxis()->SetRangeUser(0.95,1);
  corrFactor_2->Draw();

  TCanvas *canv_bJetCorr_3 = new TCanvas("canv_bJetCorr_3","bJet correction factor 3",700,700);
  canv_bJetCorr_3->cd();
  TPad *pad1_bJetCorr_3 = new TPad("pad1_bJetCorr_3","pad1_bJetCorr_3",0,0,1,1);
  pad1_bJetCorr_3->SetLeftMargin(0.2);
  pad1_bJetCorr_3->SetBottomMargin(0.2);
  pad1_bJetCorr_3->Draw();
  pad1_bJetCorr_3->cd();
  corrFactor_3->GetXaxis()->SetTitleSize(tSize);
  corrFactor_3->GetYaxis()->SetTitleSize(tSize);
  corrFactor_3->GetXaxis()->SetLabelSize(lSize);
  corrFactor_3->GetYaxis()->SetLabelSize(lSize);
  corrFactor_3->SetStats(0);
  corrFactor_3->GetYaxis()->SetTitle("Correction");
  corrFactor_3->GetXaxis()->SetTitle("#font[52]{p}_{T}^{recoJet} [GeV]");
  corrFactor_3->SetTitle("bJet Correction Factor 3");
  corrFactor_3->GetYaxis()->SetRangeUser(0.90,1);
  corrFactor_3->Draw();

  TString output_filePath = "";
  TString corrFactor_1_name = "";
  TString corrFactor_2_name = "";
  TString corrFactor_3_name = "";
  TString corrFactor_3_data_name = "";

  if(ispp){
    output_filePath = "../../rootFiles/correctionFactors/correctionFactors_pp.root";
    corrFactor_1_name = "corrFactor_1_pp";
    corrFactor_2_name = "corrFactor_2_pp";
    corrFactor_3_name = "corrFactor_3_pp";
    corrFactor_3_data_name = "corrFactor_3_data_pp";
  }
  else if(isC4){
    output_filePath = "../../rootFiles/correctionFactors/correctionFactors_C4.root";
    corrFactor_1_name = "corrFactor_1_C4";
    corrFactor_2_name = "corrFactor_2_C4";
    corrFactor_3_name = "corrFactor_3_C4";
    corrFactor_3_data_name = "corrFactor_3_data_C4";
  }
  else if(isC3){
    output_filePath = "../../rootFiles/correctionFactors/correctionFactors_C3.root";
    corrFactor_1_name = "corrFactor_1_C3";
    corrFactor_2_name = "corrFactor_2_C3";
    corrFactor_3_name = "corrFactor_3_C3";
    corrFactor_3_data_name = "corrFactor_3_data_C3";
  }
  else if(isC2){
    output_filePath = "../../rootFiles/correctionFactors/correctionFactors_C2.root";
    corrFactor_1_name = "corrFactor_1_C2";
    corrFactor_2_name = "corrFactor_2_C2";
    corrFactor_3_name = "corrFactor_3_C2";
    corrFactor_3_data_name = "corrFactor_3_data_C2";
  }
  else if(isC1){
    output_filePath = "../../rootFiles/correctionFactors/correctionFactors_C1.root";
    corrFactor_1_name = "corrFactor_1_C1";
    corrFactor_2_name = "corrFactor_2_C1";
    corrFactor_3_name = "corrFactor_3_C1";
    corrFactor_3_data_name = "corrFactor_3_data_C1";
  }

   auto wf = TFile::Open(output_filePath,"recreate");

  corrFactor_1->Write(corrFactor_1_name);
  corrFactor_2->Write(corrFactor_2_name);
  corrFactor_3->Write(corrFactor_3_name);
  corrFactor_3_data->Write(corrFactor_3_data_name);

  wf->Close();	

  return;

}

void calculateCorrectionFactors(){

  std::cout << "\n##########  GENERATING CORRECTION FACTORS ##############\n";
  std::cout << "  generating for pp...";
  correctionFactorsCalculator(1, 0, 0, 0, 0);
  std::cout << "done.\n";
  std::cout << "  generating for PbPb 50-80%...";
  correctionFactorsCalculator(0, 1, 0, 0, 0);
  std::cout << "done.\n";
  std::cout << "  generating for PbPb 30-50%...";
  correctionFactorsCalculator(0, 0, 1, 0, 0);
  std::cout << "done.\n";
  std::cout << "  generating for PbPb 10-30%...";
  correctionFactorsCalculator(0, 0, 0, 1, 0);
  std::cout << "done.\n";
  std::cout << "  generating for PbPb 0-10%...";
  correctionFactorsCalculator(0, 0, 0, 0, 1);
  std::cout << "done.\n";

  return;
  

}

void setErrorsToZero(TH1D *h){

  for(int i = 0; i < h->GetSize(); i++){
    h->SetBinError(i,0);
  }

  return;

}

void constructBJetSpectra(){

  TH1D *h_dimuonMass_pp, *h_dimuonMass_C4, *h_dimuonMass_C3, *h_dimuonMass_C2, *h_dimuonMass_C1;

  TFile *f_bPurity = TFile::Open(output_file_string.c_str());
  TH1D *bPurity_muTag_muTrig_pp, *bPurity_muTag_muTrig_C4, *bPurity_muTag_muTrig_C3, *bPurity_muTag_muTrig_C2, *bPurity_muTag_muTrig_C1;
  f_bPurity->GetObject("bFracResults_pp",bPurity_muTag_muTrig_pp);
  f_bPurity->GetObject("bFracResults_C4",bPurity_muTag_muTrig_C4);
  f_bPurity->GetObject("bFracResults_C3",bPurity_muTag_muTrig_C3);
  f_bPurity->GetObject("bFracResults_C2",bPurity_muTag_muTrig_C2);
  f_bPurity->GetObject("bFracResults_C1",bPurity_muTag_muTrig_C1);

  TFile *file_corrFactors_pp = TFile::Open("../../rootFiles/correctionFactors/correctionFactors_pp.root");
  TFile *file_corrFactors_C4 = TFile::Open("../../rootFiles/correctionFactors/correctionFactors_C4.root");
  TFile *file_corrFactors_C3 = TFile::Open("../../rootFiles/correctionFactors/correctionFactors_C3.root");
  TFile *file_corrFactors_C2 = TFile::Open("../../rootFiles/correctionFactors/correctionFactors_C2.root");
  TFile *file_corrFactors_C1 = TFile::Open("../../rootFiles/correctionFactors/correctionFactors_C1.root");

  TH1D *corrFactor_1_pp, *corrFactor_2_pp, *corrFactor_3_pp, *corrFactor_3d_pp;
  TH1D *corrFactor_1_C4, *corrFactor_2_C4, *corrFactor_3_C4, *corrFactor_3d_C4;
  TH1D *corrFactor_1_C3, *corrFactor_2_C3, *corrFactor_3_C3, *corrFactor_3d_C3;
  TH1D *corrFactor_1_C2, *corrFactor_2_C2, *corrFactor_3_C2, *corrFactor_3d_C2;
  TH1D *corrFactor_1_C1, *corrFactor_2_C1, *corrFactor_3_C1, *corrFactor_3d_C1;

  file_corrFactors_pp->GetObject("corrFactor_1_pp",corrFactor_1_pp);
  file_corrFactors_pp->GetObject("corrFactor_2_pp",corrFactor_2_pp);
  file_corrFactors_pp->GetObject("corrFactor_3_pp",corrFactor_3_pp);
  file_corrFactors_pp->GetObject("corrFactor_3_data_pp",corrFactor_3d_pp);

  file_corrFactors_C4->GetObject("corrFactor_1_C4",corrFactor_1_C4);
  file_corrFactors_C4->GetObject("corrFactor_2_C4",corrFactor_2_C4);
  file_corrFactors_C4->GetObject("corrFactor_3_C4",corrFactor_3_C4);
  file_corrFactors_C4->GetObject("corrFactor_3_data_C4",corrFactor_3d_C4);

  file_corrFactors_C3->GetObject("corrFactor_1_C3",corrFactor_1_C3);
  file_corrFactors_C3->GetObject("corrFactor_2_C3",corrFactor_2_C3);
  file_corrFactors_C3->GetObject("corrFactor_3_C3",corrFactor_3_C3);
  file_corrFactors_C3->GetObject("corrFactor_3_data_C3",corrFactor_3d_C3);

  file_corrFactors_C2->GetObject("corrFactor_1_C2",corrFactor_1_C2);
  file_corrFactors_C2->GetObject("corrFactor_2_C2",corrFactor_2_C2);
  file_corrFactors_C2->GetObject("corrFactor_3_C2",corrFactor_3_C2);
  file_corrFactors_C2->GetObject("corrFactor_3_data_C2",corrFactor_3d_C2);
  
  file_corrFactors_C1->GetObject("corrFactor_1_C1",corrFactor_1_C1);
  file_corrFactors_C1->GetObject("corrFactor_2_C1",corrFactor_2_C1);
  file_corrFactors_C1->GetObject("corrFactor_3_C1",corrFactor_3_C1);
  file_corrFactors_C1->GetObject("corrFactor_3_data_C1",corrFactor_3d_C1);

  setErrorsToZero(corrFactor_1_pp);
  setErrorsToZero(corrFactor_2_pp);
  setErrorsToZero(corrFactor_3_pp);
  setErrorsToZero(corrFactor_3d_pp);

  setErrorsToZero(corrFactor_1_C4);
  setErrorsToZero(corrFactor_2_C4);
  setErrorsToZero(corrFactor_3_C4);
  setErrorsToZero(corrFactor_3d_C4);

  setErrorsToZero(corrFactor_1_C3);
  setErrorsToZero(corrFactor_2_C3);
  setErrorsToZero(corrFactor_3_C3);
  setErrorsToZero(corrFactor_3d_C3);

  setErrorsToZero(corrFactor_1_C2);
  setErrorsToZero(corrFactor_2_C2);
  setErrorsToZero(corrFactor_3_C2);
  setErrorsToZero(corrFactor_3d_C2);

  setErrorsToZero(corrFactor_1_C1);
  setErrorsToZero(corrFactor_2_C1);
  setErrorsToZero(corrFactor_3_C1);
  setErrorsToZero(corrFactor_3d_C1);

  TH2D *H_mcJetPt_flavor_pp;
  TH1D *h_mcJetPt_pp;
  TH1D *h_mcBJetPt_pp, *h_mcBBarJetPt_pp, *h_mcBGSJetPt_pp;

  file_PYTHIA->GetObject("h_inclRecoJetPt_flavor",H_mcJetPt_flavor_pp);

  TH1D *bf = H_mcJetPt_flavor_pp->ProjectionY();
  double ss = 0.1;

  h_mcJetPt_pp = (TH1D*) H_mcJetPt_flavor_pp->ProjectionX("h_mcJetPt_pp",bf->FindBin(-5+ss),bf->FindBin(22-ss));
  h_mcBJetPt_pp = (TH1D*) H_mcJetPt_flavor_pp->ProjectionX("h_mcBJetPt_pp",bf->FindBin(5+ss),bf->FindBin(6-ss));
  h_mcBBarJetPt_pp = (TH1D*) H_mcJetPt_flavor_pp->ProjectionX("h_mcBBarJetPt_pp",bf->FindBin(-5+ss),bf->FindBin(-4-ss));
  h_mcBGSJetPt_pp = (TH1D*) H_mcJetPt_flavor_pp->ProjectionX("h_mcBGSJetPt_pp",bf->FindBin(17+ss),bf->FindBin(18-ss));
  h_mcBJetPt_pp->Add(h_mcBBarJetPt_pp);
  h_mcBJetPt_pp->Add(h_mcBGSJetPt_pp);

  TH1D *r_mcBJetPt_pp;
  r_mcBJetPt_pp = (TH1D*) h_mcBJetPt_pp->Clone("r_mcBJetPt_pp");
  r_mcBJetPt_pp->Divide(h_mcBJetPt_pp,h_mcJetPt_pp,1,1,"B");

   TH1D *h_jetPt_muTag_muTrig_pp, *h_jetPt_muTag_muTrig_C4, *h_jetPt_muTag_muTrig_C3, *h_jetPt_muTag_muTrig_C2, *h_jetPt_muTag_muTrig_C1;
  TH1D *h_jetPt_pp, *h_jetPt_C4, *h_jetPt_C3, *h_jetPt_C2, *h_jetPt_C1;

  file_pp_SingleMuon->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn",h_jetPt_muTag_muTrig_pp);
  file_PbPb_SingleMuon->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C4",h_jetPt_muTag_muTrig_C4);
  file_PbPb_SingleMuon->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C3",h_jetPt_muTag_muTrig_C3);
  file_PbPb_SingleMuon->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C2",h_jetPt_muTag_muTrig_C2);
  file_PbPb_SingleMuon->GetObject("h_inclRecoJetPt_inclRecoMuonTag_triggerOn_C1",h_jetPt_muTag_muTrig_C1);

  file_pp_MinBias->GetObject("h_inclRecoJetPt",h_jetPt_pp);
  file_PbPb_MinBias->GetObject("h_inclRecoJetPt_C4",h_jetPt_C4);
  file_PbPb_MinBias->GetObject("h_inclRecoJetPt_C3",h_jetPt_C3);
  file_PbPb_MinBias->GetObject("h_inclRecoJetPt_C2",h_jetPt_C2);
  file_PbPb_MinBias->GetObject("h_inclRecoJetPt_C1",h_jetPt_C1);

  file_pp_SingleMuon->GetObject("h_dimuonMass",h_dimuonMass_pp);
  file_PbPb_SingleMuon->GetObject("h_dimuonMass_C4",h_dimuonMass_C4);
  file_PbPb_SingleMuon->GetObject("h_dimuonMass_C3",h_dimuonMass_C3);
  file_PbPb_SingleMuon->GetObject("h_dimuonMass_C2",h_dimuonMass_C2);
  file_PbPb_SingleMuon->GetObject("h_dimuonMass_C1",h_dimuonMass_C1);

  double NZ_pp = h_dimuonMass_pp->Integral(h_dimuonMass_pp->GetXaxis()->FindBin(75),h_dimuonMass_pp->GetXaxis()->FindBin(105));
  double NZ_C4 = h_dimuonMass_C4->Integral(h_dimuonMass_C4->GetXaxis()->FindBin(75),h_dimuonMass_C4->GetXaxis()->FindBin(105));
  double NZ_C3 = h_dimuonMass_C3->Integral(h_dimuonMass_C3->GetXaxis()->FindBin(75),h_dimuonMass_C3->GetXaxis()->FindBin(105));
  double NZ_C2 = h_dimuonMass_C2->Integral(h_dimuonMass_C2->GetXaxis()->FindBin(75),h_dimuonMass_C2->GetXaxis()->FindBin(105));
  double NZ_C1 = h_dimuonMass_C1->Integral(h_dimuonMass_C1->GetXaxis()->FindBin(75),h_dimuonMass_C1->GetXaxis()->FindBin(105));

  h_jetPt_muTag_muTrig_pp = (TH1D*) h_jetPt_muTag_muTrig_pp->Rebin(N_jetPtAxisEdges-1,"h_jetPt_muTag_muTrig_pp",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_muTag_muTrig_pp);

  h_jetPt_muTag_muTrig_C4 = (TH1D*) h_jetPt_muTag_muTrig_C4->Rebin(N_jetPtAxisEdges-1,"h_jetPt_muTag_muTrig_C4",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_muTag_muTrig_C4);

  h_jetPt_muTag_muTrig_C3 = (TH1D*) h_jetPt_muTag_muTrig_C3->Rebin(N_jetPtAxisEdges-1,"h_jetPt_muTag_muTrig_C3",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_muTag_muTrig_C3);

  h_jetPt_muTag_muTrig_C2 = (TH1D*) h_jetPt_muTag_muTrig_C2->Rebin(N_jetPtAxisEdges-1,"h_jetPt_muTag_muTrig_C2",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_muTag_muTrig_C2);

  h_jetPt_muTag_muTrig_C1 = (TH1D*) h_jetPt_muTag_muTrig_C1->Rebin(N_jetPtAxisEdges-1,"h_jetPt_muTag_muTrig_C1",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_muTag_muTrig_C1);

  h_jetPt_pp = (TH1D*) h_jetPt_pp->Rebin(N_jetPtAxisEdges-1,"h_jetPt_pp",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_pp);

  h_jetPt_C4 = (TH1D*) h_jetPt_C4->Rebin(N_jetPtAxisEdges-1,"h_jetPt_C4",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_C4);

  h_jetPt_C3 = (TH1D*) h_jetPt_C3->Rebin(N_jetPtAxisEdges-1,"h_jetPt_C3",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_C3);

  h_jetPt_C2 = (TH1D*) h_jetPt_C2->Rebin(N_jetPtAxisEdges-1,"h_jetPt_C2",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_C2);

  h_jetPt_C1 = (TH1D*) h_jetPt_C1->Rebin(N_jetPtAxisEdges-1,"h_jetPt_C1",jetPtAxisEdges);
  divideByBinwidth(h_jetPt_C1);




  stylizeBJetHistograms(h_jetPt_muTag_muTrig_pp, h_jetPt_muTag_muTrig_C1, h_jetPt_muTag_muTrig_C2, h_jetPt_muTag_muTrig_C3, h_jetPt_muTag_muTrig_C4);

  
  h_jetPt_muTag_muTrig_pp->Scale(1./(NZ_pp/(0.9708*0.9708)));
  h_jetPt_muTag_muTrig_C4->Scale(1./(NZ_C4/(0.9778*0.9778)));
  h_jetPt_muTag_muTrig_C3->Scale(1./(NZ_C3/(0.9856*0.9856)));
  h_jetPt_muTag_muTrig_C2->Scale(1./(NZ_C2/(0.9069*0.9069)));
  h_jetPt_muTag_muTrig_C1->Scale(1./(NZ_C1/(0.8627*0.8627)));



  TH1D *h_bJetPt_muTag_muTrig_pp, *h_bJetPt_muTag_muTrig_C4, *h_bJetPt_muTag_muTrig_C3, *h_bJetPt_muTag_muTrig_C2, *h_bJetPt_muTag_muTrig_C1;

  h_bJetPt_muTag_muTrig_pp = (TH1D*) h_jetPt_muTag_muTrig_pp->Clone("h_bJetPt_muTag_muTrig_pp");
  h_bJetPt_muTag_muTrig_pp->Multiply(bPurity_muTag_muTrig_pp);

  h_bJetPt_muTag_muTrig_C4 = (TH1D*) h_jetPt_muTag_muTrig_C4->Clone("h_bJetPt_muTag_muTrig_C4");
  h_bJetPt_muTag_muTrig_C4->Multiply(bPurity_muTag_muTrig_C4);

  h_bJetPt_muTag_muTrig_C3 = (TH1D*) h_jetPt_muTag_muTrig_C3->Clone("h_bJetPt_muTag_muTrig_C3");
  h_bJetPt_muTag_muTrig_C3->Multiply(bPurity_muTag_muTrig_C3);

  h_bJetPt_muTag_muTrig_C2 = (TH1D*) h_jetPt_muTag_muTrig_C2->Clone("h_bJetPt_muTag_muTrig_C2");
  h_bJetPt_muTag_muTrig_C2->Multiply(bPurity_muTag_muTrig_C2);

  h_bJetPt_muTag_muTrig_C1 = (TH1D*) h_jetPt_muTag_muTrig_C1->Clone("h_bJetPt_muTag_muTrig_C1");
  h_bJetPt_muTag_muTrig_C1->Multiply(bPurity_muTag_muTrig_C1);



  TH1D *h_bJetPt_corr_pp, *h_bJetPt_corr_C4, *h_bJetPt_corr_C3, *h_bJetPt_corr_C2, *h_bJetPt_corr_C1;

  h_bJetPt_corr_pp = (TH1D*) h_bJetPt_muTag_muTrig_pp->Clone("h_bJetPt_corr_pp");
  h_bJetPt_corr_C4 = (TH1D*) h_bJetPt_muTag_muTrig_C4->Clone("h_bJetPt_corr_C4");
  h_bJetPt_corr_C3 = (TH1D*) h_bJetPt_muTag_muTrig_C3->Clone("h_bJetPt_corr_C3");  
  h_bJetPt_corr_C2 = (TH1D*) h_bJetPt_muTag_muTrig_C2->Clone("h_bJetPt_corr_C2");
  h_bJetPt_corr_C1 = (TH1D*) h_bJetPt_muTag_muTrig_C1->Clone("h_bJetPt_corr_C1");


  h_bJetPt_corr_pp->Scale(1./0.9708);
  h_bJetPt_corr_pp->Divide(corrFactor_1_pp);
  
  h_bJetPt_corr_C4->Scale(1./0.9778);
  h_bJetPt_corr_C4->Divide(corrFactor_1_C4);

  h_bJetPt_corr_C3->Scale(1./0.9856);
  h_bJetPt_corr_C3->Divide(corrFactor_1_C3);

  h_bJetPt_corr_C2->Scale(1./0.9069);
  h_bJetPt_corr_C2->Divide(corrFactor_1_C2);

  h_bJetPt_corr_C1->Scale(1./0.8627);
  h_bJetPt_corr_C1->Divide(corrFactor_1_C1);



  TH1D *bFraction_C4 = (TH1D*) h_bJetPt_corr_C4->Clone("bFraction_C4");
  bFraction_C4->Divide(bFraction_C4,h_bJetPt_corr_pp,1,1,"");

  TH1D *bFraction_C3 = (TH1D*) h_bJetPt_corr_C3->Clone("bFraction_C3");
  bFraction_C3->Divide(bFraction_C3,h_bJetPt_corr_pp,1,1,"");

  TH1D *bFraction_C2 = (TH1D*) h_bJetPt_corr_C2->Clone("bFraction_C2");
  bFraction_C2->Divide(bFraction_C2,h_bJetPt_corr_pp,1,1,"");

  TH1D *bFraction_C1 = (TH1D*) h_bJetPt_corr_C1->Clone("bFraction_C1");
  bFraction_C1->Divide(bFraction_C1,h_bJetPt_corr_pp,1,1,"");


  
  TLegend *leg = new TLegend(0.6,0.7,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  //leg->AddEntry(r_mcBJetPt_pp,"PYTHIA","l");
  //leg->AddEntry(h_jetPt_pp,"pp","p");
  leg->AddEntry(h_jetPt_C4,"PbPb 50-80%","p");  
  leg->AddEntry(h_jetPt_C3,"PbPb 30-50%","p");
  leg->AddEntry(h_jetPt_C2,"PbPb 10-30%","p");
  leg->AddEntry(h_jetPt_C1,"PbPb 0-10%","p");

  TCanvas *canv_bFraction = new TCanvas("canv_bFraction","canv_bFraction",700,700);
  canv_bFraction->cd();
  TPad *pad1_bFraction = new TPad("pad1_bFraction","pad1_bFraction",0,0,1,1);
  //pad1_bFraction->SetLogy();
  pad1_bFraction->SetLeftMargin(0.2);
  pad1_bFraction->SetBottomMargin(0.15);
  pad1_bFraction->Draw();
  pad1_bFraction->cd();
  //bFraction_C3->GetYaxis()->SetTitle("Ratio of corrected #it{b}-jet spectra (PbPb / pp)");
  bFraction_C3->GetYaxis()->SetTitleOffset(2.);
  bFraction_C3->GetYaxis()->SetTitle("#frac{1}{#it{N}_{Z}^{PbPb}} #frac{d#it{N}_{jet}^{PbPb}}{d#it{p}_{T}} #scale[3.0]{/} #frac{1}{#it{N}_{Z}^{pp}} #frac{d#it{N}_{jet}^{pp}}{d#it{p}_{T}}");
  bFraction_C3->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  bFraction_C3->GetYaxis()->SetRangeUser(0,2);
  
  bFraction_C3->GetXaxis()->SetRangeUser(80,300);
  //bFraction_C3->GetYaxis()->SetRangeUser(0,0.05);
  //bFraction_C3->GetYaxis()->SetRangeUser(0,1000);
  bFraction_C3->Draw();
  bFraction_C2->Draw("same");
  bFraction_C1->Draw("same");
  bFraction_C4->Draw("same");
  TLine *li = new TLine();
  li->SetLineStyle(7);
  li->DrawLine(80,1,300,1);
  leg->Draw();
  canv_bFraction->SaveAs("./figures/bJetPbPbToPP.pdf");
  TFile *file_BJetPbPbToPP = (TFile*) TFile::Open("./rootFiles/BJetPbPbToPP/Data/histograms_BJetPbPbToPP.root","recreate");
  bFraction_C4->Write();
  bFraction_C3->Write();
  bFraction_C2->Write();
  bFraction_C1->Write();
  file_BJetPbPbToPP->Close();

  
  return;

}


void runTemplateFitter(){

  openTemplateFitterFiles();
  calculateBPurity();
  calculateCorrectionFactors();
  constructBJetSpectra();

  return;
  
}


int main(){

  std::cout << "\n\n===== Calculate B Jets Per Z =====\n\n";

  std::cout << "Running Template Fitter...";
  runTemplateFitter();
  std::cout << "...done.\n";

  std::cout << "\n\n===== Done! =====\n\n";

  return -1;
  

}
