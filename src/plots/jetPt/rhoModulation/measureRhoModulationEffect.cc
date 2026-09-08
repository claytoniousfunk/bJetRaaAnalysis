#include <iostream>



TFile *file_noRhoModulation, *file_withRhoModulation;
TTree *tree_evt_noRhoModulation, *tree_evt_withRhoModulation;
TTree *tree_filter_noRhoModulation, *tree_filter_withRhoModulation;
TTree *tree_jet_noRhoModulation, *tree_jet_withRhoModulation, *tree_calo;

float vz_noRhoModulation = 0, vz_withRhoModulation = 0;
int hiBin_noRhoModulation = 0, hiBin_withRhoModulation = 0;

static const int NJETMAX = 9999;
int NJet_noRhoModulation = 0, NJet_withRhoModulation = 0, NCalo = 0;
float jetPt_noRhoModulation[NJETMAX], jetPt_withRhoModulation[NJETMAX], caloPt[NJETMAX];
float jetEta_noRhoModulation[NJETMAX], jetEta_withRhoModulation[NJETMAX], caloEta[NJETMAX];
float jetPhi_noRhoModulation[NJETMAX], jetPhi_withRhoModulation[NJETMAX], caloPhi[NJETMAX];

string evtBranchName = "hiEvtAnalyzer/HiTree";
string jetBranchName = "akCs4PFJetAnalyzer/t";
string caloBranchName = "akPu4CaloJetAnalyzer/t";
string filterBranchName = "skimanalysis/HltTree";

int NEvent_noRhoModulation = 0, NEvent_withRhoModulation = 0;

TH1D *h_jetPt_noRhoModulation, *h_jetPt_withRhoModulation, *h_caloPt;
TH1D *h_jetEta_noRhoModulation, *h_jetEta_withRhoModulation, *h_caloEta;
TH1D *h_jetPhi_noRhoModulation, *h_jetPhi_withRhoModulation, *h_caloPhi;
TH1D *h_vz_noRhoModulation, *h_vz_withRhoModulation;

TCanvas *canv_jetPt, *canv_jetEta, *canv_jetPhi;
TPad *pad_jetPt_upper, *pad_jetPt_lower, *pad_jetEta, *pad_jetPhi;

TLatex *la;

int hiBinCut_low = 20;
int hiBinCut_high = 60;


void initializeTH1D(){

  h_jetPt_noRhoModulation = new TH1D("h_jetPt_noRhoModulation","jet pT, no rho modulation",100,0,500);
  h_jetEta_noRhoModulation = new TH1D("h_jetEta_noRhoModulation","jet eta, no rho modulation", 100,-1.6,1.6);
  h_jetPhi_noRhoModulation = new TH1D("h_jetPhi_noRhoModulation","jet phi, no rho modulation",100,-3.14159,3.14159);

  h_jetPt_withRhoModulation = new TH1D("h_jetPt_withRhoModulation","jet pT, with rho modulation",100,0,500);
  h_jetEta_withRhoModulation = new TH1D("h_jetEta_withRhoModulation","jet eta, with rho modulation", 100,-1.6,1.6);
  h_jetPhi_withRhoModulation = new TH1D("h_jetPhi_withRhoModulation","jet phi, with rho modulation",100,-3.14159,3.14159);

  h_caloPt = new TH1D("h_caloPt","calo jet pT",100,0,500);
  h_caloEta = new TH1D("h_caloEta","calo jet eta", 100,-1.6,1.6);
  h_caloPhi = new TH1D("h_caloPhi","calo jet phi",100,-3.14159,3.14159);

  h_vz_noRhoModulation = new TH1D("h_vz_noRhoModulation","vz, no rho modulation file",60,-15,15);
  h_vz_withRhoModulation = new TH1D("h_vz_withRhoModulation","vz, with rho modulation file",60,-15,15);


}

void openFiles(){

  file_noRhoModulation = TFile::Open("~/Downloads/HiForestAOD_noRhoModulation.root");
  file_withRhoModulation = TFile::Open("~/Downloads/HiForestAOD_withRhoModulation.root");
  
}

void configureTrees(){
  
  tree_evt_noRhoModulation   = (TTree*) file_noRhoModulation->Get(evtBranchName.c_str());
  tree_evt_withRhoModulation = (TTree*) file_withRhoModulation->Get(evtBranchName.c_str());

  tree_evt_noRhoModulation->SetBranchAddress("vz",&vz_noRhoModulation);
  tree_evt_noRhoModulation->SetBranchAddress("hiBin",&hiBin_noRhoModulation);

  tree_evt_withRhoModulation->SetBranchAddress("vz",&vz_withRhoModulation);
  tree_evt_withRhoModulation->SetBranchAddress("hiBin",&hiBin_withRhoModulation);

  tree_jet_noRhoModulation   = (TTree*) file_noRhoModulation->Get(jetBranchName.c_str());
  tree_jet_withRhoModulation = (TTree*) file_withRhoModulation->Get(jetBranchName.c_str());
  tree_calo = (TTree*) file_withRhoModulation->Get(caloBranchName.c_str());

  tree_jet_noRhoModulation->SetBranchAddress("jtpt",jetPt_noRhoModulation);
  tree_jet_noRhoModulation->SetBranchAddress("jteta",jetEta_noRhoModulation);
  tree_jet_noRhoModulation->SetBranchAddress("jtphi",jetPhi_noRhoModulation);
  tree_jet_noRhoModulation->SetBranchAddress("nref",&NJet_noRhoModulation);

  tree_jet_withRhoModulation->SetBranchAddress("jtpt",jetPt_withRhoModulation);
  tree_jet_withRhoModulation->SetBranchAddress("jteta",jetEta_withRhoModulation);
  tree_jet_withRhoModulation->SetBranchAddress("jtphi",jetPhi_withRhoModulation);
  tree_jet_withRhoModulation->SetBranchAddress("nref",&NJet_withRhoModulation);

  tree_calo->SetBranchAddress("jtpt",caloPt);
  tree_calo->SetBranchAddress("jteta",caloEta);
  tree_calo->SetBranchAddress("jtphi",caloPhi);
  tree_calo->SetBranchAddress("nref",&NCalo);

}

void getNEvent(){

  NEvent_noRhoModulation = tree_evt_noRhoModulation->GetEntries();
  NEvent_withRhoModulation = tree_evt_withRhoModulation->GetEntries();
  
}

void doEventLoop_noRhoModification(){

  for(int event_i = 0; event_i < NEvent_noRhoModulation; event_i++){

    tree_evt_noRhoModulation->GetEntry(event_i);
    tree_jet_noRhoModulation->GetEntry(event_i);

    if(fabs(vz_noRhoModulation) > 15.) continue;
    if(hiBin_noRhoModulation > hiBinCut_high || hiBin_noRhoModulation < hiBinCut_low) continue;

    h_vz_noRhoModulation->Fill(vz_noRhoModulation);

    for(int jet_j = 0; jet_j < NJet_noRhoModulation; jet_j++){

      double jetPt_j = jetPt_noRhoModulation[jet_j];
      double jetEta_j = jetEta_noRhoModulation[jet_j];
      double jetPhi_j = jetPhi_noRhoModulation[jet_j];

      if(jetPt_j < 20.) continue;
      if(fabs(jetEta_j) > 1.6) continue;

      h_jetPt_noRhoModulation->Fill(jetPt_j);
      if(jetPt_j > 50.){
	h_jetEta_noRhoModulation->Fill(jetEta_j);
	h_jetPhi_noRhoModulation->Fill(jetPhi_j);
      }
      

    }
    
  }

}

void doEventLoop_withRhoModification(){

  for(int event_i = 0; event_i < NEvent_withRhoModulation; event_i++){

    tree_evt_withRhoModulation->GetEntry(event_i);
    tree_jet_withRhoModulation->GetEntry(event_i);
    tree_calo->GetEntry(event_i);

    if(fabs(vz_withRhoModulation) > 15.) continue;
    if(hiBin_withRhoModulation > hiBinCut_high || hiBin_withRhoModulation < hiBinCut_low) continue;

    h_vz_withRhoModulation->Fill(vz_withRhoModulation);

    for(int jet_j = 0; jet_j < NJet_withRhoModulation; jet_j++){

      double jetPt_j = jetPt_withRhoModulation[jet_j];
      double jetEta_j = jetEta_withRhoModulation[jet_j];
      double jetPhi_j = jetPhi_withRhoModulation[jet_j];

      if(jetPt_j < 20.) continue;
      if(fabs(jetEta_j) > 1.6) continue;

      h_jetPt_withRhoModulation->Fill(jetPt_j);
      if(jetPt_j > 50.){
	h_jetEta_withRhoModulation->Fill(jetEta_j);
	h_jetPhi_withRhoModulation->Fill(jetPhi_j);
      }

    }

    for(int jet_j = 0; jet_j < NCalo; jet_j++){

      double caloPt_j = caloPt[jet_j];
      double caloEta_j = caloEta[jet_j];
      double caloPhi_j = caloPhi[jet_j];

      if(caloPt_j < 20.) continue;
      if(fabs(caloEta_j) > 1.6) continue;

      h_caloPt->Fill(caloPt_j);
      if(caloPt_j > 50.){
	h_caloEta->Fill(caloEta_j);
	h_caloPhi->Fill(caloPhi_j);
      }

    }
    
  }

}

void measureRhoModulationEffect(){

  openFiles();
  
  configureTrees();
  
  getNEvent();

  doEventLoop_noRhoModification();

  doEventLoop_withRhoModification();

  return;

}


TH1D* divideByBinwidth(TH1D *h){

  for(int k = 0; k < h->GetSize(); k++){

    double value = h->GetBinContent(k);
    double error = h->GetBinError(k);
    double width = h->GetBinWidth(k);

    if(width != 0){

      h->SetBinContent(k,value/width);
      h->SetBinError(k,error/width);

    }

  }

  return h;

}


const int NPtEdges = 7;
double newPtAxis[NPtEdges] = {20,30,40,50,60,100,200};

void rebinPtHistograms(){
  
  h_jetPt_noRhoModulation = (TH1D*) h_jetPt_noRhoModulation->Rebin(NPtEdges-1,"h_jetPt_noRhoModulation",newPtAxis);
  h_jetPt_withRhoModulation = (TH1D*) h_jetPt_withRhoModulation->Rebin(NPtEdges-1,"h_jetPt_withRhoModulation",newPtAxis);
  h_caloPt = (TH1D*) h_caloPt->Rebin(NPtEdges-1,"h_caloPt",newPtAxis);

  return;

  divideByBinwidth(h_jetPt_noRhoModulation);
  divideByBinwidth(h_jetPt_withRhoModulation);
  divideByBinwidth(h_caloPt);

}


void normalizePtHistogramsByNEvent(){

  h_jetPt_noRhoModulation->Scale(1./h_vz_noRhoModulation->Integral());
  h_jetPt_withRhoModulation->Scale(1./h_vz_withRhoModulation->Integral());
  h_caloPt->Scale(1./h_vz_withRhoModulation->Integral());


}


void stylizeHistograms(TH1D *h1, TH1D *h2, TH1D *h3){

  h1->SetTitle("");
  h1->SetLineColor(kBlue-4);
  h2->SetLineColor(kRed-4);
  h3->SetLineColor(kGreen+2);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h3->SetLineWidth(2);

  h1->SetStats(0);
  h2->SetStats(0);
  h3->SetStats(0);

  h1->SetTitle("");
  h2->SetTitle("");
  h3->SetTitle("");

}


void drawHistograms(){

  canv_jetPt = new TCanvas("canv_jetPt","canv_jetPt",700,700);
  canv_jetPt->cd();
  pad_jetPt_upper = new TPad("pad_jetPt_upper","pad_jetPt_upper",0,0.4,1,1);
  pad_jetPt_lower = new TPad("pad_jetPt_lower","pad_jetPt_lower",0,0.,1,0.4);
  pad_jetPt_upper->SetLeftMargin(0.15);
  pad_jetPt_lower->SetLeftMargin(0.15);
  pad_jetPt_upper->SetBottomMargin(0.0);
  pad_jetPt_lower->SetBottomMargin(0.2);
  pad_jetPt_upper->SetTopMargin(0.1);
  pad_jetPt_lower->SetTopMargin(0.);
  pad_jetPt_upper->SetLogy();
  pad_jetPt_upper->Draw();
  pad_jetPt_lower->Draw();
  pad_jetPt_upper->cd();
  h_jetPt_noRhoModulation->GetYaxis()->SetRangeUser(1.e-4,2.);
  h_jetPt_noRhoModulation->GetYaxis()->SetTitleSize(0.06);
  h_jetPt_noRhoModulation->GetYaxis()->SetLabelSize(0.05);
  h_jetPt_noRhoModulation->GetYaxis()->SetTitle("1/N^{evt} dN^{jet}/dp_{T} [GeV^{-1}]");
  h_jetPt_noRhoModulation->Draw();
  h_jetPt_withRhoModulation->Draw("same");
  h_caloPt->Draw("same");
  TLegend *leg = new TLegend(0.6,0.4,0.8,0.6);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(h_jetPt_noRhoModulation,"akCs4PF, no #rho mod.");
  leg->AddEntry(h_jetPt_withRhoModulation,"akCs4PF, with #rho mod.");
  leg->AddEntry(h_caloPt,"akPu4Calo");
  leg->Draw();
  TLatex *la = new TLatex();
  la->SetTextFont(42);
  la->SetTextSize(0.05);
  la->DrawLatexNDC(0.6,0.85,"PbPb MinBias");
  la->DrawLatexNDC(0.6,0.75,Form("%i < hiBin < %i",hiBinCut_low,hiBinCut_high));
  pad_jetPt_lower->cd();
  TH1D *r1 = (TH1D*) h_jetPt_withRhoModulation->Clone("r1");
  r1->Divide(h_jetPt_noRhoModulation);
  TH1D *r2 = (TH1D*) h_caloPt->Clone("r2");
  r2->Divide(h_jetPt_noRhoModulation);
  r1->SetTitle("");
  r1->GetYaxis()->SetTitle("ratio to nom.");
  r1->GetXaxis()->SetTitle("jet p_{T} [GeV]");
  r1->GetYaxis()->SetRangeUser(0,1.5);
  r1->GetYaxis()->SetTitleSize(0.075);
  r1->GetYaxis()->SetLabelSize(0.065);
  r1->GetXaxis()->SetTitleSize(0.075);
  r1->GetXaxis()->SetLabelSize(0.065);
  r1->Draw();
  r2->Draw("same");
  TLine *li = new TLine();
  li->SetLineStyle(7);
  li->DrawLine(newPtAxis[0],1,newPtAxis[NPtEdges-1],1);

  std::cout << "<pT(no-rho-mod)> = " << h_jetPt_noRhoModulation->GetMean() << "\n";
  std::cout << "<pT(with-rho-mod)> = " << h_jetPt_withRhoModulation->GetMean() << "\n";
  std::cout << "<pT(calo)> = " << h_caloPt->GetMean() << "\n";

  canv_jetPt->SaveAs("../../../../figures/jetPt/rhoModulation/jetPt_rhoModulationEffect.pdf");
  

  canv_jetEta = new TCanvas("canv_jetEta","canv_jetEta",700,700);
  canv_jetEta->cd();
  pad_jetEta = new TPad("pad_jetEta","pad_jetEta",0,0,1,1);
  pad_jetEta->SetLeftMargin(0.15);
  pad_jetEta->SetBottomMargin(0.15);
  pad_jetEta->Draw();
  pad_jetEta->cd();
  h_jetEta_noRhoModulation->Draw();
  h_jetEta_withRhoModulation->Draw("same");
  

  delete canv_jetPt;
  return;

}




int main(){

  initializeTH1D();

  measureRhoModulationEffect();

  rebinPtHistograms();

  normalizePtHistogramsByNEvent();
  
  stylizeHistograms(h_jetPt_noRhoModulation,h_jetPt_withRhoModulation,h_caloPt);
  
  drawHistograms();

  return -1;

}
