// ROOT headers, needed because this has a main() and is built with g++ rather
// than run through the interpreter, which would supply them implicitly.
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include <iostream>



#include "../../../../headers/functions/divideByBinwidth.h"

// outDir is created if absent; all canvases are written there. Declared here
// rather than beside main() because drawRCP() saves its own figures and is
// defined earlier in the file.
const char *outDir = "../../../../figures/jetCollection/";

TFile *f_calo, *f_flow, *f_pf;
TH1D *h_calo, *h_flow, *h_pf, *r_calo, *r_flow;
TH1D *h_vz_calo, *h_vz_flow, *h_vz_pf;


int centBin = 1;

TCanvas *canv;
TPad *pad_lower, *pad_upper;

TLegend *leg;

// const int NPtEdges = 9;
// double newPtAxis[NPtEdges] = {20,30,40,50,60,80,100,150,200};

const int NPtEdges = 5;
double newPtAxis[NPtEdges] = {20,40,60,100,200};

TCanvas *canv_RCP;
TPad *pad_RCP_upper, *pad_RCP_lower;
TH1D *h_calo_C4, *h_flow_C4, *h_pf_C4;
TH1D *h_vz_calo_C4, *h_vz_flow_C4, *h_vz_pf_C4;
TH1D *RCP_calo, *RCP_flow, *RCP_pf;


double jetPtNorm = 100.;
bool normalizeByHighPt = false;
bool useRawJets = false;



void openFiles(){

  f_calo = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_noRhoModification_caloJets_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-26.root");
  f_flow = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_noRhoModification_flowJets_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-25.root");
  f_pf = TFile::Open("/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_noRhoModification_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_sameEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-25.root");

}

void fetchHistograms(){

  if(useRawJets){
    f_calo->GetObject(Form("h_inclRawJetPt_C%i",centBin),h_calo);
    f_flow->GetObject(Form("h_inclRawJetPt_C%i",centBin),h_flow);
    f_pf->GetObject(Form("h_inclRawJetPt_C%i",centBin),h_pf);
  }
  else{
    f_calo->GetObject(Form("h_inclRecoJetPt_C%i",centBin),h_calo);
    f_flow->GetObject(Form("h_inclRecoJetPt_C%i",centBin),h_flow);
    f_pf->GetObject(Form("h_inclRecoJetPt_C%i",centBin),h_pf);
  }

  f_calo->GetObject(Form("h_vz_C%i",centBin),h_vz_calo);
  f_flow->GetObject(Form("h_vz_C%i",centBin),h_vz_flow);
  f_pf->GetObject(Form("h_vz_C%i",centBin),h_vz_pf);

}

void computeRatio(){

  r_calo = (TH1D*) h_calo->Clone("r_calo");
  r_calo->Divide(h_calo, h_pf, 1, 1, "");

  r_flow = (TH1D*) h_flow->Clone("r_flow");
  r_flow->Divide(h_flow, h_pf, 1, 1, "");

}

void normalizeByNEvent(){

  h_calo->Scale(1./h_vz_calo->Integral());
  h_flow->Scale(1./h_vz_flow->Integral());
  h_pf->Scale(1./h_vz_pf->Integral());

  // The C4 histograms are the RCP denominators and must be put on the same
  // per-event footing. They were fetched along with their h_vz but never
  // scaled, so RCP was a per-event yield over a RAW yield -- wrong by the C4
  // event count (~5.8e4), which is why it read ~1e-3 instead of ~70. Worse for
  // this plot's purpose, the three files have slightly different C4 event
  // counts, so each collection was divided by a different constant and even
  // the collection-to-collection comparison was distorted.
  h_calo_C4->Scale(1./h_vz_calo_C4->Integral());
  h_flow_C4->Scale(1./h_vz_flow_C4->Integral());
  h_pf_C4->Scale(1./h_vz_pf_C4->Integral());

}



  
// void normalizeByHighPtIntegral(TH1D *h1, TH1D *h2, TH1D *h3){
  
//   h1->Scale(1./h1->Integral(h1->FindBin(jetPtNorm),h1->FindBin(500.)));
//   h2->Scale(1./h2->Integral(h2->FindBin(jetPtNorm),h2->FindBin(500.)));
//   h3->Scale(1./h3->Integral(h3->FindBin(jetPtNorm),h3->FindBin(500.)));

//   return;

// }

void normalizeByHighPtIntegral(){
  
  h_calo->Scale(1./h_calo->Integral(h_calo->FindBin(jetPtNorm),h_calo->FindBin(500.)));
  h_flow->Scale(1./h_flow->Integral(h_flow->FindBin(jetPtNorm),h_flow->FindBin(500.)));
  h_pf->Scale(1./h_pf->Integral(h_pf->FindBin(jetPtNorm),h_pf->FindBin(500.)));

  h_calo_C4->Scale(1./h_calo_C4->Integral(h_calo_C4->FindBin(jetPtNorm),h_calo_C4->FindBin(500.)));
  h_flow_C4->Scale(1./h_flow_C4->Integral(h_flow_C4->FindBin(jetPtNorm),h_flow_C4->FindBin(500.)));
  h_pf_C4->Scale(1./h_pf_C4->Integral(h_pf_C4->FindBin(jetPtNorm),h_pf_C4->FindBin(500.)));

}

// void rebin(TH1D *h1, TH1D *h2, TH1D *h3){

//   h1 = (TH1D*) h1->Rebin(NPtEdges-1,"h1",newPtAxis);
//   h2 = (TH1D*) h2->Rebin(NPtEdges-1,"h2",newPtAxis);
//   h3 = (TH1D*) h3->Rebin(NPtEdges-1,"h3",newPtAxis);

//   divideByBinwidth(h1);
//   divideByBinwidth(h2);
//   divideByBinwidth(h3);
  

//   return;
// }

void rebin(){

  h_calo = (TH1D*) h_calo->Rebin(NPtEdges-1,"h_calo",newPtAxis);
  h_flow = (TH1D*) h_flow->Rebin(NPtEdges-1,"h_flow",newPtAxis);
  h_pf = (TH1D*) h_pf->Rebin(NPtEdges-1,"h_pf",newPtAxis);

  divideByBinwidth(h_calo);
  divideByBinwidth(h_flow);
  divideByBinwidth(h_pf);


  h_calo_C4 = (TH1D*) h_calo_C4->Rebin(NPtEdges-1,"h_calo_C4",newPtAxis);
  h_flow_C4 = (TH1D*) h_flow_C4->Rebin(NPtEdges-1,"h_flow_C4",newPtAxis);
  h_pf_C4 = (TH1D*) h_pf_C4->Rebin(NPtEdges-1,"h_pf_C4",newPtAxis);

  divideByBinwidth(h_calo_C4);
  divideByBinwidth(h_flow_C4);
  divideByBinwidth(h_pf_C4);
  

}

void stylizeHistograms(){

  h_calo->SetLineColor(kGreen+2);
  h_flow->SetLineColor(kBlue-4);
  h_pf->SetLineColor(kRed-4);

  h_calo->SetLineWidth(2);
  h_flow->SetLineWidth(2);
  h_pf->SetLineWidth(2);


}

void draw(){

  canv = new TCanvas("canv","canv",700,700);
  canv->cd();
  pad_upper = new TPad("pad_upper","pad_upper",0,0.4,1,1);
  pad_lower = new TPad("pad_lower","pad_lower",0,0,1,0.4);
  pad_upper->SetLogy();
  pad_upper->SetBottomMargin(0.);
  pad_lower->SetBottomMargin(0.2);
  pad_upper->SetTopMargin(0.1);
  pad_lower->SetTopMargin(0.);
  pad_upper->SetLeftMargin(0.2);
  pad_lower->SetLeftMargin(0.2);
  pad_upper->Draw();
  pad_lower->Draw();
  pad_upper->cd();
  h_calo->SetTitle("");
  h_calo->SetStats(0);
  h_calo->GetYaxis()->SetTitleSize(0.055);
  h_calo->GetYaxis()->SetLabelSize(0.045);
  if(normalizeByHighPt) h_calo->GetYaxis()->SetTitle(Form("1/N^{jet}(p_{T} > %3.0f GeV) dN^{jet}/dp_{T} [GeV^{-1}]",jetPtNorm));
  else h_calo->GetYaxis()->SetTitle("1/N^{evt} dN^{jet}/dp_{T} [GeV^{-1}]");
  h_calo->Draw();
  h_flow->Draw("same");
  h_pf->Draw("same");
  leg = new TLegend(0.6,0.5,0.8,0.7);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(h_pf,"akCs4PF Jets");
  leg->AddEntry(h_calo,"akPu4Calo Jets");
  leg->AddEntry(h_flow,"akFlowPuCs4PF Jets");
  leg->Draw();
  pad_lower->cd();
  r_calo->SetTitle("");
  r_calo->SetStats(0);
  r_calo->GetYaxis()->SetTitleSize(0.07);
  r_calo->GetYaxis()->SetLabelSize(0.06);
  r_calo->GetYaxis()->SetTitle("Ratio to PF");
  r_calo->GetYaxis()->SetRangeUser(0,2);
  r_calo->GetXaxis()->SetTitleSize(0.07);
  r_calo->GetXaxis()->SetLabelSize(0.06);
  r_calo->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  r_calo->Draw();
  r_flow->Draw("same");
  
  

}


void defineRCPHistograms(){

  if(useRawJets){
    f_calo->GetObject("h_inclRawJetPt_C4",h_calo_C4);
    f_flow->GetObject("h_inclRawJetPt_C4",h_flow_C4);
    f_pf->GetObject("h_inclRawJetPt_C4",h_pf_C4);
  }
  else{
    f_calo->GetObject("h_inclRecoJetPt_C4",h_calo_C4);
    f_flow->GetObject("h_inclRecoJetPt_C4",h_flow_C4);
    f_pf->GetObject("h_inclRecoJetPt_C4",h_pf_C4);
  }

  f_calo->GetObject("h_vz_C4",h_vz_calo_C4);
  f_flow->GetObject("h_vz_C4",h_vz_flow_C4);
  f_pf->GetObject("h_vz_C4",h_vz_pf_C4);

}

void drawRCP(){

  RCP_calo = (TH1D*) h_calo->Clone("RCP_calo");
  RCP_flow = (TH1D*) h_flow->Clone("RCP_flow");
  RCP_pf = (TH1D*) h_pf->Clone("RCP_pf");

  RCP_calo->Divide(h_calo_C4);
  RCP_flow->Divide(h_flow_C4);
  RCP_pf->Divide(h_pf_C4);

  canv_RCP = new TCanvas("canv_RCP","canv_RCP",700,700);
  canv_RCP->cd();
  pad_RCP_upper = new TPad("pad_RCP_upper","pad_RCP_upper",0,0.,1,1);
  pad_RCP_upper->SetLeftMargin(0.19);
  pad_RCP_upper->Draw();
  pad_RCP_upper->cd();
  RCP_calo->GetYaxis()->SetTitleSize(0.045);
  RCP_calo->GetYaxis()->SetLabelSize(0.035);
  RCP_calo->GetYaxis()->SetTitle("Centrality / 50-80%");
  RCP_calo->GetXaxis()->SetTitleSize(0.045);
  RCP_calo->GetXaxis()->SetLabelSize(0.035);
  RCP_calo->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  // colours are inherited from stylizeHistograms(); without a legend the RCP
  // canvas gave no way to tell the three collections apart
  RCP_calo->SetTitle("");
  RCP_calo->GetYaxis()->SetTitleOffset(1.55);
  RCP_calo->Draw();
  RCP_flow->Draw("same");
  RCP_pf->Draw("same");

  TLegend *legRCP = new TLegend(0.55, 0.70, 0.92, 0.88);
  legRCP->SetBorderSize(0); legRCP->SetFillStyle(0); legRCP->SetTextSize(0.032);
  legRCP->AddEntry(RCP_pf,   "akCs4PF Jets",       "lp");
  legRCP->AddEntry(RCP_calo, "akPu4Calo Jets",     "lp");
  legRCP->AddEntry(RCP_flow, "akFlowPuCs4PF Jets", "lp");
  legRCP->Draw();

  TLatex laR; laR.SetNDC(); laR.SetTextFont(42); laR.SetTextSize(0.032);
  laR.DrawLatex(0.19, 0.92, Form("PbPb, C%d / 50-80%%  (no N_{coll} scaling)", centBin));

  canv_RCP->SaveAs(Form("%sjetsByCollection_RCP_C%d%s.pdf", outDir, centBin,
                        useRawJets ? "_rawPt" : ""));

  // ---- same curves, pinned to agree in the highest pT bin -----------------
  // Each RCP is divided by its OWN value there, so all three pass through 1 and
  // what remains is the pT dependence alone, with the overall offset removed.
  //
  // Note this is NOT normalizeByHighPtIntegral(): that scales the spectra, and
  // scaling numerator and denominator by their own integrals only multiplies
  // each RCP by a constant -- it does not bring the curves together.
  int bHi = RCP_calo->FindBin(jetPtNorm + 1e-6);   // the 100-200 GeV bin
  double pc = RCP_calo->GetBinContent(bHi);
  double pfl = RCP_flow->GetBinContent(bHi);
  double pp = RCP_pf->GetBinContent(bHi);
  if(pc <= 0. || pfl <= 0. || pp <= 0.){
    printf("  (pinned RCP skipped: a collection has no yield in the %.0f GeV bin)\n", jetPtNorm);
    return;
  }
  RCP_calo->Scale(1./pc);
  RCP_flow->Scale(1./pfl);
  RCP_pf->Scale(1./pp);

  TCanvas *canv_RCPpin = new TCanvas("canv_RCPpin","canv_RCPpin",700,700);
  canv_RCPpin->cd();
  TPad *padPin = new TPad("padPin","padPin",0,0.,1,1);
  padPin->SetLeftMargin(0.19);
  padPin->Draw(); padPin->cd();

  double pmax = 0.;
  for(int b = 1; b <= RCP_calo->GetNbinsX(); b++){
    pmax = TMath::Max(pmax, RCP_calo->GetBinContent(b) + RCP_calo->GetBinError(b));
    pmax = TMath::Max(pmax, RCP_flow->GetBinContent(b) + RCP_flow->GetBinError(b));
    pmax = TMath::Max(pmax, RCP_pf->GetBinContent(b)   + RCP_pf->GetBinError(b));
  }
  RCP_calo->SetTitle("");
  RCP_calo->GetYaxis()->SetTitle(Form("(C%d / 50-80%%), pinned at %.0f GeV", centBin, jetPtNorm));
  RCP_calo->GetYaxis()->SetTitleOffset(1.55);
  RCP_calo->SetMinimum(0.);
  RCP_calo->SetMaximum(pmax*1.35);
  RCP_calo->Draw();
  RCP_flow->Draw("same");
  RCP_pf->Draw("same");

  TLine *lOne = new TLine(newPtAxis[0], 1.0, newPtAxis[NPtEdges-1], 1.0);
  lOne->SetLineStyle(2); lOne->SetLineColor(kGray+2); lOne->Draw();

  legRCP->Draw();
  laR.DrawLatex(0.19, 0.92, Form("PbPb, C%d / 50-80%%, pinned at %.0f GeV", centBin, jetPtNorm));

  canv_RCPpin->SaveAs(Form("%sjetsByCollection_RCPpinned_C%d%s.pdf", outDir, centBin,
                           useRawJets ? "_rawPt" : ""));
}



// centBin comes from the command line rather than a loop: rebin() creates
// histograms by fixed name ("h_calo", ...), so a second pass in the same
// process would collide with the first. One bin per invocation sidesteps that
// entirely and lets a single bin be redone on its own.
//
//   ./plotJetsByCollection        -> C1
//   ./plotJetsByCollection 3      -> C3
int main(int argc, char **argv){

  gROOT->SetBatch(kTRUE);
  gSystem->mkdir(outDir, kTRUE);

  centBin = (argc > 1) ? atoi(argv[1]) : 1;
  if(centBin < 0 || centBin > 4){
    printf("ERROR: centBin %d out of range. This is a 4CentBins scan:\n"
           "  0 = 0-80%% inclusive, 1 = 0-10%%, 2 = 10-30%%, 3 = 30-50%%, 4 = 50-80%%\n",
           centBin);
    return 2;
  }
  printf("centBin = %d\n", centBin);
  normalizeByHighPt = false;
  useRawJets = false;
  
  openFiles();
  fetchHistograms();
  defineRCPHistograms();
  if(normalizeByHighPt) {
    normalizeByHighPtIntegral();
  }
  else normalizeByNEvent();  
  rebin();
  stylizeHistograms();
  computeRatio();
  draw();
  // RCP is this bin over C4, so for centBin 4 it is identically 1
  if(centBin != 4) drawRCP();
  else printf("  (RCP skipped: C4/C4 is identically 1)\n");

  canv->SaveAs(Form("%sjetsByCollection_C%d%s.pdf", outDir, centBin,
                    useRawJets ? "_rawPt" : ""));
  printf("wrote figures to %s\n", outDir);

  return 0;   // was -1, which reports failure to any calling script
}



