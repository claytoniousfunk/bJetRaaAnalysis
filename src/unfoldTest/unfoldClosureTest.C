//////////////////////////////////////////////////////////////////////////////////////
//
//  unfoldClosureTest.C
//
//  Split-sample closure test of the iterative Bayesian (D'Agostini) unfolding,
//  run for ONE centrality bin (or pp) at a time.
//
//    - the response matrix and the prior are built from one half of the MC
//      (even events by default)
//    - the other half (odd events) is unfolded and compared to its own
//      gen-level spectrum, so the test is statistically independent
//    - bias^2, variance, MSE and chi^2/ndf vs. N_iterations all end up on a
//      single summary plot, starting from N = 0 (the measured spectrum with no
//      unfolding applied), which is the reference the iterations improve on
//
//  The third argument selects which response is used to unfold the (always
//  pThat-weighted) spectrum; the fourth distorts the truth away from the prior so
//  that the optimal number of iterations is actually measurable.  See the comment
//  on the function itself.
//
//  usage:  root -l 'unfoldClosureTest.C("C1")'
//          root -l 'unfoldClosureTest.C("pp",12)'
//          root -l 'unfoldClosureTest.C("C1",10,"unweighted")'
//          root -l 'unfoldClosureTest.C("C1",10,"unweightedMatrix")'
//          root -l 'unfoldClosureTest.C("C1",15,"weighted","dataMC")'
//          root -l 'unfoldClosureTest.C("C1",15,"weighted","tilt:0.5")'
//          root -l 'unfoldClosureTest.C("C3",10,"weighted","none",100.,300.)'
//
//////////////////////////////////////////////////////////////////////////////////////

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif


//====================================================================================
//  sample configuration -- everything that depends on pp vs. centrality bin
//====================================================================================

struct SampleConfig {
  bool    ok    = false;
  bool    isPP  = false;
  TString label;        // short tag used in output file names
  TString title;        // pretty title drawn on the plots
  TString trainFile;    // even half, pThat-weighted   -> response + prior
  TString trainFileUnw; // even half, pThat-unweighted -> response + prior
                        // (empty when no such sample exists)
  TString testFile;     // odd half, pThat-weighted -> measured + truth, always
  TString responseName; // name of the 2D reco-vs-gen histogram

  // data, used only to derive the data/MC distortion.  Same files as the
  // original unfoldTest.C.
  TString dataMB, dataJet60, dataJet80, dataJet100;
  TString dataName;     // name of the 1D reco jet pT histogram
  double  dataMinPt;    // below this the stitched spectrum is not usable
  TString unmRecoName;  // reco jets with no gen match  -> fakes
  TString unmGenName;   // gen  jets with no reco match -> inefficiency
};

SampleConfig getSampleConfig(TString sample){

  // pp: no pThat-unweighted production exists
  const TString ppDir = "../../rootFiles/scanningOuput/PYTHIA/";
  const TString ppFile = "PYTHIA_DiJet_response_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-8-3";

  // PbPb: same selection either way, only the pThat cross-section weight differs
  const TString PHDir = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PYTHIAHYDJET/latest/response/";
  const TString PHFile = "PYTHIAHYDJET_response_DiJet_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-4-6";

  const TString PHDirUnw = "../../rootFiles/scanningOuput/PYTHIAHYDJET/";
  const TString PHFileUnw = "PYTHIAHYDJET_response_DiJet_pThat-unweighted_pThat-15_mu12_pTmu-15_tight_vzReweight_hiBinReweight_hiBinShift-10_jetTrkMaxFilter_doPThatCorrelationFilterTight_2026-7-28";

  // data, for the data/MC distortion -- the files used by the original unfoldTest.C
  const TString datDir = "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/";

  SampleConfig c;
  sample.ToLower();

  if(sample == "pp"){
    c.ok           = true;
    c.isPP         = true;
    c.label        = "pp";
    c.title        = "PYTHIA, pp";
    c.trainFile    = ppDir + ppFile + "_evenEvents.root";
    c.trainFileUnw = "";
    c.testFile     = ppDir + ppFile + "_oddEvents.root";
    c.responseName = "h_matchedRecoJetPt_genJetPt_allJets";
    c.dataMB       = datDir + "pp/latest/pp_MinBias_mu12_pTmu-14_tight_jetTrkMaxFilter_2025-10-15.root";
    c.dataJet60    = datDir + "pp/latest/pp_HighEGJet_Jet60HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_2026-3-10.root";
    c.dataJet80    = datDir + "pp/latest/pp_HighEGJet_Jet80HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_2026-3-10.root";
    c.dataJet100   = datDir + "pp/latest/pp_HighEGJet_Jet100HLT_mu12_pTmu-15to999_tight_deltaR-40_jetTrkMaxFilter_WDecayFilter_2026-3-11.root";
    c.dataName     = "h_inclRecoJetPt";
    c.dataMinPt    = 100.;   // MinBias carries the spectrum below this, poor stats
    c.unmRecoName  = "h_unmatchedRecoJetPt_allJets";
    c.unmGenName   = "h_unmatchedGenJetPt";
    return c;
  }

  TString centTitle = "";
  if     (sample == "c1") centTitle = "0-10%";
  else if(sample == "c2") centTitle = "10-30%";
  else if(sample == "c3") centTitle = "30-50%";
  else if(sample == "c4") centTitle = "50-80%";
  else{
    std::cout << "unfoldClosureTest: unknown sample \"" << sample << "\"\n"
              << "                   choose one of: pp, C1, C2, C3, C4\n";
    return c;
  }

  c.ok           = true;
  c.isPP         = false;
  c.label        = sample; c.label.ToUpper();
  c.title        = Form("PYTHIA+HYDJET %s", centTitle.Data());
  c.trainFile    = PHDir    + PHFile    + "_evenEvents.root";
  c.trainFileUnw = PHDirUnw + PHFileUnw + "_evenEvents.root";
  c.testFile     = PHDir    + PHFile    + "_oddEvents.root";
  c.responseName = Form("h_matchedRecoJetPt_genJetPt_allJets_%s", c.label.Data());
  c.dataMB       = datDir + "PbPb/latest/PbPb_MinBias_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root";
  c.dataJet60    = datDir + "PbPb/latest/PbPb_HardProbes_Jet60HLT_mu12_pTmu-15to999_tight_WDecayFilter_2026-3-11.root";
  c.dataJet80    = datDir + "PbPb/latest/PbPb_HardProbes_Jet80HLT_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root";
  c.dataJet100   = datDir + "PbPb/latest/PbPb_HardProbes_Jet100HLT_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root";
  c.dataName     = Form("h_inclRecoJetPt_%s", c.label.Data());
  // PbPb stitching sets jet60_pTmin = jet80_pTmin, so the MinBias normalisation
  // integral is an inverted range and the spectrum below 130 GeV is not usable
  c.dataMinPt    = 130.;
  c.unmRecoName  = Form("h_unmatchedRecoJetPt_allJets_%s", c.label.Data());
  c.unmGenName   = Form("h_unmatchedGenJetPt_%s",          c.label.Data());

  return c;

}


//====================================================================================
//  response loading, with the option to drop unassigned-flavour ("x") jets
//
//  x jets carry no parton flavour assignment; in the embedded PbPb sample they
//  are the background-origin category.  They are 5% of matched jets at 58 GeV in
//  C1, falling to <0.5% above 100 GeV, and are ordered by centrality
//  (C1 > C2 > C3 > C4), i.e. they behave like background rather than signal.
//
//  Rather than subtracting the xJets histogram, the response is rebuilt as the
//  sum of the five named parton flavours.  That works uniformly: in PbPb
//  allJets is exactly the sum of six categories including xJets, while pp has no
//  xJets histogram at all yet its five named flavours still fall 1.5% short of
//  allJets -- so summing removes the unassigned component in both cases.
//====================================================================================

TH2D* loadResponse(TFile *f, const TString &allName, bool removeX, const char *cloneName){

  TH2D *all = 0;
  f->GetObject(allName,all);
  if(!all) return 0;

  if(!removeX){
    TH2D *out = (TH2D*) all->Clone(cloneName);
    out->SetDirectory(0);
    return out;
  }

  const char *fl[5] = {"bJets","cJets","udJets","sJets","gJets"};
  TH2D *sum = 0;
  for(int i = 0; i < 5; i++){
    TString nm = allName; nm.ReplaceAll("allJets",fl[i]);
    TH2D *h = 0;
    f->GetObject(nm,h);
    if(!h){
      std::cout << "unfoldClosureTest: missing \"" << nm
                << "\", falling back to allJets (x jets NOT removed)\n";
      TH2D *out = (TH2D*) all->Clone(cloneName);
      out->SetDirectory(0);
      return out;
    }
    if(!sum){ sum = (TH2D*) h->Clone(cloneName); sum->SetDirectory(0); }
    else      sum->Add(h);
  }

  return sum;

}


//====================================================================================
//  projections
//
//  TH2::ProjectionX/Y default to firstbin=0, lastbin=-1, which INCLUDES the
//  under- and overflow of the other axis.  RooUnfold ignores over/underflow
//  (Overflow() is 0 by default), so the projections have to as well -- otherwise
//  the truth we compare against contains gen jets whose reco pT left the axis,
//  which the unfolded result can never reproduce.
//====================================================================================

TH1D* projX(TH2D *h, const char *name){
  return (TH1D*) h->ProjectionX(name,1,h->GetNbinsY());
}

TH1D* projY(TH2D *h, const char *name){
  return (TH1D*) h->ProjectionY(name,1,h->GetNbinsX());
}


//  Okabe-Ito qualitative palette: eight hues chosen to stay distinguishable
//  under deuteranopia, protanopia and tritanopia.  Used for the summary curves,
//  which are categorical rather than ordered.
int cbBlue()   { return TColor::GetColor("#0072B2"); }
int cbVermil() { return TColor::GetColor("#D55E00"); }
int cbGreen()  { return TColor::GetColor("#009E73"); }

//====================================================================================
//  data spectrum -> distortion function
//
//  For the distorted-truth closure we need a realistic guess at how wrong the MC
//  prior might be.  We take it from the data: stitch the trigger samples the same
//  way the original unfoldTest.C does, divide by the MC reco spectrum, and fit the
//  shape with a power law.  Only the exponent is kept -- the overall data/MC
//  normalisation is meaningless here, and a constant weight cancels in the prior.
//====================================================================================

TH1D* stitchTriggers(TH1D *mb, TH1D *j60, TH1D *j80, TH1D *j100,
                     bool isPP, const char *name){

  TH1D *out = (TH1D*) j100->Clone(name);

  const double jet80_pTmin  = 130.;
  const double jet100_pTmin = 200.;
  const double e = 0.01;
  // for PbPb the jet60 sample is not used: its window closes to zero width
  const double jet60_pTmin = isPP ? 100. : jet80_pTmin;

  double N100 = j100->Integral(j100->FindBin(jet100_pTmin+e), j100->FindBin(500.-e));
  double N80  = j80 ->Integral(j80 ->FindBin(jet100_pTmin+e), j80 ->FindBin(500.-e));

  TH1D *s80 = (TH1D*) j80->Clone(Form("%s_s80",name));
  if(N80 > 0.) s80->Scale(N100/N80);

  double N80s = s80->Integral(s80->FindBin(jet80_pTmin+e), s80->FindBin(jet100_pTmin-e));
  double N60  = j60->Integral(j60 ->FindBin(jet80_pTmin+e), j60 ->FindBin(jet100_pTmin-e));

  TH1D *s60 = (TH1D*) j60->Clone(Form("%s_s60",name));
  if(N60 > 0.) s60->Scale(N80s/N60);

  double N60s = s60->Integral(s60->FindBin(jet60_pTmin+e), s60->FindBin(jet80_pTmin-e));
  double NMB  = mb ->Integral(mb ->FindBin(jet60_pTmin+e), mb ->FindBin(jet80_pTmin-e));

  TH1D *sMB = (TH1D*) mb->Clone(Form("%s_sMB",name));
  if(NMB > 0. && N60s > 0.) sMB->Scale(N60s/NMB);

  for(int i = 0; i < out->GetSize(); i++){
    double pt = out->GetBinCenter(i);
    TH1D *src = 0;
    if     (pt < jet60_pTmin)   src = sMB;
    else if(pt < jet80_pTmin)   src = s60;
    else if(pt < jet100_pTmin)  src = s80;
    else                        src = j100;
    out->SetBinContent(i,src->GetBinContent(i));
    out->SetBinError  (i,src->GetBinError(i));
  }

  return out;

}

TH1D* loadDataSpectrum(const SampleConfig &cfg){

  TFile *fMB  = TFile::Open(cfg.dataMB);
  TFile *f60  = TFile::Open(cfg.dataJet60);
  TFile *f80  = TFile::Open(cfg.dataJet80);
  TFile *f100 = TFile::Open(cfg.dataJet100);
  if(!fMB || !f60 || !f80 || !f100 ||
     fMB->IsZombie() || f60->IsZombie() || f80->IsZombie() || f100->IsZombie()){
    std::cout << "unfoldClosureTest: could not open the data files\n";
    return 0;
  }

  TH1D *hMB, *h60, *h80, *h100;
  fMB ->GetObject(cfg.dataName,hMB);
  f60 ->GetObject(cfg.dataName,h60);
  f80 ->GetObject(cfg.dataName,h80);
  f100->GetObject(cfg.dataName,h100);
  if(!hMB || !h60 || !h80 || !h100){
    std::cout << "unfoldClosureTest: could not find \"" << cfg.dataName << "\" in the data files\n";
    return 0;
  }

  TH1D *h = stitchTriggers(hMB,h60,h80,h100,cfg.isPP,"h_data_stitched");
  h->SetDirectory(0);
  return h;

}

//  Power-law fit to (data/MC), plus the data's relative statistical error mapped
//  onto the MC binning.  The data and response histograms do not share a binning
//  (pp data is 96 bins from 20 GeV, the response 100 from 0), so everything is
//  matched by bin centre and compared as densities rather than by bin index.
//  Show where the distortion comes from: the normalised data and MC reco
//  spectra on top, their ratio and the fitted power law underneath.  The fit is
//  drawn solid inside the fit range and dashed where it is extrapolated, since
//  the weight is applied over the whole gen axis but only constrained here.
void drawDistortionDerivation(TH1D *h_data, TH1D *h_mcReco, TGraphErrors *g_ratio, TF1 *f,
                              double fitLow, double fitHigh, double pTref,
                              double expo, double expoErr, double chi2ndf,
                              const char *titleStr, const char *outPath){

  const double xLo = 40., xHi = 400.;

  TCanvas *canv = new TCanvas("canv_dist","canv_dist",750,750);
  canv->cd();
  TPad *up = new TPad("pad_dist_up","",0,0.38,1,1);
  TPad *dn = new TPad("pad_dist_dn","",0,0,1,0.38);
  up->SetLeftMargin(0.15); dn->SetLeftMargin(0.15);
  up->SetRightMargin(0.05); dn->SetRightMargin(0.05);
  up->SetBottomMargin(0.);  dn->SetBottomMargin(0.30);
  up->SetTopMargin(0.14);   dn->SetTopMargin(0.);
  up->SetLogy(); up->SetTickx(1); up->SetTicky(1);
  dn->SetTickx(1); dn->SetTicky(1);
  up->Draw(); dn->Draw();

  //---- spectra ----------------------------------------------------------------
  up->cd();
  h_data->SetStats(0); h_data->SetTitle("");
  h_data->SetLineColor(kBlack); h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerStyle(20); h_data->SetMarkerSize(0.8); h_data->SetLineWidth(2);
  h_data->GetXaxis()->SetRangeUser(xLo,xHi);
  h_data->GetYaxis()->SetTitle("normalised  (1/N) d#it{N}/d#it{p}_{T}");
  h_data->GetYaxis()->SetTitleSize(0.055); h_data->GetYaxis()->SetTitleOffset(1.25);
  h_data->GetYaxis()->SetLabelSize(0.045);
  h_data->Draw("p");
  h_mcReco->SetLineColor(cbBlue()); h_mcReco->SetLineWidth(3);
  h_mcReco->Draw("hist same");
  h_data->Draw("p same");

  TLegend *leg = new TLegend(0.58,0.62,0.93,0.84);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);
  leg->AddEntry(h_data,"data (stitched)","lp");
  leg->AddEntry(h_mcReco,"MC reco (prior)","l");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.055);
  lat.DrawLatex(0.15,0.90,titleStr);

  //---- ratio and fit ----------------------------------------------------------
  dn->cd();
  TH1F *fr = dn->DrawFrame(xLo,0.,xHi,2.4);
  fr->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
  fr->GetYaxis()->SetTitle("data / MC");
  fr->GetXaxis()->SetTitleSize(0.105); fr->GetXaxis()->SetTitleOffset(1.15);
  fr->GetYaxis()->SetTitleSize(0.095); fr->GetYaxis()->SetTitleOffset(0.70);
  fr->GetXaxis()->SetLabelSize(0.090); fr->GetYaxis()->SetLabelSize(0.090);
  fr->GetYaxis()->SetNdivisions(505);

  g_ratio->SetMarkerStyle(20); g_ratio->SetMarkerSize(0.8);
  g_ratio->SetMarkerColor(kBlack); g_ratio->SetLineColor(kBlack);
  g_ratio->Draw("P same");

  // extrapolated part, dashed
  TF1 *fx = (TF1*) f->Clone("f_dataMC_extrap");
  fx->SetRange(xLo,xHi);
  fx->SetLineColor(cbVermil()); fx->SetLineStyle(2); fx->SetLineWidth(2);
  fx->Draw("l same");
  // constrained part, solid
  TF1 *fi = (TF1*) f->Clone("f_dataMC_infit");
  fi->SetRange(fitLow,fitHigh);
  fi->SetLineColor(cbVermil()); fi->SetLineStyle(1); fi->SetLineWidth(3);
  fi->Draw("l same");

  TLine *li = new TLine(); li->SetLineStyle(7); li->SetLineColor(kGray+2);
  li->DrawLine(xLo,1.,xHi,1.);
  li->DrawLine(fitLow,0.,fitLow,2.4);
  li->DrawLine(fitHigh,0.,fitHigh,2.4);

  TLatex l2; l2.SetNDC(); l2.SetTextSize(0.075);
  l2.DrawLatex(0.19,0.88,Form("#it{w} #propto (#it{p}_{T}/%.0f)^{%.3f #pm %.3f}",pTref,expo,expoErr));
  l2.SetTextSize(0.062); l2.SetTextColor(kGray+2);
  l2.DrawLatex(0.19,0.78,Form("fit %.0f-%.0f GeV,  #chi^{2}/ndf = %.2f",fitLow,fitHigh,chi2ndf));

  canv->SaveAs(outPath);

}

bool buildDataMCDistortion(TH1D *h_data, TH1D *h_mcReco,
                           double fitLow, double fitHigh, double pTref, double dataMinPt,
                           double &expo, std::vector<double> &relErr,
                           const char *titleStr = 0, const char *outPath = 0){

  int b1 = h_data ->FindBin(fitLow), b2 = h_data ->FindBin(fitHigh);
  int m1 = h_mcReco->FindBin(fitLow), m2 = h_mcReco->FindBin(fitHigh);
  double sD = h_data ->Integral(b1,b2);
  double sM = h_mcReco->Integral(m1,m2);
  if(sD <= 0. || sM <= 0.){
    std::cout << "unfoldClosureTest: empty data or MC spectrum over the fit range\n";
    return false;
  }

  TH1D *r = (TH1D*) h_mcReco->Clone("h_dataMC_ratio");   // MC binning
  r->SetDirectory(0);
  r->Reset();

  //  normalised copies, kept only for the derivation plot
  TH1D *h_dataN = (TH1D*) h_data  ->Clone("h_data_norm");   h_dataN->SetDirectory(0);
  TH1D *h_mcN   = (TH1D*) h_mcReco->Clone("h_mcReco_norm"); h_mcN  ->SetDirectory(0);
  for(int b = 1; b <= h_dataN->GetNbinsX(); b++){
    double w = h_dataN->GetBinWidth(b);
    h_dataN->SetBinContent(b,h_dataN->GetBinContent(b)/sD/w);
    h_dataN->SetBinError  (b,h_dataN->GetBinError(b)  /sD/w);
  }
  for(int b = 1; b <= h_mcN->GetNbinsX(); b++){
    double w = h_mcN->GetBinWidth(b);
    h_mcN->SetBinContent(b,h_mcN->GetBinContent(b)/sM/w);
    h_mcN->SetBinError  (b,h_mcN->GetBinError(b)  /sM/w);
  }
  TGraphErrors *g_ratio = new TGraphErrors();

  relErr.assign(h_mcReco->GetNbinsX()+2, 0.);
  double relRef = 0.;

  for(int b = 1; b <= h_mcReco->GetNbinsX(); b++){

    double x  = h_mcReco->GetBinCenter(b);
    int    db = h_data->FindBin(x);
    if(db < 1 || db > h_data->GetNbinsX()) continue;

    double cD = h_data ->GetBinContent(db), eD = h_data ->GetBinError(db);
    double cM = h_mcReco->GetBinContent(b), eM = h_mcReco->GetBinError(b);
    if(cD <= 0. || cM <= 0.) continue;

    double rel = eD/cD;
    if(x >= dataMinPt){ relErr[b] = rel; if(relRef == 0.) relRef = rel; }

    // densities, so unequal bin widths do not bias the shape
    double dD = (cD/sD)/h_data ->GetBinWidth(db);
    double dM = (cM/sM)/h_mcReco->GetBinWidth(b);
    double ratio = dD/dM;

    double ratioErr = ratio*std::sqrt(rel*rel + (eM/cM)*(eM/cM));
    r->SetBinContent(b,ratio);
    r->SetBinError  (b,ratioErr);

    int ip = g_ratio->GetN();
    g_ratio->SetPoint(ip,x,ratio);
    g_ratio->SetPointError(ip,0.,ratioErr);

  }

  // bins with no usable data (below dataMinPt, or off the data axis) inherit the
  // precision of the lowest usable bin
  for(int b = 1; b <= h_mcReco->GetNbinsX(); b++)
    if(relErr[b] == 0.) relErr[b] = relRef;

  TF1 *f = new TF1("f_dataMC","[0]*pow(x/[2],[1])",fitLow,fitHigh);
  f->FixParameter(2,pTref);
  f->SetParameters(1.,0.3);
  TFitResultPtr res = r->Fit(f,"QNRS");
  if(!res.Get() || !res->IsValid()){
    std::cout << "unfoldClosureTest: data/MC power-law fit failed\n";
    return false;
  }

  expo = f->GetParameter(1);
  double chi2ndf = (f->GetNDF() > 0) ? f->GetChisquare()/f->GetNDF() : 0.;
  printf("  data/MC power-law fit over %.0f-%.0f GeV: exponent = %.3f +/- %.3f  (chi2/ndf = %.2f)\n",
         fitLow,fitHigh,expo,f->GetParError(1),chi2ndf);

  if(titleStr && outPath)
    drawDistortionDerivation(h_dataN,h_mcN,g_ratio,f,fitLow,fitHigh,pTref,
                             expo,f->GetParError(1),chi2ndf,titleStr,outPath);

  return true;

}

// w(pT) = (pT/pTref)^expo, frozen outside [40,500] so the extrapolation cannot
// drive a gen row to zero or blow it up
double distortionWeight(double pt, double expo, double pTref){
  double x = pt;
  if(x <  40.) x =  40.;
  if(x > 500.) x = 500.;
  return std::pow(x/pTref, expo);
}


//====================================================================================
//  closure metrics
//====================================================================================

struct ClosureMetrics {
  double chi2     = 0.;   // sum (unfolded - truth)^2 / sigma^2
  int    ndf      = 0;    // number of bins entering the sum
  double bias2    = 0.;   // sum (unfolded - truth)^2
  double variance = 0.;   // sum sigma^2
  double mse      = 0.;   // bias^2 + variance
};

ClosureMetrics computeMetrics(TH1D *h_unfold, TH1D *h_truth, double ptLow, double ptHigh){

  ClosureMetrics m;

  for(int i = 1; i <= h_unfold->GetNbinsX(); i++){

    double pt = h_unfold->GetBinCenter(i);
    if(pt < ptLow || pt > ptHigh) continue;

    double u   = h_unfold->GetBinContent(i);
    double t   = h_truth ->GetBinContent(i);
    double sig = h_unfold->GetBinError(i);

    if(t == 0.) continue;

    m.bias2    += (u - t) * (u - t);
    m.variance += sig * sig;

    if(sig > 0.){
      m.chi2 += (u - t) * (u - t) / (sig * sig);
      m.ndf++;
    }

  }

  m.mse = m.bias2 + m.variance;

  return m;

}


//====================================================================================
//  plotting helpers
//====================================================================================

//  The iterations are an ordered series, so they get a sequential palette.
//  kCividis is perceptually uniform and built for deuteranopia; the top of it is
//  a pale yellow that washes out on white, so only the lower 85% is used.
int iterColor(int iter, int N_iter_max){
  // N = 0 (the measured spectrum, no unfolding applied) is kept grey and dashed
  if(iter <= 0) return kGray+2;
  if(N_iter_max <= 1) return TColor::GetColorPalette(0);
  const double frac = 0.85;
  return TColor::GetColorPalette(int((iter-1) * frac * (TColor::GetNumberOfColors() - 1)
                                     / double(N_iter_max - 1)));
}

const char* iterLabel(int iter){
  return (iter == 0) ? "#it{N} = 0 (meas.)" : Form("#it{N} = %i",iter);
}

// spectra (top) and unfolded/truth ratio (bottom) for N = 0 ... N_iter_max
void drawClosure(int N_iter_max, TH1D **h_unfold, TH1D *h_truth,
                 double drawPtLow, double drawPtHigh,
                 double ptLow, double ptHigh,
                 double ratioLow, double ratioHigh,
                 const char *titleStr, const char *subTitleStr,
                 const char *subTitleStr2, const char *outPath){

  TLine *li = new TLine();
  li->SetLineStyle(7);

  TCanvas *canv = new TCanvas("canv_closure","canv_closure",750,750);
  canv->cd();
  TPad *pad_upper = new TPad("pad_closure_upper","",0,0.35,1,1);
  TPad *pad_lower = new TPad("pad_closure_lower","",0,0,1,0.35);
  pad_upper->SetLeftMargin(0.15);  pad_lower->SetLeftMargin(0.15);
  pad_upper->SetRightMargin(0.05); pad_lower->SetRightMargin(0.05);
  pad_upper->SetBottomMargin(0.);  pad_lower->SetBottomMargin(0.30);
  pad_upper->SetTopMargin(0.16);   pad_lower->SetTopMargin(0.);
  pad_upper->SetLogy();
  pad_upper->SetTickx(1); pad_upper->SetTicky(1);
  pad_lower->SetTickx(1); pad_lower->SetTicky(1);
  pad_upper->Draw(); pad_lower->Draw();

  TLegend *leg = new TLegend(0.60,0.32,0.93,0.84);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);
  leg->SetNColumns(2);

  //---- spectra ----------------------------------------------------------------
  pad_upper->cd();

  h_truth->SetStats(0);
  h_truth->SetTitle("");
  h_truth->SetLineColor(kBlack); h_truth->SetLineWidth(3);
  h_truth->GetXaxis()->SetRangeUser(drawPtLow,drawPtHigh);
  h_truth->GetYaxis()->SetTitle("counts");
  h_truth->GetYaxis()->SetTitleSize(0.055); h_truth->GetYaxis()->SetTitleOffset(1.3);
  h_truth->GetYaxis()->SetLabelSize(0.045);
  h_truth->Draw("hist");

  for(int n = 0; n <= N_iter_max; n++) h_unfold[n]->Draw("hist same");
  h_truth->Draw("hist same");

  leg->AddEntry(h_truth,"gen truth","l");
  for(int n = 0; n <= N_iter_max; n++) leg->AddEntry(h_unfold[n],iterLabel(n),"l");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.050);
  lat.DrawLatex(0.15,0.945,titleStr);
  lat.SetTextSize(0.030); lat.SetTextColor(kGray+2);
  lat.DrawLatex(0.15,0.885,
                (subTitleStr2 && subTitleStr2[0])
                ? Form("%s   |   %s",subTitleStr,subTitleStr2) : subTitleStr);

  //---- ratio to truth ---------------------------------------------------------
  pad_lower->cd();

  for(int n = 0; n <= N_iter_max; n++){

    TH1D *r = (TH1D*) h_unfold[n]->Clone(Form("r_closure_iter%i",n));
    r->Divide(h_unfold[n],h_truth,1,1,"");

    if(n > 0){ r->Draw("hist same"); continue; }

    r->SetStats(0);
    r->SetTitle("");
    r->GetXaxis()->SetRangeUser(drawPtLow,drawPtHigh);
    r->GetYaxis()->SetRangeUser(ratioLow,ratioHigh);
    r->GetXaxis()->SetTitle("Jet #it{p}_{T} [GeV]");
    r->GetYaxis()->SetTitle("unfolded / truth");
    r->GetXaxis()->SetTitleSize(0.100); r->GetXaxis()->SetTitleOffset(1.15);
    r->GetYaxis()->SetTitleSize(0.090); r->GetYaxis()->SetTitleOffset(0.75);
    r->GetXaxis()->SetLabelSize(0.085);
    r->GetYaxis()->SetLabelSize(0.085);
    r->GetYaxis()->SetNdivisions(505);
    r->Draw("hist");

  }

  li->DrawLine(drawPtLow,1.00,drawPtHigh,1.00);
  li->DrawLine(drawPtLow,0.95,drawPtHigh,0.95);
  li->DrawLine(drawPtLow,1.05,drawPtHigh,1.05);

  // window used for the bias/variance/chi2 metrics
  li->SetLineColor(kGray+2);
  li->DrawLine(ptLow, ratioLow,ptLow, ratioHigh);
  li->DrawLine(ptHigh,ratioLow,ptHigh,ratioHigh);
  li->SetLineColor(kBlack);

  canv->SaveAs(outPath);

}

// bias^2 / variance / MSE (left, log) and chi^2/ndf (right, log) on one canvas.
// The arrays hold N_iter_max+1 points, starting at N = 0 (no unfolding).
void drawSummary(int N_iter_max, double *x, double *bias2, double *var, double *mse,
                 double *chi2ndf, const char *titleStr, const char *subTitleStr,
                 const char *subTitleStr2, bool nOptMeaningful, const char *outPath){

  const int    N     = N_iter_max + 1;
  const double xLow  = -0.5;
  const double xHigh = N_iter_max + 0.5;
  const double dataFrac = 0.75;   // curves occupy the lower 75% of the frame,
                                  // the rest is headroom for the legend

  // ---- left axis range: only positive values survive the log scale ----------
  double yLow = 1e30, yHigh = -1e30;
  for(int i = 0; i < N; i++){
    for(double v : {bias2[i], var[i], mse[i]}){
      if(v <= 0.) continue;
      if(v < yLow)  yLow  = v;
      if(v > yHigh) yHigh = v;
    }
  }
  if(yLow > yHigh){ yLow = 1e-3; yHigh = 1.; }
  yLow  *= 0.5;
  yHigh *= 2.;
  yHigh *= std::pow(10., std::log10(yHigh/yLow) * (1./dataFrac - 1.));

  TCanvas *canv = new TCanvas("canv_summary","canv_summary",750,650);
  canv->cd();
  TPad *pad = new TPad("pad_summary","",0,0,1,1);
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.15);
  pad->SetBottomMargin(0.14); pad->SetTopMargin(0.13);
  pad->SetLogy(); pad->SetTickx(1);
  pad->Draw(); pad->cd();

  TGraph *gB = new TGraph(N,x,bias2);
  TGraph *gV = new TGraph(N,x,var);
  TGraph *gM = new TGraph(N,x,mse);
  gB->SetLineColor(cbBlue());   gB->SetMarkerColor(cbBlue());
  gV->SetLineColor(cbVermil()); gV->SetMarkerColor(cbVermil());
  gM->SetLineColor(kBlack);     gM->SetMarkerColor(kBlack);
  for(auto *g : {gB,gV,gM}){ g->SetLineWidth(2); g->SetMarkerSize(1.0); }
  // bias^2 and MSE often sit on top of each other; separate them by dash and
  // marker as well as hue
  gB->SetLineStyle(2); gB->SetMarkerStyle(24);   // open circle, dashed
  gV->SetMarkerStyle(21);                        // filled square
  gM->SetLineWidth(3); gM->SetMarkerStyle(20);   // filled circle, thick

  TH1F *frame = pad->DrawFrame(xLow,yLow,xHigh,yHigh);
  frame->GetXaxis()->SetTitle("#it{N}_{iterations}");
  frame->GetYaxis()->SetTitle("bias^{2}, variance, MSE [counts^{2}]");
  frame->GetXaxis()->SetTitleSize(0.048); frame->GetXaxis()->SetTitleOffset(1.20);
  frame->GetYaxis()->SetTitleSize(0.048); frame->GetYaxis()->SetTitleOffset(1.45);
  frame->GetXaxis()->SetLabelSize(0.042); frame->GetYaxis()->SetLabelSize(0.042);
  frame->GetXaxis()->SetNdivisions(N < 12 ? 100 + N : 510);

  // MSE is the thick line, so draw it first and let the dashed bias^2 sit on top
  // of it -- the two coincide whenever bias dominates
  gM->Draw("LP same");
  gV->Draw("LP same");
  gB->Draw("LP same");

  // ---- mark the MSE minimum (N = 0 is the no-unfolding reference, not a
  //      candidate regularisation strength, so start the search at N = 1) ----
  int iBest = 1;
  for(int i = 2; i < N; i++) if(mse[i] < mse[iBest]) iBest = i;

  TLine *li = new TLine();
  li->SetLineStyle(7); li->SetLineColor(kGray+2);
  li->DrawLine(x[iBest],yLow,x[iBest],yHigh);

  // ---- chi2/ndf on a transparent overlay pad with its own right-hand axis.
  //      Goes log if the N = 0 point drags the range over a couple of decades,
  //      which would otherwise flatten every unfolded point onto one line. ----
  double cLow = 1e30, cHigh = -1e30;
  for(int i = 0; i < N; i++){
    if(chi2ndf[i] <= 0.) continue;
    if(chi2ndf[i] < cLow)  cLow  = chi2ndf[i];
    if(chi2ndf[i] > cHigh) cHigh = chi2ndf[i];
  }
  if(cLow > cHigh){ cLow = 0.1; cHigh = 1.; }

  bool logChi2 = (cHigh / cLow) > 20.;
  if(logChi2){
    cLow  *= 0.5;
    cHigh *= 2.;
    cHigh *= std::pow(10., std::log10(cHigh/cLow) * (1./dataFrac - 1.));
  }
  else{
    cLow  -= 0.10 * (cHigh - cLow);
    cHigh  = cLow + (cHigh - cLow) / dataFrac;
  }

  canv->cd();
  TPad *overlay = new TPad("pad_summary_overlay","",0,0,1,1);
  overlay->SetFillStyle(4000); overlay->SetFrameFillStyle(0);
  overlay->SetLeftMargin(0.15); overlay->SetRightMargin(0.15);
  overlay->SetBottomMargin(0.14); overlay->SetTopMargin(0.13);
  overlay->SetLogy(logChi2);
  overlay->Draw(); overlay->cd();

  TH1F *frameOv = overlay->DrawFrame(xLow,cLow,xHigh,cHigh);
  frameOv->GetXaxis()->SetLabelOffset(999); frameOv->GetXaxis()->SetTickLength(0);
  frameOv->GetYaxis()->SetLabelOffset(999); frameOv->GetYaxis()->SetTickLength(0);

  TGraph *gC = new TGraph(N,x,chi2ndf);
  gC->SetLineColor(cbGreen()); gC->SetMarkerColor(cbGreen());
  gC->SetLineWidth(2); gC->SetLineStyle(3); gC->SetMarkerStyle(22);
  gC->SetMarkerSize(1.1);
  gC->Draw("LP same");

  TGaxis *axisC = new TGaxis(xHigh,cLow,xHigh,cHigh,cLow,cHigh,510,
                             logChi2 ? "+LG" : "+L");
  axisC->SetTitle("#it{#chi}^{2} / ndf");
  axisC->SetLineColor(cbGreen()); axisC->SetLabelColor(cbGreen()); axisC->SetTitleColor(cbGreen());
  axisC->SetTitleSize(0.048); axisC->SetTitleOffset(1.35);
  axisC->SetLabelSize(0.042);
  axisC->SetLabelFont(42);  axisC->SetTitleFont(42);
  axisC->Draw();

  TLegend *leg = new TLegend(0.19,0.64,0.86,0.80);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.033);
  leg->SetNColumns(2);
  leg->SetHeader(nOptMeaningful
                 ? Form("min MSE at #it{N} = %.0f",x[iBest])
                 : "no distortion: prior #approx truth, #it{N}_{opt} not meaningful");
  leg->AddEntry(gB,"bias^{2}","lp");
  leg->AddEntry(gV,"variance","lp");
  leg->AddEntry(gM,"MSE = bias^{2}+var","lp");
  leg->AddEntry(gC,"#it{#chi}^{2}/ndf (right)","lp");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.042);
  lat.DrawLatex(0.15,0.955,titleStr);
  lat.SetTextSize(0.031); lat.SetTextColor(kGray+2);
  lat.DrawLatex(0.15,0.905,
                (subTitleStr2 && subTitleStr2[0])
                ? Form("%s   |   %s",subTitleStr,subTitleStr2) : subTitleStr);

  canv->SaveAs(outPath);

}


//====================================================================================
//  main
//====================================================================================

void unfoldClosureTest(TString sample = "C1",              // pp, C1, C2, C3, C4
                       int    N_iter_max = 10,             // iterations to scan
                       TString responseMode = "weighted",  // see below
                       TString distortion   = "none",      // see below
                       double ptLow  =  80.,               // metric window, gen pT
                       double ptHigh = 300.){

  //  The spectrum that gets unfolded, and the truth it is compared to, are
  //  ALWAYS the pThat-weighted ones.  Only the response changes:
  //
  //   "weighted"         response matrix and prior from the pThat-weighted MC.
  //                      The standard closure test.
  //   "unweighted"       response matrix and prior from the pThat-unweighted MC.
  //                      Unfolds the weighted spectrum with an unweighted
  //                      response, exactly as it comes out of the files.
  //   "unweightedMatrix" migration probabilities P(reco|gen) from the
  //                      pThat-unweighted MC, but each gen row rescaled to the
  //                      weighted truth so the Bayesian prior is the physical
  //                      one.  Separates the effect of the matrix from the
  //                      effect of the prior, which "unweighted" conflates.
  //
  //  pThat-unweighted samples exist for PbPb only.
  //
  //  "distortion" controls what is unfolded:
  //
  //   "none"     the odd MC half is unfolded and compared to its own truth.
  //              A technical closure test.  Note that the prior is then the
  //              truth up to the even/odd split, so there is essentially no
  //              bias for the iterations to remove -- MSE is minimised at
  //              N = 1 by construction and says nothing about the optimal
  //              number of iterations.
  //   "dataMC"   the training truth is reweighted by a power law fitted to the
  //              data/MC ratio, folded through the nominal response to make a
  //              noiseless pseudo-measurement, and unfolded with the NOMINAL
  //              (unreweighted) prior.  The prior is now wrong by a realistic
  //              amount, so bias^2 genuinely falls with iterations and the MSE
  //              minimum is meaningful.
  //   "tilt"     same, but with a hand-set exponent: "tilt" (0.3) or "tilt:0.5".

  responseMode.ToLower();
  const bool useUnwMatrix     = (responseMode == "unweighted" ||
                                 responseMode == "unweightedmatrix");
  const bool reweightToTruthW = (responseMode == "unweightedmatrix");

  if(responseMode != "weighted" && !useUnwMatrix){
    std::cout << "unfoldClosureTest: unknown responseMode \"" << responseMode << "\"\n"
              << "                   choose one of: weighted, unweighted, unweightedMatrix\n";
    return;
  }

  SampleConfig cfg = getSampleConfig(sample);
  if(!cfg.ok) return;

  if(useUnwMatrix && cfg.trainFileUnw.IsNull()){
    std::cout << "unfoldClosureTest: no pThat-unweighted response exists for \""
              << cfg.label << "\" -- that production is PbPb only.\n";
    return;
  }

  const TString trainFile = useUnwMatrix ? cfg.trainFileUnw : cfg.trainFile;

  //  drop unassigned-flavour ("x") jets from the response -- see loadResponse
  const bool removeXJets = true;

  //  Fold the unmatched jets into the response.  Reco jets with no gen partner
  //  become RooUnfold "fakes"; gen jets with no reco partner become
  //  inefficiency.  Without this the response has efficiency 1 and zero fakes
  //  by construction, because measured and truth are just projections of the
  //  matched 2D -- which is self-consistent for a closure test but not what the
  //  real data contains.  RooUnfoldBayes ignores fakes unless HandleFakes(true)
  //  is set, so that is switched on below.
  const bool includeUnmatched = true;

  TString respTitle = "response: pThat-weighted";
  TString distTitle = "";
  TString outTag    = cfg.label;
  if(includeUnmatched) outTag += "_unm";
  if(responseMode == "unweighted"){
    respTitle = "response: pThat-unweighted";
    outTag    = Form("%s_unwResp",cfg.label.Data());
  }
  else if(reweightToTruthW){
    respTitle = "response: unwgt. matrix, wgt. prior";
    outTag    = Form("%s_unwMatrix",cfg.label.Data());
  }

  //---- distortion selection ---------------------------------------------------
  distortion.ToLower();
  const bool doDataMC = (distortion == "datamc");
  const bool doTilt   = distortion.BeginsWith("tilt");
  const bool doDistort = doDataMC || doTilt;

  if(distortion != "none" && !doDistort){
    std::cout << "unfoldClosureTest: unknown distortion \"" << distortion << "\"\n"
              << "                   choose one of: none, dataMC, tilt, tilt:<exponent>\n";
    return;
  }

  double tiltExpo = 0.3;
  if(doTilt && distortion.Contains(":"))
    tiltExpo = TString(distortion(distortion.Index(":")+1,distortion.Length())).Atof();

  const double drawPtLow  =  40.;   // x-range of the closure plot
  const double drawPtHigh = 400.;
  const double ratioLow   = doDistort ? 0.5 : 0.8;  // closure ratio panel range
  const double ratioHigh  = doDistort ? 1.5 : 1.2;
  const double pTref      = 150.;   // w(pTref) = 1
  const TString outDir = "../../figures/unfoldTest/";

  gStyle->SetPalette(kCividis);
  gSystem->mkdir(outDir,kTRUE);

  //---- input ------------------------------------------------------------------
  TFile *file_train = TFile::Open(trainFile);
  TFile *file_test  = TFile::Open(cfg.testFile);
  if(!file_train || file_train->IsZombie() || !file_test || file_test->IsZombie()){
    std::cout << "unfoldClosureTest: could not open\n  " << trainFile
              << "\n  " << cfg.testFile << "\n";
    return;
  }

  TH2D *h_response_train, *h_response_test;
  h_response_train = loadResponse(file_train,cfg.responseName,removeXJets,"h_resp_train");
  h_response_test  = loadResponse(file_test ,cfg.responseName,removeXJets,"h_resp_test");

  TH1D *h_unmReco_train = 0, *h_unmGen_train = 0;
  TH1D *h_unmReco_test  = 0, *h_unmGen_test  = 0;
  if(includeUnmatched){
    file_train->GetObject(cfg.unmRecoName,h_unmReco_train);
    file_train->GetObject(cfg.unmGenName ,h_unmGen_train);
    file_test ->GetObject(cfg.unmRecoName,h_unmReco_test);
    file_test ->GetObject(cfg.unmGenName ,h_unmGen_test);
    if(!h_unmReco_train || !h_unmGen_train || !h_unmReco_test || !h_unmGen_test){
      std::cout << "unfoldClosureTest: could not find the unmatched histograms ("
                << cfg.unmRecoName << ", " << cfg.unmGenName << ")\n";
      return;
    }
  }
  if(!h_response_train || !h_response_test){
    std::cout << "unfoldClosureTest: could not find \"" << cfg.responseName << "\"\n";
    return;
  }

  //  Rescale every gen row of the unweighted matrix by (weighted truth) /
  //  (unweighted truth).  P(reco|gen) is untouched -- only the gen marginal,
  //  i.e. the Bayesian prior, is swapped for the physical one.
  //  (RooUnfoldBayes::SetPriors would be the obvious route, but in this build
  //  it leaves _N0C at zero, which blows up the error propagation.)
  if(reweightToTruthW){

    TFile *file_trainW = TFile::Open(cfg.trainFile);
    TH2D  *h_response_trainW = 0;
    if(file_trainW && !file_trainW->IsZombie())
      h_response_trainW = loadResponse(file_trainW,cfg.responseName,removeXJets,"h_resp_trainW");
    if(!h_response_trainW){
      std::cout << "unfoldClosureTest: could not read the pThat-weighted response from\n  "
                << cfg.trainFile << "\n";
      return;
    }

    TH1D *h_truth_unw = projY(h_response_train ,"h_truth_unw");
    TH1D *h_truth_w   = projY(h_response_trainW,"h_truth_w");

    TH2D *h_hybrid = (TH2D*) h_response_train->Clone("h_response_hybrid");
    h_hybrid->SetDirectory(0);

    for(int iy = 1; iy <= h_hybrid->GetNbinsY(); iy++){
      double u = h_truth_unw->GetBinContent(iy);
      double s = (u > 0.) ? h_truth_w->GetBinContent(iy) / u : 0.;
      for(int ix = 1; ix <= h_hybrid->GetNbinsX(); ix++){
        h_hybrid->SetBinContent(ix,iy, s * h_hybrid->GetBinContent(ix,iy));
        h_hybrid->SetBinError  (ix,iy, s * h_hybrid->GetBinError  (ix,iy));
      }
    }

    h_response_train = h_hybrid;

  }

  // training half -> response matrix and prior
  TH1D *h_meas_train  = projX(h_response_train,"h_meas_train");
  TH1D *h_truth_train = projY(h_response_train,"h_truth_train");

  // testing half -> what gets unfolded, and what we compare against
  TH1D *h_meas_test  = projX(h_response_test,"h_meas_test");
  TH1D *h_truth_test = projY(h_response_test,"h_truth_test");

  //  unmatched reco -> fakes, unmatched gen -> inefficiency.  Added to the
  //  measured and truth spectra; RooUnfold derives both from the difference
  //  against the projections of the matched 2D.
  if(includeUnmatched){
    h_meas_train ->Add(h_unmReco_train);
    h_truth_train->Add(h_unmGen_train);
    h_meas_test  ->Add(h_unmReco_test);
    h_truth_test ->Add(h_unmGen_test);
  }

  //  The response is always built from the undistorted training half: the prior
  //  it carries is exactly what the distorted test is supposed to be wrong about.
  RooUnfoldResponse response(h_meas_train,h_truth_train,h_response_train,
                             "response",cfg.title);

  //---- distorted truth --------------------------------------------------------
  //  Scale gen row t of the training response by w(t).  That leaves P(reco|gen)
  //  untouched -- w cancels in the row normalisation -- so the only thing the
  //  unfolding has to recover from is the wrong prior.  Folding through the same
  //  matrix that will be used to unfold makes the pseudo-measurement noiseless,
  //  which keeps bias^2 free of test-sample statistical scatter.
  if(doDistort){

    double expo = tiltExpo;

    //  relative statistical precision to hand the pseudo-measurement.  Empty
    //  means "keep the MC errors"; for dataMC it is filled from the data, since
    //  otherwise the variance term reflects MC statistics and the bias/variance
    //  trade-off -- and therefore the MSE minimum -- is meaningless.
    std::vector<double> relErr;

    if(doDataMC){
      TH1D *h_data = loadDataSpectrum(cfg);
      if(!h_data) return;
      double fitLow = std::max(ptLow, cfg.dataMinPt);
      if(fitLow >= ptHigh){
        std::cout << "unfoldClosureTest: no usable data range for the fit\n";
        return;
      }
      if(!buildDataMCDistortion(h_data,h_meas_train,fitLow,ptHigh,pTref,
                                cfg.dataMinPt,expo,relErr,cfg.title,
                                Form("%sdistortion_%s.pdf",outDir.Data(),cfg.label.Data())))
        return;
    }

    TH2D *h_dist = (TH2D*) h_response_train->Clone("h_response_distorted");
    h_dist->SetDirectory(0);
    for(int iy = 1; iy <= h_dist->GetNbinsY(); iy++){
      double w = distortionWeight(h_dist->GetYaxis()->GetBinCenter(iy),expo,pTref);
      for(int ix = 1; ix <= h_dist->GetNbinsX(); ix++){
        h_dist->SetBinContent(ix,iy, w * h_dist->GetBinContent(ix,iy));
        h_dist->SetBinError  (ix,iy, w * h_dist->GetBinError  (ix,iy));
      }
    }

    h_meas_test  = projX(h_dist,"h_meas_distorted");
    h_truth_test = projY(h_dist,"h_truth_distorted");

    //  The pseudo-measurement is folded from the training half, so its
    //  unmatched jets come from there too.  Misses are a gen-level quantity and
    //  follow the reweighting; fakes are underlying-event background and do not
    //  depend on the signal truth shape, so they are left unweighted.
    if(includeUnmatched){
      TH1D *h_unmGen_w = (TH1D*) h_unmGen_train->Clone("h_unmGen_distorted");
      h_unmGen_w->SetDirectory(0);
      for(int b = 1; b <= h_unmGen_w->GetNbinsX(); b++){
        double w = distortionWeight(h_unmGen_w->GetBinCenter(b),expo,pTref);
        h_unmGen_w->SetBinContent(b, w * h_unmGen_w->GetBinContent(b));
        h_unmGen_w->SetBinError  (b, w * h_unmGen_w->GetBinError(b));
      }
      h_meas_test ->Add(h_unmReco_train);
      h_truth_test->Add(h_unmGen_w);
    }

    if(!relErr.empty())
      for(int b = 1; b <= h_meas_test->GetNbinsX() && b < (int)relErr.size(); b++)
        h_meas_test->SetBinError(b, relErr[b] * h_meas_test->GetBinContent(b));

    distTitle = Form("distorted truth: #it{w} #propto #it{p}_{T}^{%.2f}",expo);
    outTag   += doDataMC ? "_distDataMC" : Form("_distTilt%.2f",expo);

  }

  //---- scan N = 0 ... N_iter_max ----------------------------------------------
  //  N = 0 is the measured spectrum itself: the Bayesian iteration has not been
  //  applied at all, so it is the reference the unfolding has to improve on.
  const int N_points = N_iter_max + 1;

  std::vector<TH1D*>  h_unfold(N_points);
  std::vector<double> x(N_points), bias2(N_points), var(N_points),
                      mse(N_points), chi2(N_points), chi2ndf(N_points);
  std::vector<int>    ndf(N_points);

  RooUnfoldBayes unfold(&response,h_meas_test,1);
  unfold.SetVerbose(0);
  if(includeUnmatched) unfold.HandleFakes(true);   // off by default in RooUnfold

  for(int n = 0; n <= N_iter_max; n++){

    TH1D *h = 0;

    if(n == 0){
      h = (TH1D*) h_meas_test->Clone(Form("h_unfold_%s_iter0",cfg.label.Data()));
      h->SetLineStyle(2);
    }
    else{
      unfold.SetIterations(n);
      h = (TH1D*) unfold.Hunfold();
      h->SetName(Form("h_unfold_%s_iter%i",cfg.label.Data(),n));

      // errors from the diagonal of the full unfolding covariance matrix
      TMatrixD cov = unfold.Eunfold();
      for(int b = 1; b <= h->GetNbinsX(); b++){
        double v = cov(b-1,b-1);
        h->SetBinError(b, v > 0. ? std::sqrt(v) : 0.);
      }
    }

    h->SetDirectory(0);
    h->SetLineColor(iterColor(n,N_iter_max));
    h->SetLineWidth(2);
    h->SetStats(0);
    h_unfold[n] = h;

    ClosureMetrics m = computeMetrics(h,h_truth_test,ptLow,ptHigh);
    x[n]       = n;
    bias2[n]   = m.bias2;
    var[n]     = m.variance;
    mse[n]     = m.mse;
    chi2[n]    = m.chi2;
    ndf[n]     = m.ndf;
    chi2ndf[n] = (m.ndf > 0) ? m.chi2 / m.ndf : 0.;

  }

  //---- plots ------------------------------------------------------------------
  drawClosure(N_iter_max,h_unfold.data(),h_truth_test,
              drawPtLow,drawPtHigh,ptLow,ptHigh,ratioLow,ratioHigh,
              cfg.title,respTitle,distTitle,
              Form("%sclosure_%s.pdf",outDir.Data(),outTag.Data()));

  drawSummary(N_iter_max,x.data(),bias2.data(),var.data(),mse.data(),chi2ndf.data(),
              cfg.title,respTitle,distTitle,doDistort,
              Form("%sclosure_summary_%s.pdf",outDir.Data(),outTag.Data()));

  //---- table ------------------------------------------------------------------
  int iBest = 1;
  for(int n = 2; n < N_points; n++) if(mse[n] < mse[iBest]) iBest = n;

  printf("\n=== %s [%s] : %s closure, %.0f < gen pT < %.0f GeV ===\n",
         cfg.title.Data(), outTag.Data(),
         doDistort ? "distorted-truth" : "split-sample", ptLow, ptHigh);
  printf("  train: %s\n  test:  %s\n\n",trainFile.Data(),cfg.testFile.Data());
  printf("%-6s  %12s  %12s  %12s  %12s  %5s  %10s\n",
         "iter","bias^2","variance","MSE","chi2","ndf","chi2/ndf");
  for(int n = 0; n < N_points; n++){
    printf("%-6.0f  %12.4g  %12.4g  %12.4g  %12.4g  %5d  %10.3f%s\n",
           x[n],bias2[n],var[n],mse[n],chi2[n],ndf[n],chi2ndf[n],
           n == 0     ? "   (measured, no unfolding)" :
           (n == iBest && doDistort) ? "   <- min MSE" : "");
  }
  printf("\n");

}
