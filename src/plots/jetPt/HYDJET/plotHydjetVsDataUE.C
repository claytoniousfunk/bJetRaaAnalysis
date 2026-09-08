// Mean random-cone UE estimate vs centrality, HYDJET against PbPb data, both
// from mixed-event running, with the HYDJET/data ratio underneath.
//
// The UE maps are TProfile2D in (eta, phi); this collapses each to its
// entry-weighted mean over the acceptance, the same reduction used by
// plotRandConeMeanVsCent.C.
//
// The two samples are independent, so the ratio uncertainties are propagated
// uncorrelated.
//
// Usage: root -l -b -q 'plotHydjetVsDataUE.C'
// Run from: src/plots/jetPt/HYDJET/

const char *dataFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";
const char *hydjetFile = "/home/clayton/Downloads/HYDJET_scan_output_1.root";

const char *outDir  = "../../../../figures/jetPt/HYDJET/";
const char *outName = "hydjetVsData_UE.pdf";

const int    NSlice   = 16;
const double centStep = 5.0;

// Okabe-Ito, centre-symmetric markers only.
const char *serHex [2] = { "#0072B2", "#D55E00" };
const int   serMark[2] = { 20, 25 };
const double serSize[2]= { 1.0, 1.5 };

const double ratioMin = 0.9, ratioMax = 1.5;

static bool profileMean(TProfile2D *p, double &mean, double &err)
{
  double sw = 0., swy = 0., var = 0.;
  for(int ix = 1; ix <= p->GetNbinsX(); ix++){
    for(int iy = 1; iy <= p->GetNbinsY(); iy++){
      int b = p->GetBin(ix, iy);
      double n = p->GetBinEntries(b);
      if(n <= 0) continue;
      sw  += n;
      swy += n * p->GetBinContent(b);
      var += n * n * p->GetBinError(b) * p->GetBinError(b);
    }
  }
  if(sw <= 0) return false;
  mean = swy / sw;
  err  = TMath::Sqrt(var) / sw;
  return true;
}

void plotHydjetVsDataUE()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *fD = TFile::Open(dataFile);
  TFile *fH = TFile::Open(hydjetFile);
  if(!fD || fD->IsZombie()){ printf("ERROR: cannot open %s\n", dataFile);   return; }
  if(!fH || fH->IsZombie()){ printf("ERROR: cannot open %s\n", hydjetFile); return; }

  TH1D *hD = new TH1D("hUE_data",  "", NSlice, 0., NSlice*centStep);
  TH1D *hH = new TH1D("hUE_hydjet","", NSlice, 0., NSlice*centStep);
  TH1D *hR = new TH1D("hUE_ratio", "", NSlice, 0., NSlice*centStep);
  hD->SetDirectory(nullptr); hH->SetDirectory(nullptr); hR->SetDirectory(nullptr);

  printf("\n%-9s %10s %10s %9s\n", "cent [%]", "HYDJET", "data", "H/data");
  for(int si = 1; si <= NSlice; si++){
    TProfile2D *pD = nullptr, *pH = nullptr;
    fD->GetObject(Form("h_randConeEtaPhi_C%d", si), pD);
    fH->GetObject(Form("h_randConeEtaPhi_C%d", si), pH);
    double mD=0, eD=0, mH=0, eH=0;
    bool okD = pD && profileMean(pD, mD, eD);
    bool okH = pH && profileMean(pH, mH, eH);
    if(okD){ hD->SetBinContent(si, mD); hD->SetBinError(si, eD); }
    if(okH){ hH->SetBinContent(si, mH); hH->SetBinError(si, eH); }
    if(okD && okH && mD > 0. && mH > 0.){
      double r = mH/mD;
      hR->SetBinContent(si, r);
      hR->SetBinError(si, r*TMath::Sqrt((eH/mH)*(eH/mH) + (eD/mD)*(eD/mD)));
      printf("%3.0f-%-5.0f %10.2f %10.2f %9.3f\n",
             (si-1)*centStep, si*centStep, mH, mD, r);
    }
  }

  double ymax = TMath::Max(hD->GetMaximum(), hH->GetMaximum()) * 1.12;

  const double lm = 0.15, rm = 0.05, tm = 0.07, bm = 0.28, split = 0.35;
  TCanvas *c = new TCanvas("cUE", "", 700, 800);

  TPad *pUp = new TPad("pUpUE","",0,split,1,1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->Draw(); pUp->cd();

  hH->SetLineColor(TColor::GetColor(serHex[0]));
  hH->SetMarkerColor(TColor::GetColor(serHex[0]));
  hH->SetMarkerStyle(serMark[0]); hH->SetMarkerSize(serSize[0]); hH->SetLineWidth(2);
  hD->SetLineColor(TColor::GetColor(serHex[1]));
  hD->SetMarkerColor(TColor::GetColor(serHex[1]));
  hD->SetMarkerStyle(serMark[1]); hD->SetMarkerSize(serSize[1]); hD->SetLineWidth(2);

  hH->SetTitle("");
  hH->GetYaxis()->SetTitle("mean UE estimate [GeV]");
  hH->GetYaxis()->SetTitleSize(0.055);
  hH->GetYaxis()->SetLabelSize(0.050);
  hH->GetYaxis()->SetTitleOffset(1.20);
  hH->SetMinimum(0.); hH->SetMaximum(ymax);
  hH->Draw("ep");
  hD->Draw("ep same");
  hH->GetYaxis()->ChangeLabel(1, -1, 0.);

  TLegend *leg = new TLegend(0.44, 0.66, 0.94, 0.88);
  leg->SetBorderSize(0); leg->SetTextSize(0.045); leg->SetFillStyle(0);
  leg->AddEntry(hH, "HYDJET, mixed evt", "lp");
  leg->AddEntry(hD, "PbPb data, mixed evt", "lp");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.048);
  lat.DrawLatex(0.18, 0.87, "Random cone, R = 0.4");

  c->cd();
  TPad *pDn = new TPad("pDnUE","",0,0,1,split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  hR->SetLineColor(kBlack); hR->SetMarkerColor(kBlack);
  hR->SetMarkerStyle(20); hR->SetMarkerSize(1.0); hR->SetLineWidth(2);
  hR->SetTitle("");
  hR->GetXaxis()->SetTitle("centrality [%]");
  hR->GetYaxis()->SetTitle("HYDJET / data");
  hR->GetXaxis()->SetTitleSize(0.050*sc);
  hR->GetXaxis()->SetLabelSize(0.045*sc);
  hR->GetYaxis()->SetTitleSize(0.042*sc);
  hR->GetYaxis()->SetLabelSize(0.045*sc);
  hR->GetYaxis()->SetTitleOffset(1.55/sc);
  hR->GetYaxis()->SetNdivisions(505);
  hR->SetMinimum(ratioMin); hR->SetMaximum(ratioMax);
  hR->Draw("ep");

  TLine *one = new TLine(0., 1., NSlice*centStep, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
