// Check for a jet-selection-induced bias in the hiBin (centrality) distribution
// of the PbPb MinBias sample used by calculateRAA.C -- follow-up on Olga's
// question in the 2026-08 meeting about whether the fine centrality mix
// within the 10-30% class (C2 = hiBin [20,60)) is skewed toward its low-hiBin
// (i.e. more-central, 10% side), which would explain under-subtraction of
// fake jets there.
//
// h_hiBin        : inclusive events, no jet requirement
// h_hiBin_jet     : events that have >=1 inclusive reco jet
// Both are shape-normalized; the ratio isolates any bias the jet requirement
// introduces on top of the inherent hiBin distribution.
//
// Also prints the normalized yield in each 5%-wide centrality slice, so the
// four slices making up the 10-30% class (10-15, 15-20, 20-25, 25-30%) can be
// read off directly instead of eyeballing the plot.
//
// Usage: root -l -b -q 'plotHiBin_jetSelectionBias.C'
// Run from: src/plots/hiBin/

const char *inFile =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/PbPb_MinBias_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root";

const char *outDir  = "../../../figures/hiBin/";
const char *outName = "hiBin_jetSelectionBias.pdf";

// C1..C4 class boundaries in hiBin (centrality_4CentBins.h): 0-10, 10-30, 30-50, 50-80%
const int NClassEdge = 5;
const double classEdgeHiBin[NClassEdge] = {0, 20, 60, 100, 160};
const char  *classLabel[NClassEdge-1]   = {"0-10%", "10-30%", "30-50%", "50-80%"};

// Okabe-Ito, centre-symmetric markers only.
const char *serHex[2]  = {"#0072B2", "#D55E00"};
const int   serMark[2] = {20, 25};

void plotHiBin_jetSelectionBias()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TFile *f = TFile::Open(inFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", inFile); return; }

  TH1D *hIncl = nullptr, *hJet = nullptr;
  f->GetObject("h_hiBin", hIncl);
  f->GetObject("h_hiBin_jet", hJet);
  if(!hIncl || !hJet){ printf("ERROR: h_hiBin / h_hiBin_jet not found\n"); return; }
  hIncl->SetDirectory(nullptr);
  hJet->SetDirectory(nullptr);

  hIncl->Scale(1./hIncl->Integral());
  hJet->Scale(1./hJet->Integral());

  // 5%-wide slices (10 hiBin) over 0-80%, so the four 10-30% sub-slices are visible.
  const int NSlice = 16;
  const double hiBinPerSlice = 10.;
  printf("\n%-9s %14s %14s %10s\n", "cent [%]", "incl (norm)", "w/ jet (norm)", "jet/incl");
  for(int is = 0; is < NSlice; is++){
    double lo = is*hiBinPerSlice, hi = (is+1)*hiBinPerSlice;
    int b0 = hIncl->FindBin(lo+1e-6), b1 = hIncl->FindBin(hi-1e-6);
    double vIncl = hIncl->Integral(b0, b1);
    double vJet  = hJet->Integral(b0, b1);
    double ratio = (vIncl > 0) ? vJet/vIncl : 0.;
    printf("%3.0f-%-5.0f %14.4f %14.4f %10.3f\n", lo/2., hi/2., vIncl, vJet, ratio);
  }
  printf("\n");

  TH1D *hRatio = (TH1D*) hJet->Clone("hRatio");
  hRatio->Divide(hIncl);

  const double lm = 0.15, rm = 0.05, tm = 0.07, bm = 0.28, split = 0.35;
  TCanvas *c = new TCanvas("cHiBinBias", "", 700, 800);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->Draw(); pUp->cd();

  hIncl->SetLineColor(TColor::GetColor(serHex[0]));
  hIncl->SetMarkerColor(TColor::GetColor(serHex[0]));
  hIncl->SetMarkerStyle(serMark[0]); hIncl->SetLineWidth(2);
  hJet->SetLineColor(TColor::GetColor(serHex[1]));
  hJet->SetMarkerColor(TColor::GetColor(serHex[1]));
  hJet->SetMarkerStyle(serMark[1]); hJet->SetLineWidth(2);

  hIncl->SetTitle("");
  hIncl->GetXaxis()->SetRangeUser(0, 160);
  hIncl->GetYaxis()->SetTitle("1/#it{N} d#it{N}/d(hiBin), shape-normalized");
  hIncl->GetYaxis()->SetTitleSize(0.055);
  hIncl->GetYaxis()->SetLabelSize(0.050);
  hIncl->GetYaxis()->SetTitleOffset(1.20);
  hIncl->SetMinimum(0.);
  hIncl->Draw("hist");
  hJet->Draw("hist same");

  TLegend *leg = new TLegend(0.44, 0.70, 0.94, 0.88);
  leg->SetBorderSize(0); leg->SetTextSize(0.045); leg->SetFillStyle(0);
  leg->AddEntry(hIncl, "inclusive events", "l");
  leg->AddEntry(hJet, "events w/ incl. reco jet", "l");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.048);
  lat.DrawLatex(0.18, 0.90, "PbPb MinBias");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  hRatio->SetLineColor(kBlack); hRatio->SetMarkerColor(kBlack);
  hRatio->SetMarkerStyle(20); hRatio->SetLineWidth(2);
  hRatio->SetTitle("");
  hRatio->GetXaxis()->SetRangeUser(0, 160);
  hRatio->GetXaxis()->SetTitle("hiBin");
  hRatio->GetYaxis()->SetTitle("w/ jet / incl.");
  hRatio->GetXaxis()->SetTitleSize(0.050*sc);
  hRatio->GetXaxis()->SetLabelSize(0.045*sc);
  hRatio->GetYaxis()->SetTitleSize(0.042*sc);
  hRatio->GetYaxis()->SetLabelSize(0.045*sc);
  hRatio->GetYaxis()->SetTitleOffset(1.55/sc);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->SetMinimum(0.7); hRatio->SetMaximum(1.3);
  hRatio->Draw("hist");

  TLine *one = new TLine(0., 1., 160., 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  for(int ie = 1; ie < NClassEdge-1; ie++){
    TLine *cl = new TLine(classEdgeHiBin[ie], 0.7, classEdgeHiBin[ie], 1.3);
    cl->SetLineStyle(3); cl->SetLineColor(kGray+2); cl->Draw();
  }

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("Saved %s\n", out.Data());
}
