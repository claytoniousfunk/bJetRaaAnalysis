// hiBin distribution of *filled fake jets*, compared to inclusive events and
// events with a real reco jet -- follow-up on Olga's 2026-08 meeting comment
// that the fake-jet subtraction for the 10-30% class may be under-subtracting
// because the fine-hiBin mix feeding it isn't flat.
//
// Uses the actual mixed-event file that makeFakeJetFile.C reads to build
// fakeJets.root, so this is the real population behind the current
// subtraction, not a stand-in sample.
//
// h_hiBin              : inclusive events (shape used to weight makeFakeJetFile.C's
//                         per-class average, via N_events = sum of h_vz_C{si})
// h_hiBin_jet           : events with >=1 real inclusive reco jet
// h_hiBin_fakeJetFilled : hiBin, weighted by the *number of filled fake-jet
//                         entries* -- built here by summing Integral() of
//                         h_fastJetPt_PF_JEC_bkgSub_RC_C{1..16} (one entry
//                         per FastJet pseudo-jet surviving the RC-subtracted
//                         20 GeV threshold; can be 0, 1, or several per event)
//                         into the same 16 ultra-fine hiBin slices.
//
// Usage: root -l -b -q 'plotHiBin_fakeJetFilledBias.C'
// Run from: src/plots/hiBin/

const char *inFile =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PbPb/PbPb_MinBias_Part1_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_mixedEventPFClustering_pseudoJetCandPtMin-0.0_2026-8-11_ultraFineCentBins.root";

const char *outDir  = "../../../figures/hiBin/";
const char *outName = "hiBin_fakeJetFilledBias.pdf";

const int NSlice = 16;               // C1..C16, 10 hiBin wide each
const double hiBinPerSlice = 10.;

// C1..C4 class boundaries in hiBin (centrality_4CentBins.h): 0-10, 10-30, 30-50, 50-80%
const int NClassEdge = 5;
const double classEdgeHiBin[NClassEdge] = {0, 20, 60, 100, 160};

// Okabe-Ito, centre-symmetric markers only.
const char *serHex[3]  = {"#0072B2", "#D55E00", "#009E73"};
const int   serMark[3] = {20, 25, 21};

void plotHiBin_fakeJetFilledBias()
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

  // Build the fake-jet-filled hiBin distribution from the per-slice FastJet spectra.
  TH1D *hFake = new TH1D("hFake", "", NSlice, 0., NSlice*hiBinPerSlice);
  hFake->SetDirectory(nullptr);
  for(int si = 1; si <= NSlice; si++){
    TH1D *hFJ = nullptr;
    f->GetObject(Form("h_fastJetPt_PF_JEC_bkgSub_RC_C%d", si), hFJ);
    if(!hFJ){ printf("WARNING: missing h_fastJetPt_PF_JEC_bkgSub_RC_C%d\n", si); continue; }
    hFake->SetBinContent(si, hFJ->Integral());
  }

  hIncl->Scale(1./hIncl->Integral());
  hJet->Scale(1./hJet->Integral());
  hFake->Scale(1./hFake->Integral());

  // Rebin incl/jet onto the same 16 native slices as hFake, so the ratio
  // panel below can Divide() them directly.
  TH1D *hInclSlice = new TH1D("hInclSlice", "", NSlice, 0., NSlice*hiBinPerSlice);
  TH1D *hJetSlice  = new TH1D("hJetSlice",  "", NSlice, 0., NSlice*hiBinPerSlice);
  hInclSlice->SetDirectory(nullptr); hJetSlice->SetDirectory(nullptr);

  printf("\n%-9s %14s %14s %14s %10s\n", "cent [%]", "incl (norm)", "w/ jet (norm)", "fakeJet (norm)", "fake/incl");
  for(int is = 1; is <= NSlice; is++){
    double lo = (is-1)*hiBinPerSlice, hi = is*hiBinPerSlice;
    double vIncl = hIncl->GetBinContent(hIncl->FindBin(lo+1e-6));
    double vJet  = hJet->GetBinContent(hJet->FindBin(lo+1e-6));
    double vFake = hFake->GetBinContent(is);
    hInclSlice->SetBinContent(is, vIncl);
    hJetSlice->SetBinContent(is, vJet);
    double ratio = (vIncl > 0) ? vFake/vIncl : 0.;
    printf("%3.0f-%-5.0f %14.4f %14.4f %14.4f %10.3f\n", lo/2., hi/2., vIncl, vJet, vFake, ratio);
  }
  printf("\n");

  const double lm = 0.15, rm = 0.05, tm = 0.07, bm = 0.28, split = 0.35;
  TCanvas *c = new TCanvas("cHiBinFake", "", 700, 800);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->SetLogy();
  pUp->Draw(); pUp->cd();

  hIncl->SetLineColor(TColor::GetColor(serHex[0]));
  hIncl->SetLineWidth(2);
  hJet->SetLineColor(TColor::GetColor(serHex[1]));
  hJet->SetLineWidth(2);
  hFake->SetLineColor(TColor::GetColor(serHex[2]));
  hFake->SetLineWidth(2);

  hIncl->SetTitle("");
  hIncl->GetXaxis()->SetRangeUser(0, 160);
  hIncl->GetYaxis()->SetTitle("1/#it{N} d#it{N}/d(hiBin), shape-normalized");
  hIncl->GetYaxis()->SetTitleSize(0.055);
  hIncl->GetYaxis()->SetLabelSize(0.050);
  hIncl->GetYaxis()->SetTitleOffset(1.20);
  hIncl->SetMinimum(1e-4);
  hIncl->Draw("hist");
  hJet->Draw("hist same");
  hFake->Draw("hist same");

  TLegend *leg = new TLegend(0.44, 0.66, 0.94, 0.88);
  leg->SetBorderSize(0); leg->SetTextSize(0.042); leg->SetFillStyle(0);
  leg->AddEntry(hIncl, "inclusive events", "l");
  leg->AddEntry(hJet, "events w/ real reco jet", "l");
  leg->AddEntry(hFake, "filled fake jets (FastJet, RC sub > 20 GeV)", "l");
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.048);
  lat.DrawLatex(0.18, 0.90, "PbPb MinBias, mixed-event PF clustering");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  TH1D *hRatioFakeIncl = (TH1D*) hFake->Clone("hRatioFakeIncl");
  hRatioFakeIncl->Divide(hInclSlice);
  TH1D *hRatioFakeJet = (TH1D*) hFake->Clone("hRatioFakeJet");
  hRatioFakeJet->Divide(hJetSlice);

  hRatioFakeIncl->SetLineColor(TColor::GetColor(serHex[0])); hRatioFakeIncl->SetLineWidth(2);
  hRatioFakeJet->SetLineColor(TColor::GetColor(serHex[1])); hRatioFakeJet->SetLineWidth(2);

  hRatioFakeIncl->SetTitle("");
  hRatioFakeIncl->GetXaxis()->SetRangeUser(0, 160);
  hRatioFakeIncl->GetXaxis()->SetTitle("hiBin");
  hRatioFakeIncl->GetYaxis()->SetTitle("fakeJet / X");
  hRatioFakeIncl->GetXaxis()->SetTitleSize(0.050*sc);
  hRatioFakeIncl->GetXaxis()->SetLabelSize(0.045*sc);
  hRatioFakeIncl->GetYaxis()->SetTitleSize(0.042*sc);
  hRatioFakeIncl->GetYaxis()->SetLabelSize(0.045*sc);
  hRatioFakeIncl->GetYaxis()->SetTitleOffset(1.55/sc);
  hRatioFakeIncl->GetYaxis()->SetNdivisions(505);
  hRatioFakeIncl->SetMinimum(0.); hRatioFakeIncl->SetMaximum(6.);
  hRatioFakeIncl->Draw("hist");
  hRatioFakeJet->Draw("hist same");

  TLegend *leg2 = new TLegend(0.18, 0.75, 0.60, 0.93);
  leg2->SetBorderSize(0); leg2->SetTextSize(0.042*sc); leg2->SetFillStyle(0);
  leg2->AddEntry(hRatioFakeIncl, "fakeJet / inclusive", "l");
  leg2->AddEntry(hRatioFakeJet, "fakeJet / real jet", "l");
  leg2->Draw();

  TLine *one = new TLine(0., 1., 160., 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  for(int ie = 1; ie < NClassEdge-1; ie++){
    TLine *cl = new TLine(classEdgeHiBin[ie], 0., classEdgeHiBin[ie], 6.);
    cl->SetLineStyle(3); cl->SetLineColor(kGray+2); cl->Draw();
  }

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("Saved %s\n", out.Data());
}
