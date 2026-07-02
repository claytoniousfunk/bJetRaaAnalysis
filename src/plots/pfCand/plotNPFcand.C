// N PF candidates per event by centrality — PF (raw) and PFCS (constituent subtracted).
// Single canvas, two pads side by side, all centrality bins overlaid, normalised to unit area.
// Usage: root -l -q 'plotNPFcand.C'

const char *inFile =
  "/home/clayton/Downloads/PbPb_pfCandAnalyzer_output_1.root";

const char *outDir = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/";

const int NCentBins = 5;
const char *centLabel[NCentBins] = {
  "0-80% (incl.)", "0-10%", "10-30%", "30-50%", "50-80%"
};

void plotNPFcand()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TFile *f = TFile::Open(inFile);
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", inFile); return; }

  gSystem->Exec(Form("mkdir -p %s", outDir));

  // Okabe-Ito colorblind-safe palette
  const int cols[NCentBins] = {
    kBlack,
    TColor::GetColor("#E69F00"),
    TColor::GetColor("#56B4E9"),
    TColor::GetColor("#009E73"),
    TColor::GetColor("#0072B2")
  };

  TH1D *hPF[NCentBins], *hCS[NCentBins];
  for(int ci = 0; ci < NCentBins; ci++){
    f->GetObject(Form("h_nPFcand_C%d",   ci), hPF[ci]);
    f->GetObject(Form("h_nPFcandCS_C%d", ci), hCS[ci]);
    if(!hPF[ci] || !hCS[ci]){
      printf("ERROR: missing histograms for C%d\n", ci); return;
    }
    hPF[ci] = (TH1D*) hPF[ci]->Clone(Form("hPF_C%d", ci));
    hCS[ci] = (TH1D*) hCS[ci]->Clone(Form("hCS_C%d", ci));
    hPF[ci]->SetDirectory(0);
    hCS[ci]->SetDirectory(0);

    hPF[ci]->Rebin(50);
    hCS[ci]->Rebin(10);

    if(hPF[ci]->Integral() > 0) hPF[ci]->Scale(1. / hPF[ci]->Integral());
    if(hCS[ci]->Integral() > 0) hCS[ci]->Scale(1. / hCS[ci]->Integral());

    for(TH1D *h : {hPF[ci], hCS[ci]}){
      h->SetLineColor(cols[ci]);
      h->SetLineWidth(2);
      h->SetStats(0);
      h->SetTitle("");
      h->GetYaxis()->SetTitle("Fraction of events / bin");
      h->GetYaxis()->SetTitleSize(0.055);
      h->GetYaxis()->SetLabelSize(0.048);
      h->GetYaxis()->SetTitleOffset(1.4);
      h->GetXaxis()->SetTitleSize(0.055);
      h->GetXaxis()->SetLabelSize(0.048);
    }
  }

  TCanvas *c = new TCanvas("c_nPFcand", "", 1400, 640);
  c->Divide(2, 1);

  // --- left pad: raw PF candidates ---
  TPad *p1 = (TPad*) c->cd(1);
  p1->SetLogy();
  p1->SetLeftMargin(0.16);
  p1->SetRightMargin(0.04);
  p1->SetBottomMargin(0.14);
  p1->SetTopMargin(0.08);

  double ymaxPF = 0.;
  for(int ci = 0; ci < NCentBins; ci++)
    ymaxPF = TMath::Max(ymaxPF, hPF[ci]->GetMaximum());

  for(int ci = 0; ci < NCentBins; ci++){
    hPF[ci]->GetXaxis()->SetRangeUser(0., 8500.);
    hPF[ci]->GetXaxis()->SetTitle("N_{PF cand} / event");
    hPF[ci]->GetYaxis()->SetRangeUser(1e-5, ymaxPF * 1000.);
    hPF[ci]->Draw(ci == 0 ? "hist" : "hist same");
  }

  TLegend *leg1 = new TLegend(0.42, 0.55, 0.93, 0.90);
  leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.046);
  for(int ci = 0; ci < NCentBins; ci++)
    leg1->AddEntry(hPF[ci], centLabel[ci], "l");
  leg1->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC(); lat->SetTextSize(0.046);
  lat->DrawLatex(0.17, 0.935, "PbPb 5.02 TeV  MinBias  PF candidates");

  // --- right pad: PFCS candidates ---
  TPad *p2 = (TPad*) c->cd(2);
  p2->SetLogy();
  p2->SetLeftMargin(0.16);
  p2->SetRightMargin(0.04);
  p2->SetBottomMargin(0.14);
  p2->SetTopMargin(0.08);

  double ymaxCS = 0.;
  for(int ci = 0; ci < NCentBins; ci++)
    ymaxCS = TMath::Max(ymaxCS, hCS[ci]->GetMaximum());

  for(int ci = 0; ci < NCentBins; ci++){
    hCS[ci]->GetXaxis()->SetRangeUser(0., 600.);
    hCS[ci]->GetXaxis()->SetTitle("N_{PFCS cand} / event");
    hCS[ci]->GetYaxis()->SetRangeUser(1e-5, ymaxCS * 5000.);
    hCS[ci]->Draw(ci == 0 ? "hist" : "hist same");
  }

  TLegend *leg2 = new TLegend(0.42, 0.55, 0.93, 0.90);
  leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.046);
  for(int ci = 0; ci < NCentBins; ci++)
    leg2->AddEntry(hCS[ci], centLabel[ci], "l");
  leg2->Draw();

  lat->DrawLatex(0.17, 0.935, "PbPb 5.02 TeV  MinBias  PFCS candidates");

  TString outPath = Form("%snPFcand_vs_centrality.pdf", outDir);
  c->SaveAs(outPath);
  printf("Saved %s\n", outPath.Data());

  const char *gmDir = "/home/clayton/Documents/nuclear/GroupMeeting/figures/2026-07-02/";
  TString gmPath = Form("%snPFcand_vs_centrality.pdf", gmDir);
  c->SaveAs(gmPath);
  printf("Saved %s\n", gmPath.Data());

  delete c;
  for(int ci = 0; ci < NCentBins; ci++){ delete hPF[ci]; delete hCS[ci]; }
  f->Close();
}
