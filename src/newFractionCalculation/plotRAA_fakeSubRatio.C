// Ratio of RAA with fake-jet subtraction (pT > 20 GeV) to RAA without.
// Reads r_CX_r from the two saved output files and plots (with/without) per centrality.

const char *f_with_path =
  "rootFiles/JetsPerZ/histograms_withFakeSub.root";
const char *f_without_path =
  "rootFiles/JetsPerZ/histograms_noFakeSub.root";

TString outDir = "../../figures/JetsPerZ/";

const char *centLabel[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};
const int   centColor[5] = {0, kRed+1, kOrange+7, kGreen+2, kBlue+1};

void plotRAA_fakeSubRatio(){

  TFile *fWith    = TFile::Open(f_with_path);
  TFile *fWithout = TFile::Open(f_without_path);
  if(!fWith || fWith->IsZombie() || !fWithout || fWithout->IsZombie()){
    printf("ERROR: cannot open input files\n"); return;
  }

  double lm = 0.15, rm = 0.05, tm = 0.08, bm = 0.12;

  TCanvas *c = new TCanvas("c_ratio","",800,600);
  c->SetLeftMargin(lm); c->SetRightMargin(rm);
  c->SetTopMargin(tm);  c->SetBottomMargin(bm);

  TLine *unity = new TLine();
  unity->SetLineStyle(2); unity->SetLineColor(kGray+1); unity->SetLineWidth(2);

  TLegend *leg = new TLegend(0.55, 0.65, 0.92, 0.88);
  leg->SetBorderSize(0); leg->SetTextSize(0.040);

  TH1D *hFrame = nullptr;

  for(int ic = 1; ic <= 4; ic++){
    TH1D *hW, *hWo;
    fWith   ->GetObject(Form("r_C%d_r", ic), hW);
    fWithout->GetObject(Form("r_C%d_r", ic), hWo);
    if(!hW || !hWo){ printf("WARNING: missing C%d\n", ic); continue; }

    TH1D *hRatio = (TH1D*) hW->Clone(Form("ratio_C%d", ic));
    hRatio->Divide(hWo);

    hRatio->SetLineColor(centColor[ic]);
    hRatio->SetMarkerColor(centColor[ic]);
    hRatio->SetLineWidth(2);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.9);
    hRatio->SetStats(0);
    hRatio->SetTitle("");

    if(!hFrame){
      hFrame = hRatio;
      hFrame->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
      hFrame->GetYaxis()->SetTitle("R_{AA} (with sub) / R_{AA} (no sub)");
      hFrame->GetYaxis()->SetRangeUser(0.90, 1.10);
      hFrame->GetXaxis()->SetTitleSize(0.050);
      hFrame->GetXaxis()->SetLabelSize(0.045);
      hFrame->GetYaxis()->SetTitleSize(0.050);
      hFrame->GetYaxis()->SetLabelSize(0.045);
      hFrame->GetYaxis()->SetTitleOffset(1.2);
      hFrame->GetYaxis()->SetNdivisions(505);
      hFrame->Draw("ep");
    } else {
      hRatio->Draw("ep same");
    }
    leg->AddEntry(hRatio, centLabel[ic], "lp");
  }

  // draw unity line over the full x range
  double xlo = hFrame->GetXaxis()->GetXmin();
  double xhi = hFrame->GetXaxis()->GetXmax();
  unity->DrawLine(xlo, 1., xhi, 1.);

  leg->Draw();

  TLatex *lat = new TLatex();
  lat->SetNDC(); lat->SetTextSize(0.042);
  lat->DrawLatex(lm + 0.02, 1. - tm - 0.06,
    "PbPb 5.02 TeV  |#eta^{jet}| < 2.0");
  lat->DrawLatex(lm + 0.02, 1. - tm - 0.12,
    "Fake-jet sub: mixed-event pseudo-jets, #it{p}_{T}^{fake} > 20 GeV");

  TString outFile = outDir + "RAA_fakeSubRatio.pdf";
  c->SaveAs(outFile);
  printf("Saved %s\n", outFile.Data());
}
