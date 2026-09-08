// Jets-per-Z (PbPb/pp), with stitching (useMinBiasOnly=false -- the real
// calculateRAA.C chain: stitched spectra, UNFOLDED), fine pT bins, comparing
// the plain-RC and RC_geoCorr fake-jet subtraction methods, per coarse
// centrality class. See plotJetsPerZ_RCvsRCgeoCorr_noUnfold.C for the
// otherwise-identical no-unfolding twin.
//
// Top panel: the jets-per-Z ratio itself, both methods overlaid per class
// (open marker = RC sub, filled = RC_geoCorr sub).
// Bottom panel: double ratio (RC_geoCorr sub)/(RC sub) -- how much the
// subtraction METHOD choice moves the final result.
//
// Reads r_C{1..4}_fine (the fine-binned axis, {60,65,70,...,500}, not the
// coarse r_C{1..4}_r) from two calculateRAA.C runs that are identical except
// for which fake-jet estimate was subtracted:
//   histograms_stitched_RC.root        : fakeJetsPath = fakeJets.root          (plain RC)
//   histograms_stitched_RCgeoCorr.root : fakeJetsPath = fakeJets_RCgeoCorr.root (RC_geoCorr)
// Both produced with useMinBiasOnly overridden to false and
// fakeJetsPath/outputRootPath overridden per run:
//   root -l -b -q -e '.L calculateRAA.C' \
//     -e 'useMinBiasOnly=false; fakeJetsPath="../../rootFiles/fakeJets/fakeJets.root"; \
//         outputRootPath="./rootFiles/JetsPerZ/histograms_stitched_RC.root";' \
//     -e 'calculateRAA()'
// (repeat with fakeJets_RCgeoCorr.root / histograms_stitched_RCgeoCorr.root)
//
// Modeled on plotRAA_fakeSubRatio.C's two-run-diff pattern (that file compares
// doFakeJetSubtraction on/off; this one compares which estimate is used when
// it's on).

const char *f_RC_path =
  "rootFiles/JetsPerZ/histograms_stitched_RC.root";
const char *f_RCgeoCorr_path =
  "rootFiles/JetsPerZ/histograms_stitched_RCgeoCorr.root";

TString outDir = "../../figures/JetsPerZ/";

const int   NClass = 4;
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

const double plotPtMin = 60., plotPtMax = 300.;
const double topMin = 0.2, topMax = 1.05;
const double botMin = 0.65, botMax = 1.15;

// Okabe-Ito, centre-symmetric filled/open pairs per class: open = RC sub,
// filled = RC_geoCorr sub -- same convention as the MinBias-only fakeJetEffect
// plots earlier this session.
const char *classHex[NClass]   = { "#D55E00", "#0072B2", "#009E73", "#CC79A7" };
const int   markFilled[NClass] = { 20, 21, 33, 34 };
const int   markOpen[NClass]   = { 24, 25, 27, 28 };

void plotJetsPerZ_RCvsRCgeoCorr_unfold(){

  TFile *fRC  = TFile::Open(f_RC_path);
  TFile *fGeo = TFile::Open(f_RCgeoCorr_path);
  if(!fRC || fRC->IsZombie() || !fGeo || fGeo->IsZombie()){
    printf("ERROR: cannot open input files\n"); return;
  }
  gSystem->mkdir(outDir, kTRUE);
  gStyle->SetOptStat(0);

  TH1D *hRC[NClass], *hGeo[NClass], *hRatio[NClass];
  int colClass[NClass];

  printf("%-8s %14s %14s %14s %10s\n", "class", "pT [GeV]", "RC", "RC_geoCorr", "ratio");
  for(int ic = 0; ic < NClass; ic++){
    colClass[ic] = TColor::GetColor(classHex[ic]);

    fRC ->GetObject(Form("r_C%d_fine", ic+1), hRC[ic]);
    fGeo->GetObject(Form("r_C%d_fine", ic+1), hGeo[ic]);
    if(!hRC[ic] || !hGeo[ic]){ printf("WARNING: missing C%d\n", ic+1); return; }
    hRC[ic]->SetDirectory(nullptr); hGeo[ic]->SetDirectory(nullptr);

    hRatio[ic] = (TH1D*) hGeo[ic]->Clone(Form("ratio_C%d", ic+1));
    hRatio[ic]->SetDirectory(nullptr);
    hRatio[ic]->Divide(hRC[ic]);

    for(int b = 1; b <= hRatio[ic]->GetNbinsX(); b++){
      double lo = hRatio[ic]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      printf("%-8s %5.0f-%-6.0f %14.4f %14.4f %10.4f\n",
             classLabel[ic], lo, hRatio[ic]->GetXaxis()->GetBinUpEdge(b),
             hRC[ic]->GetBinContent(b), hGeo[ic]->GetBinContent(b), hRatio[ic]->GetBinContent(b));
    }
  }

  // Warn rather than silently clip if the fixed display ranges don't cover the
  // actual data in the plotted window.
  for(int ic = 0; ic < NClass; ic++){
    for(int b = 1; b <= hRatio[ic]->GetNbinsX(); b++){
      double lo = hRatio[ic]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      double vRC = hRC[ic]->GetBinContent(b), vGeo = hGeo[ic]->GetBinContent(b);
      if(vRC > topMax || vRC < topMin || vGeo > topMax || vGeo < topMin)
        printf("NOTE: %s jets-per-Z at %.0f-%.0f exceeds top-panel range [%.2f,%.2f] (RC=%.3f, RCgeo=%.3f)\n",
               classLabel[ic], lo, hRatio[ic]->GetXaxis()->GetBinUpEdge(b), topMin, topMax, vRC, vGeo);
      double r = hRatio[ic]->GetBinContent(b);
      if(r != 0. && (r > botMax || r < botMin))
        printf("NOTE: %s double ratio at %.0f-%.0f is %.3f, outside bottom-panel range [%.2f,%.2f]\n",
               classLabel[ic], lo, hRatio[ic]->GetXaxis()->GetBinUpEdge(b), r, botMin, botMax);
    }
  }

  const double lm = 0.15, rm = 0.05, tm = 0.07, bm = 0.28, split = 0.40;
  TCanvas *c = new TCanvas("cJetsPerZ_RCvsRCgeoCorr_unfold", "", 800, 800);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->Draw(); pUp->cd();

  TLegend *leg = new TLegend(0.36, 0.78, 0.94, 0.93);
  leg->SetBorderSize(0); leg->SetTextSize(0.028); leg->SetFillStyle(0); leg->SetNColumns(2);

  bool first = true;
  for(int ic = 0; ic < NClass; ic++){
    hRC[ic]->SetLineColor(colClass[ic]); hRC[ic]->SetMarkerColor(colClass[ic]);
    hRC[ic]->SetMarkerStyle(markOpen[ic]); hRC[ic]->SetMarkerSize(1.0);
    hRC[ic]->SetLineWidth(2); hRC[ic]->SetLineStyle(2);
    hRC[ic]->SetTitle("");
    hRC[ic]->GetYaxis()->SetTitle("#frac{1}{N_{Z}^{PbPb}} #frac{dN_{jet}^{PbPb}}{dp_{T}} / #frac{1}{N_{Z}^{pp}} #frac{dN_{jet}^{pp}}{dp_{T}}");
    hRC[ic]->GetYaxis()->SetTitleSize(0.048);
    hRC[ic]->GetYaxis()->SetLabelSize(0.042);
    hRC[ic]->GetYaxis()->SetTitleOffset(1.55);
    hRC[ic]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hRC[ic]->SetMinimum(topMin); hRC[ic]->SetMaximum(topMax);
    hRC[ic]->Draw(first ? "ep" : "ep same");
    first = false;

    hGeo[ic]->SetLineColor(colClass[ic]); hGeo[ic]->SetMarkerColor(colClass[ic]);
    hGeo[ic]->SetMarkerStyle(markFilled[ic]); hGeo[ic]->SetMarkerSize(1.0);
    hGeo[ic]->SetLineWidth(2);
    hGeo[ic]->Draw("ep same");

    leg->AddEntry(hRC[ic],  Form("%s, RC sub",         classLabel[ic]), "lp");
    leg->AddEntry(hGeo[ic], Form("%s, RC_{geoCorr} sub", classLabel[ic]), "lp");
  }
  leg->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.040);
  lat.DrawLatex(0.19, 0.955, "Jets per Z, with stitching, unfolded -- PbPb / pp");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  bool firstB = true;
  for(int ic = 0; ic < NClass; ic++){
    hRatio[ic]->SetLineColor(colClass[ic]); hRatio[ic]->SetMarkerColor(colClass[ic]);
    hRatio[ic]->SetMarkerStyle(markFilled[ic]); hRatio[ic]->SetMarkerSize(1.0);
    hRatio[ic]->SetLineWidth(2);
    hRatio[ic]->SetTitle("");
    hRatio[ic]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
    hRatio[ic]->GetYaxis()->SetTitle("RC_{geoCorr} / RC");
    hRatio[ic]->GetXaxis()->SetTitleSize(0.048*sc);
    hRatio[ic]->GetXaxis()->SetLabelSize(0.042*sc);
    hRatio[ic]->GetYaxis()->SetTitleSize(0.048*sc);
    hRatio[ic]->GetYaxis()->SetLabelSize(0.042*sc);
    hRatio[ic]->GetYaxis()->SetTitleOffset(1.35/sc);
    hRatio[ic]->GetYaxis()->SetNdivisions(505);
    hRatio[ic]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    hRatio[ic]->SetMinimum(botMin); hRatio[ic]->SetMaximum(botMax);
    hRatio[ic]->Draw(firstB ? "ep" : "ep same");
    firstB = false;
  }

  TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TString out = outDir + "JetsPerZ_RCvsRCgeoCorr_unfold.pdf";
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
