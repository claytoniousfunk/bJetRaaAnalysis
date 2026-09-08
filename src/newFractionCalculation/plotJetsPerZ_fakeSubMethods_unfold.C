// Jets-per-Z (PbPb/pp), with stitching (useMinBiasOnly=false -- the real
// calculateRAA.C chain), fine pT bins, UNFOLDED, comparing three fake-jet
// subtraction methods: plain RC, RC_geoCorr, and dPT (PF-PFCs). Extends
// plotJetsPerZ_RCvsRCgeoCorr_unfold.C ("for completeness") with the dPT
// method; see that file for the two-method version and the general reasoning.
// See plotJetsPerZ_fakeSubMethods_noUnfold.C for the no-unfolding twin.
//
// Top panel: jets-per-Z, all three methods overlaid per class. Colour encodes
// centrality class, marker shape encodes method (open square = RC, filled
// circle = RC_geoCorr, open cross = dPT) -- with 12 class*method combinations,
// repeating a full "class, method" legend entry per curve would be unreadable,
// so the legend is split: a class-colour key (coloured lines, no marker) and a
// separate method-marker key (black markers, no class meaning).
// Bottom panel: (RC_geoCorr sub)/(RC sub) and (dPT sub)/(RC sub), both against
// the plain-RC baseline, per class.
//
// *** dPT over-subtracts badly in 0-10%: the jets-per-Z numerator itself goes
// negative for 6 of 20 fine bins there (60-90 GeV), down to -0.34 at 60-65 GeV.
// This is not a plotting artifact -- it reproduces, at the full stitched/
// unfolded level, the same breakdown already seen in the MinBias-only shape
// study earlier this session (150% of raw yield "subtracted" in 0-10%). Every
// other class stays positive throughout. The top-panel range is widened to
// show the negative excursion rather than clip it. ***
//
// Reads r_C{1..4}_fine from three calculateRAA.C runs, identical except for
// which fake-jet estimate was subtracted:
//   histograms_stitched_RC.root        : fakeJetsPath = fakeJets.root
//   histograms_stitched_RCgeoCorr.root : fakeJetsPath = fakeJets_RCgeoCorr.root
//   histograms_stitched_dPT.root       : fakeJetsPath = fakeJets_dPT.root
// All three with useMinBiasOnly=false, e.g.:
//   root -l -b -q -e '.L calculateRAA.C' \
//     -e 'useMinBiasOnly=false; fakeJetsPath="../../rootFiles/fakeJets/fakeJets_dPT.root"; \
//         outputRootPath="./rootFiles/JetsPerZ/histograms_stitched_dPT.root";' \
//     -e 'calculateRAA()'
//
// Usage: root -l -b -q 'plotJetsPerZ_fakeSubMethods_unfold.C'
// Run from: src/newFractionCalculation/

const char *f_RC_path =
  "rootFiles/JetsPerZ/histograms_stitched_RC.root";
const char *f_RCgeoCorr_path =
  "rootFiles/JetsPerZ/histograms_stitched_RCgeoCorr.root";
const char *f_dPT_path =
  "rootFiles/JetsPerZ/histograms_stitched_dPT.root";

TString outDir = "../../figures/JetsPerZ/";

const int   NClass = 4;
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };
const char *classHex[NClass]   = { "#D55E00", "#0072B2", "#009E73", "#CC79A7" };

const int   NMethod = 3;
const char *methodLabel[NMethod] = { "RC sub", "RC_{geoCorr} sub", "dPT sub" };
// Centre-symmetric shapes only.
const int   methodMark[NMethod]  = { 25, 20, 28 };
const double methodSize[NMethod] = { 1.3, 1.0, 1.3 };

const double plotPtMin = 60., plotPtMax = 300.;
const double topMin = -0.5, topMax = 1.1;
const double botMin = -1.2, botMax = 1.6;

void plotJetsPerZ_fakeSubMethods_unfold(){

  TFile *fRC  = TFile::Open(f_RC_path);
  TFile *fGeo = TFile::Open(f_RCgeoCorr_path);
  TFile *fDpt = TFile::Open(f_dPT_path);
  if(!fRC || fRC->IsZombie() || !fGeo || fGeo->IsZombie() || !fDpt || fDpt->IsZombie()){
    printf("ERROR: cannot open input files\n"); return;
  }
  gSystem->mkdir(outDir, kTRUE);
  gStyle->SetOptStat(0);

  TH1D *hRC[NClass], *hGeo[NClass], *hDpt[NClass], *rGeo[NClass], *rDpt[NClass];
  int colClass[NClass];

  printf("%-8s %14s %14s %14s %14s\n", "class", "pT [GeV]", "RC", "RC_geoCorr", "dPT");
  for(int ic = 0; ic < NClass; ic++){
    colClass[ic] = TColor::GetColor(classHex[ic]);

    fRC ->GetObject(Form("r_C%d_fine", ic+1), hRC[ic]);
    fGeo->GetObject(Form("r_C%d_fine", ic+1), hGeo[ic]);
    fDpt->GetObject(Form("r_C%d_fine", ic+1), hDpt[ic]);
    if(!hRC[ic] || !hGeo[ic] || !hDpt[ic]){ printf("WARNING: missing C%d\n", ic+1); return; }
    hRC[ic]->SetDirectory(nullptr); hGeo[ic]->SetDirectory(nullptr); hDpt[ic]->SetDirectory(nullptr);

    rGeo[ic] = (TH1D*) hGeo[ic]->Clone(Form("rGeo_C%d", ic+1));
    rGeo[ic]->SetDirectory(nullptr); rGeo[ic]->Divide(hRC[ic]);
    rDpt[ic] = (TH1D*) hDpt[ic]->Clone(Form("rDpt_C%d", ic+1));
    rDpt[ic]->SetDirectory(nullptr); rDpt[ic]->Divide(hRC[ic]);

    for(int b = 1; b <= hRC[ic]->GetNbinsX(); b++){
      double lo = hRC[ic]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      printf("%-8s %5.0f-%-6.0f %14.4f %14.4f %14.4f\n",
             classLabel[ic], lo, hRC[ic]->GetXaxis()->GetBinUpEdge(b),
             hRC[ic]->GetBinContent(b), hGeo[ic]->GetBinContent(b), hDpt[ic]->GetBinContent(b));
    }
  }

  // Warn rather than silently clip if the fixed display ranges don't cover the
  // actual data in the plotted window.
  for(int ic = 0; ic < NClass; ic++){
    for(int b = 1; b <= hRC[ic]->GetNbinsX(); b++){
      double lo = hRC[ic]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      double vals[3] = { hRC[ic]->GetBinContent(b), hGeo[ic]->GetBinContent(b), hDpt[ic]->GetBinContent(b) };
      for(int m = 0; m < 3; m++)
        if(vals[m] > topMax || vals[m] < topMin)
          printf("NOTE: %s %s at %.0f-%.0f exceeds top-panel range [%.2f,%.2f]: %.3f\n",
                 classLabel[ic], methodLabel[m], lo, hRC[ic]->GetXaxis()->GetBinUpEdge(b), topMin, topMax, vals[m]);
      double rg = rGeo[ic]->GetBinContent(b), rd = rDpt[ic]->GetBinContent(b);
      if(rg != 0. && (rg > botMax || rg < botMin))
        printf("NOTE: %s RC_geoCorr/RC at %.0f-%.0f is %.3f, outside bottom-panel range [%.2f,%.2f]\n",
               classLabel[ic], lo, hRC[ic]->GetXaxis()->GetBinUpEdge(b), rg, botMin, botMax);
      if(rd != 0. && (rd > botMax || rd < botMin))
        printf("NOTE: %s dPT/RC at %.0f-%.0f is %.3f, outside bottom-panel range [%.2f,%.2f]\n",
               classLabel[ic], lo, hRC[ic]->GetXaxis()->GetBinUpEdge(b), rd, botMin, botMax);
    }
  }

  const double lm = 0.15, rm = 0.05, tm = 0.07, bm = 0.28, split = 0.40;
  TCanvas *c = new TCanvas("cJetsPerZ_fakeSubMethods_unfold", "", 800, 800);

  TPad *pUp = new TPad("pUp", "", 0, split, 1, 1);
  pUp->SetLeftMargin(lm); pUp->SetRightMargin(rm);
  pUp->SetTopMargin(tm);  pUp->SetBottomMargin(0.);
  pUp->Draw(); pUp->cd();

  // Class-colour key: coloured lines, no marker meaning.
  TLegend *legClass = new TLegend(0.17, 0.72, 0.45, 0.92);
  legClass->SetBorderSize(0); legClass->SetTextSize(0.032); legClass->SetFillStyle(0);
  // Method key: black markers, no colour meaning.
  TLegend *legMethod = new TLegend(0.55, 0.72, 0.94, 0.92);
  legMethod->SetBorderSize(0); legMethod->SetTextSize(0.032); legMethod->SetFillStyle(0);

  TH1D *methodKey[NMethod];
  for(int m = 0; m < NMethod; m++){
    methodKey[m] = new TH1D(Form("methodKey_%d", m), "", 1, 0, 1);
    methodKey[m]->SetDirectory(nullptr);
    methodKey[m]->SetMarkerStyle(methodMark[m]);
    methodKey[m]->SetMarkerSize(methodSize[m]);
    methodKey[m]->SetMarkerColor(kBlack);
    methodKey[m]->SetLineColor(kBlack);
    legMethod->AddEntry(methodKey[m], methodLabel[m], "p");
  }

  bool first = true;
  for(int ic = 0; ic < NClass; ic++){
    hRC[ic]->SetLineColor(colClass[ic]); hRC[ic]->SetMarkerColor(colClass[ic]);
    hRC[ic]->SetMarkerStyle(methodMark[0]); hRC[ic]->SetMarkerSize(methodSize[0]);
    hRC[ic]->SetLineWidth(2);
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
    hGeo[ic]->SetMarkerStyle(methodMark[1]); hGeo[ic]->SetMarkerSize(methodSize[1]);
    hGeo[ic]->SetLineWidth(2);
    hGeo[ic]->Draw("ep same");

    hDpt[ic]->SetLineColor(colClass[ic]); hDpt[ic]->SetMarkerColor(colClass[ic]);
    hDpt[ic]->SetMarkerStyle(methodMark[2]); hDpt[ic]->SetMarkerSize(methodSize[2]);
    hDpt[ic]->SetLineWidth(2); hDpt[ic]->SetLineStyle(2);
    hDpt[ic]->Draw("ep same");

    legClass->AddEntry(hRC[ic], classLabel[ic], "l");
  }
  legClass->Draw();
  legMethod->Draw();

  TLine *zero = new TLine(plotPtMin, 0., plotPtMax, 0.);
  zero->SetLineStyle(3); zero->SetLineColor(kGray+1); zero->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038);
  lat.DrawLatex(0.17, 0.955, "Jets per Z, with stitching, unfolded -- PbPb / pp");

  c->cd();
  TPad *pDn = new TPad("pDn", "", 0, 0, 1, split);
  pDn->SetLeftMargin(lm); pDn->SetRightMargin(rm);
  pDn->SetTopMargin(0.);  pDn->SetBottomMargin(bm);
  pDn->Draw(); pDn->cd();

  const double sc = (1.-split)/split;
  bool firstB = true;
  for(int ic = 0; ic < NClass; ic++){
    rGeo[ic]->SetLineColor(colClass[ic]); rGeo[ic]->SetMarkerColor(colClass[ic]);
    rGeo[ic]->SetMarkerStyle(methodMark[1]); rGeo[ic]->SetMarkerSize(methodSize[1]);
    rGeo[ic]->SetLineWidth(2);
    rGeo[ic]->SetTitle("");
    rGeo[ic]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
    rGeo[ic]->GetYaxis()->SetTitle("method / RC");
    rGeo[ic]->GetXaxis()->SetTitleSize(0.048*sc);
    rGeo[ic]->GetXaxis()->SetLabelSize(0.042*sc);
    rGeo[ic]->GetYaxis()->SetTitleSize(0.048*sc);
    rGeo[ic]->GetYaxis()->SetLabelSize(0.042*sc);
    rGeo[ic]->GetYaxis()->SetTitleOffset(1.35/sc);
    rGeo[ic]->GetYaxis()->SetNdivisions(505);
    rGeo[ic]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
    rGeo[ic]->SetMinimum(botMin); rGeo[ic]->SetMaximum(botMax);
    rGeo[ic]->Draw(firstB ? "ep" : "ep same");
    firstB = false;

    rDpt[ic]->SetLineColor(colClass[ic]); rDpt[ic]->SetMarkerColor(colClass[ic]);
    rDpt[ic]->SetMarkerStyle(methodMark[2]); rDpt[ic]->SetMarkerSize(methodSize[2]);
    rDpt[ic]->SetLineWidth(2); rDpt[ic]->SetLineStyle(2);
    rDpt[ic]->Draw("ep same");
  }

  TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  TString out = outDir + "JetsPerZ_fakeSubMethods_unfold.pdf";
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
