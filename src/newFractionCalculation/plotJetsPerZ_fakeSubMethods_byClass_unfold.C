// Jets-per-Z (PbPb/pp), with stitching, UNFOLDED, comparing three fake-jet
// subtraction methods (RC, RC_geoCorr, dPT) plus the unsubtracted "nominal"
// baseline -- panel-per-class layout, twin of
// plotJetsPerZ_fakeSubMethods_unfold.C's single combined panel. Splitting by
// class frees up colour to encode method directly (no split class-colour/
// method-marker legend needed) and lets each class use its own y-range
// instead of one fixed range wide enough for 0-10%'s dPT excursion alone.
//
// Reads r_C{1..4}_fine from four calculateRAA.C runs, identical except which
// fake-jet estimate (if any) was subtracted -- see plotJetsPerZ_fakeSubMethods_
// unfold.C's header for how the RC/RCgeoCorr/dPT files were produced. The
// no-sub baseline is doFakeJetSubtraction=false, e.g.:
//   root -l -b -q -e '.L calculateRAA.C' \
//     -e 'useMinBiasOnly=false; doFakeJetSubtraction=false; \
//         outputRootPath="./rootFiles/JetsPerZ/histograms_stitched_noSub.root";' \
//     -e 'calculateRAA()'
//
// *** dPT over-subtracts badly in 0-10%: the jets-per-Z numerator itself goes
// negative for 6 of 20 fine bins there (60-90 GeV), down to -0.34 at 60-65
// GeV. Not a plotting artifact -- see plotJetsPerZ_fakeSubMethods_unfold.C.
// Per-class y-ranges here are auto-computed specifically so this shows
// rather than clips. ***
//
// Bottom ratio panel is (subtracted method)/(no sub) -- i.e. what fraction of
// the raw, unsubtracted spectrum survives each subtraction method. no sub is
// the reference/denominator here, not a plotted series.
//
// Usage: root -l -b -q 'plotJetsPerZ_fakeSubMethods_byClass_unfold.C'
// Run from: src/newFractionCalculation/

TString outDir = "../../figures/JetsPerZ/";
const char *outName = "JetsPerZ_fakeSubMethods_byClass_unfold.pdf";

const int   NClass = 4;
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };

struct Method { const char *path; const char *label; const char *hex; int marker; bool inAutoRange; };
const int NMethod = 4;
Method methods[NMethod] = {
  { "rootFiles/JetsPerZ/histograms_stitched_RC.root",        "RC",           "#0072B2", 20, true  },
  { "rootFiles/JetsPerZ/histograms_stitched_RCgeoCorr.root", "RC_{geoCorr}", "#009E73", 21, true  },
  { "rootFiles/JetsPerZ/histograms_stitched_dPT.root",       "dPT",          "#D55E00", 33, true  },
  { "rootFiles/JetsPerZ/histograms_stitched_noSub.root",     "no sub",       "#000000", 34, false },
};

const double plotPtMin = 60., plotPtMax = 300.;

// Index of the reference/denominator method for the ratio panel -- "no sub",
// last in the methods[] array.
const int REF = NMethod - 1;

void plotJetsPerZ_fakeSubMethods_byClass_unfold(){

  TFile *f[NMethod];
  for(int mi = 0; mi < NMethod; mi++){
    f[mi] = TFile::Open(methods[mi].path);
    if(!f[mi] || f[mi]->IsZombie()){ printf("ERROR: cannot open %s\n", methods[mi].path); return; }
  }
  gSystem->mkdir(outDir, kTRUE);
  gStyle->SetOptStat(0);

  TH1D *hYield[NClass][NMethod];
  TH1D *hRatio[NClass][NMethod];
  double topMin[NClass], topMax[NClass], ratioMin[NClass], ratioMax[NClass];
  int col[NMethod];
  for(int mi = 0; mi < NMethod; mi++) col[mi] = TColor::GetColor(methods[mi].hex);

  printf("%-8s %14s %14s %14s %14s %14s\n", "class", "pT [GeV]", "RC", "RC_geoCorr", "dPT", "no sub");
  for(int ic = 0; ic < NClass; ic++){
    for(int mi = 0; mi < NMethod; mi++){
      f[mi]->GetObject(Form("r_C%d_fine", ic+1), hYield[ic][mi]);
      if(!hYield[ic][mi]){ printf("WARNING: missing C%d for %s\n", ic+1, methods[mi].label); return; }
      hYield[ic][mi]->SetDirectory(nullptr);
    }

    for(int mi = 0; mi < REF; mi++){
      hRatio[ic][mi] = (TH1D*) hYield[ic][mi]->Clone(Form("hRatio_C%d_m%d", ic+1, mi));
      hRatio[ic][mi]->SetDirectory(nullptr);
      hRatio[ic][mi]->Divide(hYield[ic][REF]);
    }

    topMin[ic] = 1e300; topMax[ic] = -1e300;
    ratioMin[ic] = 1e300; ratioMax[ic] = -1e300;
    for(int b = 1; b <= hYield[ic][0]->GetNbinsX(); b++){
      double lo = hYield[ic][0]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      printf("%-8s %5.0f-%-6.0f %14.4f %14.4f %14.4f %14.4f\n",
             classLabel[ic], lo, hYield[ic][0]->GetXaxis()->GetBinUpEdge(b),
             hYield[ic][0]->GetBinContent(b), hYield[ic][1]->GetBinContent(b),
             hYield[ic][2]->GetBinContent(b), hYield[ic][3]->GetBinContent(b));
      for(int mi = 0; mi < NMethod; mi++){
        if(!methods[mi].inAutoRange) continue;
        double v = hYield[ic][mi]->GetBinContent(b);
        if(v == 0.) continue;
        if(v < topMin[ic]) topMin[ic] = v;
        if(v > topMax[ic]) topMax[ic] = v;
      }
      for(int mi = 0; mi < REF; mi++){
        if(!methods[mi].inAutoRange) continue;
        double r = hRatio[ic][mi]->GetBinContent(b);
        if(r == 0.) continue;
        if(r < ratioMin[ic]) ratioMin[ic] = r;
        if(r > ratioMax[ic]) ratioMax[ic] = r;
      }
    }
    double topPad = (topMax[ic] - topMin[ic]) * 0.20;
    topMin[ic] -= topPad; topMax[ic] += topPad;
    double ratioPad = (ratioMax[ic] - ratioMin[ic]) * 0.20;
    ratioMin[ic] -= ratioPad; ratioMax[ic] += ratioPad;
  }

  TCanvas *c = new TCanvas("cJetsPerZ_fakeSubMethods_byClass_unfold", "", 1300, 1100);

  const double split = 0.35;
  const double sc = (1. - split) / split;

  for(int ic = 0; ic < NClass; ic++){
    int row = ic / 2, colIdx = ic % 2;
    double x0 = colIdx * 0.5, x1 = (colIdx + 1) * 0.5;
    double y1 = 1. - row * 0.5, y0 = 1. - (row + 1) * 0.5;

    c->cd();
    TPad *outer = new TPad(Form("outer_%d", ic), "", x0, y0, x1, y1);
    outer->SetMargin(0., 0., 0., 0.);
    outer->Draw();
    outer->cd();

    TPad *pUp = new TPad(Form("pUp_%d", ic), "", 0, split, 1, 1);
    pUp->SetLeftMargin(0.24); pUp->SetRightMargin(0.05);
    pUp->SetTopMargin(0.10);  pUp->SetBottomMargin(0.);
    pUp->Draw();

    TPad *pDn = new TPad(Form("pDn_%d", ic), "", 0, 0, 1, split);
    pDn->SetLeftMargin(0.24); pDn->SetRightMargin(0.05);
    pDn->SetTopMargin(0.);    pDn->SetBottomMargin(0.34);
    pDn->Draw();

    pUp->cd();
    bool first = true;
    for(int mi = 0; mi < NMethod; mi++){
      TH1D *h = hYield[ic][mi];
      h->SetLineColor(col[mi]); h->SetMarkerColor(col[mi]);
      h->SetMarkerStyle(methods[mi].marker); h->SetMarkerSize(0.9);
      h->SetLineWidth(2);
      h->SetTitle("");
      h->GetXaxis()->SetLabelSize(0.); h->GetXaxis()->SetTitleSize(0.);
      h->GetYaxis()->SetTitle("Jets per Z, PbPb/pp");
      h->GetYaxis()->SetTitleSize(0.058); h->GetYaxis()->SetLabelSize(0.052);
      h->GetYaxis()->SetTitleOffset(1.35);
      h->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
      h->SetMinimum(topMin[ic]); h->SetMaximum(topMax[ic]);
      h->Draw(first ? "ep" : "ep same");
      first = false;
    }
    TLine *zero = new TLine(plotPtMin, 0., plotPtMax, 0.);
    zero->SetLineStyle(3); zero->SetLineColor(kGray+1); zero->Draw();

    TLatex clat; clat.SetNDC(); clat.SetTextSize(0.062); clat.SetTextFont(62);
    clat.SetTextAlign(13);
    clat.DrawLatex(0.24 + 0.03, 1. - 0.10 - 0.03, Form("PbPb %s", classLabel[ic]));

    if(ic == 0){
      TLegend *leg = new TLegend(0.68, 0.10, 0.94, 0.40);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.044);
      for(int mi = 0; mi < NMethod; mi++)
        leg->AddEntry(hYield[ic][mi], methods[mi].label, "lp");
      leg->Draw();
    }

    pDn->cd();
    bool firstR = true;
    for(int mi = 0; mi < REF; mi++){
      TH1D *r = hRatio[ic][mi];
      r->SetLineColor(col[mi]); r->SetMarkerColor(col[mi]);
      r->SetMarkerStyle(methods[mi].marker); r->SetMarkerSize(0.9 * sc / 1.6);
      r->SetLineWidth(2);
      r->SetTitle("");
      r->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
      r->GetYaxis()->SetTitle("method / no sub");
      r->GetXaxis()->SetTitleSize(0.062 * sc); r->GetXaxis()->SetLabelSize(0.055 * sc);
      r->GetYaxis()->SetTitleSize(0.062 * sc); r->GetYaxis()->SetLabelSize(0.055 * sc);
      r->GetYaxis()->SetTitleOffset(1.25 / sc);
      r->GetYaxis()->SetNdivisions(505);
      r->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
      r->SetMinimum(ratioMin[ic]); r->SetMaximum(ratioMax[ic]);
      r->Draw(firstR ? "ep" : "ep same");
      firstR = false;
    }
    TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
    one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();
  }

  c->cd();
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.018);
  lat.DrawLatex(0.06, 0.99, "Jets per Z, with stitching, unfolded, PbPb/pp -- background-subtraction method comparison");

  TString out = outDir + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
