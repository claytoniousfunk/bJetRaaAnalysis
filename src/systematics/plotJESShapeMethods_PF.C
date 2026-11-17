// Jet energy scale of the b jets: the spectrum shape the shift is evaluated on (2026-11-16, user: shift a smooth
// fit of the measured spectrum instead of the measured spectrum itself, so its bin-to-bin fluctuations do not enter).
//   Figure 1 (fits): per system, the measured b-jet density (counts / GeV) with the global log-cubic fit of jeuShift.h,
//     and the inclusive one (the JEUdata spectrum, 5 GeV bins) with its LOCAL quadratic fits (half-width 0.3 in ln pT,
//     one per bin, drawn at the bin centres), and data / fit.
//   Figure 2 (methods): the JES assigned to the b jets per Z (pp and the four classes), to R_Z(b) and to the double
//     ratio, from the shift of the MEASURED spectrum (before 2026-11-16), of its FIT (nominal since 2026-11-16), and,
//     for the b jets per Z and R_Z(b), of the SIMULATED spectrum (jeuShiftMC_bJets_PF.C, the reco projection of the b
//     response; unquenched in PbPb, so its slope is not the data's).
// inputs  the 7-bin nominal b-jet and inclusive systematics files (fit) and their copies before the switch
//         (rootFiles/systematics/backup_2026-11-16_preJESfit/, measured), rootFiles/JEU/jeuShiftMC_bJets_PF_7bins.root
// Usage, from src/systematics/:  root -l -b -q plotJESShapeMethods_PF.C
// Output: figures/systematics/JEU/bJetJESShape_fits_PF<sfx>.pdf, bJetJESShape_methods_PF<sfx>.pdf
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include <vector>
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/gridPads.h"
#include "../../headers/functions/jeuShift.h"

void plotJESShapeMethods_PF(TString bTag = "_noFlatJEU_fakeGeoOnly_csvV2_2D7bins_nIterOpt_noRecoFloor_normSep_noStitch",
                            TString inclTag = "_noFlatJEU_fakeGeoOnly_pt7bins", TString ppVar = "_csvV2_2D7bins_nIterOpt",
                            TString mcFile = "jeuShiftMC_bJets_PF_7bins.root", TString sfx = "_7bins")
{
  initPlotStyle();
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  const TString D = Form("%s/rootFiles/systematics/", repo), Bk = D + "backup_2026-11-16_preJESfit/";
  TFile *fB = TFile::Open(D + "bJetsPerZSys_PF_MBcaloMatched_JEUppCorr" + bTag + ".root");
  TFile *fBold = TFile::Open(Bk + "bJetsPerZSys_PF_MBcaloMatched_JEUppCorr" + bTag + ".root");
  TFile *fI = TFile::Open(D + "inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr" + inclTag + ".root");
  TFile *fS = TFile::Open(Form("%s/rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_PFJets_pp%s.root", repo, ppVar.Data()));
  TFile *fM = TFile::Open(Form("%s/rootFiles/JEU/%s", repo, mcFile.Data()));
  if(!fB || !fBold || !fI || !fS || !fM){ printf("ERROR: an input did not open\n"); return; }
  const char *sys[5] = {"pp", "C4", "C3", "C2", "C1"}, *name[5] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"};
  const double xLo = 80., xHi = 200.;
  const TString outDir = Form("%s/figures/systematics/JEU/", repo); gSystem->mkdir(outDir, true);
  auto get = [&](TFile *f, TString n) -> TH1D* { TH1D *h = nullptr; f->GetObject(n, h); if(!h) printf("ERROR: %s missing in %s\n", n.Data(), f->GetName()); return h; };

  // ---- figure 1: the fits
  {
    TCanvas *cv = new TCanvas("cvFits", "", 2200, 1500);
    GridPads g = makeGridPads(cv, 5, 4, 0.07, 0.012, 0.075, 0.05);
    const char *rowLab[4] = {"#it{b} jets, d#it{N}/d#it{p}_{T} (GeV^{-1})", "data / fit", "inclusive, d#it{N}/d#it{p}_{T} (GeV^{-1})", "data / fit"};
    printf("fit chi2 / ndf (order 3): system | b jets | inclusive\n");
    double c2[2][5]; int nd[2][5];
    for(int s = 0; s < 5; s++){
      TH1D *m[2];
      if(s == 0){ m[0] = (TH1D*) get(fS, "h_bJetPt")->Clone("ppMeasC"); m[0]->SetDirectory(nullptr);
        for(int j = 1; j <= m[0]->GetNbinsX(); j++){ m[0]->SetBinContent(j, m[0]->GetBinContent(j)*m[0]->GetBinWidth(j)); m[0]->SetBinError(j, m[0]->GetBinError(j)*m[0]->GetBinWidth(j)); } }
      else m[0] = get(fB, Form("h_bJEU_measNom_%s", sys[s]));
      m[1] = get(fI, Form("h_jeuDataMeas_%s", sys[s]));
      if(!m[0] || !m[1]) return;
      for(int t = 0; t < 2; t++){
        // b jets: the global log-cubic fit; inclusive: the local quadratic (half-width 0.3 in ln pT), one per bin
        const std::vector<double> p = t == 0 ? jeuFitSpectrum(m[t]) : std::vector<double>(); const JeuFitInfo I = jeuLastFit();
        c2[t][s] = t == 0 ? I.chi2 : 0.; nd[t][s] = t == 0 ? I.ndf : 0;
        auto binFit = [&](int j){ const double a = m[t]->GetBinLowEdge(j), b = m[t]->GetBinLowEdge(j+1);
          return t == 0 ? jeuFitInt(p, a, b) : jeuFitInt(jeuLocalQuad(m[t], j, 0.3, 0.01), a, b); };
        TH1D *dens = (TH1D*) m[t]->Clone(Form("dens_%d_%d", t, s)), *rat = (TH1D*) m[t]->Clone(Form("rat_%d_%d", t, s));
        dens->SetDirectory(nullptr); rat->SetDirectory(nullptr);
        double ymin = 1e30, ymax = 0.;
        for(int j = 1; j <= dens->GetNbinsX(); j++){ const double a = dens->GetBinLowEdge(j), b = dens->GetBinLowEdge(j+1), w = b - a, c = m[t]->GetBinContent(j);
          dens->SetBinContent(j, c/w); dens->SetBinError(j, m[t]->GetBinError(j)/w);
          const double fit = c > 0. ? binFit(j) : 0.; rat->SetBinContent(j, c > 0. && fit > 0. ? c/fit : 0.); rat->SetBinError(j, c > 0. && fit > 0. ? m[t]->GetBinError(j)/fit : 0.);
          if(c > 0. && a >= 50. - 1e-6 && b <= 300. + 1e-6){ ymin = std::min(ymin, c/w); ymax = std::max(ymax, c/w); } }
        TGraph *gf = new TGraph();
        if(t == 0){ for(int k = 0; k <= 200; k++){ const double x = 50.*pow(6., k/200.); gf->SetPoint(k, x, jeuFitDens(p, x)); } }
        else { for(int j = 1; j <= m[t]->GetNbinsX(); j++){ const double x = m[t]->GetBinCenter(j); if(x < 50. || x > 300. || m[t]->GetBinContent(j) <= 0.) continue;
          gf->SetPoint(gf->GetN(), x, jeuFitDens(jeuLocalQuad(m[t], j, 0.3, 0.01), x)); } }
        // spectrum
        g.at(s, 2*t)->cd(); gPad->SetLogy(); gPad->SetLogx();
        TH1F *fr = gPad->DrawFrame(50., 0.3*ymin, 300., 4.*ymax); g.styleFrame(fr, s, 2*t, 0.075); fr->GetXaxis()->SetLabelSize(0.);
        styleLine((TH1*) gf->GetHistogram(), okabeHex[6], 4); gf->SetLineColor(TColor::GetColor(okabeHex[6])); gf->SetLineWidth(4); gf->Draw("L SAME");
        styleH(dens, "#000000", markFilledCircle, 1.6); dens->SetLineWidth(2); dens->Draw("E1 X0 SAME");
        gPad->RedrawAxis();
        if(t == 0){ TLatex nm; nm.SetNDC(); nm.SetTextFont(62); nm.SetTextSize(0.11*g.sizeScale()); nm.SetTextAlign(33); nm.DrawLatex(0.95, g.frameTop(0) - 0.05, name[s]); }   // top right: a falling spectrum leaves it empty
        // data / fit
        g.at(s, 2*t + 1)->cd(); gPad->SetLogx();
        TH1F *fr2 = gPad->DrawFrame(50., 0.7, 300., 1.3); g.styleFrame(fr2, s, 2*t + 1, 0.075);
        fr2->GetXaxis()->SetLabelSize(t == 1 ? 0.075 : 0.); fr2->GetXaxis()->SetMoreLogLabels(); fr2->GetXaxis()->SetNoExponent();
        TLine l; l.SetLineColor(kGray + 2); l.SetLineWidth(2); l.DrawLine(50., 1., 300., 1.);
        { TLine e; e.SetLineColor(kGray + 1); e.SetLineStyle(3); e.SetLineWidth(2); e.DrawLine(80., 0.7, 80., 1.3); e.DrawLine(200., 0.7, 200., 1.3); }   // reported range
        styleH(rat, "#000000", markFilledCircle, 1.6); rat->SetLineWidth(2); rat->Draw("E1 X0 SAME");
        TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.11*g.sizeScale()); if(t == 0) la.DrawLatex(g.frameLeft(s) + 0.06, 0.80, Form("#chi^{2}/ndf = %.0f/%d", I.chi2, I.ndf));
        gPad->RedrawAxis();
      }
      printf("  %-7s | %6.1f / %d | %6.1f / %d\n", name[s], c2[0][s], nd[0][s], c2[1][s], nd[1][s]);
    }
    g.at(0, 0)->cd();
    { TLegend *lg = makeLegend(g.frameLeft(0) + 0.05, 0.08, 0.98, 0.42, 0.105*g.sizeScale());
      TH1D *pd = new TH1D("lgD", "", 1, 0, 1); styleH(pd, "#000000", markFilledCircle, 1.6); TH1D *pf = new TH1D("lgF", "", 1, 0, 1); styleLine(pf, okabeHex[6], 4);
      lg->AddEntry(pd, "measured", "pe"); lg->AddEntry(pf, "smooth fit", "l"); lg->Draw(); }
    cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.028); t.SetTextAngle(90); t.SetTextAlign(21);
    const double ih = (1. - g.bm - g.tm)/4.; for(int r = 0; r < 4; r++) t.DrawLatex(0.022, 1. - g.tm - (r + 0.5)*ih, rowLab[r]);
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "", 0.032, 0.012, 0.016);
    g.cmsLabel(cv, "Internal", "measured spectra and their smooth fits,   #sqrt{#it{s}_{NN}} = 5.02 TeV", 1.0);
    savePdfTight(cv, outDir + "bJetJESShape_fits_PF" + sfx + ".pdf");
  }

  // ---- figure 2: the three methods
  {
    // MC: max over up / down; ratios per variation
    auto mcVar = [&](int s, int k) -> TH1D* { return get(fM, Form("h_jeuRatio_%s_%s", k == 0 ? "up" : "down", sys[s])); };
    auto toPct = [&](TH1D *h, const TString &n) -> TH1D* { TH1D *o = (TH1D*) h->Clone(n); o->SetDirectory(nullptr);
      for(int j = 0; j <= o->GetNbinsX() + 1; j++){ o->SetBinContent(j, 100.*fabs(h->GetBinContent(j))); o->SetBinError(j, 0.); } o->GetXaxis()->SetRangeUser(xLo + 1e-3, xHi - 1e-3); return o; };
    auto mcAssigned = [&](int s, int ref, const TString &n) -> TH1D* {   // ref < 0: the system alone
      TH1D *u = mcVar(s, 0), *d = mcVar(s, 1); TH1D *o = (TH1D*) u->Clone(n); o->SetDirectory(nullptr); o->Reset();
      for(int j = 1; j <= o->GetNbinsX(); j++){ double m = 0.;
        for(int k = 0; k < 2; k++){ const double r = mcVar(s, k)->GetBinContent(j), rr = ref < 0 ? 1. : mcVar(ref, k)->GetBinContent(j); m = std::max(m, fabs(r/rr - 1.)); }
        o->SetBinContent(j, m); }
      return o; };
    const char *hexFit = "#000000", *hexMeas = "#D55E00", *hexMC = "#0072B2";
    TCanvas *cv = new TCanvas("cvMeth", "", 2200, 1500);
    GridPads g = makeGridPads(cv, 5, 3, 0.07, 0.012, 0.075, 0.05);
    const char *rowLab[3] = {"JES on #it{b} jets per #it{N}_{Z} (%)", "JES on #it{R}_{Z}(#it{b}) (%)", "JES on #it{R}_{Z}(#it{b}) / #it{R}_{Z}(incl) (%)"};
    printf("\nassigned JES (%%) at 80-90 / 100-110 / 150-200 GeV: measured shape -> fit (MC)\n");
    // every row first, then ONE y range for the figure (2026-11-17, user: a common scale within each figure)
    std::vector<TH1D*> H[3][5]; double yMax = 0.;
    for(int r = 0; r < 3; r++) for(int s = 0; s < 5; s++){
      if(r > 0 && s == 0) continue;
      const TString q = r == 0 ? Form("b_%s", sys[s]) : r == 1 ? Form("b_PbPbOverPP_%s", sys[s]) : Form("bOverIncl_PbPbOverPP_%s", sys[s]);
      TH1D *hf = toPct(get(fB, "h_sysRel_JEU_" + q), Form("hf_%d_%d", r, s)), *hm = toPct(get(fBold, "h_sysRel_JEU_" + q), Form("hm_%d_%d", r, s));
      TH1D *hc = r < 2 ? toPct(mcAssigned(s, r == 0 ? -1 : 0, Form("mcA_%d_%d", r, s)), Form("hc_%d_%d", r, s)) : nullptr;
      H[r][s] = {hf, hm, hc};
      for(TH1D *h : H[r][s]) if(h) for(int j = 1; j <= h->GetNbinsX(); j++) if(h->GetBinLowEdge(j) >= xLo - 1e-6 && h->GetBinLowEdge(j) < xHi - 1e-6) yMax = std::max(yMax, h->GetBinContent(j));
      printf("  %-38s", (TString(rowLab[r]) + " " + name[s]).Data());
      for(double x : {85., 105., 175.}){ const int j = hf->FindBin(x); printf("  %5.1f -> %5.1f", hm->GetBinContent(j), hf->GetBinContent(j)); if(hc) printf(" (%4.1f)", hc->GetBinContent(hc->FindBin(x))); }
      printf("\n");
    }
    yMax = 5.*ceil(1.45*yMax/5.);
    // in-panel text in PIXELS (font 43 / 63): the same printed size in every pad (2026-11-17, user)
    const double k = cv->GetWw()/1800., pxText = 34.*k, pxLeg = 30.*k;
    auto hpx = [](){ return gPad->GetWh()*gPad->GetAbsHNDC(); }; auto wpx = [](){ return gPad->GetWw()*gPad->GetAbsWNDC(); };
    auto ptext = [&](double x, double y, const char *str, int font, int align, Color_t col){ TLatex t; t.SetNDC(); t.SetTextFont(font); t.SetTextSize(pxText); t.SetTextAlign(align); t.SetTextColor(col); t.DrawLatex(x, y, str); };
    for(int r = 0; r < 3; r++) for(int s = 0; s < 5; s++){
      g.at(s, r)->cd();
      TH1F *fr = gPad->DrawFrame(xLo, 0., xHi, yMax); g.styleFrame(fr, s, r, 0.075);
      fr->GetXaxis()->SetLabelSize(r == 2 ? 0.075*g.padW(0)/g.padW(s) : 0.);
      if(r > 0 && s == 0){   // pp is the reference of the ratios
        ptext(0.5*(g.frameLeft(0) + 1.), 0.5*(g.frameTop(r) + g.frameBottom(r)), "reference", 43, 22, kGray + 2); gPad->RedrawAxis(); continue; }
      if(H[r][s][2]){ styleLine(H[r][s][2], hexMC, 4); H[r][s][2]->SetLineStyle(2); H[r][s][2]->Draw("HIST SAME"); }
      styleLine(H[r][s][1], hexMeas, 4); H[r][s][1]->SetLineStyle(1); H[r][s][1]->Draw("HIST SAME");
      styleLine(H[r][s][0], hexFit, 5); H[r][s][0]->Draw("HIST SAME");
      gPad->RedrawAxis();
      if(r == 0) ptext(g.frameLeft(s) + 20./wpx(), g.frameTop(r) - 14./hpx(), name[s], 63, 13, kBlack);
      if(r == 1) ptext(g.frameLeft(s) + 20./wpx(), g.frameTop(r) - 14./hpx(), Form("%s / pp", name[s]), 63, 13, kBlack);
    }
    // legend in the pp panel of the first row (its values stay below a third of the range)
    g.at(0, 0)->cd();
    { const double top = g.frameTop(0) - (14. + 1.6*pxText)/hpx(), rowH = 1.4*pxLeg/hpx();
      TLegend *lg = makeLegend(g.frameLeft(0) + 20./wpx(), top - 3*rowH, 0.98, top); lg->SetTextFont(43); lg->SetTextSize(pxLeg);
      TH1D *a = new TH1D("lA", "", 1, 0, 1), *b = new TH1D("lB", "", 1, 0, 1), *c = new TH1D("lC", "", 1, 0, 1);
      styleLine(a, hexFit, 5); styleLine(b, hexMeas, 4); styleLine(c, hexMC, 4); c->SetLineStyle(2);
      lg->AddEntry(a, "smooth fit", "l"); lg->AddEntry(b, "measured shape", "l"); lg->AddEntry(c, "simulated shape", "l"); lg->Draw(); }
    cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.026); t.SetTextAngle(90); t.SetTextAlign(21);
    const double ih = (1. - g.bm - g.tm)/3.; for(int r = 0; r < 3; r++) t.DrawLatex(0.022, 1. - g.tm - (r + 0.5)*ih, rowLab[r]);
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "", 0.032, 0.012, 0.016);
    g.cmsLabel(cv, "Internal", "jet energy scale of the #it{b} jets: spectrum shape the shift is applied to,   #sqrt{#it{s}_{NN}} = 5.02 TeV", 0.8);
    savePdfTight(cv, outDir + "bJetJESShape_methods_PF" + sfx + ".pdf");
  }
}
