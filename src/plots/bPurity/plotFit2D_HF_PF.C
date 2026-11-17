// The NOMINAL 2D b-purity fit drawn (2026-11-16, user: make the 2D fit nominal in the AN): for one system and
// jet-pT window, the muon ptRel of all muon-tagged jets in each of the four CSVv2 categories (< 0.54, 0.54-0.84,
// 0.84-0.96, > 0.96), data against the fitted b, c and light templates, with data / fit below.
// The fit is p2dhf::fitHF of fitPtRelCSVv2_2D_ppPbPb_PF.C (HistFactory, MC tagger shapes, c / light at MC), the same
// call as systematics_CSVv2_2D_PF.C, so the purity printed here is the chain's nominal purity. The drawn model is
// N_data x [f_b b_k(j) / B + (1 - f_b) (c_k(j) + l_k(j)) / (C + L)] with the fitted f_b (of the fitted, in PbPb
// fake-subtracted, jets) and the Barlow-Beeston gammas at 1; c and light keep their MC ratio. PbPb data are the
// fake-subtracted data the fit sees.
// Outputs: figures/bPurity/PFJets_PbPb/csvV2/fit2D/fit2D_<sys>_<lo>to<hi>.pdf for every system (pp, C4..C1) and window
// 80-300 GeV. Prints f_b, the purity (b / all muon-tagged jets) and the Poisson deviance / bins per window.
// Usage, from src/plots/bPurity/:  root -l -b -q 'plotFit2D_HF_PF.C'
#include "fitPtRelCSVv2_2D_ppPbPb_PF.C"

void plotFit2D_HF_PF(TString cats = "4cat"){
  using namespace pbc; using namespace p2d;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  setvbuf(stdout, nullptr, _IOLBF, 0);
  if(!setCats(cats)){ printf("ERROR: categories 4cat or 8cat\n"); return; }
  nFit = nPR;
  setJetCollection("PFJets"); initPlotStyle();
  TFile *fDpp = TFile::Open(Form("%s/%s", repo, dataName)), *fMpp = TFile::Open(Form("%s/%s", repo, mcName)), *fDpb = TFile::Open(Form("%s/%s", repo, dataPbPb));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  if(!fDpp || !fMpp || !fDpb || !fS[0] || !fS[1] || !fS[2]){ printf("ERROR: inputs\n"); return; }
  const TString outDir = Form("%s/figures/bPurity/PFJets_PbPb/csvV2/fit2D%s/", repo, ptBinTag()); gSystem->mkdir(outDir, true);
  const int NCOL = 5; const char *colName[NCOL] = {"pp", "PbPb 50-80%", "PbPb 30-50%", "PbPb 10-30%", "PbPb 0-10%"}, *colTag[NCOL] = {"pp", "C4", "C3", "C2", "C1"};
  const int colCls[NCOL] = {-1, 3, 2, 1, 0};
  const int colF[NF] = {kRed - 4, kGreen + 2, kBlue - 4};   // b, c, light, as the CSVv2-cut fit figures
  int uid = 0;
  for(int c = 0; c < NCOL; c++){
    if(c == 0 && !buildWindows(fDpp, fMpp)) return;
    printf("\n==== %s\n  window    f_b(fit)        purity    deviance / bins\n", colName[c]);
    for(int w = 0; w < NW; w++){
      Win X;
      if(c == 0) X = W[w]; else { double nEffB = 0.; if(!buildCls(fDpb, colCls[c], w, X, nEffB)) return; }
      const p2dhf::HFResult r = p2dhf::fitHF(X, X.d, uid++);
      double nD = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) nD += X.d[k][j];
      const double kAll = X.nData > 0. ? nD/X.nData : 1., fb = r.pur;
      const double tB = X.t[0].tot, tN = X.t[1].tot + X.t[2].tot;
      double comp[NF][KMAX][NPR] = {}, dev = 0.; int nb = 0;
      for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){
        comp[0][k][j] = nD*fb*X.t[0].n[k][j]/tB;
        for(int F = 1; F < NF; F++) comp[F][k][j] = nD*(1. - fb)*X.t[F].n[k][j]/tN;
        const double mu = comp[0][k][j] + comp[1][k][j] + comp[2][k][j], d = X.d[k][j];
        if(mu > 0.){ dev += 2.*(mu - d + (d > 0. ? d*log(d/mu) : 0.)); nb++; }
      }
      printf("  %3.0f-%3.0f   %.3f +- %.3f   %.3f    %6.1f / %d\n", wLo[w], wHi[w], fb, r.epur, fb*kAll, dev, nb);

      TCanvas *cv = new TCanvas(Form("cv%d_%d", c, w), "", 400*K + 200, 900);
      GridPads g = makeGridPads(cv, K, 2, 0.08, 0.02, 0.12, 0.09);
      const double ratioH = 0.30, inner = 1. - 0.12 - 0.09, yS = 0.12 + ratioH*inner;
      // per-category y range, data drawn in jets per bin
      for(int k = 0; k < K; k++){
        TPad *pt = g.at(k, 0), *pb = g.at(k, 1);
        const double x1 = pt->GetXlowNDC(), x2 = x1 + pt->GetWNDC();
        pt->SetPad(x1, yS, x2, 1.); pt->SetTopMargin(0.09/(1. - yS)); pt->SetBottomMargin(0.);
        pb->SetPad(x1, 0., x2, yS); pb->SetBottomMargin(0.12/yS); pb->SetTopMargin(0.);
        TH1D *hd = new TH1D(Form("d%d_%d_%d", c, w, k), "", nPR, 0., nPR*ptRelBW), *ht = (TH1D*) hd->Clone(Form("t%d_%d_%d", c, w, k)), *hc[NF];
        hd->SetDirectory(nullptr); ht->SetDirectory(nullptr);
        for(int F = 0; F < NF; F++){ hc[F] = (TH1D*) hd->Clone(Form("c%d_%d_%d_%d", c, w, k, F)); hc[F]->SetDirectory(nullptr); }
        double ymax = 0.;
        for(int j = 0; j < nPR; j++){
          hd->SetBinContent(j + 1, X.d[k][j]); hd->SetBinError(j + 1, sqrt(TMath::Max(X.d[k][j], 0.)));
          double mu = 0.; for(int F = 0; F < NF; F++){ hc[F]->SetBinContent(j + 1, comp[F][k][j]); mu += comp[F][k][j]; }
          ht->SetBinContent(j + 1, mu);
          ymax = TMath::Max(ymax, TMath::Max(X.d[k][j] + sqrt(TMath::Max(X.d[k][j], 0.)), mu));
        }
        pt->cd();
        TH1F *fr = gPad->DrawFrame(0., 0., 5., 1.45*ymax); g.styleFrame(fr, k, 0, 0.05); fr->GetXaxis()->SetLabelSize(0.);
        fr->GetYaxis()->SetLabelSize(0.05); fr->GetYaxis()->SetMaxDigits(3);
        THStack *st = new THStack(Form("st%d_%d_%d", c, w, k), "");
        for(int F = NF - 1; F >= 0; F--){ hc[F]->SetFillColor(colF[F]); hc[F]->SetLineColor(colF[F]); st->Add(hc[F]); }
        st->Draw("HIST SAME");
        styleH(hd, hexData, markFilledCircle, 1.4); hd->SetLineWidth(3); hd->Draw("E1 X0 SAME");
        gPad->RedrawAxis();
        g.text(k, 0, 0.11, k == 0 ? Form("CSVv2 < %.2f", catUp[0]) : (k == K - 1 ? Form("CSVv2 > %.2f", catUp[k-1]) : Form("CSVv2 %.2f-%.2f", catUp[k-1], catUp[k])), 0.05);
        if(k == 0){
          TLegend *lg = makeLegend(0.60, 0.36, 0.97, 0.72, 0.052);
          lg->AddEntry(hd, "data", "pe"); for(int F = 0; F < NF; F++) lg->AddEntry(hc[F], Form("%s jets", flLabel[F]), "f"); lg->Draw();
        }
        if(k == K - 1){ TLatex a; a.SetNDC(); a.SetTextFont(42); a.SetTextSize(0.062); a.SetTextAlign(33);
          a.DrawLatex(0.95, 0.80, Form("b purity %.3f #pm %.3f", fb*kAll, r.epur*kAll)); }
        pb->cd();
        TH1F *fr2 = gPad->DrawFrame(0., 0., 5., 2.); g.styleFrame(fr2, k, 1, 0.11); fr2->GetYaxis()->SetNdivisions(504);
        fr2->GetXaxis()->SetLabelSize(0.11); fr2->GetYaxis()->SetLabelSize(k == 0 ? 0.11 : 0.); fr2->GetXaxis()->SetTickLength(0.07);
        TLine *u = new TLine(0., 1., 5., 1.); u->SetLineStyle(2); u->SetLineWidth(2); u->Draw();
        TH1D *rr = (TH1D*) hd->Clone(Form("r%d_%d_%d", c, w, k)); rr->Divide(ht);
        for(int j = 1; j <= nPR; j++) rr->SetBinError(j, ht->GetBinContent(j) > 0. ? hd->GetBinError(j)/ht->GetBinContent(j) : 0.);
        rr->Draw("E1 X0 SAME");
      }
      cv->cd(); TLatex a; a.SetNDC(); a.SetTextFont(42); a.SetTextSize(0.045);
      a.SetTextAlign(31); a.DrawLatex(0.98, 0.01, "muon #it{p}_{T}^{rel} (GeV)");
      a.SetTextAngle(90); a.DrawLatex(0.022, 0.98, Form("jets per %g GeV", ptRelBW)); a.DrawLatex(0.022, yS - 0.01, "data / fit");
      g.cmsLabel(cv, "Internal", Form("%s,  %.0f < #it{p}_{T}^{jet} < %.0f GeV", colName[c], wLo[w], wHi[w]));
      savePdfTight(cv, outDir + Form("fit2D_%s_%.0fto%.0f.pdf", colTag[c], wLo[w], wHi[w]));
      delete cv;
    }
  }
}
