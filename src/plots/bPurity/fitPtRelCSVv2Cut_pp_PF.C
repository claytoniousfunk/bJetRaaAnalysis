// b purity from a ptRel fit AFTER a CSVv2 cut, pp PF jets (2026-10-24, user: "rerun the purity
// measurements on CSVv2"). Copy of fitPtRelDeepCSVCut_pp_PF.C (2026-10-22) with the tagger swapped for
// CSVv2. Same inputs, templates, axis smearing and likelihood as fitPtRelCSVv2_2D_pp_PF.C (included),
// at the CMS loose / medium / tight working points (0.5426 / 0.8484 / 0.9535, on the scan's 0.02
// CSVv2 grid: 0.54 / 0.84 / 0.96) and at 0.90, the working point of the CMS b-jet shapes analyses
// (AN-2020/029).
//
// Per jet-pT window and working point:
//   P_tag    b fraction of the TAGGED jets: ptRel fit of the jets above the cut, b vs light + c at
//            PYTHIA's c share among the tagged jets, templates from the tagged PYTHIA jets
//   eps_b    PYTHIA b-tag efficiency of muon-tagged b jets (tagged / all, ptRel 0-5 GeV)
//   f_pass   tagged / all muon-tagged jets in data
//   P_all^corr = P_tag f_pass / eps_b   the purity of ALL muon-tagged jets the cut implies, i.e. the
//            b-jet count N_tag P_tag / eps_b the chain would use, over all jets
//   compared with the ptRel-only fit of all jets, P_all; their ratio P_all^corr / P_all is the
//   data / PYTHIA b-tag efficiency the two fits together imply (if both fits are right).
// Asimov closure (PYTHIA pseudo-data): as is, and with the tagger mis-modelled (b tagged x 0.85,
// c tagged x 1.2, light tagged x 1.5), against the true tagged purity and the true purity of all jets.
//
// Output: rootFiles/bPurityResults/ptRelCSVv2Cut_PF_pp_fit<max>.root
// Figures: figures/bPurity/PFJets_pp/csvV2_2D/{cutPurity,cutFits_<wp>,cutPurity_wp0p90_summary}_fit<max>.pdf
//   (the summary compares with the chain's nominal pp purity, rootFiles/systematics/bPuritySys_total_PFJets_pp.root)
// Usage, from src/plots/bPurity/:  root -l -b -q 'fitPtRelCSVv2Cut_pp_PF.C(5.0)'   (0-5 GeV nominal since 2026-10-24, user; 3.0 = the chain's edge)
#include "fitPtRelCSVv2_2D_pp_PF.C"

namespace pcut {
using namespace p2d;
const int NWP = 4;
const double wpCut[NWP] = {0.54, 0.84, 0.90, 0.96};
const char *wpName[NWP] = {"loose", "medium", "shapes", "tight"};   // "shapes" = the CMS b-jet shapes analyses
const int wpSum = 2;   // the summary figure: CSVv2 > 0.90
// the tagged jets (category 1 of the fail / pass split) as a one-category window
Win passWin(const Win &X){
  Win P; P.sigma = X.sigma;
  for(int j = 0; j < NPR; j++){
    P.d[0][j] = X.d[1][j];
    for(int F = 0; F < NF; F++){ P.t[F].n[0][j] = X.t[F].n[1][j]; P.t[F].e2[0][j] = X.t[F].e2[1][j]; }
  }
  for(int F = 0; F < NF; F++){ for(int j = 0; j < NPR; j++) P.t[F].tot += P.t[F].n[0][j]; P.t[F].eps[0] = 1.; }
  for(int j = 0; j < NPR; j++) P.nData += P.d[0][j];
  P.truthB = P.t[0].tot/(P.t[0].tot + P.t[1].tot + P.t[2].tot);
  P.rc = P.t[1].tot/(P.t[1].tot + P.t[2].tot);
  return P;
}
} // namespace pcut

void fitPtRelCSVv2Cut_pp_PF(double fitMax = 5.0){
  using namespace pcut;
  nFit = (int)(fitMax/0.5 + 0.5);
  setJetCollection("PFJets"); initPlotStyle();
  TFile *fD = TFile::Open(Form("%s/%s", repo, dataName)), *fM = TFile::Open(Form("%s/%s", repo, mcName));
  if(!fD || !fM || fD->IsZombie() || fM->IsZombie()){ printf("ERROR: inputs\n"); return; }
  const TString outDir = Form("%s/figures/bPurity/PFJets_pp/csvV2_2D/", repo), tag = Form("fit%g", fitMax);
  gSystem->mkdir(outDir, true);
  const Model &M1 = models[0];   // ptRel only, b vs light + c at PYTHIA's c share

  TH1D *hTag[NWP], *hTagMC[NWP], *hCorr[NWP], *hSF[NWP], *hAll = nullptr, *hAllMC = nullptr;
  Win Pmed[NWP][NW]; Result Rmed[NWP][NW];   // tagged windows and fits of every working point, for the fit figures
  for(int p = 0; p < NWP; p++){
    setCats(Form("cut%.2f", wpCut[p]));
    if(!buildWindows(fD, fM)) return;
    hTag[p]   = new TH1D(Form("h_purTag_%s", wpName[p]), "", NW, ptEdges + 3);
    hTagMC[p] = new TH1D(Form("h_purTagMC_%s", wpName[p]), "", NW, ptEdges + 3);
    hCorr[p]  = new TH1D(Form("h_purCorr_%s", wpName[p]), "", NW, ptEdges + 3);
    hSF[p]    = new TH1D(Form("h_sfB_%s", wpName[p]), "", NW, ptEdges + 3);
    if(!hAll){ hAll = new TH1D("h_purAll", "", NW, ptEdges + 3); hAllMC = new TH1D("h_purAllMC", "", NW, ptEdges + 3); hAll->SetDirectory(nullptr); hAllMC->SetDirectory(nullptr); }
    for(TH1D *h : {hTag[p], hTagMC[p], hCorr[p], hSF[p]}) h->SetDirectory(nullptr);
    printf("\n==== CSVv2 > %.2f (%s)  fit 0-%.1f GeV\n", wpCut[p], wpName[p], fitMax);
    printf("  window    | f_pass  eps_b(MC) eps_c eps_l | P_tag fit       P_tag PYTHIA | P_all fit (ptRel only) | P_all^corr      | implied SF_b  | dev/ndf tag\n");
    for(int w = 0; w < NW; w++){
      Win &X = W[w]; Win P = passWin(X);
      Result rA = fitWin(M1, X, X.d), rT = fitWin(M1, P, P.d);
      const double fPass = P.nData/X.nData, eB = X.t[0].eps[1];
      const double corr = rT.pur*fPass/eB, ecorr = rT.epur*fPass/eB, sf = corr/rA.pur;
      const double esf = sf*sqrt(pow(rT.epur/rT.pur, 2) + pow(rA.epur/rA.pur, 2));   // the two fits share jets: indicative only
      printf("  %3.0f-%3.0f   | %.3f   %.3f    %.3f %.3f | %.3f +- %.3f    %.3f       | %.3f +- %.3f          | %.3f +- %.3f   | %.3f        | %.1f/%d\n",
             wLo[w], wHi[w], fPass, eB, X.t[1].eps[1], X.t[2].eps[1], rT.pur, rT.epur, P.truthB, rA.pur, rA.epur, corr, ecorr, sf, rT.dev, rT.ndf);
      hTag[p]->SetBinContent(w + 1, rT.pur); hTag[p]->SetBinError(w + 1, rT.epur); hTagMC[p]->SetBinContent(w + 1, P.truthB);
      hCorr[p]->SetBinContent(w + 1, corr); hCorr[p]->SetBinError(w + 1, ecorr);
      hSF[p]->SetBinContent(w + 1, sf); hSF[p]->SetBinError(w + 1, ecorr/rA.pur);
      if(p == 0){ hAll->SetBinContent(w + 1, rA.pur); hAll->SetBinError(w + 1, rA.epur); hAllMC->SetBinContent(w + 1, X.truthB); }
      Pmed[p][w] = P; Rmed[p][w] = rT;
    }
    // Asimov closure
    printf("  Asimov closure (fit - truth): tagged purity | corrected all-jet purity\n");
    for(int s = 0; s < 2; s++){
      double mult[NF][KMAX];
      for(int F = 0; F < NF; F++){ mult[F][0] = 1.; mult[F][1] = s == 0 ? 1. : (F == 0 ? 0.85 : F == 1 ? 1.2 : 1.5); }
      printf("    %s\n", s == 0 ? "PYTHIA as is" : "tagger mis-modelled (b tagged x0.85, c x1.2, light x1.5)");
      for(int w = 0; w < NW; w++){
        Win X = W[w]; double pd[KMAX][NPR];
        const double trAll = asimov(X, mult, pd);
        for(int k = 0; k < K; k++) for(int j = 0; j < NPR; j++) X.d[k][j] = pd[k][j];
        X.nData = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < NPR; j++) X.nData += X.d[k][j];
        Win P = passWin(X);
        double tb = 0., tt = 0.; for(int j = 0; j < NPR; j++){ tb += X.t[0].n[1][j]*mult[0][1]; for(int F = 0; F < NF; F++) tt += X.t[F].n[1][j]*mult[F][1]; }
        Result rT = fitWin(M1, P, P.d);
        const double corr = rT.pur*(P.nData/X.nData)/X.t[0].eps[1];
        printf("      %3.0f-%3.0f  %+.3f (truth %.3f) | %+.3f (truth %.3f)\n", wLo[w], wHi[w], rT.pur - tb/tt, tb/tt, corr - trAll, trAll);
      }
    }
  }

  {
    TFile *fo = TFile::Open(Form("%s/rootFiles/bPurityResults/ptRelCSVv2Cut_PF_pp_%s.root", repo, tag.Data()), "RECREATE");
    for(int p = 0; p < NWP; p++){ hTag[p]->Write(); hTagMC[p]->Write(); hCorr[p]->Write(); hSF[p]->Write(); }
    hAll->Write(); hAllMC->Write();
    TNamed info("info", Form("pp PF ptRel fit after a CSVv2 cut (%.2f / %.2f / %.2f / %.2f), fit 0-%.1f GeV. data: %s. templates: %s",
                             wpCut[0], wpCut[1], wpCut[2], wpCut[3], fitMax, dataName, mcName)); info.Write();
    fo->Close();
  }

  // ---- figure 1: tagged purity (left) and corrected all-jet purity vs the ptRel-only fit (right)
  const char *hexWP[NWP] = {okabeHex[5], okabeHex[1], okabeHex[2], okabeHex[6]};
  const int markWP[NWP] = {markFilledSquare, markFilledDiamond, markFilledCircle, 34};
  {
    TCanvas *cv = new TCanvas("cvC", "", 1600, 900);
    GridPads g = makeGridPads(cv, 2, 1, 0.07, 0.02, 0.12, 0.08);
    g.at(0, 0)->cd();
    TH1F *fr = gPad->DrawFrame(80., 0., 300., 1.25); g.styleFrame(fr, 0, 0, 0.05);
    TLegend *lg = makeLegend(0.20, 0.13, 0.95, 0.45, 0.048);
    for(int p = 0; p < NWP; p++){
      styleLine(hTagMC[p], hexWP[p], 3); hTagMC[p]->SetLineStyle(7); hTagMC[p]->Draw("HIST SAME");
      styleH(hTag[p], hexWP[p], markWP[p], 2.4); hTag[p]->SetLineWidth(3); hTag[p]->Draw("E1 X0 SAME");
      lg->AddEntry(hTag[p], Form("CSVv2 > %.2f (%s)", wpCut[p], wpName[p]), "pe");
    }
    TH1D *dash = (TH1D*) hTagMC[0]->Clone("dashLeg"); dash->SetLineColor(kGray + 2); lg->AddEntry(dash, "PYTHIA", "l");
    lg->Draw();
    g.text(0, 0, 0.07, "jets passing the CSVv2 cut", 0.05, 62);
    g.at(1, 0)->cd();
    TH1F *fr2 = gPad->DrawFrame(80., 0., 300., 1.25); g.styleFrame(fr2, 1, 0, 0.05);
    styleLine(hAllMC, "#8C8C8C", 4); hAllMC->SetLineStyle(7); hAllMC->Draw("HIST SAME");
    styleH(hAll, hexData, markFilledCircle, 2.4); hAll->SetLineWidth(3); hAll->Draw("E1 X0 SAME");
    TLegend *lg2 = makeLegend(0.05, 0.56, 0.99, 0.82, 0.044); lg2->SetNColumns(2);
    lg2->AddEntry(hAll, "#it{p}_{T}^{rel} fit, all jets", "pe");
    for(int p = 0; p < NWP; p++){
      styleH(hCorr[p], hexWP[p], markWP[p], 2.4); hCorr[p]->SetLineWidth(3); hCorr[p]->Draw("E1 X0 SAME");
      lg2->AddEntry(hCorr[p], Form("%s cut / #varepsilon_{b}^{PYTHIA}", wpName[p]), "pe");
    }
    lg2->AddEntry(hAllMC, "PYTHIA", "l"); lg2->Draw();
    g.text(1, 0, 0.07, "all muon-tagged jets", 0.05, 62);
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "b purity", 0.045, 0.02, 0.02);
    g.cmsLabel(cv, "Internal", Form("pp (5.02 TeV),  fit #it{p}_{T}^{rel} < %.0f GeV", fitMax));
    savePdfTight(cv, outDir + "cutPurity_" + tag + ".pdf");
  }
  // ---- figure 2: the ptRel fits of the tagged jets, every window, one figure per working point;
  //      fit / data below. b kRed-4, c kGreen+2, light kBlue-4
  const int colF[NF] = {kRed - 4, kGreen + 2, kBlue - 4};   // b, c, light (user, 2026-10-22)
  for(int p = 0; p < NWP; p++){
    TCanvas *cv = new TCanvas(Form("cvM%d", p), "", 2200, 950);
    const double lm = 0.06, rm = 0.02, bm = 0.13, tm = 0.09;
    GridPads g = makeGridPads(cv, NW, 2, lm, rm, bm, tm);
    const double inner = 1. - bm - tm, yS = bm + 0.30*inner;   // ratio row: 30% of the frame height
    for(int w = 0; w < NW; w++){
      cur = &M1; cw = &Pmed[p][w]; cd = Pmed[p][w].d; Kf = 1;
      double mu[KMAX][NPR], va[KMAX][NPR], comp[NF][KMAX][NPR];
      const int nSave = nFit; expectAll(Rmed[p][w].par, mu, va, comp);
      double muIn = 0.; for(int j = 0; j < nSave; j++) muIn += mu[0][j];
      nFit = NPR; expectAll(Rmed[p][w].par, mu, va, comp);
      double muIn2 = 0.; for(int j = 0; j < nSave; j++) muIn2 += mu[0][j];
      nFit = nSave; const double rs = muIn/muIn2;   // normalised to the data inside the fit range, drawn over 0-5 GeV
      TH1D *hd = new TH1D(Form("md%d_%d", p, w), "", NPR, 0., 5.), *ht = new TH1D(Form("mt%d_%d", p, w), "", NPR, 0., 5.), *hc[NF];
      hd->SetDirectory(nullptr); ht->SetDirectory(nullptr);
      for(int F = 0; F < NF; F++){ hc[F] = new TH1D(Form("mc%d_%d_%d", p, w, F), "", NPR, 0., 5.); hc[F]->SetDirectory(nullptr); }
      double ymax = 0.;
      for(int j = 0; j < NPR; j++){
        const double d = Pmed[p][w].d[0][j];
        hd->SetBinContent(j + 1, 1e-3*d); hd->SetBinError(j + 1, 1e-3*sqrt(d));
        for(int F = 0; F < NF; F++) hc[F]->SetBinContent(j + 1, 1e-3*comp[F][0][j]*rs);
        ht->SetBinContent(j + 1, 1e-3*mu[0][j]*rs);
        ymax = TMath::Max(ymax, 1e-3*TMath::Max(d, mu[0][j]*rs));
      }
      TPad *pt = g.at(w, 0), *pb = g.at(w, 1);
      const double x1 = pt->GetXlowNDC(), x2 = x1 + pt->GetWNDC();
      pt->SetPad(x1, yS, x2, 1.); pt->SetTopMargin(tm/(1. - yS)); pt->SetBottomMargin(0.);
      pb->SetPad(x1, 0., x2, yS); pb->SetBottomMargin(bm/yS); pb->SetTopMargin(0.);
      pt->cd();
      TH1F *fr = gPad->DrawFrame(0., 0., 5., 1.35*ymax); g.styleFrame(fr, w, 0, 0.05);
      fr->GetXaxis()->SetLabelSize(0.); fr->GetYaxis()->SetLabelSize(w == 0 ? 0.055 : 0.);
      THStack *st = new THStack(Form("mst%d_%d", p, w), "");
      for(int F = NF - 1; F >= 0; F--){ hc[F]->SetFillColor(colF[F]); hc[F]->SetLineColor(colF[F]); st->Add(hc[F]); }
      st->Draw("HIST SAME");
      styleH(hd, hexData, markFilledCircle, 1.8); hd->SetLineWidth(3); hd->Draw("E1 X0 SAME");
      TLine *edge = new TLine(fitMax, 0., fitMax, 1.35*ymax); edge->SetLineStyle(7); edge->SetLineWidth(3); edge->SetLineColor(kGray + 2); if(fitMax < 5.) edge->Draw();
      gPad->RedrawAxis();
      g.text(w, 0, 0.10, Form("%.0f-%.0f GeV", wLo[w], wHi[w]), 0.060);
      g.text(w, 0, 0.19, Form("b purity %.2f", Rmed[p][w].pur), 0.055);
      if(w == NW - 1){
        TLegend *lg = makeLegend(0.50, 0.36, 0.97, 0.72, 0.060);
        lg->AddEntry(hd, "data", "pe"); for(int F = 0; F < NF; F++) lg->AddEntry(hc[F], Form("%s jets", flLabel[F]), "f"); lg->Draw();
      }
      pb->cd();
      TH1F *fr2 = gPad->DrawFrame(0., 0.6, 5., 1.4); g.styleFrame(fr2, w, 1, 0.12);
      fr2->GetYaxis()->SetNdivisions(202); fr2->GetXaxis()->SetTickLength(0.07);
      fr2->GetYaxis()->SetLabelSize(w == 0 ? 0.12 : 0.); fr2->GetXaxis()->SetLabelSize(0.12);
      TLine *u = new TLine(0., 1., 5., 1.); u->SetLineStyle(2); u->SetLineWidth(2); u->Draw();
      if(fitMax < 5.){ TLine *e2 = new TLine(fitMax, 0.6, fitMax, 1.4); e2->SetLineStyle(7); e2->SetLineWidth(3); e2->SetLineColor(kGray + 2); e2->Draw(); }
      TH1D *r = (TH1D*) ht->Clone(Form("mr%d_%d", p, w)); r->SetDirectory(nullptr);   // fit / data, data statistical error
      for(int j = 1; j <= NPR; j++){
        const double d = hd->GetBinContent(j);
        if(d > 0.){ r->SetBinContent(j, ht->GetBinContent(j)/d); r->SetBinError(j, ht->GetBinContent(j)/d*hd->GetBinError(j)/d); }
        else { r->SetBinContent(j, -999.); r->SetBinError(j, 0.); }
      }
      styleH(r, hexData, markFilledCircle, 1.6); r->SetLineWidth(3); r->Draw("E1 X0 SAME");
    }
    cv->cd(); TLatex a; a.SetNDC(); a.SetTextFont(42); a.SetTextSize(0.045);
    a.SetTextAlign(31); a.DrawLatex(1. - rm, 0.012, "muon #it{p}_{T}^{rel} (GeV)");
    a.SetTextAngle(90); a.DrawLatex(0.016, 1. - tm, "10^{3} jets per 0.5 GeV"); a.DrawLatex(0.016, yS - 0.01, "fit / data");
    g.cmsLabel(cv, "Internal", Form("pp (5.02 TeV),  CSVv2 > %.2f (%s)", wpCut[p], wpName[p]), 1.3);
    savePdfTight(cv, outDir + "cutFits_" + wpName[p] + "_" + tag + ".pdf");
  }
  // ---- figure 3: summary at the 0.90 (shapes) working point, against the nominal ptRel-only pp purity of the
  //      chain (bPuritySys_total_PFJets_pp.root: 2026-10-06 data, February PYTHIA, 0-3 GeV fit, total systematic)
  {
    TFile *fN = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_total_PFJets_pp.root", repo));
    TH1D *hNomSrc = fN ? (TH1D*) fN->Get("h_bPurity_sysAbs") : nullptr;
    TH1D *hNomStatSrc = fN ? (TH1D*) fN->Get("h_bPurity_nominal") : nullptr;
    if(!hNomSrc || !hNomStatSrc){ printf("ERROR: nominal pp purity file\n"); return; }
    TH1D *hNom = new TH1D("h_purNominal", "", NW, ptEdges + 3), *hNomSys = new TH1D("h_purNominalSys", "", NW, ptEdges + 3);
    hNom->SetDirectory(nullptr); hNomSys->SetDirectory(nullptr);
    for(int w = 0; w < NW; w++){
      const int b = hNomSrc->FindBin(0.5*(wLo[w] + wHi[w]));
      hNom->SetBinContent(w + 1, hNomStatSrc->GetBinContent(b)); hNom->SetBinError(w + 1, hNomStatSrc->GetBinError(b));
      hNomSys->SetBinContent(w + 1, hNomSrc->GetBinContent(b)); hNomSys->SetBinError(w + 1, hNomSrc->GetBinError(b));
    }
    printf("\nSUMMARY, CSVv2 > %.2f, fit 0-%.1f GeV\n  window    tagged     all (cut / eps_b)   nominal chain (stat, total sys)\n", wpCut[wpSum], fitMax);
    for(int w = 0; w < NW; w++) printf("  %3.0f-%3.0f   %.3f +- %.3f   %.3f +- %.3f      %.3f +- %.3f +- %.3f\n", wLo[w], wHi[w],
                                         hTag[wpSum]->GetBinContent(w + 1), hTag[wpSum]->GetBinError(w + 1), hCorr[wpSum]->GetBinContent(w + 1), hCorr[wpSum]->GetBinError(w + 1),
                                         hNom->GetBinContent(w + 1), hNom->GetBinError(w + 1), hNomSys->GetBinError(w + 1));
    TCanvas *cv = new TCanvas("cvSum", "", 1000, 900);
    GridPads g = makeGridPads(cv, 1, 1, 0.14, 0.04, 0.12, 0.07);
    g.at(0, 0)->cd();
    TH1F *fr = gPad->DrawFrame(80., 0., 300., 1.6); g.styleFrame(fr, 0, 0, 0.045);
    hNomSys->SetFillColorAlpha(kGray + 1, 0.45); hNomSys->SetLineWidth(0); hNomSys->SetMarkerSize(0); hNomSys->Draw("E2 SAME");
    // PYTHIA truth, dashed: b fraction of the tagged jets (orange) and of all muon-tagged jets (grey)
    TH1D *mcTag = (TH1D*) hTagMC[wpSum]->Clone("sumTagMC"), *mcAll = (TH1D*) hAllMC->Clone("sumAllMC");
    styleLine(mcTag, okabeHex[1], 4); mcTag->SetLineStyle(7); mcTag->Draw("HIST SAME");
    styleLine(mcAll, "#8C8C8C", 4); mcAll->SetLineStyle(7); mcAll->Draw("HIST SAME");
    styleH(hNom, hexData, markFilledCircle, 2.4); hNom->SetLineWidth(3); hNom->Draw("E1 X0 SAME");
    TH1D *tg = (TH1D*) hTag[wpSum]->Clone("sumTag"), *co = (TH1D*) hCorr[wpSum]->Clone("sumCorr");
    styleH(tg, okabeHex[1], markFilledDiamond, 2.8); tg->SetLineWidth(3); tg->Draw("E1 X0 SAME");
    styleH(co, okabeHex[5], markFilledSquare, 2.4); co->SetLineWidth(3); co->Draw("E1 X0 SAME");
    TLegend *lg = makeLegend(0.20, 0.70, 0.62, 0.90, 0.040);
    lg->AddEntry(tg, Form("CSVv2 > %.2f", wpCut[wpSum]), "pe");
    lg->AddEntry(co, Form("CSVv2 > %.2f, / #varepsilon_{b}", wpCut[wpSum]), "pe");
    lg->AddEntry(hNom, "nominal", "pe");
    lg->AddEntry(mcAll, "PYTHIA", "l");
    lg->Draw();
    g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "b purity", 0.045, 0.03, 0.02);
    g.cmsLabel(cv, "Internal", Form("pp (5.02 TeV),  CSVv2 fits 0-%.0f GeV", fitMax));
    savePdfTight(cv, outDir + "cutPurity_wp0p90_summary_" + tag + ".pdf");
  }
}
