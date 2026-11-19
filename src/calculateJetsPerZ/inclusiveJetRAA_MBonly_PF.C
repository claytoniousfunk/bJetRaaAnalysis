// Inclusive-jet R_AA from MINBIAS-TRIGGERED JETS ONLY, pp and PbPb, in the ordinary
// T_AA normalization, as a cross-check against the published CMS jet R_AA
// (2026-11-19, user: "using the MinBias-only jets from PbPb and pp, make a jet RAA
// that we can cross-check with the published from CMS").
//
//   R_AA = [ N_jet^PbPb / N_evt^PbPb ] / ( <T_AA> sigma_jet^pp )
//
// No jet triggers on either side and no N_Z: this tests the jet side of R_Z with
// the standard normalization alone.
//   PbPb  per-MinBias-event yield in each class (N_evt = the class's h_vz count in the
//         calo-matched MinBias Part1 scan), forest-CS fakes subtracted per event,
//         unfolded with the chain's calo-matched P+H response at the chain's N.
//   pp    MinBias jets, the MinBias count converted to luminosity with the Jet100-bit
//         event count (k_MB = N_evt(HighEGJet Jet100) / N_evt(MinBias & Jet100 bit),
//         the nominal stitch's MinBias normalization), divided by L_pp = 299.49 pb^-1;
//         unfolded with the chain's pp response.
//   T_AA  <Ncoll> of the 5% slices (CMS 2018 Glauber) / 67.6 mb, as inclusiveJetRAA_PF.C.
// Everything else (responses, JER scheme, reco floor 60 GeV, fakes) comes from
// inclusiveJetsPerZSys_PF.C in its MBcaloMatched mode with the event-count stitch and
// mbOnly = true. N_iter = the current nominal inclusive N (argument nItAll > 0: that N everywhere,
// outputs _nIter<N>; pp 2, 0-10% 6, 10-30% 5,
// 30-50% 3, 50-80% 2; 2026-11-19 unsmeared-JER chain).
//
// Uncertainties drawn:
//   bars   statistical (data only)
//   boxes  the jet-side systematic of the NOMINAL inclusive R_Z, BORROWED (h_sysRel_total of the
//          deck file: JEU, JER, unfolding, fakes, ...; no stitch, no N_Z sources). Not re-derived
//          for the MinBias-only chain.
//   box at 1  normalization: T_AA (Glauber, 1.8 / 2.2 / 3.3 / 5.0%), L_pp 1.9%, pp k_MB stat.
//
// Published: CMS PbPb 2.76 TeV, R = 0.4, |eta| < 2 (PRC 96 (2017) 015202, HEPData ins1487278,
// headers/publishedData/CMS_inclJetRAA_2p76TeV_PRC96.h), the only CMS R = 0.4 R_AA covering
// 80-300 GeV (the 5.02 TeV JHEP 05 (2021) 284 starts at 200 GeV, beyond the MinBias reach).
// Classes mapped to ours with <Ncoll> weights of our 5% slices (approximate):
//   0-10% = 0-5% (+) 5-10%,  10-30%, 30-50% as published,  50-80% = 50-70% (+) 70-90% x (70-80 / 70-90).
// Systematics of the merged CMS classes added linearly (correlated), stat. in quadrature.
// Ours is |eta| < 1.6 at 5.02 TeV: a different energy, so this is a sanity check of size
// and shape, not an equality test.
//
// Usage (from src/calculateJetsPerZ/):
//   root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");' inclusiveJetRAA_MBonly_PF.C+
// Output: rootFiles/systematics/inclJetRAA_MBonly_PF.root  (h_RAA_C<c> stat. errors,
//           h_sysRel_C<c> borrowed, h_xsecPP_pb, h_yieldAA_C<c>, h_cms276_C<c>)
//         figures/JetsPerZ/RAA_MBonly/inclJetRAA_MBonly_vsCMS2p76TeV.pdf
#include "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateJetsPerZ/inclusiveJetsPerZSys_PF.C"
#include "/home/clayton/Analysis/code/bJetRaaAnalysis/headers/publishedData/CMS_inclJetRAA_2p76TeV_PRC96.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"

namespace mbRAA {
const double taaRel[4] = {28./1558., 17.4/783.9, 8.9/267.0, 2.52/50.37};   // Glauber <Ncoll> per class (inclusiveJetRAA_PF.C)
const double lumiPPRel = 13./669.;
const double dEta = 3.2;                                                   // |eta| < 1.6
const int NB = 11;
const double edges[NB+1] = {80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300};   // the CMS 2.76 TeV bins from 80 GeV
double meanNcoll(int ci, int sLo = -1, int sHi = -1){
  if(sLo < 0){ sLo = coarseSliceLo[ci]; sHi = coarseSliceHi[ci]; }
  double s = 0.; int n = 0;
  for(int sl = sLo; sl <= sHi; sl++){ s += pubZ::ncollSlice[sl-1]; n++; }
  return s/n;
}
}

void inclusiveJetRAA_MBonly_PF(int nItAll = 0){   // > 0: this N for pp and every class, outputs _nIter<N>
  using namespace mbRAA;
  const char *repo = "/home/clayton/Analysis/code/bJetRaaAnalysis";
  modeMB = true; nomFloor = 60.; centMode = ""; trigAA = "plain"; prescaleK80 = true; evtCountStitch = true;
  if(!loadInputs()) return;
  int nIt[5] = {2, 6, 5, 3, 2};
  if(nItAll > 0) for(int s = 0; s < 5; s++) nIt[s] = nItAll;
  const TString nTag = nItAll > 0 ? Form("_nIter%d", nItAll) : "";
  printf("N_iter pp / 0-10 / 10-30 / 30-50 / 50-80%%: %d / %d / %d / %d / %d\n", nIt[0], nIt[1], nIt[2], nIt[3], nIt[4]);
  Opts o;
  TH1D *U[5], *Ust0 = nullptr;
  mbOnly = true;
  for(int s = 0; s < 5; s++) U[s] = unfold(resp[s], stitch(s, o, Form("MB%d", s)), nIt[s], 60., Form("U_MB%d", s));
  mbOnly = false;
  Ust0 = unfold(resp[0], stitch(0, o, "ST0"), nIt[0], 60., "U_ST0");   // pp stitched, for the MinBias-only check

  auto IE = [](TH1D *h, double a, double b, double &e){ return h->IntegralAndError(h->FindBin(a + 1e-3), h->FindBin(b - 1e-3), e); };
  TFile *fd = TFile::Open(Form("%s/rootFiles/systematics/inclJetsPerZSys_PF_MBcaloMatched_JEUppCorr_noFlatJEU_fakeGeoOnly_noRecoFloor_normSep_noStitch.root", repo));
  if(!fd){ printf("ERROR: deck file\n"); return; }

  // ---- pp cross section, pb / GeV / unit eta
  TH1D *xs = new TH1D("h_xsecPP_pb", ";jet p_{T} (GeV);d^{2}#sigma/dp_{T}d#eta (pb/GeV)", NB, edges);
  printf("\npp MinBias: k_MB = %.4g per MinBias count (+- %.2f%%), L_pp = %.2f pb^-1\n", cntKMB[0], 100.*cntKMBrel[0], pubZ::lumiPP);
  printf("  bin        sigma (pb/GeV)       MB only / stitched (pp)\n");
  for(int i = 0; i < NB; i++){
    double e, es, v = IE(U[0], edges[i], edges[i+1], e), vs = IE(Ust0, edges[i], edges[i+1], es);
    const double w = (edges[i+1] - edges[i])*dEta*pubZ::lumiPP;
    xs->SetBinContent(i+1, v/w); xs->SetBinError(i+1, e/w);
    printf("  %3.0f-%3.0f   %.4g +- %.1f%%    %.3f +- %.3f\n", edges[i], edges[i+1], v/w, 100.*e/v, v/vs, (v/vs)*e/v);
  }

  // ---- R_AA per class
  TH1D *raa[5] = {nullptr}, *sysR[5] = {nullptr}, *cms[5] = {nullptr}, *yAA[5] = {nullptr};
  double norm[5];
  for(int c = 1; c <= 4; c++){
    const int ci = c - 1;
    const double taa = meanNcoll(ci)/pubZ::sigmaNNTaa*1e-9;              // mb^-1 -> pb^-1
    norm[c] = sqrt(pow(taaRel[ci], 2) + pow(lumiPPRel, 2) + pow(cntKMBrel[0], 2));
    raa[c]  = new TH1D(Form("h_RAA_C%d", c), ";jet p_{T} (GeV);R_{AA}", NB, edges);
    sysR[c] = new TH1D(Form("h_sysRel_C%d", c), "", NB, edges);
    yAA[c]  = new TH1D(Form("h_yieldAA_C%d", c), ";jet p_{T} (GeV);1/N_{evt} d^{2}N/dp_{T}d#eta (GeV^{-1})", NB, edges);
    TH1D *hsy = nullptr; fd->GetObject(Form("h_sysRel_total_PbPbOverPP_C%d_incl", c), hsy);
    printf("\n%s: <Ncoll> %.1f, <T_AA> %.3f mb^-1, per-event via 1/k_MB = %.4g, norm. %.1f%%\n", sysLabel[c], meanNcoll(ci), meanNcoll(ci)/pubZ::sigmaNNTaa, 1./cntKMB[c], 100.*norm[c]);
    printf("  bin        R_AA (stat)        borrowed syst.\n");
    for(int i = 0; i < NB; i++){
      double e, v = IE(U[c], edges[i], edges[i+1], e);
      const double y = v/cntKMB[c]/((edges[i+1] - edges[i])*dEta);     // per MinBias event in the class
      const double ey = y*(v > 0 ? e/v : 0.);
      yAA[c]->SetBinContent(i+1, y); yAA[c]->SetBinError(i+1, ey);
      const double sp = xs->GetBinContent(i+1), esp = xs->GetBinError(i+1);
      const double r = (sp > 0 && y > 0) ? y/(taa*sp) : 0., er = r*sqrt(pow(ey/y, 2) + pow(esp/sp, 2));
      raa[c]->SetBinContent(i+1, r); raa[c]->SetBinError(i+1, r > 0 ? er : 0.);
      const double sy = hsy ? hsy->GetBinContent(hsy->FindBin(0.5*(edges[i] + edges[i+1]))) : 0.;
      sysR[c]->SetBinContent(i+1, sy);
      printf("  %3.0f-%3.0f   %.3f +- %.3f      %.1f%%\n", edges[i], edges[i+1], r, er, 100.*sy);
    }
  }

  // ---- CMS 2.76 TeV mapped onto our classes, <Ncoll>-weighted with our slices
  using namespace cmsRAA276;
  struct Piece { int k; double w; };
  const double n05 = meanNcoll(0, 1, 1), n510 = meanNcoll(0, 2, 2), n5070 = meanNcoll(0, 11, 14), n7080 = meanNcoll(0, 15, 16);
  std::vector<Piece> map[5];
  map[1] = {{0, n05}, {1, n510}};
  map[2] = {{2, 1.}};
  map[3] = {{3, 1.}};
  map[4] = {{4, 4.*n5070}, {5, 2.*n7080}};                                  // slice-count x <Ncoll>
  for(int c = 1; c <= 4; c++){
    cms[c] = new TH1D(Form("h_cms276_C%d", c), "", NB, edges);
    for(int i = 0; i < NB; i++){
      const double x = 0.5*(edges[i] + edges[i+1]);
      int b = -1; for(int j = 0; j < nMax; j++) if(x > cmsRAA276::ptEdges[j] && x < cmsRAA276::ptEdges[j+1]) b = j;
      double sw = 0., sv = 0., ss2 = 0., ssy = 0.; bool ok = b >= 0;
      for(auto &p : map[c]){ if(!ok || b >= nPt[p.k] || cmsRAA276::raa[p.k][b] <= 0.){ ok = false; break; }
        sw += p.w; sv += p.w*cmsRAA276::raa[p.k][b]; ss2 += pow(p.w*statErr[p.k][b], 2); ssy += p.w*systErr[p.k][b]; }
      if(!ok) continue;
      cms[c]->SetBinContent(i+1, sv/sw); cms[c]->SetBinError(i+1, sqrt(ss2/(sw*sw) + pow(ssy/sw, 2)));
    }
  }

  printf("\nours / CMS 2.76 TeV (ours stat. (+) borrowed syst.; CMS stat. (+) syst.; normalization not included)\n");
  printf("  bin      ");
  for(int c = 4; c >= 1; c--) printf("   %-18s", sysLabel[c]);
  printf("\n");
  for(int i = 0; i < NB; i++){
    printf("  %3.0f-%3.0f", edges[i], edges[i+1]);
    for(int c = 4; c >= 1; c--){
      const double r = raa[c]->GetBinContent(i+1), m = cms[c]->GetBinContent(i+1);
      if(r <= 0 || m <= 0){ printf("   %-18s", "-"); continue; }
      const double e = (r/m)*sqrt(pow(raa[c]->GetBinError(i+1)/r, 2) + pow(sysR[c]->GetBinContent(i+1), 2) + pow(cms[c]->GetBinError(i+1)/m, 2));
      printf("   %.3f +- %.3f     ", r/m, e);
    }
    printf("\n");
  }

  // ---- figure: 1x4 row, peripheral to central, ours / CMS below
  initPlotStyle(); gStyle->SetOptStat(0);
  const char *lab[5] = {"", "0-10%", "10-30%", "30-50%", "50-80%"};
  TCanvas *cv = new TCanvas("cv", "", 1800, 760);
  const double L = 0.08, R = 0.01, T = 0.08, B = 0.13, sp = 0.36, w = (1 - L - R)/4;
  const double xLo = 80., xHi = 239.9;   // above 240 GeV the pp MinBias-only unfolding leaves the stitched pp by 1.2-1.4
  for(int p = 0; p < 4; p++){
    const int c = 4 - p;
    const double x1 = p == 0 ? 0 : L + p*w, x2 = L + (p+1)*w + (p == 3 ? R : 0), yS = B + (1 - T - B)*sp;
    cv->cd(); TPad *t = new TPad(Form("t%d", p), "", x1, yS, x2, 1);
    t->SetLeftMargin(p == 0 ? L/(x2 - x1) : 0); t->SetRightMargin(p == 3 ? R/(x2 - x1) : 0); t->SetTopMargin(T/(1 - yS)); t->SetBottomMargin(0); t->SetTicks(1, 1); t->Draw();
    cv->cd(); TPad *bp = new TPad(Form("b%d", p), "", x1, 0, x2, yS);
    bp->SetLeftMargin(p == 0 ? L/(x2 - x1) : 0); bp->SetRightMargin(p == 3 ? R/(x2 - x1) : 0); bp->SetTopMargin(0); bp->SetBottomMargin(B/yS); bp->SetTicks(1, 1); bp->Draw();

    TGraphAsymmErrors *gs = new TGraphAsymmErrors();
    TH1D *rat = new TH1D(Form("rat%d", c), "", NB, edges);
    for(int i = 0; i < NB; i++){
      const double r = raa[c]->GetBinContent(i+1), m = cms[c]->GetBinContent(i+1);
      if(r <= 0) continue;
      const double hw = 0.5*(edges[i+1] - edges[i]);
      gs->SetPoint(gs->GetN(), edges[i] + hw, r); const double sy = r*sysR[c]->GetBinContent(i+1);
      gs->SetPointError(gs->GetN() - 1, hw, hw, sy, sy);
      if(m > 0){ rat->SetBinContent(i+1, r/m);
        rat->SetBinError(i+1, (r/m)*sqrt(pow(raa[c]->GetBinError(i+1)/r, 2) + pow(sysR[c]->GetBinContent(i+1), 2) + pow(cms[c]->GetBinError(i+1)/m, 2))); }
    }
    styleH(raa[c], "#0072B2", markFilledCircle, 1.8); raa[c]->SetLineWidth(2);
    styleH(cms[c], "#D55E00", markOpenSquare, 1.8); cms[c]->SetLineWidth(2);
    gs->SetFillColorAlpha(TColor::GetColor("#0072B2"), 0.25); gs->SetLineColor(0);

    t->cd(); TH1F *fr = t->DrawFrame(xLo, 0.001, xHi, 1.65);
    fr->GetYaxis()->SetTitle("R_{AA} (inclusive jets)"); fr->GetYaxis()->SetTitleSize(0.075); fr->GetYaxis()->SetLabelSize(0.065); fr->GetYaxis()->SetTitleOffset(0.95); fr->GetXaxis()->SetNdivisions(505);
    TLine *l1 = new TLine(xLo, 1, xHi, 1); l1->SetLineStyle(2); l1->SetLineColor(kGray+2); l1->Draw();
    TBox *nb = new TBox(xHi - 9., 1. - norm[c], xHi - 2., 1. + norm[c]); nb->SetFillColor(kGray+1); nb->Draw();
    gs->Draw("2 same"); cms[c]->Draw("E1 X0 same"); raa[c]->Draw("E1 X0 same");
    TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.075); tx.DrawLatex(p == 0 ? 0.30 : 0.08, 0.06, Form("PbPb %s", lab[c]));
    if(p == 0){
      TLegend *lg = makeLegend(0.30, 0.56, 0.97, 0.85, 0.055);
      lg->AddEntry(raa[c], "MinBias only, 5.02 TeV", "lp"); lg->AddEntry(gs, "nominal R_{Z} jet syst.", "f");
      lg->AddEntry(cms[c], "CMS 2.76 TeV (PRC 96)", "lp"); lg->AddEntry(nb, "norm. (T_{AA}, L_{pp})", "f"); lg->Draw();
    }
    bp->cd(); TH1F *fq = bp->DrawFrame(xLo, 0.45, xHi, 1.55);
    fq->GetXaxis()->SetTitle("jet p_{T} (GeV)"); fq->GetYaxis()->SetTitle("5.02 / 2.76 TeV"); fq->GetYaxis()->CenterTitle();
    fq->GetXaxis()->SetNdivisions(505); fq->GetXaxis()->SetTitleSize(0.12); fq->GetXaxis()->SetLabelSize(0.105); fq->GetYaxis()->SetTitleSize(0.095);
    fq->GetYaxis()->SetLabelSize(0.105); fq->GetYaxis()->SetTitleOffset(0.9); fq->GetYaxis()->SetNdivisions(505); fq->GetXaxis()->SetTitleOffset(0.95);
    TLine *l2 = new TLine(xLo, 1, xHi, 1); l2->SetLineStyle(2); l2->SetLineColor(kGray+2); l2->Draw();
    styleH(rat, "#000000", markFilledCircle, 1.5); rat->SetLineWidth(2); rat->Draw("E1 X0 same");
  }
  cv->cd(); TLatex tl; tl.SetNDC();
  tl.SetTextFont(62); tl.SetTextSize(0.05); tl.DrawLatex(L, 0.945, "CMS");
  tl.SetTextFont(52); tl.SetTextSize(0.042); tl.DrawLatex(L + 0.045, 0.945, "Internal");
  tl.SetTextFont(42); tl.SetTextSize(0.038); tl.DrawLatex(L + 0.14, 0.945, "anti-k_{T} R = 0.4 PF jets, |#eta| < 1.6, MinBias trigger");
  tl.SetTextAlign(31); tl.SetTextSize(0.042); tl.DrawLatex(1 - R, 0.945, "PbPb MinBias, pp MinBias (5.02 TeV)");
  gSystem->mkdir(Form("%s/figures/JetsPerZ/RAA_MBonly", repo), true);
  cv->SaveAs(Form("%s/figures/JetsPerZ/RAA_MBonly/inclJetRAA_MBonly_vsCMS2p76TeV%s.pdf", repo, nTag.Data()));

  TFile *fo = TFile::Open(Form("%s/rootFiles/systematics/inclJetRAA_MBonly_PF%s.root", repo, nTag.Data()), "recreate");
  xs->Write();
  for(int c = 1; c <= 4; c++){ raa[c]->Write(); sysR[c]->Write(); yAA[c]->Write(); cms[c]->Write(); }
  fo->Close();
}
