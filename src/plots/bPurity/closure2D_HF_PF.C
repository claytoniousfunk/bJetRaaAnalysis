// CLOSURE of the NOMINAL 2D b-purity fit (muon ptRel x CSVv2 4cat, HistFactory, MC tagger shapes; p2dhf::fitHF of
// fitPtRelCSVv2_2D_ppPbPb_PF.C) on FLUCTUATING pseudo-data (2026-11-16, user: "run closure"), pp and the four PbPb
// classes, windows 80-300 GeV, templates at the NATURAL g->bb share (the pseudo-data have it; the 1.175 factor is a
// correction to data):
//   (1) Asimov: the templates mixed at their true fractions at the data's size, no fluctuations (must close exactly).
//   (2) Poisson toys of (1) (nToy per window), templates fixed: mean bias, RMS, and pull mean / width with the fit's own
//       error, i.e. whether the HistFactory fit and its error are unbiased at the data's statistics.
//   (3) PbPb: the PYTHIA+HYDJET DiJet all-flavor map itself as pseudo-data (its own fluctuations; DiJet is also in the
//       templates), against DiJet's true b fraction, as the CSVv2-cut closure (fitPtRelCSVv2Cut_PF_PbPb.C); the pull
//       uses the fit error rescaled to DiJet's effective statistics. pp: PYTHIA is the template sample itself, so
//       there is no independent pseudo-data there (only (1), (2)).
//   (4) STRESS TEST of the tagger shapes (printed only): Asimov pseudo-data in which a fraction of the b jets is moved
//       from the three CSVv2 > 0.54 categories into CSVv2 < 0.54 (SF_b(loose) = 0.9, 0.8), fitted with the MC tagger
//       shapes: the bias a mis-modelled b-tagging efficiency of that size would leave in the nominal purity.
// Output: rootFiles/bPurityResults/closure2D_HF_PF.root (h_{asimov,toyBias,toyRMS,toyPullMean,toyPullWidth,dijet,
//         dijetErr,sf0p9,sf0p8}_<sys>, all fitted - true purity), figure figures/bPurity/PFJets_PbPb/csvV2/closure2D_HF.pdf
// Usage, from src/plots/bPurity/:  root -l -b -q 'closure2D_HF_PF.C(200)'   (-1 = redraw the figure from the stored file)
#include "fitPtRelCSVv2_2D_ppPbPb_PF.C"
#include "TRandom3.h"

void closure2D_HF_PF(int nToy = 200){
  using namespace pbc; using namespace p2d;
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  setvbuf(stdout, nullptr, _IOLBF, 0);
  setCats("4cat"); nFit = nPR;
  gsShift = 0.; gsMult = 0.;   // natural g->bb share
  setJetCollection("PFJets"); initPlotStyle();
  TFile *fDpp = TFile::Open(Form("%s/%s", repo, dataName)), *fMpp = TFile::Open(Form("%s/%s", repo, mcName)), *fDpb = TFile::Open(Form("%s/%s", repo, dataPbPb));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  if(!fDpp || !fMpp || !fDpb || !fS[0] || !fS[1] || !fS[2]){ printf("ERROR: inputs\n"); return; }
  const int NCOL = 5; const char *colName[NCOL] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"}, *colTag[NCOL] = {"pp", "C4", "C3", "C2", "C1"};
  const int colCls[NCOL] = {-1, 3, 2, 1, 0};
  enum { kAs, kTB, kTR, kPM, kPW, kDJ, kDJE, kS9, kS8, NH };
  const char *hn[NH] = {"asimov", "toyBias", "toyRMS", "toyPullMean", "toyPullWidth", "dijet", "dijetErr", "sf0p9", "sf0p8"};
  TH1D *H[NH][NCOL];
  for(int q = 0; q < NH; q++) for(int c = 0; c < NCOL; c++){ H[q][c] = new TH1D(Form("h_%s_%s", hn[q], colTag[c]), "", NW, ptEdges + 3); H[q][c]->SetDirectory(nullptr); }
  TRandom3 rng(4357);
  int uid = 0;
  // nToy < 0: no fits, redraw the figure from the stored output file
  if(nToy < 0){
    TFile *fi = TFile::Open(Form("%s/rootFiles/bPurityResults/closure2D_HF_PF%s.root", repo, ptBinTag()));
    if(!fi){ printf("ERROR: no stored closure file\n"); return; }
    for(int q = 0; q < NH; q++) for(int c = 0; c < NCOL; c++){ TH1D *h = (TH1D*) fi->Get(H[q][c]->GetName()); if(!h){ printf("ERROR: %s\n", H[q][c]->GetName()); return; } H[q][c] = (TH1D*) h->Clone(); H[q][c]->SetDirectory(nullptr); }
    fi->Close();
  }
  for(int c = 0; c < NCOL && nToy >= 0; c++){
    if(c == 0 && !buildWindows(fDpp, fMpp)) return;
    printf("\n==== %s (fitted - true purity; natural g->bb share)\n", colName[c]);
    printf("  window   truth  | Asimov   | toys: bias    RMS    fit err  pull mean  width | DiJet: fit-true  (err at DiJet stat, pull) | stress SF_b 0.9 / 0.8\n");
    for(int w = 0; w < NW; w++){
      Win X; double nEffB = 0.;
      if(c == 0) X = W[w]; else if(!buildCls(fDpb, colCls[c], w, X, nEffB)) return;
      // size = the fitted data's total (PbPb: after the uncorrelated-pair subtraction)
      double nD = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) nD += X.d[k][j];
      X.nData = nD;
      double one[NF][KMAX]; for(int F = 0; F < NF; F++) for(int k = 0; k < K; k++) one[F][k] = 1.;
      double mu[KMAX][NPR]; const double tr = asimov(X, one, mu);
      const p2dhf::HFResult ra = p2dhf::fitHF(X, mu, uid++);
      H[kAs][c]->SetBinContent(w + 1, ra.pur - tr); H[kAs][c]->SetBinError(w + 1, ra.epur);
      // (2) toys
      double sb = 0., sb2 = 0., sp = 0., sp2 = 0., se = 0.; int n = 0;
      for(int t = 0; t < nToy; t++){
        double d[KMAX][NPR];
        for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) d[k][j] = rng.Poisson(mu[k][j]);
        const p2dhf::HFResult r = p2dhf::fitHF(X, d, uid++);
        if(r.status != 0 || !(r.epur > 0.)) continue;
        const double b = r.pur - tr, p = b/r.epur;
        sb += b; sb2 += b*b; sp += p; sp2 += p*p; se += r.epur; n++;
      }
      const double mb = n ? sb/n : 0., rms = n > 1 ? sqrt(TMath::Max(0., sb2/n - mb*mb)) : 0., mp = n ? sp/n : 0., wp = n > 1 ? sqrt(TMath::Max(0., sp2/n - mp*mp)) : 0.;
      H[kTB][c]->SetBinContent(w + 1, mb); H[kTB][c]->SetBinError(w + 1, n ? rms/sqrt(n) : 0.);
      H[kTR][c]->SetBinContent(w + 1, rms); H[kPM][c]->SetBinContent(w + 1, mp); H[kPW][c]->SetBinContent(w + 1, wp);
      // (3) DiJet pseudo-data (PbPb)
      double dj = 0., djE = 0., djPull = 0.;
      if(c > 0){
        double ps[KMAX][NPR] = {}, pe[KMAX][NPR] = {};
        fillCounts(fS[0], names(h3, {"allJets"}, colCls[c]), wLo[w], wHi[w], ps, pe);
        double nAll = 0., vAll = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ nAll += ps[k][j]; vAll += pe[k][j]; }
        double d[KMAX][NPR]; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) d[k][j] = ps[k][j]*nD/nAll;   // the data's size, DiJet's shape
        const p2dhf::HFResult r = p2dhf::fitHF(X, d, uid++);
        const double nEff = vAll > 0. ? nAll*nAll/vAll : 0.;
        dj = r.pur - X.truthB; djE = nEff > 0. ? r.epur*sqrt(nD/nEff) : 0.; djPull = djE > 0. ? dj/djE : 0.;
        H[kDJ][c]->SetBinContent(w + 1, dj); H[kDJ][c]->SetBinError(w + 1, djE); H[kDJE][c]->SetBinContent(w + 1, nEff);
      }
      // (4) stress: move b jets out of the tagged categories, truth unchanged (b total conserved)
      double st[2] = {};
      for(int s = 0; s < 2; s++){
        const double sf = s == 0 ? 0.9 : 0.8;
        Win Y = X; double moved[NPR] = {};
        for(int k = 1; k < K; k++) for(int j = 0; j < nPR; j++){ const double m = (1. - sf)*Y.t[0].n[k][j]; Y.t[0].n[k][j] -= m; moved[j] += m; }
        for(int j = 0; j < nPR; j++) Y.t[0].n[0][j] += moved[j];
        double ds[KMAX][NPR]; const double ts = asimov(Y, one, ds);
        const p2dhf::HFResult r = p2dhf::fitHF(X, ds, uid++);   // fitted with the unmodified (MC) templates
        st[s] = r.pur - ts;
        H[s == 0 ? kS9 : kS8][c]->SetBinContent(w + 1, st[s]);
      }
      printf("  %3.0f-%3.0f  %.3f  | %+.4f  | %+.4f  %.4f  %.4f   %+.2f    %.2f (%d toys) |", wLo[w], wHi[w], tr, ra.pur - tr, mb, rms, n ? se/n : 0., mp, wp, n);
      if(c > 0) printf(" %+.3f  (%.3f, %+.1f; DiJet N_eff %.0f)", dj, djE, djPull, H[kDJE][c]->GetBinContent(w + 1)); else printf(" --                                        ");
      printf(" | %+.3f / %+.3f\n", st[0], st[1]);
    }
  }
  if(nToy >= 0){
    TFile *fo = TFile::Open(Form("%s/rootFiles/bPurityResults/closure2D_HF_PF%s.root", repo, ptBinTag()), "RECREATE");
    for(int q = 0; q < NH; q++) for(int c = 0; c < NCOL; c++) H[q][c]->Write();
    TNamed info("info", Form("closure of the 2D ptRel x CSVv2 (4cat) HistFactory purity fit, natural g->bb share, %d Poisson toys per window; fitted - true purity. pp %s / %s, PbPb templates %s", nToy, dataName, mcName, phFmt));
    info.Write(); fo->Close();
  }
  // ---- figure: fitted - true purity, 1x5 row; toys (mean, bars = toy RMS), Asimov, DiJet pseudo-data (bars at DiJet's statistics)
  TCanvas *cv = new TCanvas("cvCl", "", 2300, 900);
  GridPads g = makeGridPads(cv, NCOL, 1, 0.065, 0.01, 0.16, 0.10);
  const double yr = 0.30;
  for(int c = 0; c < NCOL; c++){
    g.at(c, 0)->cd();
    TH1F *fr = gPad->DrawFrame(80., -yr, 300., yr); g.styleFrame(fr, c, 0, 0.055);
    fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0);
    fr->GetXaxis()->SetLabelFont(43); fr->GetXaxis()->SetLabelSize(34);
    TLine *z = new TLine(80., 0., 300., 0.); z->SetLineStyle(2); z->SetLineWidth(2); z->SetLineColor(kGray + 2); z->Draw();
    TH1D *toy = (TH1D*) H[kTB][c]->Clone(Form("toyDraw%d", c)); toy->SetDirectory(nullptr);
    for(int w = 1; w <= NW; w++) toy->SetBinError(w, H[kTR][c]->GetBinContent(w));
    toy->SetFillColorAlpha(kGray + 1, 0.45); toy->SetLineWidth(0); toy->SetMarkerSize(0); toy->Draw("E2 SAME");
    TH1D *as = (TH1D*) H[kAs][c]->Clone(Form("asDraw%d", c)); as->SetDirectory(nullptr); for(int w = 1; w <= NW; w++) as->SetBinError(w, 0.);
    styleH(as, hexData, markOpenSquare, 2.2); as->Draw("P X0 SAME");
    if(c > 0){ styleH(H[kDJ][c], okabeHex[6], markFilledCircle, 2.2); H[kDJ][c]->SetLineWidth(3); H[kDJ][c]->Draw("E1 X0 SAME"); }
    gPad->RedrawAxis();
    g.text(c, 0, 0.10, colName[c], 0.065, 62);
    if(c == 0){
      TLegend *lg = makeLegend(0.24, 0.20, 0.98, 0.42, 0.055);
      TH1D *tl = (TH1D*) toy->Clone("tl"); lg->AddEntry(tl, "toys: mean #pm RMS", "f");
      lg->AddEntry(as, "Asimov", "p");
      TH1D *dl = new TH1D("dl", "", 1, 0., 1.); styleH(dl, okabeHex[6], markFilledCircle, 2.2); dl->SetLineWidth(3);
      lg->AddEntry(dl, "DiJet as pseudo-data", "pe"); lg->Draw();
    }
  }
  g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "fitted #minus true b purity", 0.05, 0.02, 0.02);
  g.cmsLabel(cv, "Internal", "Simulation   2D #it{p}_{T}^{rel} #times CSVv2 fit", 1.0);
  savePdfTight(cv, Form("%s/figures/bPurity/PFJets_PbPb/csvV2/closure2D_HF%s.pdf", repo, ptBinTag()));
}
