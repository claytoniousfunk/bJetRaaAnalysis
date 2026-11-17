// JET-pT BINNING of the nominal b-purity fit (2026-11-16, user: "explore optimum binning for the b purity
// measurement as to optimize our statistics. Ideally we would have more than 4 points reported").
//
// The nominal purity (2D muon ptRel x CSVv2, 4 categories, HistFactory, MC tagger shapes; p2dhf::fitHF of
// fitPtRelCSVv2_2D_ppPbPb_PF.C, included) is reported in 80-100 / 100-120 / 120-150 / 150-200 GeV. The scans'
// ptRel x jet pT x CSVv2 maps have the jet-pT edges
//     20 30 ... 110 120 | 135 150 175 200 250 300 500
// so without a rescan the finest binning of 80-200 GeV is 8 windows: 80 90 100 110 120 135 150 175 200.
//
// WHAT IS DONE. Every contiguous merge of the native windows between 60 and 300 GeV (12 native windows, 78
// intervals) is fitted on its own, for pp and the four classes, with the SAME code as the nominal: windows
// built by buildPP / buildPb of fitPtRelCSVv2Cut_lowPt_PF.C (included; identical to buildWindows / buildCls,
// explicit edges; PbPb with the per-event uncorrelated-pair subtraction, purity x (1 - fake fraction)),
// 0.2 GeV ptRel bins, fit 0-5 GeV. Per interval:
//   P, stat       the HistFactory fit to the data (what would be reported)
//   expected stat the HistFactory error on an ASIMOV data set (MC mixture at its true fractions, the data's
//                 size): the binning is RANKED on this, not on the observed error, so that the choice is not
//                 tuned on fluctuations of the data it is applied to. It contains the template statistics
//                 (one gamma per cell).
//   dev / ndf     of the hand-written likelihood of the same model (fit quality; HistFactory has no deviance)
//   N_eff         effective MC entries of the b template
// Then every partition of 80-200 GeV into native windows (2^7 = 128) is scored by its WORST expected relative
// purity error over the five systems and its bins, and the best partition for each number of bins is printed
// (the binning has to be common to pp and every class: the purity multiplies the tagged spectrum that goes
// into R_Z). A window counts as unusable if a data fit fails, has no error, or sits within 2 sigma of 0 or 1.
// Also printed: the N-weighted mean of the fine purities against the merged window's purity (does splitting
// move the answer?), and 60-80 / 200-300 GeV for the reach outside the reported range.
//
// NOT covered here: the purity systematics per window (JER / JEU template statistics grow with finer windows),
// and the rest of the chain. A finer purity binning does not by itself give more reported points: the b-jet
// unfolding's reco / truth axes and N_iter, the inclusive unfolding and every systematic are binned
// 80-100-120-150-200 as well and would have to follow.
//
// Output: rootFiles/bPurityResults/jetPtBinning_2D_PF.root (TH2D per system, x = first native window, y = last:
//         h_pur, h_err, h_errAsimov, h_truth, h_status, h_devNdf, h_nData, h_nEffB _<pp|C4..C1>)
// Figure: figures/bPurity/PFJets_PbPb/csvV2/jetPtBinning/purity_binning<outTag>.pdf: pp, 50-80 ... 0-10% in a
//         row, purity (top) and its relative stat. error (bottom) for two binnings (binsA open, binsB filled).
// Usage, from src/plots/bPurity/:
//   root -l -b -q 'studyJetPtBinning_2D_PF.C'                 fit everything (~10 min), rank, draw nominal vs finest
//   root -l -b -q 'studyJetPtBinning_2D_PF.C(false,"80,100,120,150,200","80,90,100,120,150,175,200","_6bins")'
//                                                             no fits: rank and draw from the stored file
#include "fitPtRelCSVv2Cut_lowPt_PF.C"
#include "fitPtRelCSVv2_2D_ppPbPb_PF.C"

namespace jpb {
using namespace plow;
const int NA = 12;                                  // native windows 60-300 GeV
const double aEdge[NA + 1] = {60., 70., 80., 90., 100., 110., 120., 135., 150., 175., 200., 250., 300.};
const int iLo = 2, iHi = 10;                        // 80 and 200 GeV: the reported range
const int NSYS = 5; const char *sysTag[NSYS] = {"pp", "C4", "C3", "C2", "C1"}, *sysName[NSYS] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"};
const int sysCls[NSYS] = {-1, 3, 2, 1, 0};
enum { kP, kE, kEA, kTR, kST, kDV, kND, kNB, NQ };
const char *qName[NQ] = {"h_pur", "h_err", "h_errAsimov", "h_truth", "h_status", "h_devNdf", "h_nData", "h_nEffB"};
TH2D *H[NQ][NSYS];
inline double val(int q, int s, int i, int j){ return H[q][s]->GetBinContent(i + 1, j); }   // interval = native windows i .. j-1
// a data fit that can be reported: converged, with an error, away from the limits
inline bool usable(int s, int i, int j){
  const double p = val(kP, s, i, j), e = val(kE, s, i, j);
  return val(kST, s, i, j) == 0. && e > 0. && p - 2.*e > 0. && p + 2.*e < 1.;
}
inline double relExp(int s, int i, int j){ const double t = val(kTR, s, i, j); return t > 0. ? val(kEA, s, i, j)/t : 9.; }
inline double relObs(int s, int i, int j){ const double p = val(kP, s, i, j); return p > 0. && val(kE, s, i, j) > 0. ? val(kE, s, i, j)/p : 9.; }
std::vector<int> parseEdges(const TString &str){
  std::vector<int> v; std::unique_ptr<TObjArray> a(str.Tokenize(","));
  for(int i = 0; i < a->GetEntries(); i++){
    const double x = ((TObjString*) a->At(i))->String().Atof(); int k = -1;
    for(int e = 0; e <= NA; e++) if(fabs(aEdge[e] - x) < 1e-6) k = e;
    if(k < 0){ printf("ERROR: edge %g is not on the scans' jet-pT grid\n", x); return {}; }
    v.push_back(k);
  }
  return v;
}
TString edgeStr(const std::vector<int> &v){ TString s; for(size_t i = 0; i < v.size(); i++) s += Form("%s%.0f", i ? " " : "", aEdge[v[i]]); return s; }
} // namespace jpb

void studyJetPtBinning_2D_PF(bool refit = true, TString binsA = "80,100,120,150,200", TString binsB = "80,90,100,110,120,135,150,175,200", TString outTag = ""){
  using namespace jpb;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  setvbuf(stdout, nullptr, _IOLBF, 0);
  if(!setCats("4cat")) return;
  nFit = nPR;
  setJetCollection("PFJets"); initPlotStyle();
  const TString outFile = Form("%s/rootFiles/bPurityResults/jetPtBinning_2D_PF.root", repo);
  const TString outDir = Form("%s/figures/bPurity/PFJets_PbPb/csvV2/jetPtBinning/", repo); gSystem->mkdir(outDir, true);
  for(int q = 0; q < NQ; q++) for(int s = 0; s < NSYS; s++){ H[q][s] = new TH2D(Form("%s_%s", qName[q], sysTag[s]), ";first native window;last native window", NA, 0., NA, NA, 0., NA); H[q][s]->SetDirectory(nullptr); }

  if(refit){
    TFile *fDpp = TFile::Open(Form("%s/%s", repo, dataName)), *fMpp = TFile::Open(Form("%s/%s", repo, mcName)), *fDpb = TFile::Open(Form("%s/%s", repo, dataPbPb));
    for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
    if(!fDpp || !fMpp || !fDpb || !fS[0] || !fS[1] || !fS[2]){ printf("ERROR: inputs\n"); return; }
    int uid = 0;
    for(int s = 0; s < NSYS; s++){
      printf("\n==== %s: 2D ptRel x CSVv2 (4cat, HistFactory), every merge of the native jet-pT windows\n", sysName[s]);
      printf("  window    | N data  | b N_eff | MC truth | P +- stat (rel)         st | expected stat (rel) | dev/ndf\n");
      for(int i = 0; i < NA; i++) for(int j = i + 1; j <= NA; j++){
        Win X; double fake = 0., nEffB = 0., nEffL = 0.;
        if(s == 0){
          if(!buildPP(fDpp, fMpp, aEdge[i], aEdge[j], X)) return;
          double sn = 0., sv = 0.; for(int k = 0; k < K; k++) for(int b = 0; b < nPR; b++){ sn += X.t[0].n[k][b]; sv += X.t[0].e2[k][b]; }
          nEffB = sv > 0. ? sn*sn/sv : 0.;
        } else { if(!buildPb(fDpb, sysCls[s], aEdge[i], aEdge[j], X, nEffB, nEffL)) return; fake = lastFakeFrac; }
        const p2dhf::HFResult r = p2dhf::fitHF(X, X.d, uid++);
        // Asimov at the fitted data's size (PbPb: after the pair subtraction)
        double one[NF][KMAX]; for(int F = 0; F < NF; F++) for(int k = 0; k < K; k++) one[F][k] = 1.;
        double pd[KMAX][NPR], nD = 0.; for(int k = 0; k < K; k++) for(int b = 0; b < nPR; b++) nD += X.d[k][b];
        Win Y = X; Y.nData = nD; const double tr = asimov(Y, one, pd);
        const p2dhf::HFResult ra = p2dhf::fitHF(X, pd, uid++);
        Result rm = fitWin(models[1], X, X.d);
        const double sc = 1. - fake, p = r.pur*sc, e = r.epur*sc, ea = ra.epur*sc, t = tr*sc;
        const double q[NQ] = {p, e, ea, t, (double) r.status, rm.ndf > 0 ? rm.dev/rm.ndf : 0., X.nData, nEffB};
        for(int k = 0; k < NQ; k++) H[k][s]->SetBinContent(i + 1, j, q[k]);
        printf("  %3.0f-%3.0f   | %7.0f | %7.0f | %.3f    | %.3f +- %.3f (%5.1f%%)  %d | %.3f (%5.1f%%)       | %.2f\n", aEdge[i], aEdge[j], X.nData, nEffB, t, p, e, p > 0. ? 100.*e/p : 0., r.status, ea, t > 0. ? 100.*ea/t : 0., q[kDV]);
      }
    }
    TFile *fo = TFile::Open(outFile, "RECREATE");
    for(int q = 0; q < NQ; q++) for(int s = 0; s < NSYS; s++) H[q][s]->Write();
    TH1D hE("h_nativeEdges", "", NA, aEdge); hE.Write();
    TNamed info("info", Form("b purity (2D ptRel x CSVv2 4cat, HistFactory, MC tagger shapes) for every contiguous merge of the native jet-pT windows; TH2 bin (i+1, j) = native windows i .. j-1 of h_nativeEdges. pp %s / %s; PbPb %s / %s", dataName, mcName, dataPbPb, phFmt));
    info.Write(); fo->Close();
  } else {
    TFile *fi = TFile::Open(outFile); if(!fi){ printf("ERROR: no stored file, run with refit = true\n"); return; }
    for(int q = 0; q < NQ; q++) for(int s = 0; s < NSYS; s++){ TH2D *h = (TH2D*) fi->Get(H[q][s]->GetName()); if(!h){ printf("ERROR: %s\n", H[q][s]->GetName()); return; } H[q][s] = (TH2D*) h->Clone(); H[q][s]->SetDirectory(nullptr); }
    fi->Close();
  }

  // ---- every partition of 80-200 GeV, scored by the worst expected relative error
  struct Part { std::vector<int> e; double worstExp = 0., worstObs = 0., meanExp = 0.; int nBad = 0; int ws = 0, wi = 0, wj = 0; };
  auto score = [&](const std::vector<int> &e){
    Part P; P.e = e; int n = 0;
    for(size_t b = 0; b + 1 < e.size(); b++) for(int s = 0; s < NSYS; s++){
      const double x = relExp(s, e[b], e[b+1]), o = relObs(s, e[b], e[b+1]);
      if(x > P.worstExp){ P.worstExp = x; P.ws = s; P.wi = e[b]; P.wj = e[b+1]; }
      P.worstObs = TMath::Max(P.worstObs, o); P.meanExp += x; n++;
      if(!usable(s, e[b], e[b+1])) P.nBad++;
    }
    P.meanExp /= n; return P;
  };
  const int nIn = iHi - iLo - 1;   // inner edges
  std::vector<Part> all;
  for(int m = 0; m < (1 << nIn); m++){
    std::vector<int> e = {iLo}; for(int b = 0; b < nIn; b++) if(m & (1 << b)) e.push_back(iLo + 1 + b);
    e.push_back(iHi); all.push_back(score(e));
  }
  printf("\n==== partitions of 80-200 GeV into the native windows, ranked on the WORST expected relative stat. error of the purity\n");
  printf("     (over pp + four classes and every bin; 'unusable' = data fits failed / without error / within 2 sigma of a limit)\n");
  for(int nb = 1; nb <= nIn + 1; nb++){
    std::vector<Part> v; for(const Part &P : all) if((int) P.e.size() - 1 == nb) v.push_back(P);
    std::sort(v.begin(), v.end(), [](const Part &a, const Part &b){ return a.nBad != b.nBad ? a.nBad < b.nBad : a.worstExp < b.worstExp; });
    printf("  %d bins (%zu partitions)\n", nb, v.size());
    for(size_t r = 0; r < TMath::Min(v.size(), (size_t) 3); r++)
      printf("    %-36s worst expected %5.1f%% (%s %.0f-%.0f), mean %4.1f%%, worst observed %5.1f%%, unusable fits %d\n", edgeStr(v[r].e).Data(), 100.*v[r].worstExp, sysName[v[r].ws], aEdge[v[r].wi], aEdge[v[r].wj], 100.*v[r].meanExp, 100.*v[r].worstObs, v[r].nBad);
  }

  // ---- the two binnings asked for: table, split-vs-merged check, figure
  std::vector<int> eA = parseEdges(binsA), eB = parseEdges(binsB);
  if(eA.size() < 2 || eB.size() < 2) return;
  for(const std::vector<int> *e : {&eA, &eB}){
    const Part P = score(*e);
    printf("\n==== binning %s: purity +- stat (expected stat), per system\n", edgeStr(*e).Data());
    printf("  window    |"); for(int s = 0; s < NSYS; s++) printf(" %-27s|", sysName[s]); printf("\n");
    for(size_t b = 0; b + 1 < e->size(); b++){
      printf("  %3.0f-%3.0f   |", aEdge[(*e)[b]], aEdge[(*e)[b+1]]);
      for(int s = 0; s < NSYS; s++) printf(" %.3f +- %.3f (%.3f)%s      |", val(kP, s, (*e)[b], (*e)[b+1]), val(kE, s, (*e)[b], (*e)[b+1]), val(kEA, s, (*e)[b], (*e)[b+1]), usable(s, (*e)[b], (*e)[b+1]) ? " " : "*");
      printf("\n");
    }
    printf("  worst expected %.1f%%, worst observed %.1f%%, unusable (*) %d\n", 100.*P.worstExp, 100.*P.worstObs, P.nBad);
  }
  // does splitting move the answer? N-weighted mean of the native-window purities against each nominal window's own fit
  printf("\n==== N-weighted mean purity of the native windows minus the merged window's fit, in units of the merged fit's stat. error\n");
  { const std::vector<int> nom = parseEdges("80,100,120,150,200");
    printf("  window    |"); for(int s = 0; s < NSYS; s++) printf(" %-22s|", sysName[s]); printf("\n");
    for(size_t b = 0; b + 1 < nom.size(); b++){
      printf("  %3.0f-%3.0f   |", aEdge[nom[b]], aEdge[nom[b+1]]);
      for(int s = 0; s < NSYS; s++){
        double sw = 0., sp = 0.; for(int i = nom[b]; i < nom[b+1]; i++){ const double n = val(kND, s, i, i + 1); sw += n; sp += n*val(kP, s, i, i + 1); }
        const double m = sp/sw, p = val(kP, s, nom[b], nom[b+1]), e = val(kE, s, nom[b], nom[b+1]);
        printf(" %.3f vs %.3f (%+.1f)   |", m, p, e > 0. ? (m - p)/e : 0.);
      }
      printf("\n");
    } }

  // ---- figure: purity (top) and relative stat. error (bottom), binsA open circles, binsB filled squares
  auto mk = [&](const std::vector<int> &e, int s, int mode, const char *nm){
    std::vector<double> x; for(int k : e) x.push_back(aEdge[k]);
    TH1D *h = new TH1D(Form("%s_%s", nm, sysTag[s]), "", (int) x.size() - 1, x.data()); h->SetDirectory(nullptr);
    for(size_t b = 0; b + 1 < e.size(); b++){
      const double p = val(kP, s, e[b], e[b+1]), er = val(kE, s, e[b], e[b+1]);
      if(mode == 0){ h->SetBinContent(b + 1, p); h->SetBinError(b + 1, er); }
      else if(mode == 1) h->SetBinContent(b + 1, p > 0. ? 100.*er/p : 0.);
      else h->SetBinContent(b + 1, 100.*relExp(s, e[b], e[b+1]));
    }
    return h;
  };
  const double ptLo = aEdge[TMath::Min(eA.front(), eB.front())], ptHi = aEdge[TMath::Max(eA.back(), eB.back())];
  double rMax = 0.;
  for(int s = 0; s < NSYS; s++) for(const std::vector<int> *e : {&eA, &eB}) for(size_t b = 0; b + 1 < e->size(); b++){
    const double o = 100.*relObs(s, (*e)[b], (*e)[b+1]), x = 100.*relExp(s, (*e)[b], (*e)[b+1]); if(o < 100.) rMax = TMath::Max(rMax, o); if(x < 100.) rMax = TMath::Max(rMax, x); }
  rMax = TMath::Min(60., 1.25*rMax);
  TCanvas *cv = new TCanvas("cvBin", "", 2300, 1100);
  GridPads g = makeGridPads(cv, NSYS, 2, 0.060, 0.01, 0.11, 0.075);
  const double split = g.bm + 0.40*(1. - g.bm - g.tm), fs = 0.055;
  for(int c = 0; c < NSYS; c++){
    TPad *t = g.at(c, 0), *b = g.at(c, 1); double x1, y1, x2, y2;
    t->GetPadPar(x1, y1, x2, y2); t->SetPad(x1, split, x2, y2); t->SetTopMargin(g.tm/(1. - split));
    b->GetPadPar(x1, y1, x2, y2); b->SetPad(x1, 0., x2, split); b->SetBottomMargin(g.bm/split);
  }
  for(int c = 0; c < NSYS; c++){
    g.at(c, 0)->cd();
    TH1F *fr = gPad->DrawFrame(ptLo, 0., ptHi, 1.0); g.styleFrame(fr, c, 0, fs);
    fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0); fr->GetYaxis()->ChangeLabel(1, -1, 0.); fr->GetYaxis()->SetNdivisions(505);
    TH1D *hA = mk(eA, c, 0, "pA"), *hB = mk(eB, c, 0, "pB");
    styleH(hB, okabeHex[5], markFilledSquare, 2.2); hB->SetLineWidth(3); hB->Draw("E1 SAME");
    styleH(hA, hexData, markOpenCircle, 3.0); hA->SetLineWidth(3); hA->Draw("E1 SAME");
    gPad->RedrawAxis();
    g.text(c, 0, 0.10, sysName[c], 0.065, 62);
    if(c == 0){
      const double fL = gPad->GetLeftMargin(), fW = 1. - fL - gPad->GetRightMargin(), fT = 1. - gPad->GetTopMargin();
      TLegend *lg = new TLegend(fL + 0.04*fW, fT - 0.36, fL + fW, fT - 0.15);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(43); lg->SetTextSize(30);
      lg->AddEntry(hA, Form("%d bins", (int) eA.size() - 1), "pe"); lg->AddEntry(hB, Form("%d bins", (int) eB.size() - 1), "pe"); lg->Draw();
    }
    g.at(c, 1)->cd();
    TH1F *fb = gPad->DrawFrame(ptLo, 0., ptHi, rMax); g.styleFrame(fb, c, 1, fs);
    fb->GetXaxis()->SetLabelFont(43); fb->GetXaxis()->SetLabelSize(34); fb->GetYaxis()->SetLabelFont(43); fb->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0);
    fb->GetYaxis()->SetNdivisions(505); fb->GetXaxis()->SetTickLength(0.03*0.60/0.40);
    TH1D *xA = mk(eA, c, 2, "xA"), *xB = mk(eB, c, 2, "xB"), *oA = mk(eA, c, 1, "oA"), *oB = mk(eB, c, 1, "oB");
    styleLine(xA, hexData, 3); xA->SetLineStyle(7); xA->Draw("HIST SAME");
    styleLine(xB, okabeHex[5], 3); xB->SetLineStyle(7); xB->Draw("HIST SAME");
    styleH(oB, okabeHex[5], markFilledSquare, 2.2); oB->Draw("P SAME");
    styleH(oA, hexData, markOpenCircle, 3.0); oA->Draw("P SAME");
    gPad->RedrawAxis();
    if(c == 0){
      const double fL = gPad->GetLeftMargin(), fW = 1. - fL - gPad->GetRightMargin();
      TLegend *lg = new TLegend(fL + 0.04*fW, 0.80, fL + fW, 0.97);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(43); lg->SetTextSize(28);
      lg->AddEntry(oA, "data", "p"); lg->AddEntry(xA, "expected", "l"); lg->Draw();
    }
  }
  cv->cd();
  { TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.040); t.SetTextAlign(21);
    t.DrawLatex(g.lm + (1. - g.lm - g.rm)/2., 0.016, "#it{p}_{T}^{jet} (GeV)");
    t.SetTextAngle(90);
    t.DrawLatex(0.018, split + 0.5*(1. - g.tm - split), "b purity of muon-tagged jets");
    t.SetTextSize(0.032); t.DrawLatex(0.018, g.bm + 0.5*(split - g.bm), "stat. error (%)"); }
  g.cmsLabel(cv, "Internal", "pp, PbPb (5.02 TeV)   #it{p}_{T}^{rel} #times CSVv2 fit", 1.0);
  savePdfTight(cv, outDir + "purity_binning" + outTag + ".pdf");
}
