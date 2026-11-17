// PURITY SYSTEMATICS of the nominal 2D fit in a FINER JET-pT BINNING (2026-11-16, user: after
// studyJetPtBinning_2D_PF.C found that 7 bins over 80-200 GeV keep the worst statistical error, "Yes, 7 bins":
// do the purity systematics hold up there?).
//
// The sources of systematics_CSVv2_2D_PF.C (included: its scaleFl, pcorr and the template-ratio tRatio, same
// HistFactory fit), evaluated in an arbitrary list of jet-pT edges on the scans' native grid, windows built by
// buildPP / buildPb of fitPtRelCSVv2Cut_lowPt_PF.C (included; explicit edges):
//   bGS        g->bb share x 0.94, x 1.41 (nominal natural x 1.175)
//   cMult      c yield x 0.8, 0.9, 1.1, 1.2
//   JEU, JER   template ptRel x T_var / T_0 per flavor (T2, T3 / T1)
//   fakeMuJet  PbPb: without the per-event uncorrelated muon-jet pair subtraction
// each max |P_var / P_nom - 1| over its variations, total in quadrature. Two binnings are run: binsA (default
// the nominal 80-100-120-150-200, which must reproduce the chain's bPuritySys_*_csvV2_2D.root: checked and
// printed) and binsB (default the 7 bins 80 90 100 110 120 135 150 200).
// JEU and JER come from the T maps of the same MC, so in a narrow window part of what they give is the
// statistical scatter of those maps rather than a shift; the table prints the stat. error beside the total.
//
// Output: rootFiles/bPurityResults/jetPtBinning_sys_2D_PF<outTag>.root (per binning X = A, B and system:
//         h_pur<X>, h_sys<X>_<source>, h_sys<X>_total _<pp|C4..C1>, relative)
// Figure: figures/bPurity/PFJets_PbPb/csvV2/jetPtBinning/puritySys_binning<outTag>.pdf: pp, 50-80 ... 0-10%;
//         binsB total (band) and its sources (lines), binsA total (black dashed), binsB stat. error (points)
// Usage, from src/plots/bPurity/:  root -l -b -q 'studyJetPtBinning_sys_2D_PF.C'      (~5 min)
#include "fitPtRelCSVv2Cut_lowPt_PF.C"
#include "../../systematics/systematics_CSVv2_2D_PF.C"

namespace jps {
using namespace csys2d; using namespace plow;
const int NSYS = 5; const char *sysTag[NSYS] = {"pp", "C4", "C3", "C2", "C1"}, *sysName[NSYS] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"};
const int sysCls[NSYS] = {-1, 3, 2, 1, 0};          // pbc class; csys2d system index = class + 1
const int NQ = 5; const char *qKey[NQ] = {"bGS", "cMult", "JEU", "JER", "fakeMuJet"};
const char *qLab[NQ] = {"g#rightarrowb#bar{b} share", "c fraction", "JES", "JER", "uncorrelated #mu-jet pairs"};
bool buildE(int s, double lo, double hi, Win &X){
  if(s == 0) return buildPP(fDpp, fMpp, lo, hi, X);
  double a = 0., b = 0.; return buildPb(fDAA, sysCls[s], lo, hi, X, a, b);
}
bool shiftedE(int s, double lo, double hi, const Win &X, int T, double &p){
  Win Y = X;
  for(int F = 0; F < NF; F++){
    double R[NPR]; if(!tRatio(s == 0 ? 0 : sysCls[s] + 1, lo, hi, F, T, R)) return false;
    Y.t[F].tot = 0.;
    for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ Y.t[F].n[k][j] *= R[j]; Y.t[F].e2[k][j] *= R[j]*R[j]; Y.t[F].tot += Y.t[F].n[k][j]; }
  }
  Y.truthB = Y.t[0].tot/(Y.t[0].tot + Y.t[1].tot + Y.t[2].tot); Y.rc = Y.t[1].tot/(Y.t[1].tot + Y.t[2].tot);
  p = pcorr(Y).p;
  return true;
}
std::vector<double> parse(const TString &str){
  std::vector<double> v; std::unique_ptr<TObjArray> a(str.Tokenize(","));
  for(int i = 0; i < a->GetEntries(); i++) v.push_back(((TObjString*) a->At(i))->String().Atof());
  return v;
}
} // namespace jps

// jesWin (2026-11-16, user: "Yes, merged JES/JER"): a coarser edge list for the JEU and JER sources of binsB. Each binsB
// window then takes the JEU / JER of the jesWin window that contains it (evaluated there: nominal and shifted fits of
// the merged window), instead of its own; the other sources stay per window. "" = every source per window.
// e.g. ("80,100,120,150,200", "80,90,100,110,120,135,150,200", "_mergedJESJER", "80,100,120,150,200")
void studyJetPtBinning_sys_2D_PF(TString binsA = "80,100,120,150,200", TString binsB = "80,90,100,110,120,135,150,200", TString outTag = "", TString jesWin = "", TString writeTag = ""){
  // writeTag (2026-11-16, user: "run the b-jet chain in the 7 bins as a variant"): also write the chain's purity files for
  // binning B, rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2<writeTag>.root and bPuritySys_PFJets_PbPb_csvV2<writeTag>.root,
  // in the layout of systematics_CSVv2_2D_PF.C on the axis 50 60 70 | binsB | 300 500 (binsB must span 80-200 or 80-300;
  // every bin outside binsB is copied from the nominal _csvV2_2D files, which have the same edges there). Read by the
  // chain with BJET_PURITY_VARIANT=csvV2_<writeTag without the leading underscore>, e.g. "_2D7bins" <-> csvV2__2D7bins.
  using namespace jps;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  setvbuf(stdout, nullptr, _IOLBF, 0);
  if(!setCats("4cat")) return;
  nFit = nPR; gsMult = 1.175;
  setJetCollection("PFJets"); initPlotStyle();
  fDpp = TFile::Open(Form("%s/%s", repo, dataName)); fMpp = TFile::Open(Form("%s/%s", repo, mcName)); fDAA = TFile::Open(Form("%s/%s", repo, dataPbPb));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  TFile *fCpp = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2_2D.root", repo)), *fCAA = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_PFJets_PbPb_csvV2_2D.root", repo));
  if(!fDpp || !fMpp || !fDAA || !fS[0] || !fS[1] || !fS[2] || !fCpp || !fCAA){ printf("ERROR: inputs\n"); return; }
  const std::vector<double> edges[2] = {parse(binsA), parse(binsB)};
  const char *bn[2] = {"A", "B"};
  const std::vector<double> jesE = parse(jesWin);
  if(jesE.size() >= 2) printf("JEU / JER of binning B taken from the windows %s\n", jesWin.Data());
  TH1D *hP[2][NSYS], *hS[2][NSYS][NQ], *hT[2][NSYS], *hStat[2][NSYS];
  std::vector<double> sgn[NSYS][6];   // binning B: signed P_var / P_nom - 1 of bGS (2) and cMult (4), per window
  for(int b = 0; b < 2; b++){
    const int nw = (int) edges[b].size() - 1;
    for(int s = 0; s < NSYS; s++){
      auto bk = [&](const TString &n){ TH1D *h = new TH1D(n + "_" + sysTag[s], "", nw, edges[b].data()); h->SetDirectory(nullptr); return h; };
      hP[b][s] = bk(Form("h_pur%s", bn[b])); hT[b][s] = bk(Form("h_sys%s_total", bn[b])); hStat[b][s] = bk(Form("h_stat%s", bn[b]));
      for(int q = 0; q < NQ; q++) hS[b][s][q] = bk(Form("h_sys%s_%s", bn[b], qKey[q]));
      printf("\n==== %s, binning %s: relative systematic of the 2D purity (%%)\n", sysName[s], b == 0 ? binsA.Data() : binsB.Data());
      printf("  window    P +- stat (rel %%)        "); for(int q = 0; q < NQ; q++) printf(" %10s", qKey[q]); printf("   total%s\n", b == 0 ? "   chain file total" : "");
      for(int w = 0; w < nw; w++){
        const double lo = edges[b][w], hi = edges[b][w+1];
        gsMult = 1.175; fakeSubPerEvt = true;
        Win X; if(!buildE(s, lo, hi, X)) return;
        const Pc n = pcorr(X);
        auto dev = [&](double v){ return fabs(v/n.p - 1.); };
        double S[NQ] = {};
        { int k = 0; for(double gm : {0.94, 1.41}){ gsMult = gm; Win Y; if(!buildE(s, lo, hi, Y)) return; const double p = pcorr(Y).p; S[0] = TMath::Max(S[0], dev(p)); if(b == 1) sgn[s][k].push_back(p/n.p - 1.); k++; } }
        gsMult = 1.175;
        { int k = 2; for(double m : {0.8, 0.9, 1.1, 1.2}){ Win Y = X; scaleFl(Y, 1, m); const double p = pcorr(Y).p; S[1] = TMath::Max(S[1], dev(p)); if(b == 1) sgn[s][k].push_back(p/n.p - 1.); k++; } }
        int jw = -1; if(b == 1) for(size_t k = 0; k + 1 < jesE.size(); k++) if(0.5*(lo + hi) > jesE[k] && 0.5*(lo + hi) < jesE[k+1]) jw = (int) k;
        if(jw >= 0){   // JEU, JER of the merged window that contains this one
          if(lo < jesE[jw] - 1e-6 || hi > jesE[jw+1] + 1e-6){ printf("ERROR: window %.0f-%.0f straddles the JES/JER windows\n", lo, hi); return; }
          Win Z; if(!buildE(s, jesE[jw], jesE[jw+1], Z)) return;
          const double pz = pcorr(Z).p;
          for(int T : {2, 3}){ double p; if(!shiftedE(s, jesE[jw], jesE[jw+1], Z, T, p)) return; S[2] = TMath::Max(S[2], fabs(p/pz - 1.)); }
          { double p; if(!shiftedE(s, jesE[jw], jesE[jw+1], Z, 1, p)) return; S[3] = fabs(p/pz - 1.); }
        } else {
          for(int T : {2, 3}){ double p; if(!shiftedE(s, lo, hi, X, T, p)) return; S[2] = TMath::Max(S[2], dev(p)); }
          { double p; if(!shiftedE(s, lo, hi, X, 1, p)) return; S[3] = dev(p); }
        }
        if(s > 0){ fakeSubPerEvt = false; Win Y; const bool ok = buildE(s, lo, hi, Y); fakeSubPerEvt = true; if(!ok) return; S[4] = dev(pcorr(Y).p); }
        double t = 0.; for(int q = 0; q < NQ; q++){ t += S[q]*S[q]; hS[b][s][q]->SetBinContent(w + 1, 100.*S[q]); }
        hP[b][s]->SetBinContent(w + 1, n.p); hP[b][s]->SetBinError(w + 1, n.e); hT[b][s]->SetBinContent(w + 1, 100.*sqrt(t)); hStat[b][s]->SetBinContent(w + 1, n.p > 0. ? 100.*n.e/n.p : 0.);
        printf("  %3.0f-%3.0f   %.3f +- %.3f (%5.1f)  ", lo, hi, n.p, n.e, n.p > 0. ? 100.*n.e/n.p : 0.); for(int q = 0; q < NQ; q++) printf(" %10.1f", 100.*S[q]); printf("   %5.1f", 100.*sqrt(t));
        if(b == 0){   // against the chain's purity file (same code path through other window builders)
          TH1 *cp = s == 0 ? (TH1*) fCpp->Get("h_bPurity_nominal") : (TH1*) fCAA->Get(Form("h_bPurity_nominal_%s", sysTag[s]));
          TH1 *ct = s == 0 ? (TH1*) fCpp->Get("h_sysRel_total") : (TH1*) fCAA->Get(Form("h_sysRel_total_%s", sysTag[s]));
          if(cp && ct){ const int ib = cp->FindBin(0.5*(lo + hi)); printf("   %5.1f (P %.3f)%s", 100.*ct->GetBinContent(ib), cp->GetBinContent(ib), fabs(cp->GetBinContent(ib) - n.p) > 2e-3 || fabs(ct->GetBinContent(ib) - sqrt(t)) > 2e-3 ? "  <-- DIFFERS" : ""); }
        }
        printf("\n");
      }
    }
  }
  {
    TFile *fo = TFile::Open(Form("%s/rootFiles/bPurityResults/jetPtBinning_sys_2D_PF%s.root", repo, outTag.Data()), "RECREATE");
    for(int b = 0; b < 2; b++) for(int s = 0; s < NSYS; s++){ hP[b][s]->Write(); hT[b][s]->Write(); hStat[b][s]->Write(); for(int q = 0; q < NQ; q++) hS[b][s][q]->Write(); }
    TNamed info("info", Form("relative systematics (%%) of the 2D ptRel x CSVv2 purity in two jet-pT binnings, A = %s, B = %s; sources as systematics_CSVv2_2D_PF.C", binsA.Data(), binsB.Data()));
    info.Write(); fo->Close();
  }

  // ---- the chain's purity files for binning B
  if(writeTag != ""){
    std::vector<double> ax = {50., 60., 70.}; for(double x : edges[1]) ax.push_back(x);
    if(fabs(edges[1].front() - 80.) > 1e-6 || (fabs(edges[1].back() - 200.) > 1e-6 && fabs(edges[1].back() - 300.) > 1e-6)){ printf("ERROR: writeTag needs binsB over 80-200 or 80-300 GeV\n"); return; }
    if(fabs(edges[1].back() - 200.) < 1e-6) ax.push_back(300.);
    ax.push_back(500.);
    const int nAx = (int) ax.size() - 1, nw = (int) edges[1].size() - 1;
    const char *srcKey7[7] = {"bGS", "cMult", "fitRangeHigh", "tplStat", "JEU", "JER", "fakeMuJet"};
    const int qOf[7] = {0, 1, -1, -1, 2, 3, 4};
    const TString base = Form("%s/rootFiles/systematics/", repo);
    bool firstAA = true;
    for(int s = 0; s < NSYS; s++){
      const bool pp = s == 0; const TString sfx = pp ? "" : Form("_%s", sysTag[s]);
      TFile *fN = pp ? fCpp : fCAA;
      auto old = [&](const TString &n){ TH1 *h = (TH1*) fN->Get(n + sfx); if(!h) printf("ERROR: %s%s missing in the nominal 2D purity file\n", n.Data(), sfx.Data()); return h; };
      auto mkH = [&](const TString &n){ TH1D *h = new TH1D(n + sfx, "", nAx, ax.data()); h->SetDirectory(nullptr); return h; };
      std::vector<TH1D*> out;
      TH1D *p = mkH("h_bPurity_nominal"), *tot = mkH("h_sysRel_total"), *unc = mkH("h_sysRel_uncorr"), *abs = mkH("h_sysAbs_total"), *pa = mkH("h_bPurity_sysAbs");
      TH1D *src[7], *vr[6]; for(int q = 0; q < 7; q++) src[q] = mkH(Form("h_sysRel_%s", srcKey7[q]));
      for(int k = 0; k < 6; k++) vr[k] = mkH(k < 2 ? Form("h_varRel_bGS_%d", k) : Form("h_varRel_cMult_%d", k - 2));
      TH1 *oP = old("h_bPurity_nominal"), *oT = old("h_sysRel_total"), *oU = old("h_sysRel_uncorr"), *oS[7], *oV[6];
      for(int q = 0; q < 7; q++) oS[q] = old(Form("h_sysRel_%s", srcKey7[q]));
      for(int k = 0; k < 6; k++) oV[k] = old(k < 2 ? Form("h_varRel_bGS_%d", k) : Form("h_varRel_cMult_%d", k - 2));
      if(!oP || !oT || !oU) return;
      for(int i = 1; i <= nAx; i++){
        const double c = p->GetBinCenter(i);
        if(c > edges[1].front() && c < edges[1].back()){
          const int w = hP[1][s]->FindBin(c) - 1;
          p->SetBinContent(i, hP[1][s]->GetBinContent(w + 1)); p->SetBinError(i, hP[1][s]->GetBinError(w + 1));
          double t = 0., u = 0.;
          for(int q = 0; q < 7; q++){ const double v = qOf[q] < 0 ? 0. : 0.01*hS[1][s][qOf[q]]->GetBinContent(w + 1); src[q]->SetBinContent(i, v); t += v*v; if(q >= 2) u += v*v; }
          tot->SetBinContent(i, sqrt(t)); unc->SetBinContent(i, sqrt(u));
          for(int k = 0; k < 6; k++) vr[k]->SetBinContent(i, sgn[s][k][w]);
        } else {   // buffer bins: the nominal 2D files (same edges there)
          const int j = oP->FindBin(c);
          if(fabs(oP->GetBinLowEdge(j) - p->GetBinLowEdge(i)) > 1e-6 || fabs(oP->GetBinLowEdge(j + 1) - p->GetBinLowEdge(i + 1)) > 1e-6){ printf("ERROR: buffer bin %g-%g is not a bin of the nominal file\n", p->GetBinLowEdge(i), p->GetBinLowEdge(i + 1)); return; }
          p->SetBinContent(i, oP->GetBinContent(j)); p->SetBinError(i, oP->GetBinError(j)); tot->SetBinContent(i, oT->GetBinContent(j)); unc->SetBinContent(i, oU->GetBinContent(j));
          for(int q = 0; q < 7; q++) if(oS[q]) src[q]->SetBinContent(i, oS[q]->GetBinContent(j));
          for(int k = 0; k < 6; k++) if(oV[k]) vr[k]->SetBinContent(i, oV[k]->GetBinContent(j));
        }
        const double a = tot->GetBinContent(i)*p->GetBinContent(i);
        abs->SetBinContent(i, a); pa->SetBinContent(i, p->GetBinContent(i)); pa->SetBinError(i, a);
      }
      const TString outName = base + (pp ? Form("bPuritySys_total_PFJets_pp_csvV2%s.root", writeTag.Data()) : Form("bPuritySys_PFJets_PbPb_csvV2%s.root", writeTag.Data()));
      TFile *fo = TFile::Open(outName, pp || firstAA ? "RECREATE" : "UPDATE");
      p->Write(); for(int q = 0; q < 7; q++) src[q]->Write(); tot->Write(); abs->Write(); pa->Write(); unc->Write(); for(int k = 0; k < 6; k++) vr[k]->Write();
      if(pp || firstAA){ TNamed info("info", Form("%s PF b purity from the 2D muon ptRel x CSVv2 fit (4cat, HistFactory) in the jet-pT windows %s (studyJetPtBinning_sys_2D_PF.C); JEU / JER %s; fitRangeHigh, tplStat zero; bins outside those windows copied from the nominal _csvV2_2D file. VARIANT, not nominal.", pp ? "pp" : "PbPb (C1..C4)", binsB.Data(), jesWin == "" ? "per window" : Form("from the merged windows %s", jesWin.Data()))); info.Write(); }
      fo->Close(); if(!pp) firstAA = false;
      printf("wrote %s (%s)\n", outName.Data(), sysTag[s]);
    }
  }

  // ---- figure
  const double ptLo = edges[1].front(), ptHi = edges[1].back();
  double ymax = 0.; for(int b = 0; b < 2; b++) for(int s = 0; s < NSYS; s++) ymax = TMath::Max(ymax, TMath::Max(hT[b][s]->GetMaximum(), hStat[1][s]->GetMaximum()));
  const double yCap = 30.;   // totals above it are written in the panel instead of setting the scale
  ymax = TMath::Min(1.30*ymax, yCap);
  TCanvas *cv = new TCanvas("cvJS", "", 2300, 750);
  GridPads g = makeGridPads(cv, NSYS, 1, 0.055, 0.01, 0.16, 0.10);
  const char *hexQ[NQ] = {okabeHex[1], okabeHex[2], okabeHex[5], okabeHex[7], okabeHex[6]};
  const int styQ[NQ] = {1, 1, 7, 3, 1};
  TLegend *lg = nullptr;
  for(int c = 0; c < NSYS; c++){
    g.at(c, 0)->cd();
    TH1F *fr = gPad->DrawFrame(ptLo, 0., ptHi, ymax); g.styleFrame(fr, c, 0, 0.055);
    hT[1][c]->SetFillColorAlpha(kGray, 0.6); hT[1][c]->SetLineColor(kGray + 2); hT[1][c]->SetLineWidth(3); hT[1][c]->Draw("HIST SAME");
    if(c == 0){ lg = makeLegend(0.30, 0.38, 0.98, 0.86, 0.048); lg->AddEntry(hT[1][c], Form("total, %d bins", hT[1][c]->GetNbinsX()), "f"); }
    for(int q = 0; q < NQ; q++){
      if(q == 4 && c == 0){ TH1D *d = new TH1D("jsDummy", "", 1, 0., 1.); styleLine(d, hexQ[q], 4); lg->AddEntry(d, qLab[q], "l"); continue; }
      styleLine(hS[1][c][q], hexQ[q], 4); hS[1][c][q]->SetLineStyle(styQ[q]); hS[1][c][q]->Draw("HIST SAME");
      if(c == 0) lg->AddEntry(hS[1][c][q], qLab[q], "l");
    }
    styleLine(hT[0][c], hexData, 4); hT[0][c]->SetLineStyle(2); hT[0][c]->Draw("HIST SAME");
    styleH(hStat[1][c], hexData, markOpenCircle, 2.4); hStat[1][c]->Draw("P SAME");
    if(c == 0){ lg->AddEntry(hT[0][c], Form("total, %d bins", hT[0][c]->GetNbinsX()), "l"); lg->AddEntry(hStat[1][c], Form("stat., %d bins", hStat[1][c]->GetNbinsX()), "p"); }
    gPad->RedrawAxis();
    g.text(c, 0, 0.08, sysName[c], 0.06, 62);
    if(hT[1][c]->GetMaximum() > ymax) g.text(c, 0, 0.16, Form("off scale: up to %.0f%%", hT[1][c]->GetMaximum()), 0.045);
  }
  if(lg){ g.at(0, 0)->cd(); lg->Draw(); }
  g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "relative uncertainty of the b purity (%)", 0.05, 0.02, 0.02);
  g.cmsLabel(cv, "Internal", "pp, PbPb (5.02 TeV)   #it{p}_{T}^{rel} #times CSVv2 fit");
  savePdfTight(cv, Form("%s/figures/bPurity/PFJets_PbPb/csvV2/jetPtBinning/puritySys_binning%s.pdf", repo, outTag.Data()));
}
