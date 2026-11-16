// The SIMPLE SURFACE FIT of the b purity and its systematics ONE SOURCE AT A TIME, PF pp and the four PbPb classes
// (2026-11-16, user: "the simplest possible fit to a 2d surface from templates ... fine bins, 2d minimization. Then I
// want to redo the systematics one by one"; choices the same day: b fraction the only free parameter with c / light
// fixed at MC, plain Poisson likelihood, 25 uniform CSVv2 bins, a STUDY (nominal files, chain and AN untouched), and
// the c-yield source replaced by LIGHT x 0.8-1.2).
// Fit, per system and jet-pT window 80-300 GeV: the muon ptRel x CSVv2 map of all muon-tagged jets, 25 CSVv2 bins of
// 0.04 over [0, 1] (the first also holds the jets without a discriminant value) x 25 ptRel bins of 0.2 GeV over 0-5 GeV,
//   mu_kj = N [ f_b B_kj + (1 - f_b) (r_c C_kj + (1 - r_c) L_kj) ],   B, C, L unit-normalized MC templates,
//   r_c = c / (c + light) from the MC yields, N = the data in the fitted cells, f_b free,
// binned Poisson likelihood WITH the template-statistics term (Barlow-Beeston lite, one factor per cell solved in closed
// form, no extra fit parameters; RESTORED 2026-11-16, user, after the plain-Poisson version first asked for moved the
// purity by up to +0.16 in 0-10% below 120 GeV, where the model fits badly, a shift no toy source covered). The
// template statistics are therefore in the fit's statistical error and "tplStat" is not a source of this fit.
// BJET_SS_LIKELIHOOD=plain = the plain Poisson likelihood (templates exact), outputs tagged _plain; its template
// statistics are the "tplStat" toys. Cells with no simulated jet of any flavor are left out (data share printed).
// Minimizer: Minuit2 through p2d::fitWinMinuit, the existing likelihood code. Templates, g->bb share (natural x 1.175), PbPb fake-pair subtraction: the nominal 2D fit's (csys2d::build).
// P = f_b x (fitted / all muon-tagged jets), i.e. b / ALL muon-tagged jets as every other purity here.
//
// Steps (argument), each prints a table and writes its own root file and figure:
//   "nominal"    P and its statistical error, sum pull^2 / cells (cells with fit >= 5 jets; the pull's sigma includes the
//                template variance, so it is a chi2 / ndf-like number), Asimov closure, against the 4-category
//                HistFactory nominal (bPuritySys_*_csvV2_2D.root)
//   "bGS"        g->bb share of the b template x 0.94, x 1.41 (nominal x 1.175)
//   "lightMult"  light yield x 0.8, 0.9, 1.1, 1.2 (changes r_c; replaces cMult)
//   "JEU"        template ptRel scaled by T2 / T0, T3 / T0 per flavor (as systematics_CSVv2_2D_PF.C)
//   "JER"        the same with T1 / T0
//   "fakeMuJet"  PbPb: without the per-event uncorrelated muon-jet pair subtraction
//   "files"      the chain's purity files from the stored steps, tag _csvV2_surf (BJET_PURITY_VARIANT=csvV2__surf)
//   "tplStat"    (plain likelihood only) toys: every template cell fluctuated within its MC error (Gaussian, floored at 0), RMS of P
// A source step needs the nominal step's file. Each source figure shows the signed P_var / P_nom - 1 per variation and,
// dashed, +- the same source of the 4-category nominal where it exists.
// Outputs: rootFiles/systematics/simpleSurface/simpleSurface_PF_<step>.root
//          figures/systematics/bPurity/simpleSurface/simpleSurface_<step>.pdf
// Usage, from src/systematics/:  root -l -b -q 'systematics_simpleSurface_PF.C("nominal")'   then one source at a time
#define P2D_KMAX 25
#include "systematics_CSVv2_2D_PF.C"   // csys2d::build, scaleFl, tRatio; includes the 2D fit macros
#include "TRandom3.h"

namespace ssurf {
using namespace csys2d;
// CSVv2 bins: 25 nominal; BJET_SS_NCSV=<n> (a divisor of 50 up to 25: 5, 10) for the binning check, outputs tagged _csv<n>
inline int initNCSV(){ const char *e = std::getenv("BJET_SS_NCSV"); const int n = e ? atoi(e) : 25; if(n < 2 || n > 25 || 50 % n != 0){ printf("ERROR: BJET_SS_NCSV=%s\n", e); gSystem->Exit(1); } return n; }
const int NCSV = initNCSV();
inline bool initPlain(){ const char *e = std::getenv("BJET_SS_LIKELIHOOD"); const TString v = e ? e : ""; if(v != "" && v != "plain"){ printf("ERROR: BJET_SS_LIKELIHOOD=%s (plain)\n", e); gSystem->Exit(1); } return v == "plain"; }
const bool plainL = initPlain();
inline TString csvTag(){ return (NCSV == 25 ? TString("") : TString::Format("_csv%d", NCSV)) + (plainL ? "_plain" : ""); }
const int order[NSY] = {0, 4, 3, 2, 1};   // display: pp, 50-80 ... 0-10%
inline void setUniform(){ K = NCSV; for(int k = 0; k < K; k++) catUp[k] = (k + 1.)/NCSV; catUp[K - 1] = 1.01; }
struct Fit { double p = 0., e = 0., fb = 0., noMC = 0., chi2 = 0.; int nChi = 0, status = -1; };
// the simple fit of window X (a copy): data in cells without MC dropped; plainL: template errors zeroed -> plain Poisson
Fit psimple(Win X, bool withChi2 = false){
  static Win Y; Y = X; Fit o;
  double sd = 0., kept = 0.;
  for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){
    sd += Y.d[k][j];
    if(Y.t[0].n[k][j] + Y.t[1].n[k][j] + Y.t[2].n[k][j] <= 0.){ o.noMC += Y.d[k][j]; Y.d[k][j] = 0.; }
    kept += Y.d[k][j];
    if(plainL) for(int F = 0; F < NF; F++) Y.t[F].e2[k][j] = 0.;
  }
  o.noMC = sd > 0. ? o.noMC/sd : 0.;
  const double kAll = X.nData > 0. ? sd/X.nData : 1.;
  Result r = fitWinMinuit(models[1], Y, Y.d);
  o.fb = r.pur; o.p = r.pur*kAll; o.e = r.epur*kAll; o.status = r.status;
  if(withChi2){   // pulls with the template variance, from the templates' own errors
    static Win Z; Z = Y; for(int F = 0; F < NF; F++) for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) Z.t[F].e2[k][j] = X.t[F].e2[k][j];
    static double mu[KMAX][NPR], va[KMAX][NPR];
    cur = &models[1]; cw = &Z; cd = Z.d; Kf = K; expectAll(r.par, mu, va);
    for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) if(mu[k][j] >= 5.){ o.chi2 += pow(Z.d[k][j] - mu[k][j], 2)/(mu[k][j] + va[k][j]); o.nChi++; }
  }
  return o;
}
// template ptRel scaled by T_var / T_0 per flavor (csys2d::shifted with the simple fit)
bool shiftedSimple(int sy, int w, const Win &X, int T, double &p){
  static Win Y; Y = X;
  for(int F = 0; F < NF; F++){
    double R[NPR]; if(!tRatio(sy, w, F, T, R)) return false;
    Y.t[F].tot = 0.;
    for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ Y.t[F].n[k][j] *= R[j]; Y.t[F].e2[k][j] *= R[j]*R[j]; Y.t[F].tot += Y.t[F].n[k][j]; }
  }
  Y.truthB = Y.t[0].tot/(Y.t[0].tot + Y.t[1].tot + Y.t[2].tot); Y.rc = Y.t[1].tot/(Y.t[1].tot + Y.t[2].tot);
  p = psimple(Y).p;
  return true;
}
TH1D* book(const TString &n){ TH1D *h = new TH1D(n, "", NW, ptEdges + 3); h->SetDirectory(nullptr); return h; }
// a histogram of the 4-category nominal purity files at our windows (chain binning there)
TH1D* fromChain(TFile *fpp, TFile *fAA, const char *key, int sy, const TString &n, bool withErr){
  TH1 *h = sy == 0 ? (TH1*) fpp->Get(key) : (TH1*) fAA->Get(Form("%s_%s", key, syTag[sy]));
  if(!h) return nullptr;
  TH1D *o = book(n);
  for(int w = 0; w < NW; w++){ const int b = h->FindBin(0.5*(wLo[w] + wHi[w])); o->SetBinContent(w + 1, h->GetBinContent(b)); o->SetBinError(w + 1, withErr ? h->GetBinError(b) : 0.); }
  return o;
}
} // namespace ssurf

void systematics_simpleSurface_PF(TString step = "nominal", int nToy = 200){
  using namespace ssurf;
  setvbuf(stdout, nullptr, _IOLBF, 0);
  setUniform(); nFit = nPR; gsMult = 1.175;
  setJetCollection("PFJets"); initPlotStyle();
  fDpp = TFile::Open(Form("%s/%s", repo, dataName)); fMpp = TFile::Open(Form("%s/%s", repo, mcName)); fDAA = TFile::Open(Form("%s/%s", repo, dataPbPb));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  TFile *f4pp = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2_2D.root", repo)), *f4AA = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_PFJets_PbPb_csvV2_2D.root", repo));
  if(!fDpp || !fMpp || !fDAA || !fS[0] || !fS[1] || !fS[2] || !f4pp || !f4AA){ printf("ERROR: inputs\n"); return; }
  const TString rootDir = Form("%s/rootFiles/systematics/simpleSurface/", repo), figDir = Form("%s/figures/systematics/bPurity/simpleSurface/", repo);
  gSystem->mkdir(rootDir, true); gSystem->mkdir(figDir, true);
  const TString outRoot = rootDir + "simpleSurface_PF_" + step + csvTag() + ".root", outFig = figDir + "simpleSurface_" + step + csvTag() + ".pdf";

  // =================================================================== the chain's purity files
  // "files": bPuritySys_total_PFJets_pp_csvV2_surf.root, bPuritySys_PFJets_PbPb_csvV2_surf.root from the stored steps, in
  // the layout of the 4-category files (cloned from them: chain binning), read with BJET_PURITY_VARIANT=csvV2__surf.
  // 80-300 GeV: this fit's purity and statistical error; buffer bins: the 4-category file's value x the ratio of the two
  // fits in the adjacent window. The light-yield source is stored under the name cMult (h_sysRel_cMult, h_varRel_cMult_<k>),
  // which bJetsPerZSys_PF_PbPb.C reads as the pp-PbPb correlated c / light source; fitRangeHigh and tplStat zero.
  if(step == "files"){
    if(NCSV != 25 || plainL){ printf("ERROR: the purity files are written for the nominal configuration only\n"); return; }
    const char *stp[5] = {"bGS", "lightMult", "JEU", "JER", "fakeMuJet"}, *key[5] = {"bGS", "cMult", "JEU", "JER", "fakeMuJet"};
    const int nVar[5] = {2, 4, 0, 0, 0};
    TFile *fN = TFile::Open(rootDir + "simpleSurface_PF_nominal.root"), *fSt[5];
    for(int q = 0; q < 5; q++){ fSt[q] = TFile::Open(rootDir + "simpleSurface_PF_" + stp[q] + ".root"); if(!fSt[q] || !fN){ printf("ERROR: step %s (or nominal) not run\n", stp[q]); return; } }
    const TString base = Form("%s/rootFiles/systematics/", repo);
    for(int sy = 0; sy < NSY; sy++){
      const bool pp = sy == 0; const TString sfx = pp ? "" : Form("_%s", syTag[sy]); TFile *f4 = pp ? f4pp : f4AA;
      TH1D *p4 = (TH1D*) f4->Get("h_bPurity_nominal" + sfx), *hN = (TH1D*) fN->Get(Form("h_bPurity_%s", syTag[sy]));
      if(!p4 || !hN){ printf("ERROR: inputs %s\n", syTag[sy]); return; }
      auto win = [&](double c){ int k = 0; while(k < NW - 1 && c > wHi[k]) k++; return c < wLo[0] ? 0 : k; };
      auto mk = [&](const TString &n, const char *title){ TH1D *h = (TH1D*) p4->Clone(n + sfx); h->SetDirectory(nullptr); h->Reset(); h->SetTitle(title); return h; };
      TH1D *p = (TH1D*) p4->Clone("h_bPurity_nominal" + sfx); p->SetDirectory(nullptr); p->SetTitle("b purity, simple surface fit");
      TH1D *tot = mk("h_sysRel_total", "total"), *unc = mk("h_sysRel_uncorr", "uncorrelated purity sources (JEU, JER, fakeMuJet)");
      TH1D *zero1 = mk("h_sysRel_fitRangeHigh", "fit range (zero)"), *zero2 = mk("h_sysRel_tplStat", "template statistics (in the fit error)");
      std::vector<TH1D*> out = {zero1, zero2};
      for(int i = 1; i <= p->GetNbinsX(); i++){
        const double c = p->GetBinCenter(i); const int w = win(c); const bool inside = c > wLo[0] && c < wHi[NW-1];
        const double r = hN->GetBinContent(w + 1)/p4->GetBinContent(p4->FindBin(0.5*(wLo[w] + wHi[w])));
        if(inside){ p->SetBinContent(i, hN->GetBinContent(w + 1)); p->SetBinError(i, hN->GetBinError(w + 1)); }
        else { p->SetBinContent(i, p4->GetBinContent(i)*r); p->SetBinError(i, p4->GetBinError(i)*r); }
      }
      for(int q = 0; q < 5; q++){
        TH1D *hs = (TH1D*) fSt[q]->Get(Form("h_sysRel_%s_%s", stp[q], syTag[sy])); if(!hs){ printf("ERROR: %s %s\n", stp[q], syTag[sy]); return; }
        TH1D *o = mk(Form("h_sysRel_%s", key[q]), q == 1 ? "light yield (stored as cMult)" : stp[q]); out.push_back(o);
        for(int i = 1; i <= p->GetNbinsX(); i++){
          const int w = win(p->GetBinCenter(i)); const double v = 0.01*hs->GetBinContent(w + 1);
          o->SetBinContent(i, v); tot->SetBinContent(i, tot->GetBinContent(i) + v*v); if(q >= 2) unc->SetBinContent(i, unc->GetBinContent(i) + v*v);
        }
        for(int k = 0; k < nVar[q]; k++){
          TH1D *hv = (TH1D*) fSt[q]->Get(Form("h_varRel_%s_%d_%s", stp[q], k, syTag[sy])); if(!hv){ printf("ERROR: var %s %d %s\n", stp[q], k, syTag[sy]); return; }
          TH1D *ov = mk(Form("h_varRel_%s_%d", key[q], k), Form("%s variation %d: P_var / P_nom - 1 (signed)", stp[q], k)); out.push_back(ov);
          for(int i = 1; i <= p->GetNbinsX(); i++) ov->SetBinContent(i, 0.01*hv->GetBinContent(win(p->GetBinCenter(i)) + 1));
        }
      }
      for(int i = 1; i <= p->GetNbinsX(); i++){ tot->SetBinContent(i, sqrt(tot->GetBinContent(i))); unc->SetBinContent(i, sqrt(unc->GetBinContent(i))); }
      TH1D *abs = (TH1D*) tot->Clone("h_sysAbs_total" + sfx), *pa = (TH1D*) p->Clone("h_bPurity_sysAbs" + sfx);
      for(int i = 1; i <= p->GetNbinsX(); i++){ const double a = tot->GetBinContent(i)*p->GetBinContent(i); abs->SetBinContent(i, a); pa->SetBinError(i, a); }
      const TString outName = pp ? base + "bPuritySys_total_PFJets_pp_csvV2_surf.root" : base + "bPuritySys_PFJets_PbPb_csvV2_surf.root";
      TFile *fo = TFile::Open(outName, pp || sy == 1 ? "RECREATE" : "UPDATE");
      p->Write(); for(TH1D *h : out) h->Write(); tot->Write(); abs->Write(); pa->Write(); unc->Write();
      if(pp || sy == 1){ TNamed info("info", "PF b purity from the SIMPLE SURFACE FIT (25 CSVv2 x 25 ptRel cells, f_b free, c / light at MC, Poisson + Barlow-Beeston lite), 80-300 GeV; buffer bins = the 4-category file x the adjacent-window ratio. Sources from systematics_simpleSurface_PF.C: bGS, light yield x 0.8-1.2 (stored as cMult), JEU, JER, fakeMuJet; fitRangeHigh and tplStat zero. STUDY, not nominal"); info.Write(); }
      fo->Close();
      printf("wrote %s (%s):", outName.Data(), syTag[sy]); for(int i = 1; i <= p->GetNbinsX(); i++) printf("  %.0f-%.0f %.3f (%.1f%%)", p->GetBinLowEdge(i), p->GetBinLowEdge(i + 1), p->GetBinContent(i), 100.*tot->GetBinContent(i)); printf("\n");
    }
    return;
  }

  // =================================================================== nominal
  if(step == "nominal"){
    TH1D *hP[NSY], *hC[NSY], *h4[NSY];
    for(int sy = 0; sy < NSY; sy++){
      hP[sy] = book(Form("h_bPurity_%s", syTag[sy])); hC[sy] = book(Form("h_pull2_%s", syTag[sy]));
      h4[sy] = fromChain(f4pp, f4AA, "h_bPurity_nominal", sy, Form("h_bPurity4cat_%s", syTag[sy]), true);
      if(!h4[sy]){ printf("ERROR: 4-category nominal %s\n", syTag[sy]); return; }
    }
    for(int ix = 0; ix < NSY; ix++){
      const int sy = order[ix];
      printf("\n==== %s: simple surface fit (%d x %d cells, f_b free, c / light at MC, %s)\n", syName[sy], K, nPR, plainL ? "plain Poisson" : "Poisson + template statistics");
      printf("  window   | P (stat)        st | 4-cat HistFactory | sum pull^2 / cells (fit >= 5) | data in cells without MC | c / (c + light) MC | Asimov fit - truth\n");
      for(int w = 0; w < NW; w++){
        static Win X; if(!build(sy, w, X)) return;
        const Fit f = psimple(X, true);
        hP[sy]->SetBinContent(w + 1, f.p); hP[sy]->SetBinError(w + 1, f.e); hC[sy]->SetBinContent(w + 1, f.nChi > 0 ? f.chi2/f.nChi : 0.);
        // Asimov: the templates at their true mixture and the data's size
        double one[NF][KMAX]; for(int F = 0; F < NF; F++) for(int k = 0; k < K; k++) one[F][k] = 1.;
        static Win A; A = X; double sd = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) sd += X.d[k][j];
        A.nData = sd; const double tr = asimov(A, one, A.d); const Fit fa = psimple(A);
        printf("  %3.0f-%3.0f  | %.3f +- %.3f (%d) | %.3f +- %.3f    | %7.1f / %-4d = %.2f         | %.2f%%                    | %.3f              | %+.4f\n", wLo[w], wHi[w], f.p, f.e, f.status,
               h4[sy]->GetBinContent(w + 1), h4[sy]->GetBinError(w + 1), f.chi2, f.nChi, f.nChi > 0 ? f.chi2/f.nChi : 0., 100.*f.noMC, X.rc, fa.fb - tr);
      }
    }
    TFile *fo = TFile::Open(outRoot, "RECREATE"); for(int sy = 0; sy < NSY; sy++){ hP[sy]->Write(); hC[sy]->Write(); h4[sy]->Write(); }
    TNamed info("info", Form("simple surface fit of the b purity: %d CSVv2 x %d ptRel cells, f_b free, c / light at MC, %s; pp, C1..C4 = 0-10 ... 50-80%%", K, nPR, plainL ? "plain Poisson" : "Poisson + Barlow-Beeston lite")); info.Write(); fo->Close();

    TCanvas *cv = new TCanvas("cvN", "", 2300, 1000);
    GridPads g = makeGridPads(cv, NSY, 2, 0.055, 0.01, 0.12, 0.075);
    const double split = g.bm + 0.36*(1. - g.bm - g.tm);
    for(int c = 0; c < NSY; c++){
      TPad *t = g.at(c, 0), *b = g.at(c, 1); double x1, y1, x2, y2;
      t->GetPadPar(x1, y1, x2, y2); t->SetPad(x1, split, x2, y2); t->SetTopMargin(g.tm/(1. - split));
      b->GetPadPar(x1, y1, x2, y2); b->SetPad(x1, 0., x2, split); b->SetBottomMargin(g.bm/split);
    }
    double cMax = 1.5; for(int sy = 0; sy < NSY; sy++) cMax = TMath::Max(cMax, hC[sy]->GetMaximum());
    for(int ix = 0; ix < NSY; ix++){
      const int sy = order[ix];
      g.at(ix, 0)->cd();
      TH1F *fr = gPad->DrawFrame(80., 0.2, 300., 0.8); g.styleFrame(fr, ix, 0, 0.055);
      fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(ix == 0 ? 34 : 0); fr->GetYaxis()->ChangeLabel(1, -1, 0.); fr->GetYaxis()->SetNdivisions(507);
      styleH(h4[sy], hexData, markOpenCircle, 3.4); h4[sy]->SetLineWidth(3); h4[sy]->Draw("E1 X0 SAME");
      styleH(hP[sy], okabeHex[5], markFilledCircle, 2.0); hP[sy]->SetLineWidth(3); hP[sy]->Draw("E1 X0 SAME");
      gPad->RedrawAxis();
      g.text(ix, 0, 0.10, syName[sy], 0.065, 62);
      if(ix == 0){
        const double fL = gPad->GetLeftMargin(), fW = 1. - fL - gPad->GetRightMargin(), fT = 1. - gPad->GetTopMargin();
        TLegend *lg = new TLegend(fL + 0.04*fW, fT - 0.36, fL + fW, fT - 0.14);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(43); lg->SetTextSize(30);
        lg->AddEntry(h4[sy], "4 categories (nominal)", "pe"); lg->AddEntry(hP[sy], "simple surface fit", "pe"); lg->Draw();
      }
      g.at(ix, 1)->cd();
      TH1F *fb = gPad->DrawFrame(80., 0., 300., 1.15*cMax); g.styleFrame(fb, ix, 1, 0.055);
      fb->GetXaxis()->SetLabelFont(43); fb->GetXaxis()->SetLabelSize(34); fb->GetYaxis()->SetLabelFont(43); fb->GetYaxis()->SetLabelSize(ix == 0 ? 34 : 0);
      fb->GetYaxis()->SetNdivisions(505); fb->GetXaxis()->SetTickLength(0.03*0.64/0.36);
      TLine *one = new TLine(80., 1., 300., 1.); one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->SetLineWidth(2); one->Draw();
      styleH(hC[sy], okabeHex[5], markFilledCircle, 2.0); hC[sy]->Draw("P X0 SAME");
      gPad->RedrawAxis();
    }
    cv->cd();
    { TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.040); t.SetTextAlign(21);
      t.DrawLatex(g.lm + (1. - g.lm - g.rm)/2., 0.018, "#it{p}_{T}^{jet} (GeV)");
      t.SetTextAngle(90); t.DrawLatex(0.016, split + 0.5*(1. - g.tm - split), "b purity of muon-tagged jets");
      t.SetTextSize(0.030); t.DrawLatex(0.016, g.bm + 0.5*(split - g.bm), "#Sigma pull^{2} / cells"); }
    g.cmsLabel(cv, "Internal", Form("pp, PbPb (5.02 TeV)   #it{p}_{T}^{rel} #times CSVv2, %d #times %d cells", K, nPR), 1.0);
    savePdfTight(cv, outFig);
    return;
  }

  // =================================================================== one source
  std::vector<TString> vLab; const char *oldKey = nullptr; TString title; bool toys = false;
  if(step == "bGS"){ vLab = {"g#rightarrowb#bar{b} share #times 0.94", "#times 1.41"}; oldKey = "h_sysRel_bGS"; title = "g#rightarrowb#bar{b} share of the b template"; }
  else if(step == "lightMult"){ vLab = {"light #times 0.8", "#times 0.9", "#times 1.1", "#times 1.2"}; oldKey = "h_sysRel_cMult"; title = "light yield"; }
  else if(step == "JEU"){ vLab = {"JES up", "JES down"}; oldKey = "h_sysRel_JEU"; title = "jet energy scale"; }
  else if(step == "JER"){ vLab = {"JER smeared"}; oldKey = "h_sysRel_JER"; title = "jet energy resolution"; }
  else if(step == "fakeMuJet"){ vLab = {"no pair subtraction"}; oldKey = "h_sysRel_fakeMuJet"; title = "uncorrelated #mu-jet pairs"; }
  else if(step == "tplStat"){ if(!plainL){ printf("tplStat: the template statistics are in the fit error of the nominal likelihood; this step is for BJET_SS_LIKELIHOOD=plain\n"); return; } vLab = {"toy RMS"}; title = "template statistics"; toys = true; }
  else { printf("ERROR: step %s unknown (nominal, bGS, lightMult, JEU, JER, fakeMuJet, tplStat)\n", step.Data()); return; }
  TFile *fNom = TFile::Open(rootDir + "simpleSurface_PF_nominal" + csvTag() + ".root");
  if(!fNom){ printf("ERROR: run the nominal step first\n"); return; }
  const int NV = (int) vLab.size();
  TH1D *hV[4][NSY], *hS[NSY], *hOld[NSY];
  TRandom3 rng(20261116);
  for(int ix = 0; ix < NSY; ix++){
    const int sy = order[ix];
    for(int v = 0; v < NV; v++) hV[v][sy] = book(Form("h_varRel_%s_%d_%s", step.Data(), v, syTag[sy]));
    hS[sy] = book(Form("h_sysRel_%s_%s", step.Data(), syTag[sy]));
    hOld[sy] = oldKey ? fromChain(f4pp, f4AA, oldKey, sy, Form("h_sysRel4cat_%s_%s", step.Data(), syTag[sy]), false) : nullptr;
    TH1D *hN = (TH1D*) fNom->Get(Form("h_bPurity_%s", syTag[sy])); if(!hN){ printf("ERROR: nominal %s\n", syTag[sy]); return; }
    printf("\n==== %s: %s, P_var / P_nom - 1 (%%)\n  window   | P_nom  |", syName[sy], step.Data());
    for(int v = 0; v < NV; v++) printf(" %-22s|", TString(vLab[v]).ReplaceAll("#rightarrow", "->").ReplaceAll("#bar{b}", "b").ReplaceAll("#times", "x").Data());
    printf(" assigned | 4-category nominal\n");
    for(int w = 0; w < NW; w++){
      static Win X; gsMult = 1.175; if(!build(sy, w, X)) return;
      const double p0 = psimple(X).p;
      if(fabs(p0 - hN->GetBinContent(w + 1)) > 1e-6) printf("  WARNING: nominal %.4f differs from the stored %.4f\n", p0, hN->GetBinContent(w + 1));
      double pv[4] = {p0, p0, p0, p0}; bool skip = false;
      if(step == "bGS"){ int v = 0; for(double gm : {0.94, 1.41}){ gsMult = gm; static Win Y; if(!build(sy, w, Y)) return; pv[v++] = psimple(Y).p; } gsMult = 1.175; }
      else if(step == "lightMult"){ int v = 0; for(double m : {0.8, 0.9, 1.1, 1.2}){ static Win Y; Y = X; scaleFl(Y, 2, m); pv[v++] = psimple(Y).p; } }
      else if(step == "JEU"){ int v = 0; for(int T : {2, 3}) if(!shiftedSimple(sy, w, X, T, pv[v++])) return; }
      else if(step == "JER"){ if(!shiftedSimple(sy, w, X, 1, pv[0])) return; }
      else if(step == "fakeMuJet"){
        if(sy == 0) skip = true;
        else { fakeSubPerEvt = false; static Win Y; const bool ok = build(sy, w, Y); fakeSubPerEvt = true; if(!ok) return; pv[0] = psimple(Y).p; }
      }
      double rms = 0.;
      if(toys){
        double s1 = 0., s2 = 0.;
        for(int t = 0; t < nToy; t++){
          static Win Y; Y = X;
          for(int F = 0; F < NF; F++){
            Y.t[F].tot = 0.;
            for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++){ Y.t[F].n[k][j] = TMath::Max(0., rng.Gaus(X.t[F].n[k][j], sqrt(X.t[F].e2[k][j]))); Y.t[F].tot += Y.t[F].n[k][j]; }
          }
          Y.truthB = Y.t[0].tot/(Y.t[0].tot + Y.t[1].tot + Y.t[2].tot); Y.rc = Y.t[1].tot/(Y.t[1].tot + Y.t[2].tot);
          const double p = psimple(Y).p/p0 - 1.; s1 += p; s2 += p*p;
        }
        const double mean = s1/nToy; rms = sqrt(TMath::Max(0., s2/nToy - mean*mean));
        pv[0] = p0*(1. + rms);
        printf("  %3.0f-%3.0f  | %.3f  | RMS %5.2f, mean %+5.2f     |", wLo[w], wHi[w], p0, 100.*rms, 100.*mean);
      } else {
        printf("  %3.0f-%3.0f  | %.3f  |", wLo[w], wHi[w], p0);
        for(int v = 0; v < NV; v++) printf(skip ? "          -            |" : " %+6.2f                |", 100.*(pv[v]/p0 - 1.));
      }
      double a = 0.; for(int v = 0; v < NV; v++){ const double r = skip ? 0. : pv[v]/p0 - 1.; hV[v][sy]->SetBinContent(w + 1, 100.*r); a = TMath::Max(a, fabs(r)); }
      hS[sy]->SetBinContent(w + 1, 100.*a);
      printf(" %5.2f    | %s\n", 100.*a, hOld[sy] ? Form("%.2f", 100.*hOld[sy]->GetBinContent(w + 1)) : "-");
    }
  }
  TFile *fo = TFile::Open(outRoot, "RECREATE");
  for(int sy = 0; sy < NSY; sy++){ for(int v = 0; v < NV; v++) hV[v][sy]->Write(); hS[sy]->Write(); if(hOld[sy]) hOld[sy]->Write(); }
  fo->Close();

  // ---- figure: signed variations per system, the assigned magnitude as a band, the 4-category source dashed
  double yM = 1.; for(int sy = 0; sy < NSY; sy++) for(int w = 1; w <= NW; w++){ yM = TMath::Max(yM, hS[sy]->GetBinContent(w)); }   // this fit's own variations set the scale; a larger 4-category value is written in the panel
  const int lgCols = NV > 2 ? 2 : 1, lgRows = (NV + 2 + lgCols - 1)/lgCols; const double lgH = 0.055*lgRows;
  yM /= TMath::Max(0.25, 1. - 2.*(lgH + 0.05)/0.74);   // room below the points for the legend
  const char *vHex[4] = {okabeHex[5], okabeHex[6], okabeHex[1], okabeHex[3]}; const int vMark[4] = {markFilledCircle, markFilledSquare, markFilledDiamond, markOpenCircle}; const double vSize[4] = {2.2, 2.2, 3.0, 2.4};
  TCanvas *cv = new TCanvas("cvV", "", 2300, 760);
  GridPads g = makeGridPads(cv, NSY, 1, 0.06, 0.01, 0.16, 0.10);
  for(int ix = 0; ix < NSY; ix++){
    const int sy = order[ix];
    g.at(ix, 0)->cd();
    TH1F *fr = gPad->DrawFrame(80., -yM, 300., yM); g.styleFrame(fr, ix, 0, 0.055);
    fr->GetXaxis()->SetLabelFont(43); fr->GetXaxis()->SetLabelSize(34); fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(ix == 0 ? 34 : 0);
    TH1D *up = (TH1D*) hS[sy]->Clone(Form("bandUp%d", sy)), *dn = (TH1D*) hS[sy]->Clone(Form("bandDn%d", sy)); dn->Scale(-1.);
    for(TH1D *h : {up, dn}){ h->SetDirectory(nullptr); h->SetFillColorAlpha(kGray, 0.55); h->SetLineWidth(0); h->Draw("HIST SAME"); }
    TLine *z = new TLine(80., 0., 300., 0.); z->SetLineColor(kGray + 2); z->SetLineWidth(2); z->Draw();
    TH1D *oU = nullptr;
    if(hOld[sy]){
      oU = (TH1D*) hOld[sy]->Clone(Form("oldUp%d", sy)); oU->Scale(100.); TH1D *oD = (TH1D*) oU->Clone(Form("oldDn%d", sy)); oD->Scale(-1.);
      for(TH1D *h : {oU, oD}){ h->SetDirectory(nullptr); styleLine(h, "#555555", 3); h->SetLineStyle(7); h->Draw("HIST SAME"); }
    }
    if(!toys) for(int v = 0; v < NV; v++){ styleH(hV[v][sy], vHex[v], vMark[v], vSize[v]); hV[v][sy]->SetLineWidth(3); hV[v][sy]->Draw("P X0 SAME"); }
    gPad->RedrawAxis();
    g.text(ix, 0, 0.09, syName[sy], 0.065, 62);
    if(hOld[sy] && 100.*hOld[sy]->GetMaximum() > yM) g.text(ix, 0, 0.17, Form("4-category: up to %.1f%%", 100.*hOld[sy]->GetMaximum()), 0.045);
    if(ix == 0){
      const double fL = gPad->GetLeftMargin(), fW = 1. - fL - gPad->GetRightMargin(), fB = gPad->GetBottomMargin();
      TLegend *lg = new TLegend(fL + 0.04*fW, fB + 0.03, fL + fW, fB + 0.03 + lgH); lg->SetNColumns(lgCols); lg->SetColumnSeparation(0.02);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(43); lg->SetTextSize(lgCols > 1 ? 22 : 28);
      if(!toys) for(int v = 0; v < NV; v++) lg->AddEntry(hV[v][sy], vLab[v], "p");
      lg->AddEntry(up, toys ? "toy RMS" : "assigned", "f");
      if(oU) lg->AddEntry(oU, "4-category nominal", "l");
      lg->Draw();
    }
  }
  g.titles(cv, "#it{p}_{T}^{jet} (GeV)", "change of the b purity (%)", 0.05, 0.018, 0.02);
  g.cmsLabel(cv, "Internal", Form("simple surface fit:  %s", title.Data()), 1.4);
  savePdfTight(cv, outFig);
}
