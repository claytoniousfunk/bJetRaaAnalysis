// Optimal Bayes iterations of the PF b-jet unfolding with the MC PRIOR, the procedure of the
// inclusive jets (2026-10-21, user: Bayes nominal for the b jets again, "like we do for the
// inclusive jets"). Copy of optimalIterations_bJets_PF.C (data-driven prior, data-like closure on
// even/odd halves), changed:
//   prior    none: the response is used as it is (MC prior), as the chain now does
//            (bUseDataDrivenPrior = false, bJetUnfoldExt.h)
//   closure  same-sample, as the inclusive PbPb closure: pseudo-data = the chain's own response
//            with its gen rows tilted by (g/100)^-k, response = the untilted matrix; the pseudo
//            measurement carries the real data's relative errors, so the variance is the data's
//   N        the minimum of MSE(N) = mean over 80-200 GeV of (U/T - 1)^2 + (sigma_U/T)^2 at the
//            inclusive closure's tilt, k = 0.3
// It also prints, at that N, the closure for a DATA-LIKE tilt -- PbPb: the data-driven prior's
// fitted k of the chain file (priorTilt_<cls>, about 1 in 0-10%); pp: the same fit in
// unfoldBJetSpectrum_caloJets_pp.C, entered below -- because with the MC prior a steep data
// spectrum is where the Bayes unfolding at small N is biased (dataDrivenPrior.h header).
// Responses: PbPb the chain's (phRespPath, PFJets: sibling 2026-02-18 muTaggedJetsNoTrigger, nominal
// classes), pp the chain's (ppRespFiles), outliers clipped and (since 2026-11-05) JER-smeared x1.10 as
// in the chain (jerSmear.h; BJET_JER_SCHEME=old for the unsmeared response).
//
// Output: rootFiles/unfoldTest/optimalIterations_bJets_PF_MCprior.root (g_{mse,bias2,var}_<sys>, optimalN,
//         h_uOverT_<sys> = unfolded / truth at N_opt on the reported bins, since 2026-11-06),
//         figures/unfoldTest/optimalIterations_bJets_PF_MCprior.pdf
// Usage, from src/unfoldTest/:
//   root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");
//                     gInterpreter->AddIncludePath("/home/clayton/Programs/RooUnfold/build");' optimalIterations_bJets_PF_MCprior.C
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TVectorD.h"
#include "../../headers/functions/jetCollectionInputs.h"
#include "../../headers/plotting/coarseCent.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/functions/dataDrivenPrior.h"
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "../../headers/functions/bJetUnfoldExt.h"
#include "../../headers/functions/clipResponseOutliers.h"
#include "../../headers/functions/jerSmear.h"   // JER nominal: the response smeared x1.10, as the chain (2026-11-05)
#include "../../headers/functions/reverseUnfoldReweightCore.h"   // rvu: edges, recoOnEdges, genOnEdges

namespace oimc {
const char *chainPath = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/BJetPerZ_PFJets/backup_2026-10-21_preLeakFit/bJetChain_PFJets_PbPb_dedicatedPurity.root";   // data-driven-prior chain: bJetPt_meas_<C>, priorTilt_<C>
const char *ppSpecPath= "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_PFJets_pp.root";
const char *outPath   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/unfoldTest/optimalIterations_bJets_PF_MCprior.root";
const char *figPath   = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/unfoldTest/optimalIterations_bJets_PF_MCprior.pdf";
int NMax = 10;   // scan range; 3rd argument of the macro (2026-11-16)
const double floorB = bRecoFloorPF(), kClosure = 0.3;   // 60 since 2026-11-05 (was 70)
double kDataPP = 0.;   // set from the pp chain log if known (data-driven prior fit); 0 -> skipped

TH1D* pseudoMeas(const TH2D *P, TH1D *relSrc, const char *name){
  TH1D *pd = rvu::recoOnEdges(P, name);
  for(int b = 1; b < rvu::NEdge; b++){
    int j = relSrc->FindBin(pd->GetBinCenter(b)); double m = relSrc->GetBinContent(j);
    pd->SetBinError(b, m > 0. ? pd->GetBinContent(b)*relSrc->GetBinError(j)/m : 0.);
  }
  return pd;
}
struct Row { double mse = 0., bias2 = 0., var = 0.; std::vector<double> r, er, lo, hi; };
// unfold pd with the UNTILTED response R (MC prior), against truth tt
Row closeAt(const TH2D *R, TH1D *pd, TH1D *tt, int N, const char *tag){
  const std::vector<double> extE = bExtEdges(rvu::NEdge, rvu::edges);
  TH2D *Re = rebin2DTo(const_cast<TH2D*>(R), extE, "oiRe");
  TH1D *U = unfoldOnExt<RooUnfoldResponse, RooUnfoldBayes>(Re, pd, rvu::NEdge, rvu::edges, floorB, N, Form("oiU_%s_%d", tag, N));
  Row o; int nb = 0;
  for(int b = 1; b < rvu::NEdge; b++){
    if(!bReportedBin(rvu::edges[b-1], rvu::edges[b])) continue;
    const double t = tt->GetBinContent(b), u = U->GetBinContent(b), e = U->GetBinError(b);
    if(t <= 0.) continue;
    const double r = u/t, er = e/t;
    o.r.push_back(r); o.er.push_back(er); o.lo.push_back(rvu::edges[b-1]); o.hi.push_back(rvu::edges[b]); o.bias2 += (r - 1.)*(r - 1.); o.var += er*er; nb++;
  }
  if(nb){ o.bias2 /= nb; o.var /= nb; } o.mse = o.bias2 + o.var;
  delete Re; delete U;
  return o;
}
} // namespace oimc
using namespace oimc;

// inTag (2026-11-16): "" = the inputs above (the nominal derivation). Otherwise the data's relative errors come from
// the chain run of that tag, bJetChain_PFJets_PbPb_dedicatedPurity<inTag>.root and bJetSpectrum_PFJets_pp<inTag>.root,
// and the outputs carry <inTag>. With BJET_PT_7BINS defined (7-bin axis, reverseUnfoldReweightCore.h edges) the tag
// must name it, e.g. "_csvV2_2D7bins"; "_csvV2_2D" without the symbol is the same-inputs 4-bin reference.
// nMax (2026-11-16, user: "exploratory: extend the unfolding out to 50 iterations"): the scan range, default 10. Another
// value appends "_N<nMax>" to the output names (with an empty inTag too), so the nominal 10-iteration files stay as they are.
void optimalIterations_bJets_PF_MCprior(double kPP = 0., TString inTag = "", int nMax = 10){
  kDataPP = kPP; NMax = nMax;
  if(!ptBinGuard(inTag)) return;
  if(inTag != ""){
    const char *b = "/home/clayton/Analysis/code/bJetRaaAnalysis";
    chainPath = strdup(Form("%s/rootFiles/BJetPerZ_PFJets/bJetChain_PFJets_PbPb_dedicatedPurity%s.root", b, inTag.Data()));
    ppSpecPath = strdup(Form("%s/rootFiles/MuTagMuTrigBJetSpectra/Data/bJetSpectrum_PFJets_pp%s.root", b, inTag.Data()));
    outPath = strdup(Form("%s/rootFiles/unfoldTest/optimalIterations_bJets_PF_MCprior%s.root", b, inTag.Data()));
    figPath = strdup(Form("%s/figures/unfoldTest/optimalIterations_bJets_PF_MCprior%s.pdf", b, inTag.Data()));
  }
  if(nMax != 10){ TString o = outPath, f = figPath; o.ReplaceAll(".root", Form("_N%d.root", nMax)); f.ReplaceAll(".pdf", Form("_N%d.pdf", nMax)); outPath = strdup(o.Data()); figPath = strdup(f.Data()); }
  setJetCollection("PFJets");
  initPlotStyle();
  TFile *fR = TFile::Open(phRespPath), *fC = TFile::Open(chainPath), *fS = TFile::Open(ppSpecPath);
  for(TFile *f : {fR, fC, fS}) if(!f || f->IsZombie()){ printf("ERROR: an input did not open\n"); return; }
  gSystem->mkdir(gSystem->DirName(outPath), true); gSystem->mkdir(gSystem->DirName(figPath), true);
  TFile *out = TFile::Open(outPath, "RECREATE");

  const char *sysName[5] = {"pp", "0-10%", "10-30%", "30-50%", "50-80%"};
  const char *sysTag[5]  = {"pp", "C1", "C2", "C3", "C4"};
  TGraph *gM[5], *gB[5], *gV[5]; int best[5];
  for(int s = 0; s < 5; s++){
    TH2D *R = nullptr; TH1D *rel = nullptr; double kData = 0.;
    if(s == 0){
      for(const std::string &fn : ppRespFiles){ TFile *f = TFile::Open(Form("%s%s", ppRespDir, fn.c_str())); TH2D *r = nullptr;
        if(f && !f->IsZombie()) f->GetObject(ppRespName, r);
        if(!r){ printf("ERROR: %s missing\n", ppRespName); return; }
        if(!R){ R = (TH2D*) r->Clone("Rpp"); R->SetDirectory(nullptr); } else R->Add(r); }
      fS->GetObject("h_bJetPt", rel);
      kData = kDataPP;
    } else {
      TH2D *r = nullptr; fR->GetObject(Form(phRespMuTagFmt, s), r);
      if(!r){ printf("ERROR: response missing for %s\n", sysTag[s]); return; }
      R = (TH2D*) r->Clone(Form("R_%s", sysTag[s])); R->SetDirectory(nullptr);
      fC->GetObject(Form("bJetPt_meas_%s", sysTag[s]), rel);
      TVectorD *kv = nullptr; fC->GetObject(Form("priorTilt_%s", sysTag[s]), kv); if(kv) kData = (*kv)[0];
    }
    if(!R || !rel){ printf("ERROR: inputs missing for %s\n", sysTag[s]); return; }
    clipResponseOutliers(R, 100., false, sysTag[s]);
    { TH2D *sm = jerSmear::nominal(R, Form("%s_jerNom", R->GetName())); delete R; R = sm; }   // as the chain (jerSmear.h)
    // tilt 0.3 pseudo-data: choose N
    TH2D *P = closureTiltRows(R, kClosure, Form("P_%s", sysTag[s]));   // clamped, harder (2026-11-05, user)
    TH1D *pd = pseudoMeas(P, rel, Form("pd_%s", sysTag[s])), *tt = rvu::genOnEdges(P, Form("tt_%s", sysTag[s]));
    gM[s] = new TGraph(); gB[s] = new TGraph(); gV[s] = new TGraph();
    printf("\n===== %s =====  tilt %.1f closure, unfolded / true in the reported bins 80-200 GeV\n", sysName[s], kClosure);
    double mMin = 1e9; best[s] = 1; std::vector<double> mseV; int nMin = 1;
    std::vector<Row> rows;   // every N, for the within-x%-of-truth criterion below
    for(int N = 1; N <= NMax; N++){
      Row r = closeAt(R, pd, tt, N, sysTag[s]);
      gM[s]->SetPoint(N-1, N, r.mse); gB[s]->SetPoint(N-1, N, r.bias2); gV[s]->SetPoint(N-1, N, r.var);
      if(N <= 10 || N % 5 == 0){ printf("  N = %2d  MSE %.5f (bias^2 %.5f, var %.5f) |", N, r.mse, r.bias2, r.var);
      for(size_t i = 0; i < r.r.size(); i++) printf("  %.3f+-%.3f", r.r[i], r.er[i]);
      printf("\n"); }
      if(r.mse < mMin){ mMin = r.mse; nMin = N; }
      mseV.push_back(r.mse); rows.push_back(r);
    }
    best[s] = plateauN(mseV);   // 5% plateau rule (2026-11-05, user)
    printf("  -> optimal N = %d (5%% plateau; absolute minimum at N = %d)\n", best[s], nMin);
    // BIAS CRITERION (2026-11-16, user: "what if we add a criteria that the unfolding has to get to within 2% of truth"):
    // the smallest N at which EVERY reported bin has |unfolded / true - 1| <= x, for x = 1, 2, 3, 5%, with what it costs:
    // the largest per-bin relative stat. error there against the one at the MSE optimum. Also the smallest largest
    // deviation reached anywhere in the scan. Stored: g_maxDev_<sys> (largest |u/t - 1| vs N), g_maxErr_<sys>.
    { TGraph *gD = new TGraph(), *gE = new TGraph(); double dBest = 1e9; int nBest = 0;
      auto mx = [](const Row &r, bool err){ double m = 0.; for(size_t i = 0; i < r.r.size(); i++) m = TMath::Max(m, err ? r.er[i] : fabs(r.r[i] - 1.)); return m; };
      for(size_t i = 0; i < rows.size(); i++){ const double d = mx(rows[i], false); gD->SetPoint(i, i + 1, d); gE->SetPoint(i, i + 1, mx(rows[i], true)); if(d < dBest){ dBest = d; nBest = i + 1; } }
      printf("  largest |u/t - 1| over the reported bins: %.1f%% at the MSE optimum (stat. up to %.1f%%); smallest in the scan %.1f%% at N = %d (stat. up to %.1f%%)\n",
             100.*mx(rows[best[s]-1], false), 100.*mx(rows[best[s]-1], true), 100.*dBest, nBest, 100.*mx(rows[nBest-1], true));
      for(double x : {0.01, 0.02, 0.03, 0.05}){ int n = 0; for(size_t i = 0; i < rows.size() && !n; i++) if(mx(rows[i], false) <= x) n = i + 1;
        if(n) printf("    within %.0f%% of truth in every bin: first at N = %2d, stat. up to %.1f%% (x %.2f of the MSE optimum's), MSE x %.2f\n", 100.*x, n, 100.*mx(rows[n-1], true), mx(rows[n-1], true)/mx(rows[best[s]-1], true), rows[n-1].mse/rows[best[s]-1].mse);
        else  printf("    within %.0f%% of truth in every bin: NOT reached up to N = %d\n", 100.*x, NMax); }
      out->cd(); gD->Write(Form("g_maxDev_%s", sysTag[s])); gE->Write(Form("g_maxErr_%s", sysTag[s])); }
    { // unfolded / truth on the reported bins at N_opt, for the AN closure figure (2026-11-06)
      Row r = closeAt(R, pd, tt, best[s], Form("opt%s", sysTag[s]));
      std::vector<double> e(r.lo); if(!r.hi.empty()) e.push_back(r.hi.back());
      TH1D *h = new TH1D(Form("h_uOverT_%s", sysTag[s]), Form("%s;p_{T};unfolded / truth, N_opt = %d", sysName[s], best[s]), (int) r.r.size(), e.data());
      for(size_t i = 0; i < r.r.size(); i++){ h->SetBinContent(i+1, r.r[i]); h->SetBinError(i+1, r.er[i]); }
      out->cd(); h->Write(); delete h;
    }
    // data-like tilt at the optimum and its neighbours
    if(kData != 0.){
      TH2D *Pd = tiltGenRows(R, kData, Form("Pd_%s", sysTag[s]));
      TH1D *pdd = pseudoMeas(Pd, rel, Form("pdd_%s", sysTag[s])), *ttd = rvu::genOnEdges(Pd, Form("ttd_%s", sysTag[s]));
      for(int N : {TMath::Max(1, best[s] - 1), best[s], best[s] + 1, best[s] + 3}){
        Row r = closeAt(R, pdd, ttd, N, Form("d%s", sysTag[s]));
        printf("  data-like tilt k = %.2f, N = %2d: unfolded / true", kData, N);
        for(size_t i = 0; i < r.r.size(); i++) printf("  %.3f", r.r[i]);
        printf("\n");
      }
      delete Pd; delete pdd; delete ttd;
    }
    out->cd(); gM[s]->Write(Form("g_mse_%s", sysTag[s])); gB[s]->Write(Form("g_bias2_%s", sysTag[s])); gV[s]->Write(Form("g_var_%s", sysTag[s]));
    delete R; delete P; delete pd; delete tt;
  }
  { TVectorD v(5); for(int s = 0; s < 5; s++) v[s] = best[s]; out->cd(); v.Write("optimalN"); }
  out->Close();
  printf("\noptimal N (MC prior, tilt %.1f): pp %d, PbPb %d / %d / %d / %d\n", kClosure, best[0], best[1], best[2], best[3], best[4]);

  // figure: MSE, bias^2, variance vs N, peripheral to central
  const int ord[5] = {0, 4, 3, 2, 1};
  TCanvas *cv = new TCanvas("coimc", "", 5*380, 380); cv->Divide(5, 1, 0.001, 0.001);
  for(int p = 0; p < 5; p++){
    const int s = ord[p];
    cv->cd(p+1); gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.10);
    double ymax = 0.; for(int i = 0; i < gM[s]->GetN(); i++) ymax = std::max(ymax, gM[s]->GetY()[i]);
    TH1F *fr = gPad->DrawFrame(0.5, 0., NMax + 0.5, 1.3*ymax);
    fr->GetXaxis()->SetTitle("Bayes iterations #it{N}"); fr->GetYaxis()->SetTitle("mean over 80-200 GeV");
    fr->GetXaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleSize(0.055); fr->GetYaxis()->SetTitleOffset(1.5);
    fr->GetXaxis()->SetTickLength(0.015); fr->GetYaxis()->SetTickLength(0.012);
    auto st = [](TGraph *g, const char *hex, int m){ int c = TColor::GetColor(hex); g->SetLineColor(c); g->SetMarkerColor(c); g->SetMarkerStyle(m); g->SetLineWidth(2); };
    st(gM[s], hexData, markFilledCircle); st(gB[s], okabeHex[5], markOpenSquare); st(gV[s], okabeHex[1], markOpenDiamond);
    gB[s]->Draw("LP same"); gV[s]->Draw("LP same"); gM[s]->Draw("LP same");
    if(p == 0){ TLegend *lg = makeLegend(0.45, 0.66, 0.95, 0.86, 0.050); lg->AddEntry(gM[s], "MSE", "lp"); lg->AddEntry(gB[s], "bias^{2}", "lp"); lg->AddEntry(gV[s], "variance", "lp"); lg->Draw(); }
    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.055);
    la.DrawLatex(0.22, 0.92, Form("%s, optimal #it{N} = %d", sysName[s], best[s]));
  }
  savePdfTight(cv, figPath); delete cv;
  printf("wrote %s and %s\n", outPath, figPath);
}
