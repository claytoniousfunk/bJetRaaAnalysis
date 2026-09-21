// b-jet fraction template fits in pp with calo jets, with and without the
// tagging muon's 4-momentum added to the jet.
//
// FIT: the pp branch of templateFitter() in calculateBJetsPerZ.cc, reproduced
// step for step in fitOneBin(), with calculateBPurity()'s settings: 2-template
// fit (b vs light+c) with c folded into light at its MC-truth fraction,
// gluon-splitting enhancement bGS_multiplier = 1.175, ptRel rebinned to
// 0-5 GeV in 0.1 GeV steps, RooFit fb over 0-5 GeV with SumW2Error.
//
// DELIBERATE DIFFERENCES from templateFitter():
//  1. The jet pT window is found on EACH histogram's own axis. templateFitter()
//     finds it on the data axis and reuses the bin numbers on the MC templates,
//     which is only right when both axes are identical. They are for the
//     2026-02 golden files (96 bins, 20-500 GeV on both sides), but the scans in
//     this repo book 100 bins from 0 GeV: pairing a new file with an old one
//     would shift the MC window by 20 GeV with no error.
//  2. A fit is refused, not reported, when (a) the data has entries in a ptRel
//     bin inside the fit range where both templates are empty -- the likelihood
//     is then -inf for every fb -- or (b) the minimizer does not converge.
//     templateFitter() takes whatever the minimizer leaves, which after a failed
//     fit is fb near its 0.5 start with an error of ~0.5: a number that means
//     nothing but reads like one. (a) happens below 80 GeV, where the MC thins
//     out at high ptRel; templateFitter() never fits that low.
//  3. The closure's true b-fraction error adds the b and bGS counts' errors in
//     quadrature as absolute errors; templateFitter() added their relative
//     errors, which overstates it.
//  4. The hand-computed chi2/ndof skips bins with no error instead of dividing
//     by zero, and normalizes by the bins it used.
//
// VARIANTS are (pp data, PYTHIA) file pairs. The muon-added scans can write the
// usual histogram names into separate files, or new names into the same files;
// either works -- set the paths and, if needed, the names in the table.
//
// OUTPUT: per-window fit figures and the b fraction vs jet pT in
// figures/templateFits/pp_caloJets/, the fractions as TH1Ds in
// rootFiles/bPurityResults/bFraction_pp_caloJets_<variant>.root. Failed fits
// are stored as -999 (kAbsent), so they are left out of plots rather than
// drawn as a fraction of zero.
//
// Usage, from src/calculateBJetsPerZ/:
//   root -l -b -q 'fitBFraction_pp_caloJets.C'
//   root -l -b -q 'fitBFraction_pp_caloJets.C("dataNoMu.root","mcNoMu.root","dataMu.root","mcMu.root")'
// Non-empty arguments override the paths in the table. A variant with an empty
// path is skipped, so either can be run on its own. A fifth argument is appended
// to every output name, to keep a test run from overwriting real results.

#include <cmath>
#include <memory>
#include "TH1D.h"
#include "../../headers/functions/divideByBinwidth.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/plotting/ratioPanel.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TNamed.h"
#include "TSystem.h"

// ---------------------------------------------------------------- inputs ----

struct Variant {
  const char *tag;        // output file / directory tag
  const char *label;      // legend text
  TString     dataFile;   // pp SingleMuon, calo jets
  TString     mcFile;     // PYTHIA DiJet, calo jets
  const char *dataHist;   // ptRel (x) vs jet pT (y), tagged jets, trigger on
  const char *mcHistFmt;  // same, per flavour; %s = allJets, bJets, cJets, ...
};

Variant variants[] = {
  {"noMuAdded", "calo jet",
   "",   // TODO pp_SingleMuon_caloJets_...root
   "",   // TODO PYTHIA_DiJet_caloJets_...root
   "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn",
   "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_%s_T0"},
  {"muAdded", "calo jet + #mu p^{4}",
   "",   // TODO
   "",   // TODO
   "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn",
   "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn_%s_T0"},
};
const int NVariants = sizeof(variants)/sizeof(Variant);

// PLACEHOLDER until the calo binning is settled: calculateBPurity()'s pp edges
// with the calo reach below 80 GeV prepended. Every edge must be a bin boundary
// of both files' jet pT axes (5 GeV in current scans); this is checked per fit.
const double jetPtEdges[] = {50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500};
const int    NPt = sizeof(jetPtEdges)/sizeof(double) - 1;

// ------------------------------------------------------------ fit settings ---

const int M = 51;
double muRelPtAxis[M] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0};

const double low_x = 0.0, high_x = 5.0;    // RooFit "fit_region"
const double c_multiplier   = 1.0;
const double bGS_multiplier = 1.175;
// Also fit PYTHIA's inclusive sample and compare with its true b fraction. With
// bGS_multiplier != 1 the b template no longer matches the sample it is fitted
// to, so the fit comes out LOW by construction: -0.03 to -0.07 at 1.175 on the
// 2026-02 PYTHIA file, against agreement to 0.002 with the multiplier at 1.0.
const bool   doMCClosure    = true;

const char *hexB = "#D55E00";   // b-jet template
const char *hexL = "#0072B2";   // light+c template

// ----------------------------------------------------------------- helpers ---

struct FitResult {
  bool    ok = false;
  double  fb = -1., efb = 0.;          // RooFit b fraction
  double  fbTrue = -1., efbTrue = 0.;  // MC closure only
  double  chi2ndf = -1.;
  double  nJets = 0.;                  // weighted count in the window
  int     status = -1;                 // RooFit minimizer status, 0 = converged
  bool    gsApplied = false;           // bGS enhancement actually applied (bGS template non-empty)
  double  xFrac = -1.;                 // MC closure only: unmatched-flavour (xJets) share of the fitted sample
  TString why;                         // reason, when !ok
};

// Jet pT window on h's OWN y axis. Fails if an edge is not a bin boundary,
// since the projection would then silently widen or narrow.
static TH1D* projectWindow(TH2D *H, double lo, double hi, const char *name)
{
  TAxis *a = H->GetYaxis();
  int b1 = a->FindBin(lo + 1e-6), b2 = a->FindBin(hi - 1e-6);
  if(fabs(a->GetBinLowEdge(b1) - lo) > 1e-6 || fabs(a->GetBinUpEdge(b2) - hi) > 1e-6){
    printf("ERROR: %s: %.0f-%.0f GeV is not on its jet pT bin boundaries\n", H->GetName(), lo, hi);
    return nullptr;
  }
  TH1D *p = H->ProjectionX(name, b1, b2);
  p->SetDirectory(nullptr);
  return p;
}

static TH1D* cloneH(TH1D *h, const char *name)
{
  TH1D *c = (TH1D*) h->Clone(name);
  c->SetDirectory(nullptr);
  return c;
}

// One jet pT window. isData = false fits PYTHIA's own allJets histogram against
// its own templates (templateFitter()'s !isData closure).
static FitResult fitOneBin(TFile *fD, TFile *fM, const Variant &v, bool isData,
                           double lo, double hi, const TString &figDir)
{
  FitResult r;
  TString sfx = Form("%s_%s_%.0f", v.tag, isData ? "data" : "mc", lo);

  TH2D *H_data = nullptr;
  if(isData) fD->GetObject(v.dataHist, H_data);
  else       fM->GetObject(Form(v.mcHistFmt, "allJets"), H_data);
  if(!H_data){ r.why = isData ? "data histogram missing" : "MC allJets histogram missing"; return r; }

  const int NFl = 7;
  const char *flavour[NFl] = {"bJets","bGSJets","cJets","uJets","dJets","sJets","gJets"};
  TH1D *h[NFl] = {nullptr};
  for(int i = 0; i < NFl; i++){
    TH2D *H = nullptr;
    fM->GetObject(Form(v.mcHistFmt, flavour[i]), H);
    if(!H){ r.why = Form("%s missing", Form(v.mcHistFmt, flavour[i])); return r; }
    h[i] = projectWindow(H, lo, hi, Form("h_%s_%s", flavour[i], sfx.Data()));
    if(!h[i]){ r.why = "jet pT edges off the MC axis"; return r; }
  }
  TH1D *h_data = projectWindow(H_data, lo, hi, Form("h_data_%s", sfx.Data()));
  if(!h_data){ r.why = "jet pT edges off the data axis"; return r; }

  TH1D *h_b = h[0], *h_bGS = h[1], *h_c = h[2];
  TH1D *h_l = cloneH(h[3], Form("h_l_%s", sfx.Data()));
  h_l->Add(h[4]); h_l->Add(h[5]); h_l->Add(h[6]);

  // MC-truth flavour fractions over the whole ptRel axis, before any rescaling
  double nAll = h_l->Integral() + h_b->Integral() + h_c->Integral() + h_bGS->Integral();
  r.nJets = h_data->Integral();
  if(nAll <= 0. || r.nJets <= 0. || h_b->Integral() <= 0.){ r.why = "empty window"; return r; }
  double c_truth = h_c->Integral() / nAll;
  double l_truth = h_l->Integral() / nAll;

  // closure truth: b + bGS over the fit range, as templateFitter() counts it
  if(!isData){
    const double eps = 0.001;
    int x1 = h_data->FindBin(low_x + eps), x2 = h_data->FindBin(high_x - eps);
    double en_tot = 0., en_b = 0., en_bGS = 0.;
    double n_tot = h_data->IntegralAndError(x1, x2, en_tot);
    double n_b   = h_b  ->IntegralAndError(x1, x2, en_b);
    double n_bGS = h_bGS->IntegralAndError(x1, x2, en_bGS);
    double n_bt = n_b + n_bGS, en_bt = sqrt(en_b*en_b + en_bGS*en_bGS);
    if(n_tot > 0. && n_bt > 0.){
      r.fbTrue  = n_bt / n_tot;
      // b is a subset of the total, so treating them as independent makes this
      // an upper bound
      r.efbTrue = r.fbTrue * sqrt(pow(en_bt/n_bt, 2) + pow(en_tot/n_tot, 2));
    }
    // Jets with no flavour match are in the fitted sample but in neither
    // template, so the fit has to put them somewhere -- and their ptRel is
    // b-like (<ptRel> 1.35 against 1.52 for b, ~0.9 for c/u/d/s). Their share
    // is what limits the closure: removed from the sample, the fit reproduces
    // the truth to 0.002 (checked 2026-09-15 on the refparton_flavor calo
    // PYTHIA, where they are 18% of tagged jets above 50 GeV).
    TH2D *Hx = nullptr; fM->GetObject(Form(v.mcHistFmt, "xJets"), Hx);
    TH1D *hx = Hx ? projectWindow(Hx, lo, hi, Form("h_xJets_%s", sfx.Data())) : nullptr;
    if(hx && n_tot > 0.){ r.xFrac = hx->Integral(x1, x2) / n_tot; delete hx; }
  }

  // Gluon-splitting enhancement. As in templateFitter(), bGS_multiplier - 1 is
  // ADDED to the gluon-splitting fraction of b jets, not multiplied into it.
  double N_b = h_b->Integral(), N_bGS = h_bGS->Integral();
  if(N_bGS > 0.){
    double f_bGS = N_bGS / (N_b + N_bGS);
    double target = f_bGS + (bGS_multiplier - 1.);
    if(target >= 1.){ r.why = Form("bGS fraction %.3f + shift reaches 1", f_bGS); return r; }
    h_bGS->Scale(N_b * target / (N_bGS * (1. - target)));
    r.gsApplied = (bGS_multiplier != 1.);
  }
  h_b->Add(h_bGS);

  // 2-template fit: c folded into light at its MC-truth fraction
  if(h_c->Integral() > 0.) h_c->Scale(c_multiplier * c_truth / h_c->Integral());
  if(h_l->Integral() > 0.) h_l->Scale(l_truth / h_l->Integral());
  h_l->Add(h_c);
  if(h_l->Integral() <= 0.){ r.why = "no light+c template"; return r; }
  h_l->Scale(1./h_l->Integral());
  h_b->Scale(1./h_b->Integral());

  // h_roo keeps the data's raw counts for the likelihood; h_data is unit-normalized for drawing
  TH1D *h_roo = cloneH(h_data, Form("h_roo_%s", sfx.Data()));
  h_data->Scale(1./h_data->Integral());

  TH1D *R_roo  = rebinTo(h_roo,  M-1, muRelPtAxis, Form("R_roo_%s",  sfx.Data()));
  TH1D *R_data = rebinTo(h_data, M-1, muRelPtAxis, Form("R_data_%s", sfx.Data()));
  TH1D *R_b    = rebinTo(h_b,    M-1, muRelPtAxis, Form("R_b_%s",    sfx.Data()));
  TH1D *R_l    = rebinTo(h_l,    M-1, muRelPtAxis, Form("R_l_%s",    sfx.Data()));
  divideByBinwidth(R_roo); divideByBinwidth(R_data);
  divideByBinwidth(R_b);   divideByBinwidth(R_l);

  // Difference 2a: data where both templates are empty makes the likelihood
  // -inf at every fb, so no fit result there means anything.
  { TString where; int nEmpty = 0;
    for(int b = R_roo->FindBin(low_x + 1e-6); b <= R_roo->FindBin(high_x - 1e-6); b++){
      if(R_roo->GetBinContent(b) > 0. && R_b->GetBinContent(b) <= 0. && R_l->GetBinContent(b) <= 0.){
        if(nEmpty < 4) where += Form(" %.1f", R_roo->GetBinLowEdge(b));
        nEmpty++;
      }
    }
    if(nEmpty > 0){
      r.why = Form("%d ptRel bin(s) in the fit range with data but no template, at%s%s GeV",
                   nEmpty, where.Data(), nEmpty > 4 ? " ..." : "");
      return r;
    }
  }

  RooRealVar x(Form("x_%s", sfx.Data()), "muon ptRel", muRelPtAxis[0], muRelPtAxis[M-1]);
  x.setRange("fit_region", low_x, high_x);
  RooDataHist dataHist(Form("dataHist_%s", sfx.Data()), "data", x, RooFit::Import(*R_roo));
  RooDataHist bHist(Form("bHist_%s", sfx.Data()), "b template", x, RooFit::Import(*R_b));
  RooDataHist lHist(Form("lHist_%s", sfx.Data()), "light+c template", x, RooFit::Import(*R_l));
  RooHistPdf  bPdf(Form("bPdf_%s", sfx.Data()), "b PDF", x, bHist, 0);
  RooHistPdf  lPdf(Form("lPdf_%s", sfx.Data()), "light+c PDF", x, lHist, 0);
  RooRealVar  fb(Form("fb_%s", sfx.Data()), "b fraction", 0.5, 0.0, 1.0);
  RooAddPdf   model(Form("model_%s", sfx.Data()), "b + light+c", RooArgList(bPdf, lPdf), RooArgList(fb));
  std::unique_ptr<RooFitResult> fitRes{model.fitTo(dataHist,
                                                   RooFit::Minimizer("Minuit","migrad"),
                                                   RooFit::SumW2Error(true),
                                                   RooFit::Range("fit_region"),
                                                   RooFit::PrintLevel(-1),
                                                   RooFit::Save())};
  r.fb     = fb.getVal();
  r.efb    = fb.getError();
  r.status = fitRes ? fitRes->status() : -1;
  // difference 2b
  if(r.status != 0 || !std::isfinite(r.fb) || !std::isfinite(r.efb)){
    r.why = Form("fit did not converge (status %d)", r.status);
    return r;
  }

  // fitted composition, and the hand chi2 templateFitter() prints
  TH1D *b_sc = cloneH(R_b, Form("b_sc_%s", sfx.Data())); b_sc->Scale(r.fb);
  TH1D *l_sc = cloneH(R_l, Form("l_sc_%s", sfx.Data())); l_sc->Scale(1. - r.fb);
  TH1D *fitSum = cloneH(b_sc, Form("fitSum_%s", sfx.Data())); fitSum->Add(l_sc);
  // both carry statistical errors (data, and finite MC templates)
  TH1D *ratio = makeRatio(R_data, fitSum, Form("ratio_%s", sfx.Data()), RatioErr::kBoth);
  double chi2 = 0.; int nUsed = 0;
  for(int b = 1; b <= ratio->GetNbinsX(); b++){
    double val = ratio->GetBinContent(b), err = ratio->GetBinError(b);
    if(val == kAbsent || err <= 0.) continue;
    chi2 += (val - 1.)*(val - 1.)/(err*err); nUsed++;
  }
  r.chi2ndf = nUsed > 0 ? chi2/nUsed : -1.;
  r.ok = true;

  if(figDir != ""){
    TCanvas *c = new TCanvas(Form("c_%s", sfx.Data()), "", 600, 700);
    TPad *pT = nullptr, *pB = nullptr;
    splitPads(pT, pB);

    pT->cd();
    // Stacked by hand: the total in the light+c color, b drawn over its lower
    // part. A THStack here segfaulted inside SaveAs under ROOT 6.40 -- removing
    // it, and only it, stopped the crash.
    int cb = TColor::GetColor(hexB), cl = TColor::GetColor(hexL);
    TH1D *stackTop = cloneH(fitSum, Form("stackTop_%s", sfx.Data()));
    stackTop->SetFillColorAlpha(cl, 0.6); stackTop->SetLineColor(cl); stackTop->SetMarkerSize(0);
    b_sc->SetFillColorAlpha(cb, 0.6);     b_sc->SetLineColor(cb);     b_sc->SetMarkerSize(0);
    styleH(R_data, hexData, markFilledCircle);

    TH1D *frame = cloneH(R_data, Form("frame_%s", sfx.Data()));
    frame->Reset();
    frame->SetMinimum(0.);
    frame->SetMaximum(1.55 * TMath::Max(R_data->GetMaximum(), fitSum->GetMaximum()));
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("normalized entries / GeV");
    frame->GetYaxis()->SetTitleSize(0.050); frame->GetYaxis()->SetTitleOffset(1.45);
    frame->GetYaxis()->SetLabelSize(0.042);
    frame->Draw("AXIS");
    stackTop->Draw("hist same");
    b_sc->Draw("hist same");
    R_data->Draw("E same");
    frame->Draw("AXIS same");

    TLegend *leg = makeLegend(0.62, 0.62, 0.93, 0.83, 0.040);
    leg->AddEntry(R_data, isData ? "pp data" : "PYTHIA", "lp");
    leg->AddEntry(b_sc, "#it{b} jets", "f");
    leg->AddEntry(stackTop, "light + #it{c} jets", "f");
    leg->Draw();

    TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.040);
    la.DrawLatex(0.21, 0.885, Form("%s, anti-#it{k}_{T} #it{R} = 0.4 %s",
                                   isData ? "pp 5.02 TeV" : "PYTHIA 5.02 TeV", v.label));
    la.DrawLatex(0.21, 0.835, Form("%.0f < #it{p}_{T}^{jet} < %.0f GeV", lo, hi));
    la.DrawLatex(0.21, 0.780, Form("#it{f}_{#it{b}} = %.3f #pm %.3f", r.fb, r.efb));
    if(!isData) la.DrawLatex(0.21, 0.720, Form("#it{f}_{#it{b}}^{true} = %.3f #pm %.3f", r.fbTrue, r.efbTrue));
    la.DrawLatex(0.21, isData ? 0.725 : 0.660, Form("#chi^{2}/ndof = %.2f", r.chi2ndf));

    pB->cd();
    styleH(ratio, hexData, markFilledCircle);
    styleRatioAxes(ratio, "muon #it{p}_{T}^{rel} [GeV]", isData ? "data / fit" : "MC / fit");
    ratio->SetMinimum(0.5);
    ratio->SetMaximum(TMath::Min(2.0, ratioMax(ratio, fitSum, low_x, high_x, 1.5)));
    ratio->Draw("E");
    unityLine(muRelPtAxis[0], muRelPtAxis[M-1])->Draw();
    ratio->Draw("E same");

    c->SaveAs(Form("%s/templateFit_%s_%.0f-%.0f.pdf", figDir.Data(), isData ? "data" : "mcClosure", lo, hi));
    delete c; delete stackTop; delete frame;
  }

  for(int i = 0; i < NFl; i++) delete h[i];
  delete h_data; delete h_l; delete h_roo;
  delete R_roo; delete R_data; delete R_b; delete R_l;
  delete b_sc; delete l_sc; delete fitSum; delete ratio;
  return r;
}

// Bins start at kAbsent, so a failed fit is left out of every plot rather than
// drawn as a measured fraction of zero. The written histograms carry it too.
static TH1D* bookResult(const char *name)
{
  TH1D *h = new TH1D(name, "", NPt, jetPtEdges);
  h->SetDirectory(nullptr);
  for(int b = 1; b <= NPt; b++) h->SetBinContent(b, kAbsent);
  return h;
}

// makeRatio() only knows about empty denominators; carry the sentinel through
static void maskAbsent(TH1D *r, TH1D *a, TH1D *b)
{
  for(int i = 1; i <= r->GetNbinsX(); i++)
    if(a->GetBinContent(i) == kAbsent || b->GetBinContent(i) == kAbsent){
      r->SetBinContent(i, kAbsent); r->SetBinError(i, 0.);
    }
}

// ------------------------------------------------------------------ driver ---

void fitBFraction_pp_caloJets(const char *dataNoMu = "", const char *mcNoMu = "",
                              const char *dataMu = "",   const char *mcMu = "",
                              const char *outTag = "")
{
  initPlotStyle();
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  if(dataNoMu[0]) variants[0].dataFile = dataNoMu;
  if(mcNoMu[0])   variants[0].mcFile   = mcNoMu;
  if(dataMu[0])   variants[1].dataFile = dataMu;
  if(mcMu[0])     variants[1].mcFile   = mcMu;

  TString figBase = "../../figures/templateFits/pp_caloJets";
  TString resDir  = "../../rootFiles/bPurityResults";
  gSystem->mkdir(resDir, kTRUE);

  TH1D *hData[NVariants], *hMCFit[NVariants], *hMCTrue[NVariants];
  bool   gsAny[NVariants];      // any fit in the variant applied the bGS enhancement
  double xFracMax[NVariants];   // largest unmatched-flavour share seen in the closure
  for(int iv = 0; iv < NVariants; iv++){
    hData[iv] = hMCFit[iv] = hMCTrue[iv] = nullptr;
    gsAny[iv] = false; xFracMax[iv] = -1.;
  }

  for(int iv = 0; iv < NVariants; iv++){
    const Variant &v = variants[iv];
    printf("\n=== %s (%s) ===\n", v.tag, v.label);
    if(v.dataFile == "" || v.mcFile == ""){ printf("  no input paths set, skipped\n"); continue; }

    TFile *fD = TFile::Open(v.dataFile), *fM = TFile::Open(v.mcFile);
    if(!fD || fD->IsZombie()){ printf("ERROR: cannot open %s\n", v.dataFile.Data()); continue; }
    if(!fM || fM->IsZombie()){ printf("ERROR: cannot open %s\n", v.mcFile.Data()); continue; }
    printf("  data: %s\n  MC:   %s\n", v.dataFile.Data(), v.mcFile.Data());
    for(TFile *f : {fD, fM}){
      TNamed *p = nullptr; f->GetObject("provenance", p);
      if(p) printf("  provenance: %s\n", p->GetTitle());
    }

    // report the two jet pT axes: projectWindow() handles a mismatch, but it is
    // worth knowing the files come from different booking conventions
    { TH2D *a = nullptr, *b = nullptr;
      fD->GetObject(v.dataHist, a); fM->GetObject(Form(v.mcHistFmt, "bJets"), b);
      if(a && b){
        printf("  jet pT axis  data %d bins [%g,%g]   MC %d bins [%g,%g]\n",
               a->GetNbinsY(), a->GetYaxis()->GetXmin(), a->GetYaxis()->GetXmax(),
               b->GetNbinsY(), b->GetYaxis()->GetXmin(), b->GetYaxis()->GetXmax());
      } }

    TString figDir = Form("%s/%s%s", figBase.Data(), v.tag, outTag);
    gSystem->mkdir(figDir, kTRUE);

    hData[iv]   = bookResult(Form("bFrac_data_%s", v.tag));
    hMCFit[iv]  = bookResult(Form("bFrac_mcFit_%s", v.tag));
    hMCTrue[iv] = bookResult(Form("bFrac_mcTrue_%s", v.tag));

    printf("  %-9s %9s %16s %8s", "jet pT", "N data", "f_b data", "chi2/n");
    if(doMCClosure) printf(" %16s %16s %8s %7s", "f_b MC fit", "f_b MC true", "fit-true", "x frac");
    printf("\n");

    for(int i = 0; i < NPt; i++){
      double lo = jetPtEdges[i], hi = jetPtEdges[i+1];
      FitResult rd = fitOneBin(fD, fM, v, true, lo, hi, figDir);
      FitResult rm;
      if(doMCClosure) rm = fitOneBin(fD, fM, v, false, lo, hi, figDir);

      if((rd.ok && rd.gsApplied) || (rm.ok && rm.gsApplied)) gsAny[iv] = true;
      if(rm.ok && rm.xFrac > xFracMax[iv]) xFracMax[iv] = rm.xFrac;

      if(rd.ok){ hData[iv]->SetBinContent(i+1, rd.fb); hData[iv]->SetBinError(i+1, rd.efb); }
      if(rm.ok){
        hMCFit[iv] ->SetBinContent(i+1, rm.fb);     hMCFit[iv] ->SetBinError(i+1, rm.efb);
        hMCTrue[iv]->SetBinContent(i+1, rm.fbTrue); hMCTrue[iv]->SetBinError(i+1, rm.efbTrue);
      }

      printf("  %3.0f-%-5.0f %9.0f", lo, hi, rd.nJets);
      if(rd.ok) printf(" %8.4f+-%.4f %8.2f", rd.fb, rd.efb, rd.chi2ndf);
      else      printf(" %25s", "FAILED");
      if(doMCClosure){
        if(rm.ok) printf(" %8.4f+-%.4f %8.4f+-%.4f %+8.4f %7.3f", rm.fb, rm.efb, rm.fbTrue, rm.efbTrue, rm.fb - rm.fbTrue, rm.xFrac);
        else      printf(" %42s", "FAILED");
      }
      printf("\n");
      if(!rd.ok) printf("              data:       %s\n", rd.why.Data());
      if(doMCClosure && !rm.ok) printf("              MC closure: %s\n", rm.why.Data());
    }

    TString resName = Form("%s/bFraction_pp_caloJets_%s%s.root", resDir.Data(), v.tag, outTag);
    TFile *wf = TFile::Open(resName, "recreate");
    hData[iv]->Write();
    if(doMCClosure){ hMCFit[iv]->Write(); hMCTrue[iv]->Write(); }
    TNamed("inputs", Form("data=%s mc=%s dataHist=%s mcHistFmt=%s",
                          v.dataFile.Data(), v.mcFile.Data(), v.dataHist, v.mcHistFmt)).Write();
    TNamed("settings", Form("fit %.1f-%.1f GeV, 2-template, c_multiplier=%.3f, bGS_multiplier=%.3f, failed fits = %g",
                            low_x, high_x, c_multiplier, bGS_multiplier, kAbsent)).Write();
    wf->Close();
    printf("  results: %s\n", resName.Data());

    fD->Close(); fM->Close();
  }

  // ---- summary: data b fraction vs jet pT, variants overlaid ----------------
  const char *hexV[2] = {hexMC, hexCorrected};
  const int   markV[2] = {markFilledCircle, markFilledSquare};
  int nDone = 0;
  for(int iv = 0; iv < NVariants; iv++) if(hData[iv]) nDone++;
  if(nDone == 0){ printf("\nno variant was fitted\n"); return; }

  bool both = (NVariants >= 2 && hData[0] && hData[1]);
  // Same 700x800 either way: a 700x600 canvas came out shifted and squashed
  // in the PDF, where every 700x800 figure renders cleanly.
  TCanvas *cs = new TCanvas("c_summary", "", 700, 800);
  TPad *pT = nullptr, *pB = nullptr;
  if(both) splitPads(pT, pB);
  else{
    pT = new TPad("padTop", "", 0, 0, 1, 1);
    // room above the frame for the header text
    pT->SetLeftMargin(0.17); pT->SetBottomMargin(0.14); pT->SetTopMargin(0.20);
    pT->Draw();
  }

  pT->cd();
  // PYTHIA truth (b + bGS share of the tagged jets, as in the closure) sits with
  // each variant's data in the same color, open marker. Truth is per variant
  // because adding the muon moves jets between pT bins.
  const int markTrueV[2] = {markOpenCircle, markOpenSquare};
  int nEntries = 0;
  for(int iv = 0; iv < NVariants; iv++)
    if(hData[iv]) nEntries += (doMCClosure && hMCTrue[iv]) ? 2 : 1;
  TLegend *leg = makeLegend(0.55, 0.76 - 0.065*nEntries, 0.93, 0.76, 0.040);
  bool first = true;
  for(int iv = 0; iv < NVariants; iv++){
    if(!hData[iv]) continue;
    styleH(hData[iv], hexV[iv % 2], markV[iv % 2]);
    if(first){
      hData[iv]->SetMinimum(0.); hData[iv]->SetMaximum(1.);
      hData[iv]->GetYaxis()->SetTitle("#it{b}-jet fraction");
      hData[iv]->GetYaxis()->SetTitleSize(0.050); hData[iv]->GetYaxis()->SetTitleOffset(1.35);
      hData[iv]->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
      if(both) hData[iv]->GetXaxis()->SetLabelSize(0);
      else{
        hData[iv]->GetXaxis()->SetTitleSize(0.045); hData[iv]->GetXaxis()->SetTitleOffset(1.2);
        hData[iv]->GetXaxis()->SetLabelSize(0.040);
      }
      hData[iv]->Draw("E1");
      first = false;
    }
    else hData[iv]->Draw("E1 same");
    // the variant label only distinguishes anything when both are drawn
    TString suffix = both ? Form(", %s", variants[iv].label) : "";
    leg->AddEntry(hData[iv], "pp data" + suffix, "lp");
    if(doMCClosure && hMCTrue[iv]){
      // a clone: the closure canvas below restyles hMCTrue
      TH1D *t = cloneH(hMCTrue[iv], Form("bFrac_mcTrue_summary_%s", variants[iv].tag));
      styleH(t, hexV[iv % 2], markTrueV[iv % 2]);
      t->Draw("E1 same");
      leg->AddEntry(t, "PYTHIA truth" + suffix, "lp");
    }
  }
  leg->Draw();
  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.040);
  la.DrawLatex(0.21, 0.885, "pp 5.02 TeV, anti-#it{k}_{T} #it{R} = 0.4 calo jets");
  la.DrawLatex(0.21, 0.835, Form("2-template fit, %.0f < #it{p}_{T}^{rel} < %.0f GeV", low_x, high_x));

  if(both){
    pB->cd();
    // The two variants fit the same events, binned by a different jet pT, so
    // they are neither independent nor the same jets per bin. The propagated
    // error is an upper bound.
    TH1D *r = makeRatio(hData[1], hData[0], "r_variants", RatioErr::kBoth);
    maskAbsent(r, hData[1], hData[0]);
    styleH(r, hexCorrected, markFilledSquare);
    styleRatioAxes(r, "#it{p}_{T}^{jet} [GeV]", "#mu added / not");
    r->SetMinimum(0.5); r->SetMaximum(1.5);
    r->Draw("E1");
    unityLine(jetPtEdges[0], jetPtEdges[NPt])->Draw();
    r->Draw("E1 same");
  }
  gSystem->mkdir(figBase, kTRUE);
  cs->SaveAs(Form("%s/bFraction_vs_jetPt%s.pdf", figBase.Data(), outTag));
  delete cs;

  // ---- per-variant MC closure: fit vs truth ---------------------------------
  if(doMCClosure){
    for(int iv = 0; iv < NVariants; iv++){
      if(!hMCFit[iv]) continue;
      TCanvas *cc = new TCanvas(Form("c_closure_%d", iv), "", 700, 800);
      TPad *qT = nullptr, *qB = nullptr;
      splitPads(qT, qB);
      qT->cd();
      styleH(hMCTrue[iv], hexData, markOpenCircle);
      styleH(hMCFit[iv],  hexMC,   markFilledSquare);
      hMCTrue[iv]->SetMinimum(0.); hMCTrue[iv]->SetMaximum(1.);
      hMCTrue[iv]->GetYaxis()->SetTitle("#it{b}-jet fraction");
      hMCTrue[iv]->GetYaxis()->SetTitleSize(0.050); hMCTrue[iv]->GetYaxis()->SetTitleOffset(1.35);
      hMCTrue[iv]->GetXaxis()->SetLabelSize(0);
      hMCTrue[iv]->Draw("E1");
      hMCFit[iv]->Draw("E1 same");
      TLegend *lc = makeLegend(0.55, 0.62, 0.93, 0.76, 0.040);
      lc->AddEntry(hMCTrue[iv], "PYTHIA truth", "lp");
      lc->AddEntry(hMCFit[iv],  "template fit", "lp");
      lc->Draw();
      la.DrawLatex(0.21, 0.885, Form("PYTHIA closure, %s", variants[iv].label));
      // Name what keeps the fit from the truth, but only what actually applies:
      // this note used to print whenever bGS_multiplier != 1, including on a
      // sample with an empty bGS template where no enhancement was made.
      double yNote = 0.835;
      if(gsAny[iv]){
        la.DrawLatex(0.21, yNote, Form("#it{b} template #it{g}#rightarrow#it{b}#bar{#it{b}} fraction %+.3f: fit #neq truth by construction",
                                       bGS_multiplier - 1.));
        yNote -= 0.045;
      }
      if(xFracMax[iv] > 0.02)
        la.DrawLatex(0.21, yNote, Form("unmatched-flavor jets (#leq %.0f%%) in no template", 100.*xFracMax[iv]));

      qB->cd();
      TH1D *rc = makeRatio(hMCFit[iv], hMCTrue[iv], Form("r_closure_%d", iv), RatioErr::kNumerator);
      maskAbsent(rc, hMCFit[iv], hMCTrue[iv]);
      styleH(rc, hexMC, markFilledSquare);
      styleRatioAxes(rc, "#it{p}_{T}^{jet} [GeV]", "fit / truth");
      rc->SetMinimum(0.5); rc->SetMaximum(1.5);
      rc->Draw("E1");
      unityLine(jetPtEdges[0], jetPtEdges[NPt])->Draw();
      rc->Draw("E1 same");
      cc->SaveAs(Form("%s/%s%s/mcClosure_vs_jetPt.pdf", figBase.Data(), variants[iv].tag, outTag));
      delete cc;
    }
  }

  printf("\nfigures in %s\n", figBase.Data());
}
