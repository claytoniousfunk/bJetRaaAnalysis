// Shared machinery for the pp CALO-jet b-purity template fits: inputs, template
// construction and the RooFit fit. Used by
//   src/plots/bPurity/plotBPurity_caloJets_pp.C        (fit figures)
//   src/systematics/systematics_bGS_caloJets_pp.C      (bGS-share systematic)
// so both see exactly the same templates and fit.
//
// Calo data with calo templates only -- never mix with PF.
#pragma once
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TBox.h"
#include "TLatex.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"

const char *dataPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/pp/"
  "pp_SingleMuon_caloJets_manualJEC_mu12_pTmu-15to999_tight_deltaR-40_"
  "mu12TriggerEfficiencyCorrection_jetTrkMaxFilter_WDecayFilter_2026-9-29.root";
const char *mcPath =
  "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIA/"
  "PYTHIA_DiJet_caloJets_PFflavor_PFbHadNum_manualJEC_pThat-15_mu12_pTmu-15to999_"
  "tight_vzReweight_jetTrkMaxFilter_removeHYDJETjet0p35CutOnGen_WDecayFilter_2026-10-4.root";
const char *hName = "h_muptrel_recoJetPt_inclRecoMuonTag_triggerOn";

// Jet-pT windows, matched to the unfolding binning. 50-60, 60-70 and 70-80 are
// UNDERFLOW bins kept for the unfolding and not reported; the first reported
// bin starts at ptReportMin. Edges sit on the 5 GeV bins of the input
// histograms.
//
// The pp SingleMuon calo data has NO jets below 70 GeV: that forest stores
// calo jets only above ~35 GeV raw (MinBias: ~20), which the calo JEC maps to
// ~65-70 GeV at |eta| < 1.6. 50-60 and 60-70 therefore cannot be fitted from
// this data (see fittable()), and 70-80 sits on the threshold turn-on.
const int NPt = 9;
const double ptEdges[NPt+1] = {50, 60, 70, 80, 100, 120, 150, 200, 300, 500};
const double ptReportMin = 80.;

// Minimum data entries in 0-5 GeV ptRel for a window to be fitted at all.
// Below it the window is skipped and left empty in every output, rather than
// fitted on nothing.
const double minFitEntries = 50.;

const double fitLo = 0.0, fitHi = 5.0;

int uidP = 0;

// ptRel projection in a jet-pT window. Window edges are on 5 GeV bin
// boundaries; FindBin(hi - eps) keeps the window from picking up the next bin.
TH1D* projWin(TFile *f, const char *name, double lo, double hi)
{
  TH2D *H = nullptr; f->GetObject(name, H);
  if(!H){ printf("ERROR: missing %s\n", name); return nullptr; }
  TH1D *p = H->ProjectionX(Form("pw%d", uidP++), H->GetYaxis()->FindBin(lo + 1e-6),
                           H->GetYaxis()->FindBin(hi - 1e-6));
  p->SetDirectory(nullptr);
  // weighted-MC bins can come out slightly negative; a template cannot
  for(int i = 0; i <= p->GetNbinsX()+1; i++) if(p->GetBinContent(i) < 0.) p->SetBinContent(i, 0.);
  return p;
}

// copy of the 0-5 GeV part, 0.1 GeV bins, as the fit sees it
TH1D* range05(TH1D *h, const char *name)
{
  TH1D *o = new TH1D(name, "", 50, 0., 5.); o->SetDirectory(nullptr); o->Sumw2();
  for(int i = 1; i <= 50; i++){ o->SetBinContent(i, h->GetBinContent(i)); o->SetBinError(i, h->GetBinError(i)); }
  return o;
}

inline bool fittable(TH1D *data05){ return data05 && data05->Integral() >= minFitEntries; }

// Grey band over the underflow bins (ptEdges[0] .. ptReportMin), drawn in the
// current pad's user coordinates between ymin and ymax, with a small label
// starting at labelY (pass labelY below ymin for no label).
inline void drawUnderflowBand(double ymin, double ymax, double labelY, double textSize = 0.030)
{
  TBox *b = new TBox(ptEdges[0], ymin, ptReportMin, ymax);
  b->SetFillColorAlpha(kGray, 0.45); b->SetLineWidth(0); b->Draw();
  if(labelY > ymin){
    // vertical: the band is too narrow for the word written horizontally
    TLatex *t = new TLatex(0.5*(ptEdges[0] + ptReportMin), labelY, "underflow");
    t->SetTextAlign(12); t->SetTextAngle(90); t->SetTextFont(42); t->SetTextSize(textSize);
    t->SetTextColor(kGray+2); t->Draw();
  }
}

// For drawing only: bins whose nominal purity is empty (window not fitted) are
// pushed far below any frame, so they are absent rather than drawn at zero.
inline void blankUnfitted(TH1D *h, TH1D *nominal)
{
  for(int i = 1; i <= h->GetNbinsX(); i++)
    if(nominal->GetBinContent(i) <= 0.){ h->SetBinContent(i, -999.); h->SetBinError(i, 0.); }
}

struct Tpl {
  TH1D *bFC, *bGS, *lc;     // unit-normalised over 0-5 GeV
  double fGSnat;            // PYTHIA's bGS share of all b in this window
  double truthB;            // MC-truth b fraction (b + bGS) over 0-5 GeV
};

// cMult scales the c weight inside the light+c template relative to its MC
// truth, as c_multiplier does in templateFitter() (1.0 = MC truth).
// tplIdx picks the MC jet-pT variant the templates are binned in
// (templateIndexNames: 0 nominal, 1 JER smear, 2 JEU up, 3 JEU down).
Tpl buildTemplates(TFile *fM, double lo, double hi, double cMult = 1.0, int tplIdx = 0)
{
  auto T = [&](const char *fl){ return projWin(fM, Form("%s_%sJets_T%d", hName, fl, tplIdx), lo, hi); };
  TH1D *b = T("b"), *g = T("bGS"), *c = T("c"), *l = T("u");
  l->Add(T("d")); l->Add(T("s")); l->Add(T("g"));

  Tpl t;
  double nb = b->Integral(1,50), ng = g->Integral(1,50), nc = c->Integral(1,50), nl = l->Integral(1,50);
  t.truthB = (nb + ng)/(nb + ng + nc + nl);
  t.fGSnat = ng/(nb + ng);

  // c folded into light at cMult x its MC-truth weight relative to light
  TH1D *lc = range05(l, Form("lc%d", uidP++)); lc->Scale(nl/lc->Integral());
  TH1D *cc = range05(c, Form("cc%d", uidP++)); cc->Scale(cMult*nc/cc->Integral());
  lc->Add(cc);

  t.bFC = range05(b, Form("bfc%d", uidP++));
  t.bGS = range05(g, Form("bgs%d", uidP++));
  t.lc  = lc;
  t.bFC->Scale(1./t.bFC->Integral()); t.bGS->Scale(1./t.bGS->Integral()); t.lc->Scale(1./t.lc->Integral());
  return t;
}

// b template with the bGS share set to s (share of all b)
TH1D* bWithGS(const Tpl &t, double s)
{
  TH1D *b = (TH1D*) t.bFC->Clone(Form("bmix%d", uidP++)); b->SetDirectory(nullptr);
  b->Scale(1. - s); b->Add(t.bGS, s);
  return b;
}

struct FitRes { double fb, efb, fGS, efGS; double chi2; int ndf; TH1D *compB, *compBGS, *compL, *total; };

// Fit data with N templates (2 or 3) over [lo, hi] (default the nominal range).
// The fractions are fractions of the jets INSIDE the fit range, as in
// templateFitter(). Components come back as counts per 0.1 GeV bin, normalised
// so that inside the fit range they add up to the data.
FitRes doFit(TH1D *data05, std::vector<TH1D*> tpl, double lo = fitLo, double hi = fitHi)
{
  FitRes r{}; const int N = (int) tpl.size();
  RooRealVar x(Form("x%d", uidP), "p_{T}^{rel}", 0., 5.); x.setBins(50);
  x.setRange("fitRange", lo, hi);
  RooDataHist dh(Form("dh%d", uidP), "", x, RooFit::Import(*data05));
  std::vector<RooDataHist*> th; std::vector<RooHistPdf*> tp;
  for(int i = 0; i < N; i++){
    th.push_back(new RooDataHist(Form("th%d_%d", uidP, i), "", x, RooFit::Import(*tpl[i])));
    tp.push_back(new RooHistPdf(Form("tp%d_%d", uidP, i), "", x, *th.back(), 0));
  }
  RooRealVar f1(Form("f1_%d", uidP), "", 0.5, 0., 1.), f2(Form("f2_%d", uidP), "", 0.1, 0., 1.);
  RooArgList pdfs, coefs;
  for(auto p : tp) pdfs.add(*p);
  coefs.add(f1); if(N == 3) coefs.add(f2);
  RooAddPdf model(Form("model%d", uidP++), "", pdfs, coefs);
  model.fitTo(dh, RooFit::SumW2Error(true), RooFit::Range("fitRange"), RooFit::PrintLevel(-1));

  r.fb = f1.getVal(); r.efb = f1.getError();
  if(N == 3){ r.fGS = f2.getVal(); r.efGS = f2.getError(); }

  // components on the data's normalisation inside the fit range
  const int b1 = data05->FindBin(lo + 1e-6), b2 = data05->FindBin(hi - 1e-6);
  double nD = data05->Integral(b1, b2);
  double w[3] = {f1.getVal(), N == 3 ? f2.getVal() : 0., 0.};
  w[N-1] = 1. - f1.getVal() - (N == 3 ? f2.getVal() : 0.);
  r.compB = (TH1D*) tpl[0]->Clone(Form("cb%d", uidP++)); r.compB->SetDirectory(nullptr); r.compB->Scale(w[0]*nD/tpl[0]->Integral(b1, b2));
  if(N == 3){ r.compBGS = (TH1D*) tpl[1]->Clone(Form("cg%d", uidP++)); r.compBGS->SetDirectory(nullptr); r.compBGS->Scale(w[1]*nD/tpl[1]->Integral(b1, b2)); }
  r.compL = (TH1D*) tpl[N-1]->Clone(Form("cl%d", uidP++)); r.compL->SetDirectory(nullptr); r.compL->Scale(w[N-1]*nD/tpl[N-1]->Integral(b1, b2));
  r.total = (TH1D*) r.compB->Clone(Form("ct%d", uidP++)); r.total->SetDirectory(nullptr);
  if(N == 3) r.total->Add(r.compBGS);
  r.total->Add(r.compL);

  // chi2 for display, inside the fit range: data errors plus the templates' MC
  // statistical errors
  r.chi2 = 0.; r.ndf = 0;
  for(int i = b1; i <= b2; i++){
    double d = data05->GetBinContent(i), ed = data05->GetBinError(i), m = r.total->GetBinContent(i);
    double e2 = ed*ed;
    for(int k = 0; k < N; k++) e2 += pow(w[k]*nD/tpl[k]->Integral(b1, b2)*tpl[k]->GetBinError(i), 2);
    if(e2 > 0.){ r.chi2 += (d - m)*(d - m)/e2; r.ndf++; }
  }
  r.ndf -= (N - 1);

  for(auto p : tp) delete p;
  for(auto h : th) delete h;
  return r;
}

