// Two-dimensional b-purity fit in muon ptRel x CSVv2 for pp AND the four PbPb classes, against the
// nominal CSVv2-cut purity (2026-11-12, user: "try doing a 2d fit of the pTrel / csvv2 distributions to
// extract the b-purity. Make a plot of the b-purity extracted this way as compared to the c/light-fixed
// nominal (with a csvv2 cut) method"). pp + PbPb counterpart of fitPtRelCSVv2_2D_pp_PF.C (included
// through fitPtRelCSVv2Cut_PF_PbPb.C): same models, likelihood (Poisson + Barlow-Beeston lite),
// categories and ptRel binning (0.2 GeV nominal, BJET_PTREL_BIN), fit 0-5 GeV.
//
// Per system and jet-pT window (80-300 GeV), ONE fit of the ptRel x CSVv2-category map of ALL
// muon-tagged jets; f_b = b fraction of all muon-tagged jets = the purity directly (no eps_b division):
//   "MC"     every flavor's split across the CSVv2 categories from MC (PYTHIA / PYTHIA+HYDJET), c / light
//            at the MC ratio: the plain 2D template fit, the same MC assumptions as the nominal
//   "bFree"  the b tagger shape free (data / MC b-tag efficiency measured by the fit), c and light MC
//   "full"   b free, c and light tagger shapes with Gaussian priors (sigma 0.3 on the log weights), c share free
// Templates and data are built exactly as the nominal fits build them: pp by buildWindows, PbPb by
// buildCls (DiJet / BJet / MuJet shapes, DiJet fractions, natural g->bb x 1.175, uncorrelated muon-jet
// pairs subtracted per event; f_b is then rescaled by (1 - fake fraction), so it is b / ALL muon-tagged
// jets like the nominal).
// NOMINAL (reference): the chain's purity files, CSVv2 loose cut + tagged-jet ptRel fit with c / light at
// MC, P = P_tag f_pass / eps_b(MC): bPuritySys_total_PFJets_pp_csvV2loose.root,
// bPuritySys_PFJets_PbPb_csvV2loose.root (h_bPurity_nominal, total systematic h_bPurity_sysAbs).
// An Asimov check of each 2D model (pseudo-data = the templates at the data's size) is printed.
//
// Output: rootFiles/bPurityResults/ptRelCSVv2_2D_PF_ppPbPb_<cats>_fit5.root
// Figure: figures/bPurity/PFJets_PbPb/csvV2/purity_2DvsCut_<cats>.pdf (pp, 50-80 ... 0-10% in a row,
//         purity on top, 2D / nominal below). Only the "MC" model is drawn (2026-11-12, user: "just use the
//         MC tagger shapes"); bFree and full are still fitted, printed and written to the output file.
//         NOMINAL 2D FIT = the HISTFACTORY fit since 2026-11-12 (user: "make the histfactory the nominal method"):
//         it is what purity_2DvsCut_<cats>.pdf draws (the 2026-11-12 deck, slide 5); the hand-written minuit fit
//         of the same model is the cross-check, purity_2DvsCut_minuit_<cats>.pdf. In the output file the nominal
//         is h_bPurity2D_MC_HistFactory_<sys>, minuit h_bPurity2D_MC_<sys>. This is the nominal of the 2D-fit
//         study only; the chain's purity is still the CSVv2-cut fit.
// Usage, from src/plots/bPurity/:  root -l -b -q 'fitPtRelCSVv2_2D_ppPbPb_PF.C("4cat")'
#include "fitPtRelCSVv2Cut_PF_PbPb.C"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"
#include "RooStats/ModelConfig.h"
#include "RooWorkspace.h"
#include "TMatrixDSym.h"

// HISTFACTORY fit of the "MC" model (2026-11-12, user: "add the HistFactory version as a cross-check"; the NOMINAL
// 2D fit since the same day, user: "make the histfactory the nominal method").
// The same templates and data, but the likelihood built by RooStats::HistFactory instead of by hand:
//   one channel per CSVv2 category, nFit ptRel bins each; two samples, b and non-b = c + light summed
//   (that IS c / light at the MC ratio), each scaled to the data total x its MC truth fraction, so the
//   MC tagger shapes enter through the per-channel content; free NormFactors mu_b, mu_nb in [0, 10];
//   ActivateStatError on both samples with StatErrorConfig(0, Gaussian): one constrained gamma per bin
//   on the total MC, width = its relative MC error at the MC-truth mix (no 5% threshold).
//   f_b = mu_b t_b / (mu_b t_b + mu_nb (1 - t_b)), t_b = MC b fraction; error from the (mu_b, mu_nb) covariance.
//   Cells with no MC in either template are left out (floor, no gamma, data emptied), as the hand-written
//   fit effectively does; the share of data in them is printed.
// Differences from the hand-written likelihood, both expected: the total normalization floats (same
// Poisson minimum), and Barlow-Beeston uses one gamma per bin with its width fixed at the MC mix,
// where the hand-written fit profiles each bin analytically at every f_b (cf. the 1D RooFit engine,
// within 0.012 of minuit except the thin 0-10% light template).
namespace p2dhf {
using namespace p2d;
struct HFResult { double pur = 0., epur = 0., dropped = 0.; int status = -1; };
HFResult fitHF(const Win &X, const double (*data)[NPR], int id){
  HFResult r;
  double tot[2] = {}, nD = 0.;
  for(int k = 0; k < K; k++) for(int j = 0; j < nFit; j++) nD += data[k][j];
  tot[0] = X.t[0].tot; tot[1] = X.t[1].tot + X.t[2].tot;
  const double tb = tot[0]/(tot[0] + tot[1]);
  if(tot[0] <= 0. || tot[1] <= 0. || nD <= 0.) return r;
  RooStats::HistFactory::Measurement meas(Form("m%d", id), "");
  meas.SetPOI("mu_b"); meas.SetLumi(1.); meas.SetLumiRelErr(0.01); meas.AddConstantParam("Lumi");
  for(int k = 0; k < K; k++){
    TH1D *hd = new TH1D(Form("hfD%d_%d", id, k), "", nFit, 0., nFit*ptRelBW), *hb = (TH1D*) hd->Clone(Form("hfB%d_%d", id, k)), *hn = (TH1D*) hd->Clone(Form("hfN%d_%d", id, k));
    for(TH1D *h : {hd, hb, hn}) h->SetDirectory(nullptr);   // OWNED by HistFactory (Sample::SetHisto, Channel::SetData): never delete them here
    const double sb = nD*tb/tot[0], sn = nD*(1. - tb)/tot[1];
    const double floorV = 1e-9*nD/(K*nFit);
    for(int j = 0; j < nFit; j++){
      const double tB = X.t[0].n[k][j], tN = X.t[1].n[k][j] + X.t[2].n[k][j];
      if(tB + tN <= 0.){
        // no MC in this cell: the hand-written fit clamps its expectation to a constant, i.e. the cell carries no
        // information on f_b. Here: a negligible floor with no stat gamma and the data cell emptied (a floor alone
        // lets data there drag the total normalization; Asimov toy: mu_nb 1.37 instead of 1.02)
        hb->SetBinContent(j + 1, floorV); hn->SetBinContent(j + 1, floorV); hd->SetBinContent(j + 1, 0.);
        r.dropped += data[k][j]; continue;
      }
      hd->SetBinContent(j + 1, data[k][j]);
      hb->SetBinContent(j + 1, sb*tB); hb->SetBinError(j + 1, sb*sqrt(X.t[0].e2[k][j]));
      hn->SetBinContent(j + 1, sn*tN); hn->SetBinError(j + 1, sn*sqrt(X.t[1].e2[k][j] + X.t[2].e2[k][j]));
    }
    RooStats::HistFactory::Channel ch(Form("csv%d", k));
    ch.SetData(hd);
    ch.SetStatErrorConfig(1e-6, "Gaussian");   // every bin with MC error gets a gamma; the floor cells (error 0) none
    RooStats::HistFactory::Sample b("b"), nb("nonb");
    b.SetHisto(hb); b.ActivateStatError(); b.AddNormFactor("mu_b", 1., 0., 10.);
    nb.SetHisto(hn); nb.ActivateStatError(); nb.AddNormFactor("mu_nb", 1., 0., 10.);
    ch.AddSample(b); ch.AddSample(nb);
    meas.AddChannel(ch);
  }
  RooWorkspace *ws = RooStats::HistFactory::HistoToWorkspaceFactoryFast::MakeCombinedModel(meas);
  if(!ws) return r;
  RooStats::ModelConfig *mc = (RooStats::ModelConfig*) ws->obj("ModelConfig");
  RooAbsData *dat = ws->data("obsData");
  RooRealVar *mb = ws->var("mu_b"), *mn = ws->var("mu_nb");
  if(mc && dat && mb && mn){
    std::unique_ptr<RooFitResult> fr(mc->GetPdf()->fitTo(*dat, RooFit::Save(true), RooFit::PrintLevel(-1), RooFit::Strategy(1),
                                     RooFit::GlobalObservables(*mc->GetGlobalObservables()), RooFit::Minimizer("Minuit2", "Migrad"), RooFit::Offset(true)));
    if(fr){
      r.status = fr->status();
      const double a = mb->getVal()*tb, c = mn->getVal()*(1. - tb), s = a + c;
      r.pur = a/s;
      if(fr->covQual() < 0){ r.epur = -1.; delete ws; return r; }   // no covariance: error not available
      TMatrixDSym V = fr->reducedCovarianceMatrix(RooArgList(*mb, *mn));
      const double ga = tb*c/(s*s), gc = -(1. - tb)*a/(s*s);   // d f_b / d mu_b, d f_b / d mu_nb
      r.epur = sqrt(TMath::Max(0., ga*ga*V(0, 0) + gc*gc*V(1, 1) + 2.*ga*gc*V(0, 1)));
    }
  }
  delete ws;
  return r;
}
} // namespace p2dhf

void fitPtRelCSVv2_2D_ppPbPb_PF(TString cats = "4cat", TString nomWP = "loose"){
  using namespace pbc;
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  setvbuf(stdout, nullptr, _IOLBF, 0);   // line-buffered, so progress shows in a log file
  if(!setCats(cats)){ printf("ERROR: categories are 4cat or 8cat\n"); return; }
  nFit = nPR;   // fit 0-5 GeV
  setJetCollection("PFJets"); initPlotStyle();
  TFile *fDpp = TFile::Open(Form("%s/%s", repo, dataName)), *fMpp = TFile::Open(Form("%s/%s", repo, mcName));
  TFile *fDpb = TFile::Open(Form("%s/%s", repo, dataPbPb));
  for(int s = 0; s < NS; s++) fS[s] = TFile::Open(Form("%s/%s", repo, Form(phFmt, sName[s])));
  TFile *fNpp = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_total_PFJets_pp_csvV2%s.root", repo, nomWP.Data()));
  TFile *fNpb = TFile::Open(Form("%s/rootFiles/systematics/bPuritySys_PFJets_PbPb_csvV2%s.root", repo, nomWP.Data()));
  if(!fDpp || !fMpp || !fDpb || !fS[0] || !fS[1] || !fS[2] || !fNpp || !fNpb){ printf("ERROR: inputs\n"); return; }
  const TString tag = cats + "_fit5" + fitEngineTag();

  const int NCOL = 5; const char *colName[NCOL] = {"pp", "50-80%", "30-50%", "10-30%", "0-10%"}, *colTag[NCOL] = {"pp", "C4", "C3", "C2", "C1"};
  const int colCls[NCOL] = {-1, 3, 2, 1, 0};   // pbc class index (0 = 0-10%)
  const int NM = 3; const int mIdx[NM] = {1, 2, 3};   // models[]: MC, bFree, full
  TH1D *hHF[NCOL]; int hfUid = 0;
  TH1D *hP[NM][NCOL], *hNom[NCOL], *hNomSys[NCOL], *hTruth[NCOL];
  double dev[NM][NCOL][NW] = {}; int ndf[NM][NCOL][NW] = {};
  for(int c = 0; c < NCOL; c++){
    auto bk = [&](const TString &n){ TH1D *h = new TH1D(n + "_" + colTag[c], "", NW, ptEdges + 3); h->SetDirectory(nullptr); return h; };
    for(int m = 0; m < NM; m++) hP[m][c] = bk(Form("h_bPurity2D_%s", models[mIdx[m]].name));
    hHF[c] = bk("h_bPurity2D_MC_HistFactory");
    hNom[c] = bk("h_bPurityNominal"); hNomSys[c] = bk("h_bPurityNominalSys"); hTruth[c] = bk("h_bFractionMC");
    TH1 *nP = c == 0 ? (TH1*) fNpp->Get("h_bPurity_nominal") : (TH1*) fNpb->Get(Form("h_bPurity_nominal_%s", colTag[c]));
    TH1 *nS = c == 0 ? (TH1*) fNpp->Get("h_bPurity_sysAbs") : (TH1*) fNpb->Get(Form("h_bPurity_sysAbs_%s", colTag[c]));
    if(!nP || !nS){ printf("ERROR: nominal purity %s missing\n", colTag[c]); return; }
    if(c == 0 && !buildWindows(fDpp, fMpp)) return;
    printf("\n==== %s, 2D ptRel x CSVv2 (%s), fit 0-5 GeV, %g GeV bins\n", colName[c], cats.Data(), ptRelBW);
    printf("  window   | nominal (%s cut)   | MC truth | 2D MC tagger       dev/ndf | MC, HistFactory (st) - minuit | 2D b free          dev/ndf | 2D all floated     dev/ndf | b-tag data/MC per cat (b free) | Asimov fit - truth: MC / bFree / full\n", nomWP.Data());
    for(int w = 0; w < NW; w++){
      Win X; double fakeFrac = 0.;
      if(c == 0) X = W[w];
      else { double nEffB = 0.; if(!buildCls(fDpb, colCls[c], w, X, nEffB)) return; fakeFrac = lastFakeFrac; }
      const int nb = nP->FindBin(0.5*(wLo[w] + wHi[w]));
      if(fabs(nP->GetBinLowEdge(nb) - wLo[w]) > 1e-6 || fabs(nP->GetBinLowEdge(nb + 1) - wHi[w]) > 1e-6){ printf("ERROR: nominal binning %s\n", colTag[c]); return; }
      hNom[c]->SetBinContent(w + 1, nP->GetBinContent(nb)); hNom[c]->SetBinError(w + 1, nP->GetBinError(nb));
      hNomSys[c]->SetBinContent(w + 1, nS->GetBinContent(nb)); hNomSys[c]->SetBinError(w + 1, nS->GetBinError(nb));
      hTruth[c]->SetBinContent(w + 1, X.truthB);
      printf("  %3.0f-%3.0f  | %.3f +- %.3f (sys %.3f) | %.3f    |", wLo[w], wHi[w], nP->GetBinContent(nb), nP->GetBinError(nb), nS->GetBinError(nb), X.truthB);
      Result rr[NM];
      for(int m = 0; m < NM; m++){
        rr[m] = fitWin(models[mIdx[m]], X, X.d);
        const double p = rr[m].pur*(1. - fakeFrac), e = rr[m].epur*(1. - fakeFrac);
        hP[m][c]->SetBinContent(w + 1, p); hP[m][c]->SetBinError(w + 1, e);
        dev[m][c][w] = rr[m].dev; ndf[m][c][w] = rr[m].ndf;
        printf(" %.3f +- %.3f (%d%d) %5.1f/%-3d |", p, e, rr[m].status, rr[m].hesse, rr[m].dev, rr[m].ndf);
        if(m == 0){   // HistFactory cross-check of the MC model
          const p2dhf::HFResult h = p2dhf::fitHF(X, X.d, hfUid++);
          const double ph = h.pur*(1. - fakeFrac), eh = h.epur*(1. - fakeFrac);
          hHF[c]->SetBinContent(w + 1, ph); hHF[c]->SetBinError(w + 1, eh);
          printf(" %.3f +- %.3f (%d) %+.3f%s |", ph, eh, h.status, ph - p, h.dropped > 0. ? Form(" [%.1f%% data in no-MC cells]", 100.*h.dropped/X.nData) : "");
        }
      }
      for(int k = 0; k < K; k++) printf(" %.2f", rr[1].sf[0][k]);
      // Asimov: the templates at the data's size, every model against the truth
      double one[NF][KMAX]; for(int F = 0; F < NF; F++) for(int k = 0; k < K; k++) one[F][k] = 1.;
      double pd[KMAX][NPR]; const double nSave = X.nData; double sd = 0.; for(int k = 0; k < K; k++) for(int j = 0; j < nPR; j++) sd += X.d[k][j];
      X.nData = sd; const double tr = asimov(X, one, pd); X.nData = nSave;
      printf("  |");
      for(int m = 0; m < NM; m++){ Result ra = fitWin(models[mIdx[m]], X, pd); printf(" %+.3f", ra.pur - tr); }
      { const p2dhf::HFResult ha = p2dhf::fitHF(X, pd, hfUid++); printf("  HistFactory %+.3f", ha.pur - tr); }
      printf("\n");
    }
  }
  {
    TFile *fo = TFile::Open(Form("%s/rootFiles/bPurityResults/ptRelCSVv2_2D_PF_ppPbPb_%s%s.root", repo, tag.Data(), ptBinTag()), "RECREATE");
    for(int c = 0; c < NCOL; c++){ for(int m = 0; m < NM; m++) hP[m][c]->Write(); hHF[c]->Write(); hNom[c]->Write(); hNomSys[c]->Write(); hTruth[c]->Write(); }
    TNamed info("info", Form("2D ptRel x CSVv2 (%s) b purity of all muon-tagged jets, fit 0-5 GeV, %g GeV ptRel bins; pp, C4..C1 = 50-80 ... 0-10%%. Nominal = CSVv2 %s cut purity files. pp data %s, templates %s; PbPb data %s, templates %s",
                             cats.Data(), ptRelBW, nomWP.Data(), dataName, mcName, dataPbPb, phFmt));
    info.Write(); fo->Close();
  }

  // ---- figure: 1x5 row, purity on top, 2D / nominal below
  const char *mLeg[NM]  = {"2D fit", "2D, b tagger shape free", "2D, all tagger shapes free"};
  // only the MC-tagger-shape model is drawn (2026-11-12, user); the others are fitted, printed and written
  const char *mHex[NM]  = {okabeHex[5], okabeHex[6], okabeHex[3]};
  const int   mMark[NM] = {markFilledSquare, markFilledDiamond, 34};
  const double mSize[NM] = {2.2, 2.9, 2.6};
  const double ptLo = 80., ptHi = 300., yLo = 0., yHi = 1.2, rLo = 0.5, rHi = 1.5;
  // one figure: the nominal (CSVv2 cut, syst. band), the MC b fraction and ONE 2D series h2[c]; outName "" or "_HistFactory"
  auto drawFig = [&](TH1D **h2, const char *leg2, const char *hex2, int mark2, double size2, const TString &outName){
    TCanvas *cv = new TCanvas("cv2D_" + outName, "", 2300, 1000);
    GridPads g = makeGridPads(cv, NCOL, 2, 0.055, 0.01, 0.12, 0.075);
    const double split = g.bm + 0.33*(1. - g.bm - g.tm);
    for(int c = 0; c < NCOL; c++){
      TPad *t = g.at(c, 0), *b = g.at(c, 1); double x1, y1, x2, y2;
      t->GetPadPar(x1, y1, x2, y2); t->SetPad(x1, split, x2, y2); t->SetTopMargin(g.tm/(1. - split));
      b->GetPadPar(x1, y1, x2, y2); b->SetPad(x1, 0., x2, split); b->SetBottomMargin(g.bm/split);
    }
    const double fs = 0.055;
    for(int c = 0; c < NCOL; c++){
      g.at(c, 0)->cd();
      TH1F *fr = gPad->DrawFrame(ptLo, yLo, ptHi, yHi); g.styleFrame(fr, c, 0, fs);
      fr->GetYaxis()->SetLabelFont(43); fr->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0);
      fr->GetYaxis()->ChangeLabel(1, -1, 0.); fr->GetYaxis()->SetNdivisions(507);
      hNomSys[c]->SetFillColorAlpha(kGray + 1, 0.45); hNomSys[c]->SetLineWidth(0); hNomSys[c]->SetMarkerSize(0); hNomSys[c]->Draw("E2 SAME");
      styleLine(hTruth[c], "#8C8C8C", 4); hTruth[c]->SetLineStyle(7); hTruth[c]->Draw("HIST SAME");
      styleH(h2[c], hex2, mark2, size2); h2[c]->SetLineWidth(3); h2[c]->Draw("E1 X0 SAME");
      styleH(hNom[c], hexData, markOpenCircle, 3.0); hNom[c]->SetLineWidth(3); hNom[c]->Draw("E1 X0 SAME");
      gPad->RedrawAxis();
      g.text(c, 0, 0.10, colName[c], 0.065, 62);
      if(c == 0){   // in the empty top of the pp panel
        const double fL = gPad->GetLeftMargin(), fW = 1. - fL - gPad->GetRightMargin(), fT = 1. - gPad->GetTopMargin();
        TLegend *lg = new TLegend(fL + 0.04*fW, fT - 0.38, fL + fW, fT - 0.14);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextFont(43); lg->SetTextSize(30);
        if(c == 0){ lg->AddEntry(hNom[c], Form("CSVv2 > %s", nomWP == "loose" ? "0.54" : "0.84"), "pe"); lg->AddEntry(hTruth[c], "MC", "l"); }
        lg->AddEntry(h2[c], leg2, "pe");
        lg->Draw();
      }
      g.at(c, 1)->cd();
      TH1F *fb = gPad->DrawFrame(ptLo, rLo, ptHi, rHi); g.styleFrame(fb, c, 1, fs);
      fb->GetXaxis()->SetLabelFont(43); fb->GetXaxis()->SetLabelSize(34);
      fb->GetYaxis()->SetLabelFont(43); fb->GetYaxis()->SetLabelSize(c == 0 ? 34 : 0);
      fb->GetYaxis()->SetNdivisions(504); fb->GetXaxis()->SetTickLength(0.03*0.67/0.33);
      TLine *one = new TLine(ptLo, 1., ptHi, 1.); one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->SetLineWidth(2); one->Draw();
      {
        TH1D *r = (TH1D*) h2[c]->Clone("r2D_" + outName + "_" + colTag[c]); r->SetDirectory(nullptr);
        for(int w = 1; w <= NW; w++){
          const double d = hNom[c]->GetBinContent(w), v = d > 0. ? h2[c]->GetBinContent(w)/d : -999.;
          r->SetBinContent(w, v); r->SetBinError(w, 0.);
          if(d > 0. && (v < rLo || v > rHi)) printf("NOTE: %s %s %.0f-%.0f GeV 2D / nominal %.2f is off scale (not drawn)\n", leg2, colName[c], wLo[w-1], wHi[w-1], v);
        }
        r->Draw("P X0 SAME");
      }
      gPad->RedrawAxis();
    }
    cv->cd();
    { TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.040); t.SetTextAlign(21);
      t.DrawLatex(g.lm + (1. - g.lm - g.rm)/2., 0.018, "#it{p}_{T}^{jet} (GeV)");
      t.SetTextAngle(90);
      t.DrawLatex(0.016, split + 0.5*(1. - g.tm - split), "b purity of muon-tagged jets");
      t.SetTextSize(0.028); t.DrawLatex(0.016, g.bm + 0.5*(split - g.bm), "2D / CSVv2 cut"); }
    g.cmsLabel(cv, "Internal", Form("pp, PbPb (5.02 TeV)   #it{p}_{T}^{rel} #times CSVv2, %d categories", K), 1.0);
    const TString out = Form("%s/figures/bPurity/PFJets_PbPb/csvV2/purity_2DvsCut%s_%s%s.pdf", repo, outName.Data(), tag.Data(), ptBinTag());
    savePdfTight(cv, out);
  };
  // NOMINAL 2D fit = HistFactory since 2026-11-12 (user: "make the histfactory the nominal method"); the 2026-11-12 deck, slide 5
  drawFig(hHF, mLeg[0], mHex[0], mMark[0], mSize[0], "");
  // the hand-written minuit likelihood of the same model, cross-check figure
  drawFig(hP[0], "2D fit (minuit)", okabeHex[6], markFilledDiamond, 2.9, "_minuit");
}
