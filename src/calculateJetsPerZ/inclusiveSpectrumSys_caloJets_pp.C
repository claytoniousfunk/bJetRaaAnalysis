// Spectrum systematics of the pp CALO-jet INCLUSIVE jet spectrum (jets per Z),
// standalone: the same three spectrum-level sources as the b-jet chain
// (unfoldBJetSpectrum_caloJets_pp.C) --
//   iterations     max(|U(N-1) - U(N)|, |U(N+1) - U(N)|) / U(N)
//   response stat  the response matrix's MC statistics, by toys: every cell
//                  fluctuated as Poisson(neff) x sigma^2/c, neff = (c/sigma)^2,
//                  prior / reco floor / measured projection / RooUnfoldResponse
//                  rebuilt from it, data unfolded; RMS over toys / nominal
//   spectrum JEU   the PYTHIA calo reco spectra with the jet pT shifted up/down
//                  by the JEU (stored by src/unfoldTest/unfoldPYTHIAJetPt_JEUShift.C)
//                  unfolded with THIS response, floor and iteration count;
//                  max(|up/nom - 1|, |down/nom - 1|) at gen level
// -- added in quadrature.
//
// The unfolding reproduces calculateJetsPerZ.cc's pp calo path without
// touching that file:
//   measured  jetsPerZ_fine_pp from its NO-unfolding output
//             (jetsPerZSpectra_caloJets_coarseMuEff.root): trigger-stitched,
//             per-Z normalised, reco floor already applied, native 5 GeV bins.
//             Bayesian unfolding is linear in the input scale, so unfolding the
//             per-Z spectrum gives the unfolded spectrum per Z.
//   response  h_matchedRecoJetPt_genJetPt_allJets, 2026-09-22 calo response
//             (the file calculateJetsPerZ.cc opens), truth from the untruncated
//             matrix, reco bins with upper edge <= 60 GeV zeroed.
//   Unfolded in the native 5 GeV bins, then rebinned for the outputs. The N = 1
//   unfolding is checked against calculateJetsPerZ.cc's own N = 1 output
//   (jetsPerZ_fine_pp in ..._unfold1iter.root) and the largest difference is
//   printed.
//
// TWO NOMINALS, both reported (the choice is still open):
//   A  N = 8   unfoldClosureTest.C("pp",15,"weighted","tilt:0.5",80.,500.,true,
//              1,"","","_inclMetric80to500",-1,"","",false,60.) minimum MSE
//              (2026-10-04: MSE flat to ~6% for N = 5-15); varied to 7 and 9.
//   B  N = 1   calculateJetsPerZ.cc's current value, which that closure shows is
//              under-converged (MSE 21x the minimum); varied to 2 only -- N = 0
//              is no unfolding, so the iteration systematic is one-sided.
//
// OUTPUT BINNINGS: "incl" = the inclusive spectrum's own axis (calculateJetsPerZ
// newAxis, 60-500 GeV), and "bjet" = the b-jet chain's axis (ptEdges in
// caloJetBPurityFit_pp.h, 50-500 GeV) for the b / inclusive ratio. Every
// variation is unfolded in native bins and rebinned before the relative
// systematic is formed.
//
// CAVEATS: the spectrum-JEU reco spectra come from the 2026-10-02 PYTHIA scan,
// not the response sample, and use the AK4PF uncertainty file (no pp AK4Calo
// one exists). Stat errors are rebinned in quadrature, ignoring bin-to-bin
// correlations, as calculateJetsPerZ.cc does.
//
// Output  rootFiles/systematics/inclJetSpectrumSys_caloJets_pp.root
//           h_incl_unfolded_<cfg>_<bins>          unfolded jets per Z per 5 GeV-bin
//                                                 content rebinned (counts per Z), stat
//           h_incl_statRel_<cfg>_<bins>
//           h_incl_sysRel_{iterations,responseStat,spectrumJEU,total}_<cfg>_<bins>
//           h_incl_sysRel_spectrumJEU_{up,down}_<cfg>_<bins>  signed, for correlated ratios
//           h_incl_respToys_meanShift_<cfg>_<bins>  toys' mean / nominal - 1 (check)
//           <cfg> = N8, N1, N2;  <bins> = incl, bjet
//         figures/systematics/inclJetSpectrumSys_caloJets_pp_N8.pdf, _N1.pdf, _N2.pdf
//         figures/systematics/inclJetSpectrum_N8overN1_caloJets_pp.pdf
//         figures/systematics/inclJetSpectrum_caloJets_pp.pdf   the N = 8 spectrum, d2N/dpT deta (not per Z)
//
// Usage: root -l -b -q -e 'gSystem->Load("/home/clayton/Programs/RooUnfold/build/libRooUnfold.so");' \
//                         inclusiveSpectrumSys_caloJets_pp.C
// Run from: src/calculateJetsPerZ/

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "TFile.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TPad.h"
#include "TBox.h"
#include "TColor.h"
#include "../../headers/plotting/plotStyle.h"
#include "../../headers/functions/caloJetBPurityFit_pp.h"   // ptEdges (b-jet binning), ptReportMin

const char *measPath   = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateJetsPerZ/rootFiles/JetsPerZ/jetsPerZSpectra_caloJets_coarseMuEff.root";
const char *checkPath  = "/home/clayton/Analysis/code/bJetRaaAnalysis/src/calculateJetsPerZ/rootFiles/JetsPerZ/jetsPerZSpectra_caloJets_coarseMuEff_unfold1iter.root";
const char *respPath   = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/scanningOuput/PYTHIA/PYTHIA_DiJet_response_caloJets_manualJEC_pThat-15_mu12_pTmu-15_tight_jetTrkMaxFilter_2026-9-22.root";
const char *respName   = "h_matchedRecoJetPt_genJetPt_allJets";
const char *jeuPath    = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/JEU/unfold_PYTHIA_caloJets_JEUShift.root";
const double recoFloorIncl = 60.;

const char *outRoot = "/home/clayton/Analysis/code/bJetRaaAnalysis/rootFiles/systematics/inclJetSpectrumSys_caloJets_pp.root";
const char *figDir  = "/home/clayton/Analysis/code/bJetRaaAnalysis/figures/systematics";

const int nToysIncl = 1000;
const int toySeedIncl = 20261005;

// the inclusive spectrum's own axis (calculateJetsPerZ.cc newAxis)
const int NInclE = 14;
const double inclEdges[NInclE] = {60,70,80,90,100,110,130,150,180,200,240,280,350,500};

struct Cfg { const char *tag; int nNom; int nLo; int nHi; };   // nLo/nHi < 1 => not used
const int NCfg = 3;
const Cfg cfgs[NCfg] = { {"N8", 8, 7, 9}, {"N1", 1, 0, 2}, {"N2", 2, 1, 3} };   // N2: same N as the b-jet unfolding, for the b-jet fraction

struct Binning { const char *tag; int n; const double *edges; };

TH2D* floorResponse(TH2D *r, const char *name)
{
  TH2D *t = (TH2D*) r->Clone(name); t->SetDirectory(nullptr);
  for(int ix = 0; ix <= t->GetNbinsX()+1; ix++){
    if(t->GetXaxis()->GetBinUpEdge(ix) > recoFloorIncl + 1e-6) continue;
    for(int iy = 0; iy <= t->GetNbinsY()+1; iy++){ t->SetBinContent(ix, iy, 0.); t->SetBinError(ix, iy, 0.); }
  }
  return t;
}

// counts summed into the edges of a binning (native bins are 5 GeV, all edges on that grid)
TH1D* rebinCounts(TH1D *h, const Binning &b, const char *name)
{
  TH1D *o = new TH1D(name, h->GetTitle(), b.n - 1, b.edges); o->SetDirectory(nullptr); o->Sumw2();
  for(int i = 1; i <= h->GetNbinsX(); i++){
    int j = o->FindBin(h->GetBinCenter(i)); if(j < 1 || j > b.n - 1) continue;
    o->SetBinContent(j, o->GetBinContent(j) + h->GetBinContent(i));
    o->SetBinError(j, sqrt(pow(o->GetBinError(j), 2) + pow(h->GetBinError(i), 2)));
  }
  return o;
}

void inclusiveSpectrumSys_caloJets_pp()
{
  initPlotStyle();
  gSystem->mkdir(gSystem->DirName(outRoot), kTRUE);
  gSystem->mkdir(figDir, kTRUE);

  // ---- inputs -----------------------------------------------------------------
  TFile *fM = TFile::Open(measPath), *fR = TFile::Open(respPath), *fJ = TFile::Open(jeuPath), *fC = TFile::Open(checkPath);
  if(!fM || !fR || !fJ){ printf("ERROR: cannot open inputs\n"); return; }
  TH1D *hMeas = nullptr, *jNom = nullptr, *jUp = nullptr, *jDn = nullptr, *hCheck = nullptr;
  TH2D *hResp = nullptr;
  fM->GetObject("jetsPerZ_fine_pp", hMeas);
  fR->GetObject(respName, hResp);
  fJ->GetObject("reco_nominal", jNom); fJ->GetObject("reco_up", jUp); fJ->GetObject("reco_down", jDn);
  if(fC) fC->GetObject("jetsPerZ_fine_pp", hCheck);
  if(!hMeas || !hResp || !jNom || !jUp || !jDn){ printf("ERROR: missing input histograms\n"); return; }
  hMeas = (TH1D*) hMeas->Clone("hMeas"); hMeas->SetDirectory(nullptr);
  // the floor is already applied in the file; apply it again so the check is exact either way
  for(int i = 0; i <= hMeas->GetNbinsX()+1; i++)
    if(hMeas->GetXaxis()->GetBinUpEdge(i) <= recoFloorIncl + 1e-6){ hMeas->SetBinContent(i, 0.); hMeas->SetBinError(i, 0.); }

  TH1D *truthMC = hResp->ProjectionY("truthMC", 1, hResp->GetNbinsX()); truthMC->SetDirectory(nullptr);
  TH2D *respF = floorResponse(hResp, "respF");
  TH1D *measMC = respF->ProjectionX("measMC", 1, respF->GetNbinsY()); measMC->SetDirectory(nullptr);
  RooUnfoldResponse resp(measMC, truthMC, respF, "respIncl", "pp calo inclusive response");

  auto unfold = [&](RooUnfoldResponse &r, TH1D *meas, int n, const char *name){
    RooUnfoldBayes u(&r, meas, n); u.SetVerbose(0);
    TH1D *h = (TH1D*) u.Hunfold(); h->SetName(name); h->SetDirectory(nullptr); return h;
  };

  // ---- check against calculateJetsPerZ.cc's N = 1 ------------------------------
  {
    TH1D *u1 = unfold(resp, hMeas, 1, "check_N1");
    if(hCheck){
      double maxRel = 0.;
      for(int i = 1; i <= u1->GetNbinsX(); i++){
        double a = u1->GetBinContent(i), b = hCheck->GetBinContent(i);
        if(u1->GetBinLowEdge(i) >= 60. && b != 0.) maxRel = TMath::Max(maxRel, fabs(a/b - 1.));
      }
      printf("\ncheck: N = 1 unfolding vs calculateJetsPerZ.cc's stored N = 1 spectrum, max |rel. diff| above 60 GeV = %.2e\n", maxRel);
    } else printf("\ncheck: %s not found, N = 1 reproduction not verified\n", checkPath);
  }

  const Binning bins[2] = { {"incl", NInclE, inclEdges}, {"bjet", NPt + 1, ptEdges} };

  TFile *fo = TFile::Open(outRoot, "recreate");
  TString infoTxt;

  for(int c = 0; c < NCfg; c++){
    const Cfg &cf = cfgs[c];
    TH1D *uNom = unfold(resp, hMeas, cf.nNom, Form("uNom_%s", cf.tag));
    TH1D *uLo  = cf.nLo >= 1 ? unfold(resp, hMeas, cf.nLo, Form("uLo_%s", cf.tag)) : nullptr;
    TH1D *uHi  = unfold(resp, hMeas, cf.nHi, Form("uHi_%s", cf.tag));
    // spectrum JEU: unfold the shifted PYTHIA reco spectra with this response and N
    TH1D *mJ[3] = {jNom, jUp, jDn}, *uJ[3];
    for(int k = 0; k < 3; k++){
      TH1D *m = (TH1D*) mJ[k]->Clone(Form("mJ%d_%s", k, cf.tag)); m->SetDirectory(nullptr);
      for(int i = 0; i <= m->GetNbinsX()+1; i++) if(m->GetXaxis()->GetBinUpEdge(i) <= recoFloorIncl + 1e-6){ m->SetBinContent(i, 0.); m->SetBinError(i, 0.); }
      uJ[k] = unfold(resp, m, cf.nNom, Form("uJ%d_%s", k, cf.tag));
    }
    // response toys, native bins
    std::vector<TH1D*> toys;
    {
      TRandom3 rnd(toySeedIncl + c);
      for(int t = 0; t < nToysIncl; t++){
        TH2D *r2 = (TH2D*) hResp->Clone(Form("toyR_%s_%d", cf.tag, t)); r2->SetDirectory(nullptr);
        for(int ix = 1; ix <= r2->GetNbinsX(); ix++) for(int iy = 1; iy <= r2->GetNbinsY(); iy++){
          double cc = hResp->GetBinContent(ix, iy), e = hResp->GetBinError(ix, iy);
          if(cc <= 0. || e <= 0.) continue;
          double neff = (cc/e)*(cc/e);
          r2->SetBinContent(ix, iy, rnd.Poisson(neff)*e*e/cc);
        }
        TH1D *tTruth = r2->ProjectionY(Form("toyT_%s_%d", cf.tag, t), 1, r2->GetNbinsX()); tTruth->SetDirectory(nullptr);
        TH2D *rF = floorResponse(r2, Form("toyRF_%s_%d", cf.tag, t));
        TH1D *tMeas = rF->ProjectionX(Form("toyM_%s_%d", cf.tag, t), 1, rF->GetNbinsY()); tMeas->SetDirectory(nullptr);
        RooUnfoldResponse rt(tMeas, tTruth, rF, Form("toyRoo_%s_%d", cf.tag, t), "");
        toys.push_back(unfold(rt, hMeas, cf.nNom, Form("toyU_%s_%d", cf.tag, t)));
        delete r2; delete rF; delete tTruth; delete tMeas;
      }
    }

    for(int b = 0; b < 2; b++){
      const Binning &bn = bins[b];
      TString sfx = Form("%s_%s", cf.tag, bn.tag);
      TH1D *N = rebinCounts(uNom, bn, Form("h_incl_unfolded_%s", sfx.Data()));
      N->SetTitle(Form("pp calo inclusive jets per Z, unfolded N = %d;#it{p}_{T}^{jet} [GeV];jets per Z per bin", cf.nNom));
      TH1D *L = uLo ? rebinCounts(uLo, bn, Form("rL_%s", sfx.Data())) : nullptr, *H = rebinCounts(uHi, bn, Form("rH_%s", sfx.Data()));
      TH1D *J[3]; for(int k = 0; k < 3; k++) J[k] = rebinCounts(uJ[k], bn, Form("rJ%d_%s", k, sfx.Data()));
      std::vector<TH1D*> T; for(size_t t = 0; t < toys.size(); t++) T.push_back(rebinCounts(toys[t], bn, Form("rT%zu_%s", t, sfx.Data())));

      const int nb = bn.n - 1;
      auto book = [&](const char *what, const char *ttl){
        TH1D *h = new TH1D(Form("h_incl_%s_%s", what, sfx.Data()), Form("%s;#it{p}_{T}^{jet} [GeV];relative uncertainty", ttl), nb, bn.edges);
        h->SetDirectory(nullptr); return h; };
      TH1D *sIt = book("sysRel_iterations",   Form("iterations, N = %d%s", cf.nNom, uLo ? " +/- 1" : " vs N+1 (one-sided)"));
      TH1D *sRs = book("sysRel_responseStat", Form("response-matrix MC stat. (RMS of %d toys)", nToysIncl));
      TH1D *sJE = book("sysRel_spectrumJEU",  "spectrum JEU (gen level)");
      TH1D *sJu = book("sysRel_spectrumJEU_up",   "spectrum JEU up: shifted/nominal - 1 (gen level, signed)");
      TH1D *sJd = book("sysRel_spectrumJEU_down", "spectrum JEU down: shifted/nominal - 1 (gen level, signed)");
      TH1D *sTo = book("sysRel_total",        "total (quadrature)");
      TH1D *sSt = book("statRel",             "stat. (propagated through the unfolding)");
      TH1D *sMS = book("respToys_meanShift",  "response toys: mean / nominal - 1 (check)");

      printf("\n=== %s, %s binning ===\n%-9s %11s %8s | %8s %8s %8s | %8s  %s\n", cf.tag, bn.tag, "jet pT", "jets/Z", "stat",
             "iter", "respMC", "specJEU", "TOTAL", "(toy mean shift)");
      for(int i = 1; i <= nb; i++){
        double v = N->GetBinContent(i);
        if(v <= 0.){ printf("%3.0f-%-5.0f  empty (below the reco floor)\n", bn.edges[i-1], bn.edges[i]); continue; }
        double it = fabs(H->GetBinContent(i)/v - 1.); if(L) it = TMath::Max(it, fabs(L->GetBinContent(i)/v - 1.));
        double jn = J[0]->GetBinContent(i);
        double ju = jn > 0. ? J[1]->GetBinContent(i)/jn - 1. : 0., jd = jn > 0. ? J[2]->GetBinContent(i)/jn - 1. : 0.;
        double je = TMath::Max(fabs(ju), fabs(jd));
        sJu->SetBinContent(i, ju); sJd->SetBinContent(i, jd);
        double s1 = 0., s2 = 0.; for(auto t : T){ double x = t->GetBinContent(i)/v; s1 += x; s2 += x*x; }
        double mean = s1/T.size(), rs = sqrt(TMath::Max(0., s2/T.size() - mean*mean));
        double to = sqrt(it*it + je*je + rs*rs);
        sIt->SetBinContent(i, it); sRs->SetBinContent(i, rs); sJE->SetBinContent(i, je); sTo->SetBinContent(i, to);
        sSt->SetBinContent(i, N->GetBinError(i)/v); sMS->SetBinContent(i, mean - 1.);
        printf("%3.0f-%-5.0f %11.4g %7.2f%% | %7.2f%% %7.2f%% %7.2f%% | %7.2f%%  (%+.2f%%)\n", bn.edges[i-1], bn.edges[i], v,
               100*sSt->GetBinContent(i), 100*it, 100*rs, 100*je, 100*to, 100*(mean - 1.));
      }
      fo->cd();
      N->Write(); sSt->Write(); sIt->Write(); sRs->Write(); sJE->Write(); sJu->Write(); sJd->Write(); sTo->Write(); sMS->Write();

      // breakdown figure, inclusive binning only
      if(b == 0){
        TCanvas *cv = new TCanvas(Form("cv_%s", cf.tag), "", 800, 700);
        cv->SetLeftMargin(0.13); cv->SetBottomMargin(0.13); cv->SetRightMargin(0.04); cv->SetTopMargin(0.06);
        TH1D *lT = (TH1D*) sTo->Clone(Form("lT_%s", cf.tag)), *lI = (TH1D*) sIt->Clone(Form("lI_%s", cf.tag));
        TH1D *lR = (TH1D*) sRs->Clone(Form("lR_%s", cf.tag)), *lJ = (TH1D*) sJE->Clone(Form("lJ_%s", cf.tag)), *lS = (TH1D*) sSt->Clone(Form("lS_%s", cf.tag));
        for(TH1D *h : {lT, lI, lR, lJ, lS}){ h->Scale(100.); for(int i = 0; i <= nb+1; i++) h->SetBinError(i, 0.); }
        double ymax = TMath::Max(lT->GetMaximum(), lS->GetMaximum());
        TH1F *fr = cv->DrawFrame(bn.edges[0], 0., bn.edges[nb], 1.9*ymax);   // headroom for header + legend
        fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]");
        fr->GetYaxis()->SetTitle("relative uncertainty [%]");
        fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
        fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
        styleLine(lI, "#E69F00", 3); lI->SetLineStyle(7);
        styleLine(lR, "#CC79A7", 3); lR->SetLineStyle(9);
        styleLine(lJ, "#0072B2", 3); lJ->SetLineStyle(3);
        styleLine(lS, "#009E73", 3); lS->SetLineStyle(2);
        styleLine(lT, "#000000", 4);
        // JEU drawn last: it is nearly the whole total and would otherwise hide under it
        lS->Draw("HIST same"); lI->Draw("HIST same"); lR->Draw("HIST same"); lT->Draw("HIST same"); lJ->SetLineWidth(4); lJ->Draw("HIST same");
        TLegend *leg = makeLegend(0.40, 0.58, 0.95, 0.85, 0.032);
        leg->AddEntry(lJ, "spectrum JEU (gen level)", "l");
        leg->AddEntry(lI, uLo ? Form("iterations (#it{N} = %d #pm 1)", cf.nNom) : Form("iterations (#it{N} = %d vs %d, one-sided)", cf.nNom, cf.nHi), "l");
        leg->AddEntry(lR, "response-matrix MC stat.", "l");
        leg->AddEntry(lT, "total systematic (quadrature)", "l");
        leg->AddEntry(lS, "stat. (for reference)", "l");
        leg->Draw();
        TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
        la.DrawLatex(0.17, 0.885, Form("pp 5.02 TeV, calo jets: inclusive jet spectrum, #it{N} = %d", cf.nNom));
        savePdfTight(cv, Form("%s/inclJetSpectrumSys_caloJets_pp_%s.pdf", figDir, cf.tag));
        delete cv;
      }
    }
    infoTxt += Form("%s: N = %d (variations %s%d); ", cf.tag, cf.nNom, cf.nLo >= 1 ? Form("%d, ", cf.nLo) : "", cf.nHi);
    for(auto t : toys) delete t;
  }

  // ---- figure: N = 8 / N = 1 -----------------------------------------------------
  {
    TH1D *a = (TH1D*) fo->Get("h_incl_unfolded_N8_incl"), *b = (TH1D*) fo->Get("h_incl_unfolded_N1_incl");
    if(a && b){
      TH1D *r = (TH1D*) a->Clone("ratio81"); r->SetDirectory(nullptr); r->Divide(b);
      for(int i = 1; i <= r->GetNbinsX(); i++) r->SetBinError(i, 0.);
      TCanvas *cv = new TCanvas("cv81", "", 800, 600);
      cv->SetLeftMargin(0.13); cv->SetBottomMargin(0.13); cv->SetRightMargin(0.04); cv->SetTopMargin(0.06);
      double lo = 1e9, hi = -1e9; for(int i = 1; i <= r->GetNbinsX(); i++) if(r->GetBinContent(i) > 0.){ lo = TMath::Min(lo, r->GetBinContent(i)); hi = TMath::Max(hi, r->GetBinContent(i)); }
      double pad = 0.2*(hi - lo + 0.02);
      TH1F *fr = cv->DrawFrame(inclEdges[0], TMath::Min(lo - pad, 0.98), inclEdges[NInclE-1], TMath::Max(hi + 2.5*pad, 1.02));
      fr->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); fr->GetYaxis()->SetTitle("unfolded #it{N} = 8 / #it{N} = 1");
      fr->GetXaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetTitleSize(0.048);
      fr->GetXaxis()->SetLabelSize(0.040); fr->GetYaxis()->SetLabelSize(0.040);
      TLine *one = new TLine(inclEdges[0], 1., inclEdges[NInclE-1], 1.); one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();
      styleH(r, "#000000", markFilledCircle, 1.2);
      r->Draw("P same");
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.038);
      la.DrawLatex(0.17, 0.88, "pp 5.02 TeV, calo jets: inclusive jet spectrum");
      la.SetTextSize(0.032);
      la.DrawLatex(0.17, 0.83, "closure-test optimum (#it{N} = 8) vs calculateJetsPerZ.cc (#it{N} = 1)");
      savePdfTight(cv, Form("%s/inclJetSpectrum_N8overN1_caloJets_pp.pdf", figDir));
      delete cv;
    }
  }

  // ---- figure: the spectrum itself ---------------------------------------------
  // d^2N / dpT deta, NOT per Z: the per-Z input is multiplied back by the pp Z
  // count calculateJetsPerZ.cc divided it by (NZ[0] in its output file), then
  // divided by the pT bin width and by 3.2 for |eta| < 1.6. The stitched spectrum
  // is in Jet100-sample units (calculateJetsPerZ.cc, stitchSamples). Nominal
  // N = 8 only, reported range only (from ptReportMin = 80 GeV, as the b-jet
  // figures; the 60-80 GeV bins stay in the output file).
  {
    TH1D *c8 = (TH1D*) fo->Get("h_incl_unfolded_N8_incl");
    TH1D *s8 = (TH1D*) fo->Get("h_incl_sysRel_total_N8_incl");
    TH1D *t8 = (TH1D*) fo->Get("h_incl_statRel_N8_incl");
    TVectorD *nzv = nullptr; fM->GetObject("NZ", nzv);
    if(c8 && s8 && t8 && nzv){
      const double nzPP = (*nzv)[0];
      printf("\nspectrum figure: per-Z spectrum x N_Z(pp) = %.1f (undoes calculateJetsPerZ.cc's normalisation)\n", nzPP);
      TH1D *d8 = (TH1D*) c8->Clone("spec8"); d8->SetDirectory(nullptr);
      for(int i = 1; i <= d8->GetNbinsX(); i++){
        double w = d8->GetBinWidth(i)*3.2;
        d8->SetBinContent(i, d8->GetBinContent(i)*nzPP/w); d8->SetBinError(i, d8->GetBinError(i)*nzPP/w);
      }
      TCanvas *cv = new TCanvas("cvSpec", "", 700, 800);
      TPad *pTop = new TPad("pTopS", "", 0, 0.34, 1, 1), *pBot = new TPad("pBotS", "", 0, 0, 1, 0.34);
      pTop->SetBottomMargin(0.02); pTop->SetLeftMargin(0.17); pTop->SetTopMargin(0.07);
      pBot->SetTopMargin(0.02);    pBot->SetLeftMargin(0.17); pBot->SetBottomMargin(0.32);
      pTop->Draw(); pBot->Draw();
      pTop->cd(); pTop->SetLogy();   // log: the spectrum falls by ~4 orders of magnitude, no negative bins
      auto shown = [&](int i){ return d8->GetBinLowEdge(i) >= ptReportMin - 1e-6; };
      double ymax = 0., ymin = 1e30;
      for(int i = 1; i <= d8->GetNbinsX(); i++) if(shown(i) && d8->GetBinContent(i) > 0.){ ymax = TMath::Max(ymax, d8->GetBinContent(i)); ymin = TMath::Min(ymin, d8->GetBinContent(i)); }
      TH1D *fr = (TH1D*) d8->Clone("frSpec"); fr->Reset(); fr->SetTitle("");
      fr->SetMinimum(0.3*ymin); fr->SetMaximum(40.*ymax);
      fr->GetXaxis()->SetLabelSize(0);
      fr->GetXaxis()->SetRangeUser(ptReportMin, inclEdges[NInclE-1]);
      fr->GetYaxis()->SetTitle("d^{2}#it{N}_{jet} / d#it{p}_{T} d#eta [GeV^{-1}]");
      fr->GetYaxis()->SetTitleSize(0.050); fr->GetYaxis()->SetTitleOffset(1.45); fr->GetYaxis()->SetLabelSize(0.045);
      fr->Draw("AXIS");
      int colSys = TColor::GetColor("#D55E00");
      for(int i = 1; i <= d8->GetNbinsX(); i++){
        double v = d8->GetBinContent(i); if(v <= 0. || !shown(i)) continue;
        double e = s8->GetBinContent(i)*v;
        TBox *bx = new TBox(d8->GetBinLowEdge(i), v - e, d8->GetBinLowEdge(i+1), v + e);
        bx->SetFillColorAlpha(colSys, 0.30); bx->SetLineColor(colSys); bx->Draw("l"); bx->Draw();
      }
      styleH(d8, "#000000", markFilledCircle, 1.2);
      for(int i = 1; i <= d8->GetNbinsX(); i++) if(!shown(i)){ d8->SetBinContent(i, -999.); d8->SetBinError(i, 0.); }
      d8->Draw("E1 X0 same");
      TLegend *leg = makeLegend(0.22, 0.05, 0.56, 0.18, 0.036);   // lower left is empty: the spectrum falls
      leg->AddEntry(d8, "unfolded, #it{N} = 8 (stat.)", "lp");
      TBox *lb = new TBox(); lb->SetFillColorAlpha(colSys, 0.30); lb->SetLineColor(colSys);
      leg->AddEntry(lb, "total systematic", "f");
      leg->Draw();
      TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.046);
      la.DrawLatex(0.22, 0.84, "pp 5.02 TeV, calo jets: inclusive jets");
      la.SetTextSize(0.036);
      la.DrawLatex(0.22, 0.78, "anti-#it{k}_{T} #it{R} = 0.4, |#eta| < 1.6");
      la.DrawLatex(0.22, 0.73, "Jet100-sample-equivalent yield");

      pBot->cd();
      TH1D *rS = (TH1D*) s8->Clone("rSspec"), *rT = (TH1D*) t8->Clone("rTspec");
      for(TH1D *h : {rS, rT}){ h->SetDirectory(nullptr); h->Scale(100.); for(int i = 0; i <= h->GetNbinsX()+1; i++) h->SetBinError(i, 0.); }
      styleLine(rS, "#D55E00", 3); styleLine(rT, "#0072B2", 3); rT->SetLineStyle(2);
      double m = 0.;
      for(int i = 1; i <= rS->GetNbinsX(); i++) if(shown(i)) m = TMath::Max(m, TMath::Max(rS->GetBinContent(i), rT->GetBinContent(i)));
      TH1D *rF = (TH1D*) rS->Clone("rFspec"); rF->Reset(); rF->SetTitle("");
      rF->SetMinimum(0.); rF->SetMaximum(1.7*m);   // headroom for the legend row
      rF->GetXaxis()->SetRangeUser(ptReportMin, inclEdges[NInclE-1]);
      rS->GetXaxis()->SetRangeUser(ptReportMin, inclEdges[NInclE-1]); rT->GetXaxis()->SetRangeUser(ptReportMin, inclEdges[NInclE-1]);
      rF->GetXaxis()->SetTitle("#it{p}_{T}^{jet} [GeV]"); rF->GetYaxis()->SetTitle("rel. unc. [%]");
      rF->GetXaxis()->SetTitleSize(0.105); rF->GetXaxis()->SetTitleOffset(1.25); rF->GetXaxis()->SetLabelSize(0.090);
      rF->GetYaxis()->SetTitleSize(0.085); rF->GetYaxis()->SetTitleOffset(0.85); rF->GetYaxis()->SetLabelSize(0.085);
      rF->GetYaxis()->SetNdivisions(505);
      rF->Draw("AXIS");
      rS->Draw("HIST same"); rT->Draw("HIST same");
      TLegend *legB = makeLegend(0.20, 0.78, 0.95, 0.97, 0.075);
      legB->SetNColumns(2);
      legB->AddEntry(rS, "total syst.", "l"); legB->AddEntry(rT, "stat.", "l");
      legB->Draw();
      cv->SaveAs(Form("%s/inclJetSpectrum_caloJets_pp.pdf", figDir));
      delete cv;
    } else printf("WARNING: spectrum figure skipped (missing histograms or NZ in %s)\n", measPath);
  }

  fo->cd();
  TNamed info("info", Form("pp calo inclusive jets per Z: spectrum systematics. Measured = jetsPerZ_fine_pp (calculateJetsPerZ.cc, "
                            "no unfolding); response %s from %s, reco floor %.0f GeV; %s iterations / response MC stat "
                            "(%d Poisson toys) / spectrum JEU (shifted PYTHIA reco spectra unfolded with this response), "
                            "quadrature. Binnings: incl = calculateJetsPerZ newAxis, bjet = b-jet chain ptEdges.",
                            respName, gSystem->BaseName(respPath), recoFloorIncl, infoTxt.Data(), nToysIncl));
  info.Write();
  fo->Close();
  printf("\nwritten %s\nfigures in %s\n", outRoot, figDir);
}
