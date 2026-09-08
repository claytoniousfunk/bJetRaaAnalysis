// Task #5 follow-up: is the jetTrkMax filter's rejection rate
// centrality-dependent, and does it explain the 0-10% ~ 10-30% degeneracy?
//
// passesJetTrkMaxFilter (headers/functions/jet_filter/) keeps a jet only if
//   0.01 < jetTrkMax/jetPt < 0.98
// and is applied as a `continue` BEFORE any histogram is filled, with no
// trkMax diagnostic histograms written out.  Rejected jets therefore leave no
// trace in a scan's output, so the rejection rate cannot be recovered from a
// single file -- it needs a matched pair of scans differing only in the flag.
//
// The filename encodes the flag reliably: configureOutputDatasetName_PbPb.h
// appends "_jetTrkMaxFilter" only when doJetTrkMaxFilter is true.
//
// Matched pair used here (identical selections, differ only in the flag):
//   ON  : PbPb_HardProbes_<trig>_..._tight_jetTrkMaxFilter_WDecayFilter_2026-5-5
//   OFF : PbPb_HardProbes_<trig>_..._tight_WDecayFilter_2026-4-8
// Both are native 4-class files, both pre-2026-08-12, so the centrality
// convention matches on both sides.
//
// *** WHY COARSE CLASSES AND NOT 5% SLICES ***
// No filter-OFF scan exists with the current 16-slice ultra-fine scheme.  The
// one ultra-fine MinBias file with the filter ON (2026-5-7) predates commit
// c19858f2 ("extend number of ultra-cent bins out to 160", 2026-05-12) and
// carries only 8 slices on a different set of edges, so it cannot be merged
// against the 4-class filter-OFF file.  Coarse classes are therefore the
// finest granularity currently measurable; a fine-slice version needs a
// dedicated filter-OFF rescan.
//
// HardProbes rather than MinBias because that is the only matched pair
// available.  The jet population is trigger-biased toward real jets, so this
// likely UNDERSTATES the rejection rate for the combinatorial-rich MinBias
// sample that dominates the low-pT region of the stitch.
//
// Usage: root -l -b -q 'plotJetTrkMaxFilter_passFraction_vsCent.C'
// Run from: src/plots/jetPt/

const char *dirM =
  "/home/clayton/Analysis/code/bJetMuonTaggingAnalysis/rootFiles/scanningOutput/PbPb/latest/";

const char *outDir  = "../../../figures/jetPt/";
const char *outName = "jetTrkMaxFilter_passFraction_vsCent.pdf";

const int   NTrig = 2;
const char *trig[NTrig]      = { "Jet80HLT", "Jet100HLT" };
const char *trigLabel[NTrig] = { "HardProbes Jet80", "HardProbes Jet100" };
const int   trigMark[NTrig]  = { 20, 24 };

const int   NClass = 4;
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };
const char *classHex[NClass]   = { "#D55E00", "#0072B2", "#009E73", "#CC79A7" };

const double plotPtMin = 60., plotPtMax = 260.;

const int    NEdge = 11;
double       ptEdge[NEdge] = { 60,70,80,90,100,120,140,170,200,230,260 };

static TH1D* perEvent(TFile *f, int c, const char *name)
{
  TH1D *h = nullptr, *v = nullptr;
  f->GetObject(Form("h_inclRecoJetPt_C%d", c), h);
  f->GetObject(Form("h_vz_C%d", c), v);
  if(!h || !v) return nullptr;
  TH1D *r = (TH1D*) h->Rebin(NEdge-1, name, ptEdge);
  r->SetDirectory(nullptr);
  r->Scale(1./v->Integral());
  return r;
}

void plotJetTrkMaxFilter_passFraction_vsCent()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  TH1D *hPass[NTrig][NClass];
  int col[NClass];
  for(int c = 0; c < NClass; c++) col[c] = TColor::GetColor(classHex[c]);

  for(int t = 0; t < NTrig; t++){
    TFile *fon  = TFile::Open(Form("%sPbPb_HardProbes_%s_mu12_pTmu-15to999_tight_jetTrkMaxFilter_WDecayFilter_2026-5-5.root", dirM, trig[t]));
    TFile *foff = TFile::Open(Form("%sPbPb_HardProbes_%s_mu12_pTmu-15to999_tight_WDecayFilter_2026-4-8.root", dirM, trig[t]));
    if(!fon || fon->IsZombie() || !foff || foff->IsZombie()){
      printf("ERROR: cannot open matched pair for %s\n", trig[t]); return;
    }

    printf("\n=== %s : trkMax-filter pass fraction (per event, ON/OFF) ===\n", trigLabel[t]);
    printf("%-9s", "pT [GeV]");
    for(int c = 0; c < NClass; c++) printf(" %10s", classLabel[c]);
    printf("\n");

    for(int c = 0; c < NClass; c++){
      TH1D *on  = perEvent(fon,  c+1, Form("on_%d_%d",  t, c));
      TH1D *off = perEvent(foff, c+1, Form("off_%d_%d", t, c));
      if(!on || !off){ printf("ERROR: missing class %d\n", c+1); return; }
      hPass[t][c] = (TH1D*) on->Clone(Form("pass_%d_%d", t, c));
      hPass[t][c]->Divide(off);
    }
    for(int b = 1; b <= hPass[t][0]->GetNbinsX(); b++){
      double lo = hPass[t][0]->GetXaxis()->GetBinLowEdge(b);
      if(lo < plotPtMin || lo >= plotPtMax) continue;
      printf("%4.0f-%-4.0f", lo, hPass[t][0]->GetXaxis()->GetBinUpEdge(b));
      for(int c = 0; c < NClass; c++) printf(" %10.4f", hPass[t][c]->GetBinContent(b));
      printf("\n");
    }
  }

  TCanvas *c = new TCanvas("cTrkMaxPass", "", 900, 700);
  c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
  c->SetTopMargin(0.09);  c->SetBottomMargin(0.14);

  TH1F *fr = c->DrawFrame(plotPtMin, 0.84, plotPtMax, 1.08);
  fr->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  fr->GetYaxis()->SetTitle("fraction of jets passing trkMax filter");
  fr->GetXaxis()->SetTitleSize(0.048); fr->GetXaxis()->SetLabelSize(0.042);
  fr->GetYaxis()->SetTitleSize(0.048); fr->GetYaxis()->SetLabelSize(0.042);
  fr->GetYaxis()->SetTitleOffset(1.35);

  TLine *one = new TLine(plotPtMin, 1., plotPtMax, 1.);
  one->SetLineStyle(2); one->SetLineColor(kGray+1); one->Draw();

  for(int t = 0; t < NTrig; t++){
    for(int cc = 0; cc < NClass; cc++){
      hPass[t][cc]->SetLineColor(col[cc]); hPass[t][cc]->SetMarkerColor(col[cc]);
      hPass[t][cc]->SetMarkerStyle(trigMark[t]); hPass[t][cc]->SetMarkerSize(1.1);
      hPass[t][cc]->SetLineWidth(2);
      if(t == 1) hPass[t][cc]->SetLineStyle(2);
      hPass[t][cc]->GetXaxis()->SetRangeUser(plotPtMin, plotPtMax);
      hPass[t][cc]->Draw("ep same");
    }
  }

  TLegend *legC = new TLegend(0.47, 0.17, 0.68, 0.41);
  legC->SetBorderSize(0); legC->SetFillStyle(0); legC->SetTextSize(0.038);
  for(int cc = 0; cc < NClass; cc++) legC->AddEntry(hPass[0][cc], classLabel[cc], "lp");
  legC->Draw();

  TLegend *legT = new TLegend(0.66, 0.17, 0.95, 0.31);
  legT->SetBorderSize(0); legT->SetFillStyle(0); legT->SetTextSize(0.036);
  TH1D *keyT[NTrig];
  for(int t = 0; t < NTrig; t++){
    keyT[t] = new TH1D(Form("keyT_%d", t), "", 1, 0, 1);
    keyT[t]->SetDirectory(nullptr);
    keyT[t]->SetMarkerStyle(trigMark[t]); keyT[t]->SetMarkerColor(kBlack);
    keyT[t]->SetLineColor(kBlack); keyT[t]->SetMarkerSize(1.1);
    if(t == 1) keyT[t]->SetLineStyle(2);
    legT->AddEntry(keyT[t], trigLabel[t], "lp");
  }
  legT->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038);
  lat.DrawLatex(0.14, 0.945, "jetTrkMax filter: PbPb rejection rate is centrality-dependent");

  TString out = TString(outDir) + outName;
  c->SaveAs(out);
  printf("\nSaved %s\n", out.Data());
}
