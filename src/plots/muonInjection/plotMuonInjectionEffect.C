// Effect of clustering a mixed-event muon INTO the donor event, measured
// against the same muon matched to jets clustered without it.
//
// Reads the histograms added by the muon-injection study in
// PbPb_pfCandAnalyzer.C (commit e278b25c):
//
//   noInject       donor PF candidates clustered alone, muon matched afterwards
//                  -- what the existing mixed-event templates do
//   inject         muon added to the donor candidates and everything
//                  reclustered, so it contributes pT and pulls the axis
//   injectConstit  as inject, but tagged by the jet anti-kT actually assigned
//                  the muon to, rather than the nearest jet within dR < 0.4
//
// Three observables, each with an inject/noInject ratio panel:
//   jet pT   inclusive donor jet spectrum. Filled per injection in both
//            collections, so they share a denominator and the ratio is the
//            migration caused by adding the muon's pT.
//   dR       muon to jet axis. Should shrink under injection: the axis is
//            pulled toward the muon.
//   ptRel    muon transverse to that axis. Should shrink for the same reason,
//            and this is the quantity the b-purity fit uses.
//
// NORMALISATION: per event, dividing by h_vz. The injection fills are NOT
// trigger-gated (unlike h_mixedMuonPtRel_recoJetPt), so h_vz is the matching
// denominator here, not h_vz_triggerOn.
//
// Note the noInject collection here is not identical to
// h_mixedMuonPtRel_recoJetPt: that one matches to real akCs4PF reco jets,
// while this matches to fastJets reclustered from the donor's PF candidates.
// The fastJet system is used because real reco jets cannot be reclustered with
// an added muon. Reading the inject/noInject ratio as the correction to apply
// to the reco-jet template assumes the two jet definitions respond alike.
//
// Usage: root -l -b -q 'plotMuonInjectionEffect.C("<scan file>")'
// Run from: src/plots/muonInjection/

#include "../../../headers/functions/divideByBinwidth.h"

const char *inFileDefault = "";   // no default yet; pass the scan file explicitly

const char *outDir = "../../../figures/muonInjection/";

const int NClass = 4;
const int sliceLo[NClass] = { 1,  3,  7, 11};
const int sliceHi[NClass] = { 2,  6, 10, 16};
const char *classLabel[NClass] = { "0-10%", "10-30%", "30-50%", "50-80%" };
const char *classTag[NClass]   = { "0to10pct", "10to30pct", "30to50pct", "50to80pct" };

const int    NPt = 3;
const double ptLo[NPt] = {50, 60,  80};
const double ptHi[NPt] = {60, 80, 120};

// Okabe-Ito. noInject is the reference, inject the thing being measured.
const char *hexNo     = "#0072B2";   // blue
const char *hexInj    = "#D55E00";   // vermillion
const char *hexConst  = "#009E73";   // green
const int   markNo    = 20;          // circle
const int   markInj   = 21;          // square
const int   markConst = 33;          // diamond   -- no triangles

// binning, matching the companion plots so figures can be read together
const int    NEdgeDR = 11;
double       edgeDR[NEdgeDR] = {0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50};
const double edgePtRel[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,
                            2.0,2.4,2.8,3.4,4.2,5.0};
const int    NEdgePtRel = (int)(sizeof(edgePtRel)/sizeof(double));
const int    rebinJetPt = 4;

double classEvents(TFile *f, int ci){
  double n = 0.;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *h = nullptr;
    f->GetObject(Form("h_vz_C%d", si), h);
    if(!h) return -1.;
    n += h->Integral();
  }
  return n;
}

// sum a TH2 over a coarse class
TH2D* sum2(TFile *f, const char *base, int ci, const char *tag){
  TH2D *s = nullptr;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH2D *h = nullptr;
    f->GetObject(Form("%s_C%d", base, si), h);
    if(!h){ if(s) delete s; return nullptr; }
    if(!s){ s = (TH2D*) h->Clone(Form("s2_%s_%s_%d", base, tag, ci)); s->SetDirectory(nullptr); }
    else s->Add(h);
  }
  return s;
}
TH1D* sum1(TFile *f, const char *base, int ci, const char *tag){
  TH1D *s = nullptr;
  for(int si = sliceLo[ci]; si <= sliceHi[ci]; si++){
    TH1D *h = nullptr;
    f->GetObject(Form("%s_C%d", base, si), h);
    if(!h){ if(s) delete s; return nullptr; }
    if(!s){ s = (TH1D*) h->Clone(Form("s1_%s_%s_%d", base, tag, ci)); s->SetDirectory(nullptr); }
    else s->Add(h);
  }
  return s;
}

void styleH(TH1 *h, const char *hex, int mark){
  int c = TColor::GetColor(hex);
  h->SetLineColor(c); h->SetMarkerColor(c);
  h->SetMarkerStyle(mark); h->SetMarkerSize(1.0); h->SetLineWidth(2);
  h->SetStats(0);
}

// ratio with the denominator treated as exact-enough: both come from the same
// injections, so they are strongly correlated and independent-error propagation
// would overstate the uncertainty. Numerator error only.
TH1D* ratioOf(TH1D *num, TH1D *den, const char *name){
  TH1D *r = (TH1D*) num->Clone(name); r->SetDirectory(nullptr);
  for(int b = 1; b <= r->GetNbinsX(); b++){
    double d = den->GetBinContent(b);
    if(d == 0.){ r->SetBinContent(b, -999); r->SetBinError(b, 0); continue; }
    r->SetBinContent(b, num->GetBinContent(b)/d);
    r->SetBinError  (b, num->GetBinError(b)/d);
  }
  return r;
}

// one canvas: overlay + inject/noInject ratio panel
void drawPair(TH1D *hNo, TH1D *hInj, TH1D *hConst,
              const char *xTitle, const char *yTitle,
              const char *headLine, const char *subLine,
              double xMin, double xMax, const char *outName)
{
  TCanvas *c = new TCanvas(Form("c_%s", outName), "", 700, 800);
  TPad *pT = new TPad("pT","",0,0.34,1,1);
  TPad *pB = new TPad("pB","",0,0,1,0.34);
  pT->SetBottomMargin(0.02); pT->SetLeftMargin(0.17); pT->SetTopMargin(0.07);
  pB->SetTopMargin(0.02);    pB->SetLeftMargin(0.17); pB->SetBottomMargin(0.32);
  pT->Draw(); pB->Draw();

  pT->cd();
  styleH(hNo, hexNo, markNo); styleH(hInj, hexInj, markInj);
  if(hConst) styleH(hConst, hexConst, markConst);

  double ymax = TMath::Max(hNo->GetMaximum(), hInj->GetMaximum());
  if(hConst) ymax = TMath::Max(ymax, hConst->GetMaximum());
  if(ymax <= 0.){ printf("  %s: empty, skipped\n", outName); delete c; return; }

  hNo->GetXaxis()->SetRangeUser(xMin, xMax);
  hNo->GetXaxis()->SetLabelSize(0);
  hNo->GetYaxis()->SetTitle(yTitle);
  hNo->GetYaxis()->SetTitleSize(0.055); hNo->GetYaxis()->SetTitleOffset(1.45);
  hNo->GetYaxis()->SetLabelSize(0.045);
  hNo->SetTitle(""); hNo->SetMinimum(0.); hNo->SetMaximum(ymax*1.75);
  hNo->Draw("E");
  hInj->Draw("E same");
  if(hConst) hConst->Draw("E same");

  // sits below the three header lines, which run down to 0.785
  TLegend *leg = new TLegend(0.40,0.52,0.95,0.75);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.038);
  leg->AddEntry(hNo,  Form("no inject  (#mu=%.3f)", hNo->GetMean()), "lp");
  leg->AddEntry(hInj, Form("inject  (#mu=%.3f)",    hInj->GetMean()), "lp");
  if(hConst) leg->AddEntry(hConst, Form("inject, constituent  (#mu=%.3f)", hConst->GetMean()), "lp");
  leg->Draw();

  TLatex la; la.SetNDC(); la.SetTextFont(42); la.SetTextSize(0.042);
  la.DrawLatex(0.21, 0.885, headLine);
  la.DrawLatex(0.21, 0.835, subLine);
  // the headline number: how far the mean moved
  la.SetTextFont(62);
  la.DrawLatex(0.21, 0.785, Form("#LT%s#GT: %.3f #rightarrow %.3f  (%+.1f%%)",
                                 xTitle, hNo->GetMean(), hInj->GetMean(),
                                 hNo->GetMean() != 0. ?
                                 100.*(hInj->GetMean()/hNo->GetMean() - 1.) : 0.));

  pB->cd();
  TH1D *rI = ratioOf(hInj, hNo, Form("rI_%s", outName));
  TH1D *rC = hConst ? ratioOf(hConst, hNo, Form("rC_%s", outName)) : nullptr;
  styleH(rI, hexInj, markInj);
  if(rC) styleH(rC, hexConst, markConst);

  // Range from bins where the DENOMINATOR is genuinely populated (>2% of its
  // own peak). An earlier version capped the scan at 6, which was tuned for
  // ratios near unity -- but injection multiplies the tagged yield by ~9, so
  // that ceiling pushed almost every real point off-scale and left the panel
  // showing one bin. The blow-ups worth excluding come from an empty noInject
  // denominator, not from a large honest ratio, so gate on the denominator
  // instead of on the ratio value.
  double rmax = 1.2, noPeak = hNo->GetMaximum();
  for(int b = 1; b <= rI->GetNbinsX(); b++){
    double x = rI->GetBinCenter(b);
    if(x < xMin || x > xMax) continue;
    if(hNo->GetBinContent(b) < 0.02*noPeak) continue;
    double v = rI->GetBinContent(b) + rI->GetBinError(b);
    if(v > rmax) rmax = v;
  }
  rI->GetXaxis()->SetRangeUser(xMin, xMax);
  rI->SetTitle("");
  rI->GetXaxis()->SetTitle(xTitle);
  rI->GetXaxis()->SetTitleSize(0.105); rI->GetXaxis()->SetTitleOffset(1.25);
  rI->GetXaxis()->SetLabelSize(0.090);
  rI->GetYaxis()->SetTitle("inject / no inject");
  rI->GetYaxis()->SetTitleSize(0.085); rI->GetYaxis()->SetTitleOffset(0.80);
  rI->GetYaxis()->SetLabelSize(0.085); rI->GetYaxis()->SetNdivisions(505);
  rI->SetMinimum(0.); rI->SetMaximum(rmax*1.25);
  rI->Draw("E");
  TLine *l1 = new TLine(xMin, 1.0, xMax, 1.0);
  l1->SetLineStyle(2); l1->SetLineColor(kGray+2); l1->Draw();
  rI->Draw("E same");
  if(rC) rC->Draw("E same");

  c->SaveAs(Form("%s%s.pdf", outDir, outName));
  delete c;
}

void plotMuonInjectionEffect(const char *scanFile = nullptr)
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  gSystem->mkdir(outDir, kTRUE);

  const char *usePath = (scanFile && scanFile[0]) ? scanFile : inFileDefault;
  if(!usePath || !usePath[0]){
    printf("ERROR: pass the scan file explicitly, e.g.\n"
           "  root -l -b -q 'plotMuonInjectionEffect.C(\"/path/to/scan.root\")'\n");
    return;
  }
  TFile *f = TFile::Open(gSystem->ExpandPathName(usePath));
  if(!f || f->IsZombie()){ printf("ERROR: cannot open %s\n", usePath); return; }
  printf("scan file: %s\n\n", usePath);

  if(!f->GetListOfKeys()->FindObject("h_donorJetPt_noInject_C1")){
    printf("ERROR: this scan has no muon-injection histograms.\n"
           "       They were added to PbPb_pfCandAnalyzer.C on 2026-09-08 (e278b25c)\n"
           "       and require a rescan to populate.\n");
    return;
  }

  // Entry counts matter as much as the means here: the injected muon makes its
  // own jet, so it is tagged far more often than when matched to jets built
  // without it. That change in population IS the under-representation this
  // study is about, and a mean alone would hide it.
  printf("%-8s %-8s %-6s %9s %9s %9s %9s %9s   %s\n",
         "class","jetpT","obs","N(noInj)","N(inj)","<noInj>","<inject>","<constit>","shift");

  for(int ci = 0; ci < NClass; ci++){

    double nEvt = classEvents(f, ci);
    if(nEvt <= 0.){ printf("  %s: no h_vz, skipped\n", classLabel[ci]); continue; }

    // ---------------- inclusive donor jet pT ----------------
    TH1D *jNo  = sum1(f, "h_donorJetPt_noInject", ci, "jNo");
    TH1D *jInj = sum1(f, "h_donorJetPt_inject",   ci, "jInj");
    if(jNo && jInj){
      jNo->Rebin(rebinJetPt); jInj->Rebin(rebinJetPt);
      jNo->Scale(1./nEvt);  divideByBinwidth(jNo);
      jInj->Scale(1./nEvt); divideByBinwidth(jInj);
      printf("%-8s %-8s %-6s %9.0f %9.0f %9.3f %9.3f %9s   %s\n",
             classLabel[ci], "incl", "jetpT",
             jNo->GetEntries(), jInj->GetEntries(),
             jNo->GetMean(), jInj->GetMean(), "-",
             jNo->GetMean() > 0. ? Form("%+.1f%%", 100.*(jInj->GetMean()/jNo->GetMean()-1.)) : "-");
      drawPair(jNo, jInj, nullptr,
               "#it{p}_{T}^{jet} [GeV]", "d#it{N}/d#it{p}_{T}^{jet} per event",
               Form("PbPb SingleMuon, %s", classLabel[ci]),
               "donor fastJets, inclusive",
               20., 150., Form("muonInjection_jetPt_%s", classTag[ci]));
    }
    else printf("  %s: donor jet pT histograms missing\n", classLabel[ci]);

    // ---------------- dR and ptRel, per jet pT window ----------------
    struct Obs { const char *base; const char *tag; const char *xT;
                 const double *edges; int nEdge; double xMax; };
    TH2D *dNo = sum2(f, "h_injMuonDR_donorJetPt_noInject",       ci, "dNo");
    TH2D *dIn = sum2(f, "h_injMuonDR_donorJetPt_inject",         ci, "dIn");
    TH2D *dCo = sum2(f, "h_injMuonDR_donorJetPt_injectConstit",  ci, "dCo");
    TH2D *pNo = sum2(f, "h_injMuonPtRel_donorJetPt_noInject",      ci, "pNo");
    TH2D *pIn = sum2(f, "h_injMuonPtRel_donorJetPt_inject",        ci, "pIn");
    TH2D *pCo = sum2(f, "h_injMuonPtRel_donorJetPt_injectConstit", ci, "pCo");

    for(int p = 0; p < NPt; p++){
      for(int which = 0; which < 2; which++){
        TH2D *H[3];
        if(which == 0){ H[0]=dNo; H[1]=dIn; H[2]=dCo; }
        else          { H[0]=pNo; H[1]=pIn; H[2]=pCo; }
        if(!H[0] || !H[1] || !H[2]) continue;

        TH1D *h[3];
        for(int k = 0; k < 3; k++){
          int b1 = H[k]->GetYaxis()->FindBin(ptLo[p] + 1e-6);
          int b2 = H[k]->GetYaxis()->FindBin(ptHi[p] - 1e-6);
          TH1D *raw = H[k]->ProjectionX(Form("px_%d_%d_%d_%d", ci, p, which, k), b1, b2);
          raw->SetDirectory(nullptr);
          h[k] = (which == 0)
               ? (TH1D*) raw->Rebin(NEdgeDR-1,    Form("rb_%d_%d_%d_%d", ci,p,which,k), edgeDR)
               : (TH1D*) raw->Rebin(NEdgePtRel-1, Form("rb_%d_%d_%d_%d", ci,p,which,k), edgePtRel);
          h[k]->SetDirectory(nullptr);
          h[k]->Scale(1./nEvt);
          divideByBinwidth(h[k]);
          delete raw;
        }

        printf("%-8s %3.0f-%-4.0f %-6s %9.0f %9.0f %9.3f %9.3f %9.3f   %s\n",
               classLabel[ci], ptLo[p], ptHi[p], which == 0 ? "dR" : "ptRel",
               h[0]->GetEntries(), h[1]->GetEntries(),
               h[0]->GetMean(), h[1]->GetMean(), h[2]->GetMean(),
               h[0]->GetMean() > 0. ? Form("%+.1f%%", 100.*(h[1]->GetMean()/h[0]->GetMean()-1.)) : "-");

        drawPair(h[0], h[1], h[2],
                 which == 0 ? "#it{#Delta}#it{R}(#mu,jet)" : "#it{p}_{T}^{rel} [GeV]",
                 which == 0 ? "d#it{N}/d(#it{#Delta}#it{R}) per event"
                            : "d#it{N}/d#it{p}_{T}^{rel} per event",
                 Form("PbPb SingleMuon, %s", classLabel[ci]),
                 Form("%.0f < #it{p}_{T}^{jet} < %.0f GeV", ptLo[p], ptHi[p]),
                 0., which == 0 ? 0.5 : 5.0,
                 Form("muonInjection_%s_%s_pt%.0fto%.0f",
                      which == 0 ? "dR" : "ptRel", classTag[ci], ptLo[p], ptHi[p]));
      }
    }
  }

  printf("\nfigures written to %s\n", outDir);
}
