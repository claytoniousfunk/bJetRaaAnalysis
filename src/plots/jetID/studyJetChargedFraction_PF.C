// Charged / photon / neutral-hadron energy fractions of PF jets, PbPb data vs PYTHIA+HYDJET
// and pp data vs PYTHIA, straight from the forests (2026-11-19, user: "run the charged energy
// fraction check"). Motivation: the inclusive R_AA is low in every class against ATLAS and
// CMS, at reco level already, with the configs, triggers, N_Z and unfolding ruled out. A PbPb
// DATA calorimeter (photon / neutral-hadron) scale that the simulation does not have would
// raise the data's charged fraction relative to P+H in PbPb but not in pp; a 7% jet energy
// scale offset carried only by the ~40% neutral part moves the charged fraction by ~+0.04.
//
// Jets: akCs4PF (PbPb) / ak4PF (pp), |eta| < 1.6, pT = the scans' manual JEC on rawpt
// (Autumn18_HI_V8 / Spring18_ppRef5TeV_V6, data with L2L3Residual), the jetTrkMax filter
// (0.01 < trackMax / pT < 0.98). Events: pprimaryVertexFilter, pclusterCompatibilityFilter
// (PbPb, as the MinBias / HardProbes scans), hiBin < 180, MC hiBin - 10 (config_PYTHIAHYDJET.h).
// MC: genuine jets only (refpt > 0), forest generator weight. No trigger requirement:
// a jet-triggered data sample (HardProbes Jet80 / Jet100, HighEGJet) is fine for a per-jet
// fraction.
// Fractions per jet: X / S with S = chargedSum + photonSum + neutralSum + eSum + muSum
// (the HiJetID sums of the jet's constituents, CS-subtracted in PbPb); S / rawpt is kept as a
// sanity check. Jets with a muon or electron fraction above 0.1 are dropped in every sample (the
// local pp data forest is SingleMuon, whose jets often hold the trigger muon).
//
// NOTE (2026-11-19, local forests): in PbPb the HiJetID sums are NOT background-subtracted:
// S / rawpt (p_sumRaw) is 1.03-1.11 in 50-80% (P+H, 40-350 GeV) but 1.3-2.2 in 10-30%, so in
// central classes the fractions mix the jet with the underlying event and test HYDJET's UE
// composition as much as the jet energy scale. 50-80% is the clean PbPb test; S / rawpt - 1
// (the in-cone UE) is itself a data / MC comparison. pp data (SingleMuon, lepton veto) and
// PYTHIA agree: charged fraction data - MC within -0.02 / +0.02 at 50-200 GeV.
//
// Usage (from src/plots/jetID/):
//   fill:  root -l -b -q 'studyJetChargedFraction_PF.C+("fill","PbPb","<file.root | list.txt>",<maxFiles>,"<outTag>")'
//          systems "PbPb", "PH" (PYTHIA+HYDJET), "pp", "PYTHIA"
//          -> rootFiles/jetID/chargedFraction_<system><outTag>.root
//   plot:  root -l -b -q 'studyJetChargedFraction_PF.C+("plot","<PbPb tag>","<PH tag>",0,"<pp tag>|<PYTHIA tag>")'
//          -> figures/jetID/chargedFraction_dataVsMC.pdf
// The PbPb data forests are on EOS only (run on lxplus), e.g. the first 300 files of
// fileNames/fileNames_HIHardProbes_withCaloAndFlowJets_fresh_partial.txt.
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include <fstream>
#include <vector>
#include <string>
#include "../../../JetEnergyCorrections/JetCorrector.h"
#include "../../../headers/plotting/plotStyle.h"

namespace jcf {
const char *repo = "../../..";   // run from src/plots/jetID/ (works on lxplus too)
const int NB = 10;
const double edges[NB+1] = {40, 50, 60, 70, 80, 100, 120, 150, 200, 250, 350};
const int NC = 5;                                       // 0 = pp, 1-4 = 0-10, 10-30, 30-50, 50-80%
const char *clsLab[NC] = {"pp", "0-10%", "10-30%", "30-50%", "50-80%"};
const int hiLo[NC] = {-1, 0, 20, 60, 100}, hiHi[NC] = {-1, 20, 60, 100, 160};
const char *frac[4] = {"chf", "phf", "nhf", "sumRaw"};
const char *fracLab[4] = {"charged-hadron fraction", "photon fraction", "neutral-hadron fraction", "#Sigma constituents / raw p_{T}"};
const int NJ = 500;

void addFiles(TChain &c, const std::string &in, int maxFiles){
  if(in.size() > 5 && in.substr(in.size() - 5) == ".root"){ c.Add(in.c_str()); return; }
  std::ifstream f(in); std::string l; int n = 0;
  while(std::getline(f, l) && (maxFiles <= 0 || n < maxFiles)){
    if(l.empty() || l[0] == '#') continue;
    if(l.rfind("/eos/cms", 0) == 0) l = "root://eoscms.cern.ch/" + l;
    c.Add(l.c_str()); n++;
  }
}

void fill(TString sys, std::string in, int maxFiles, TString outTag){
  const bool isPbPb = (sys == "PbPb" || sys == "PH"), isMC = (sys == "PH" || sys == "PYTHIA");
  const char *jt = isPbPb ? "akCs4PFJetAnalyzer/t" : "ak4PFJetAnalyzer/t";
  TChain cj(jt), ce("hiEvtAnalyzer/HiTree"), cs("skimanalysis/HltTree");
  addFiles(cj, in, maxFiles); addFiles(ce, in, maxFiles); addFiles(cs, in, maxFiles);
  cj.AddFriend(&ce); cj.AddFriend(&cs);
  printf("%s: %d files, %lld events\n", sys.Data(), cj.GetNtrees(), cj.GetEntries());

  std::vector<std::string> F;
  const std::string jd = std::string(repo) + "/JetEnergyCorrections/";
  if(isPbPb){ F.push_back(jd + "Autumn18_HI_V8_" + (isMC ? "MC" : "DATA") + "_L2Relative_AK4PF.txt");
              if(!isMC) F.push_back(jd + "Autumn18_HI_V8_DATA_L2L3Residual_AK4PF.txt"); }
  else      { F.push_back(jd + (isMC ? "Spring18_ppRef5TeV_V6_MC_L2Relative_AK4PF.txt" : "Spring18_ppRef5TeV_V6_DATA_L2Relative_AK4PF.txt"));
              if(!isMC) F.push_back(jd + "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt"); }
  for(auto &s : F) if(gSystem->AccessPathName(s.c_str())){ printf("ERROR: JEC file %s missing\n", s.c_str()); return; }
  JetCorrector JEC(F);

  int nref, hiBin = 0, pv = 1, cc = 1; float weight = 1.;
  float rawpt[NJ], jteta[NJ], jtphi[NJ], refpt[NJ], trackMax[NJ], chS[NJ], phS[NJ], neS[NJ], eS[NJ], muS[NJ];
  cj.SetBranchStatus("*", 0);
  for(const char *b : {"nref", "rawpt", "jteta", "jtphi", "trackMax", "chargedSum", "photonSum", "neutralSum", "eSum", "muSum"}) cj.SetBranchStatus(b, 1);
  cj.SetBranchAddress("nref", &nref); cj.SetBranchAddress("rawpt", rawpt); cj.SetBranchAddress("jteta", jteta); cj.SetBranchAddress("jtphi", jtphi);
  cj.SetBranchAddress("trackMax", trackMax); cj.SetBranchAddress("chargedSum", chS); cj.SetBranchAddress("photonSum", phS);
  cj.SetBranchAddress("neutralSum", neS); cj.SetBranchAddress("eSum", eS); cj.SetBranchAddress("muSum", muS);
  if(isMC){ cj.SetBranchStatus("refpt", 1); cj.SetBranchAddress("refpt", refpt); cj.SetBranchStatus("weight", 1); cj.SetBranchAddress("weight", &weight); }
  if(isPbPb){ cj.SetBranchStatus("hiBin", 1); cj.SetBranchAddress("hiBin", &hiBin);
              cj.SetBranchStatus("pprimaryVertexFilter", 1); cj.SetBranchAddress("pprimaryVertexFilter", &pv);
              cj.SetBranchStatus("pclusterCompatibilityFilter", 1); cj.SetBranchAddress("pclusterCompatibilityFilter", &cc); }

  TProfile *p[NC][4]; TH1D *hn[NC];
  for(int c = 0; c < NC; c++){
    for(int k = 0; k < 4; k++) p[c][k] = new TProfile(Form("p_%s_C%d", frac[k], c), Form(";jet p_{T} (GeV);%s", fracLab[k]), NB, edges);
    hn[c] = new TH1D(Form("h_nJets_C%d", c), ";jet p_{T} (GeV);jets (unweighted)", NB, edges);
  }
  const Long64_t N = cj.GetEntries();
  for(Long64_t i = 0; i < N; i++){
    if(cj.GetEntry(i) <= 0) continue;
    if(i % 200000 == 0) printf("  %lld / %lld\n", i, N);
    int c = 0;
    if(isPbPb){
      if(!pv || !cc) continue;
      const int hb = isMC ? hiBin - 10 : hiBin;
      if(hb < 0 || hb >= 160) continue;
      c = -1; for(int k = 1; k < NC; k++) if(hb >= hiLo[k] && hb < hiHi[k]) c = k;
      if(c < 0) continue;
    }
    const double w = isMC ? weight : 1.;
    for(int j = 0; j < nref && j < NJ; j++){
      if(fabs(jteta[j]) >= 1.6) continue;
      if(isMC && refpt[j] <= 0.) continue;
      JEC.SetJetPT(rawpt[j]); JEC.SetJetEta(jteta[j]); JEC.SetJetPhi(jtphi[j]);
      const double pt = JEC.GetCorrectedPT();
      if(pt < edges[0] || pt >= edges[NB]) continue;
      const double tf = trackMax[j]/pt; if(tf <= 0.01 || tf >= 0.98) continue;
      const double S = chS[j] + phS[j] + neS[j] + eS[j] + muS[j];
      if(S <= 0.) continue;
      if(muS[j]/S > 0.1 || eS[j]/S > 0.1) continue;   // lepton veto, every sample alike (muon-triggered data)
      p[c][0]->Fill(pt, chS[j]/S, w); p[c][1]->Fill(pt, phS[j]/S, w); p[c][2]->Fill(pt, neS[j]/S, w); p[c][3]->Fill(pt, S/rawpt[j], w);
      hn[c]->Fill(pt);
    }
  }
  gSystem->mkdir(Form("%s/rootFiles/jetID", repo), true);
  TFile *fo = TFile::Open(Form("%s/rootFiles/jetID/chargedFraction_%s%s.root", repo, sys.Data(), outTag.Data()), "recreate");
  for(int c = 0; c < NC; c++){ for(int k = 0; k < 4; k++) p[c][k]->Write(); hn[c]->Write(); }
  fo->Close();
  printf("wrote rootFiles/jetID/chargedFraction_%s%s.root\n", sys.Data(), outTag.Data());
}

TProfile* getP(TFile *f, int k, int c){ return f ? (TProfile*) f->Get(Form("p_%s_C%d", frac[k], c)) : nullptr; }

// tags: PbPb data, PH, pp data | PYTHIA ("<pp>|<PYTHIA>")
void plot(TString tAA, TString tPH, TString tPP){
  TString tpp = tPP(0, tPP.Index("|")), tPY = tPP(tPP.Index("|") + 1, tPP.Length());
  auto op = [&](const char *s, TString t){ return TFile::Open(Form("%s/rootFiles/jetID/chargedFraction_%s%s.root", repo, s, t.Data())); };
  TFile *fAA = op("PbPb", tAA), *fPH = op("PH", tPH), *fpp = op("pp", tpp), *fPY = op("PYTHIA", tPY);
  initPlotStyle(); gStyle->SetOptStat(0);
  // columns pp, 50-80, 30-50, 10-30, 0-10%; rows: chf, data - MC of chf, data - MC of phf, data - MC of nhf
  const int order[NC] = {0, 4, 3, 2, 1};
  TCanvas *cv = new TCanvas("cv", "", 2000, 1300);
  const double L = 0.07, R = 0.01, T = 0.06, B = 0.08, w = (1 - L - R)/NC, h = (1 - T - B)/4;
  printf("\n%-8s %-9s  %-22s %-22s  data - MC: chf     phf     nhf     sum/raw\n", "system", "pT", "chf data", "chf MC");
  for(int col = 0; col < NC; col++){
    const int c = order[col];
    TFile *fd = c == 0 ? fpp : fAA, *fm = c == 0 ? fPY : fPH;
    for(int row = 0; row < 4; row++){
      const double x1 = col == 0 ? 0 : L + col*w, x2 = L + (col+1)*w + (col == NC-1 ? R : 0);
      const double y2 = 1 - T - row*h + (row == 0 ? T : 0), y1 = 1 - T - (row+1)*h - (row == 3 ? B : 0);
      cv->cd(); TPad *pd = new TPad(Form("p%d_%d", row, col), "", x1, y1 < 0 ? 0 : y1, x2, y2);
      pd->SetLeftMargin(col == 0 ? L/(x2 - x1) : 0); pd->SetRightMargin(col == NC-1 ? R/(x2 - x1) : 0);
      pd->SetTopMargin(row == 0 ? T/(y2 - y1) : 0); pd->SetBottomMargin(row == 3 ? B/(y2 - y1) : 0); pd->SetTicks(1, 1); pd->Draw(); pd->cd();
      const bool top = (row == 0);
      TH1F *fr = pd->DrawFrame(40, top ? 0.45 : -0.10, 349, top ? 0.85 : 0.10);
      fr->GetXaxis()->SetNdivisions(505); fr->GetYaxis()->SetNdivisions(505);
      const double ts = 0.075*(h/(y2 - y1));
      fr->GetXaxis()->SetLabelSize(row == 3 ? 1.1*ts : 0); fr->GetXaxis()->SetTitleSize(1.2*ts); fr->GetXaxis()->SetTitle(row == 3 ? "jet p_{T} (GeV)" : "");
      fr->GetYaxis()->SetLabelSize(col == 0 ? ts : 0); fr->GetYaxis()->SetTitleSize(col == 0 ? ts : 0); fr->GetYaxis()->SetTitleOffset(1.0);
      const char *yt[4] = {"charged fraction", "#Delta charged (data - MC)", "#Delta photon (data - MC)", "#Delta neutral had. (data - MC)"};
      fr->GetYaxis()->SetTitle(yt[row]);
      if(!top){ TLine *l = new TLine(40, 0, 349, 0); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw(); }
      if(top){
        TProfile *d = getP(fd, 0, c), *m = getP(fm, 0, c);
        if(m){ m = (TProfile*) m->Clone(); styleH(m, "#D55E00", markOpenSquare, 1.6); m->SetLineWidth(2); m->Draw("E1 X0 same"); }
        if(d){ d = (TProfile*) d->Clone(); styleH(d, "#0072B2", markFilledCircle, 1.6); d->SetLineWidth(2); d->Draw("E1 X0 same"); }
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(1.1*ts); tx.DrawLatex(col == 0 ? 0.30 : 0.06, 0.06/(1.0) + 0.04, c == 0 ? "pp" : Form("PbPb %s", clsLab[c]));
        if(col == 0 && d && m){ TLegend *lg = makeLegend(0.30, 0.68, 0.95, 0.90, 0.9*ts); lg->AddEntry(d, "data", "lp"); lg->AddEntry(m, c == 0 ? "PYTHIA / P+H" : "MC", "lp"); lg->Draw(); }
      } else {
        const int k = row - 1;
        TProfile *d = getP(fd, k, c), *m = getP(fm, k, c);
        if(d && m){
          TH1D *df = new TH1D(Form("df%d_%d", row, col), "", NB, edges);
          for(int b = 1; b <= NB; b++){ if(d->GetBinEntries(b) <= 0 || m->GetBinEntries(b) <= 0) continue;
            df->SetBinContent(b, d->GetBinContent(b) - m->GetBinContent(b)); df->SetBinError(b, hypot(d->GetBinError(b), m->GetBinError(b))); }
          styleH(df, "#000000", markFilledCircle, 1.4); df->SetLineWidth(2); df->Draw("E1 X0 same");
        }
      }
    }
    for(int b = 1; b <= NB; b++){
      TProfile *d0 = getP(fd, 0, c), *m0 = getP(fm, 0, c);
      if(!d0 || !m0 || d0->GetBinEntries(b) <= 0 || m0->GetBinEntries(b) <= 0) continue;
      printf("%-8s %3.0f-%3.0f   %.3f +- %.3f (%5.0f)   %.3f +- %.3f            ", clsLab[c], edges[b-1], edges[b], d0->GetBinContent(b), d0->GetBinError(b),
             ((TH1D*) fd->Get(Form("h_nJets_C%d", c)))->GetBinContent(b), m0->GetBinContent(b), m0->GetBinError(b));
      for(int k = 0; k < 4; k++){ TProfile *d = getP(fd, k, c), *m = getP(fm, k, c); printf("%+.3f  ", d->GetBinContent(b) - m->GetBinContent(b)); }
      printf("\n");
    }
  }
  cv->cd(); TLatex t; t.SetNDC(); t.SetTextFont(62); t.SetTextSize(0.03); t.DrawLatex(L, 1 - 0.7*T, "CMS");
  t.SetTextFont(52); t.SetTextSize(0.026); t.DrawLatex(L + 0.045, 1 - 0.7*T, "Internal");
  t.SetTextFont(42); t.SetTextAlign(31); t.DrawLatex(1 - R, 1 - 0.7*T, "PF jets, |#eta| < 1.6, 5.02 TeV");
  gSystem->mkdir(Form("%s/figures/jetID", repo), true);
  cv->SaveAs(Form("%s/figures/jetID/chargedFraction_dataVsMC.pdf", repo));
}
}

void studyJetChargedFraction_PF(TString mode = "fill", TString a = "PH", TString b = "", int maxFiles = 0, TString c = ""){
  if(mode == "fill") jcf::fill(a, b.Data(), maxFiles, c);
  else jcf::plot(a, b, c);
}
