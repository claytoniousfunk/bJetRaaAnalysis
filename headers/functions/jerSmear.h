// Jet energy resolution scheme of the PF unfolding (2026-11-05, user), after CMS-HIN-18-014
// (arXiv:2102.13080, Sec. 5, item 2b): the resolution in data is 10-15% worse than in the
// simulation, so the NOMINAL response is smeared to a 10% worse resolution, and the JER
// systematic is a FURTHER 10% on top of that nominal (1.1 x 1.1 = 1.21 x the simulation).
// Until 2026-11-05 the nominal was the unsmeared response and the variation +20%
// (extra width sqrt(1.2^2 - 1) = 0.663 sigma).
//
// The smearing: every generator-level column of the response gets an extra Gaussian smearing of
// its reconstructed axis, of width extra x the column's own relative resolution x pT^gen
// (sigma from the reco / gen spread of the column), bin-integrated from each reco bin centre.
// Weight smeared past the reco range goes to its under / overflow, so each column's total over
// the gen bins is kept and the gen projection (prior, tag frequency) is unchanged. Applying it to an already smeared
// response measures that response's own sigma, so variation(nominal(R)) is 1.1 x 1.1 x sigma.
//
// INPUT of the variation (2026-11-05, user): the JER source is measured on the MC's own reco
// spectrum (the nominal response's reco projection, i.e. the MC truth folded with the nominal
// response), unfolded with the nominal and with the varied response at the chain's N:
// var / nom - 1 of those two. On data the same ratio also carries the stitching, the fake-jet
// subtraction and the data's fluctuations; the MC isolates the added smearing. The relative shift
// is then applied to the data's nominal. BJET_JER_INPUT=data = the earlier data-based variation.
//
// BJET_JER_SCHEME=old  restores the pre-2026-11-05 scheme (nominal unsmeared, variation +20%),
//                      to reproduce the old outputs.
//
// Used by inclusiveJetsPerZSys_PF.C, unfoldBJetSpectrum_caloJets_pp.C (PFJets),
// bJetChain_caloJets_PbPb.C (PFJets) and bJetsPerZSys_PF_PbPb.C (PFJets). The calo chain is
// unchanged. bJetsPerZSys_PF_PbPb.C rebuilds the chain's unfolding and stops if it is not
// reproduced, so the chain and the systematics must run with the same scheme.
#ifndef JERSMEAR_H
#define JERSMEAR_H

#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

namespace jerSmear {
// relative resolution of a response per generator-level column, sigma(pT^reco / pT^gen), computed exactly as smear()
// measures the width it widens (reco bin centres, gen bin centre); for the extraction figures (2026-11-17)
inline TH1D* resolution(TH2D *r, const char *name){
  TH1D *o = r->ProjectionY(name, 1, r->GetNbinsX()); o->Reset(); o->SetDirectory(nullptr);
  for(int iy = 1; iy <= r->GetNbinsY(); iy++){
    const double g = r->GetYaxis()->GetBinCenter(iy); double sw = 0., s1 = 0., s2 = 0.;
    for(int ix = 1; ix <= r->GetNbinsX(); ix++){ const double w = r->GetBinContent(ix, iy); if(w <= 0.) continue; const double x = r->GetXaxis()->GetBinCenter(ix)/g; sw += w; s1 += w*x; s2 += w*x*x; }
    if(sw <= 0.) continue;
    const double m = s1/sw; o->SetBinContent(iy, std::sqrt(TMath::Max(0., s2/sw - m*m))); o->SetBinError(iy, 0.);
  }
  return o;
}

inline bool oldScheme(){
  static const bool o = std::getenv("BJET_JER_SCHEME") && TString(std::getenv("BJET_JER_SCHEME")) == "old";
  return o;
}
inline bool onData(){
  static const bool d = std::getenv("BJET_JER_INPUT") && TString(std::getenv("BJET_JER_INPUT")) == "data";
  return d;
}
inline const char* inputText(){ return onData() ? "JER variation measured on the data" : "JER variation measured on the MC reco spectrum"; }
inline double nomScale(){ return oldScheme() ? 1.0 : 1.10; }   // nominal resolution / simulation
inline double varScale(){ return oldScheme() ? 1.2 : 1.10; }   // variation resolution / nominal
inline double extraOf(double scale){ return std::sqrt(TMath::Max(0., scale*scale - 1.)); }
inline double nomExtra(){ return extraOf(nomScale()); }        // 0.458
inline double varExtra(){ return extraOf(varScale()); }        // 0.458 (old: 0.663)
inline const char* varLabel(){ return oldScheme() ? "jet energy resolution (+20%)" : "jet energy resolution (+10% on the 10%-smeared nominal)"; }
inline const char* tagText(){ return oldScheme() ? "JER nominal unsmeared, variation +20%" : "JER nominal x1.10 (HIN-18-014), variation x1.10 on it"; }

// the smearing: the former smearResp / smearRespJER of the PF macros, verbatim apart from the
// width argument (so BJET_JER_SCHEME=old reproduces the old variation bin for bin). Gen under /
// overflow rows and reco under / overflow INPUT cells are not carried (the old behaviour; the
// unfolding reads neither); smeared weight that leaves the reco range lands in its under / overflow.
inline TH2D* smear(TH2D *r, const char *name, double extra){
  TH2D *o = (TH2D*) r->Clone(name); o->Reset(); o->SetDirectory(nullptr);
  TAxis *ax = r->GetXaxis();
  for(int iy = 1; iy <= r->GetNbinsY(); iy++){
    double g = r->GetYaxis()->GetBinCenter(iy), sw = 0., s1 = 0., s2 = 0.;
    for(int ix = 1; ix <= r->GetNbinsX(); ix++){
      double w = r->GetBinContent(ix, iy); if(w <= 0.) continue;
      double x = ax->GetBinCenter(ix)/g; sw += w; s1 += w*x; s2 += w*x*x;
    }
    if(sw <= 0.) continue;
    double mean = s1/sw, sig = extra*std::sqrt(TMath::Max(0., s2/sw - mean*mean))*g;
    for(int ix = 1; ix <= r->GetNbinsX(); ix++){
      double w = r->GetBinContent(ix, iy), e2 = std::pow(r->GetBinError(ix, iy), 2); if(w == 0.) continue;
      double x0 = ax->GetBinCenter(ix);
      if(sig <= 0.){ o->SetBinContent(ix, iy, o->GetBinContent(ix, iy) + w); continue; }
      for(int jx = 0; jx <= r->GetNbinsX()+1; jx++){
        double lo = (jx == 0) ? -1e9 : ax->GetBinLowEdge(jx), hi = (jx == r->GetNbinsX()+1) ? 1e9 : ax->GetBinUpEdge(jx);
        double p = 0.5*(TMath::Erf((hi - x0)/(std::sqrt(2.)*sig)) - TMath::Erf((lo - x0)/(std::sqrt(2.)*sig)));
        if(p < 1e-9) continue;
        o->SetBinContent(jx, iy, o->GetBinContent(jx, iy) + p*w);
        o->SetBinError(jx, iy, std::sqrt(std::pow(o->GetBinError(jx, iy), 2) + p*p*e2));
      }
    }
  }
  return o;
}

// the nominal response: smeared by nomExtra(); returns a new histogram (the caller owns both)
inline TH2D* nominal(TH2D *r, const char *name){
  if(nomExtra() <= 0.){ TH2D *o = (TH2D*) r->Clone(name); o->SetDirectory(nullptr); return o; }
  return smear(r, name, nomExtra());
}
// the JER variation of a NOMINAL response
inline TH2D* variation(TH2D *rNom, const char *name){ return smear(rNom, name, varExtra()); }

inline void announce(const char *who){
  printf("%s: %s (extra widths nominal %.3f, variation %.3f x sigma); %s\n", who, tagText(), nomExtra(), varExtra(), inputText());
}

}   // namespace jerSmear
#endif
