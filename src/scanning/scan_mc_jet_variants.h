#pragma once
// Computes all NTemplateIndices jet pT variants (nominal, JER-smear, JEU-shift-up, JEU-shift-down)
// used in PYTHIA and PYTHIAHYDJET MC scans.
// Caller must call JEU.SetJetPT/Eta/Phi before calling.

// PYTHIA: single JER function
inline void computeJetPtVariants(
    double recoJetPt_i,
    JetUncertainty& JEU,
    TRandom* randomGenerator,
    TF1* JER_fxn,
    double& recoJetPt_JEUShiftUp_i,
    double& recoJetPt_JEUShiftDown_i,
    double& recoJetPt_JERSmear_i
) {
    recoJetPt_JEUShiftUp_i   = recoJetPt_i * (1 + JEU.GetUncertainty().second);
    recoJetPt_JEUShiftDown_i = recoJetPt_i * (1 - JEU.GetUncertainty().first);
    double sigma = 0.663 * JER_fxn->Eval(recoJetPt_i);
    recoJetPt_JERSmear_i = recoJetPt_i * randomGenerator->Gaus(1.0, sigma);
}

// PYTHIAHYDJET: centrality-indexed JER function array
// sigma defaults to 0.2 when CentralityIndex > 4 (preserves original behavior).
inline void computeJetPtVariants(
    double recoJetPt_i,
    JetUncertainty& JEU,
    TRandom* randomGenerator,
    TF1** JER_fxn,
    int CentralityIndex,
    double& recoJetPt_JEUShiftUp_i,
    double& recoJetPt_JEUShiftDown_i,
    double& recoJetPt_JERSmear_i
) {
    recoJetPt_JEUShiftUp_i   = recoJetPt_i * (1 + JEU.GetUncertainty().second);
    recoJetPt_JEUShiftDown_i = recoJetPt_i * (1 - JEU.GetUncertainty().first);
    double sigma = 0.2;
    if(CentralityIndex <= 4) sigma = 0.663 * JER_fxn[CentralityIndex]->Eval(recoJetPt_i);
    recoJetPt_JERSmear_i = recoJetPt_i * randomGenerator->Gaus(1.0, sigma);
}
