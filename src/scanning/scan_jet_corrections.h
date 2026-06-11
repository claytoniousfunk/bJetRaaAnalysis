#pragma once
// Applies JEU shift, JER smear, BJet neutrino energy shift, and JERCorrection
// to a JEC-corrected jet pT.
// Requires config globals: apply_JEU_shift_up, apply_JEU_shift_down,
//   apply_JER_smear, doBJetNeutrinoEnergyShift, doJERCorrection.
// Returns updated pT, or -1.0 if the jet should be skipped (neutrino dice-roll failed).

inline double applyJEU_JER(
    double x,
    JetUncertainty& JEU,
    TF1* JER_fxn,
    TRandom* randomGenerator,
    TF1* fitFxn_PYTHIA_JERCorrection,
    TH1D* neutrino_tag_fraction,
    TH2D* neutrino_energy_map
) {
    double correctedPt_up = 1.0;
    double correctedPt_down = 1.0;
    if(apply_JEU_shift_up){
        correctedPt_up = x * (1 + JEU.GetUncertainty().second);
        x = correctedPt_up;
    } else if(apply_JEU_shift_down){
        correctedPt_down = x * (1 - JEU.GetUncertainty().first);
        x = correctedPt_down;
    }

    double mu = 1.0;
    double sigma = 0.2;
    double smear = 0.0;
    if(apply_JER_smear){
        sigma = 0.663*JER_fxn->Eval(x);
        smear = randomGenerator->Gaus(mu, sigma);
        x = x * smear;
    }

    double skipDoBJetNeutrinoEnergyShift_diceRoll = 0.0;
    double smear_doBJetNeutrinoEnergyShift = 0.0;
    if(doBJetNeutrinoEnergyShift){
        skipDoBJetNeutrinoEnergyShift_diceRoll = randomGenerator->Rndm();
        if(skipDoBJetNeutrinoEnergyShift_diceRoll > neutrino_tag_fraction->GetBinContent(neutrino_tag_fraction->FindBin(x))) return -1.0;
        TH1D* neutrino_energy_map_proj = (TH1D*) neutrino_energy_map->ProjectionX(
            "neutrino_energy_map_proj",
            neutrino_energy_map->GetYaxis()->FindBin(x),
            neutrino_energy_map->GetYaxis()->FindBin(x)+1);
        smear_doBJetNeutrinoEnergyShift = neutrino_energy_map_proj->GetRandom();
        x += smear_doBJetNeutrinoEnergyShift;
    }

    double mu_JERCorrection = 1.0;
    double sigma_JERCorrection = 0.2;
    double smear_JERCorrection = 0.0;
    double k_JERCorrection = 0.0;
    if(doJERCorrection){
        k_JERCorrection = TMath::Sqrt(fitFxn_PYTHIA_JERCorrection->Eval(x)*fitFxn_PYTHIA_JERCorrection->Eval(x) - 1.);
        sigma_JERCorrection = k_JERCorrection*JER_fxn->Eval(x);
        smear_JERCorrection = randomGenerator->Gaus(mu_JERCorrection, sigma_JERCorrection);
        x = x * smear_JERCorrection;
    }

    return x;
}
