#pragma once
// Shared reco-muon selection helpers for PYTHIA and PYTHIAHYDJET MC scans.
// Requires config globals: muPtCut, muPtMaxCut, fillMu12, fillMu5, fillMu7,
//   skipGenParticles, epsilon, epsilon_mm.

// Returns true if muon m passes matchFlagR check, optional kinematic cuts,
// and muon quality cuts (tight or hybridSoft per fillMu* flags).
// checkKinematic=false skips pT range and eta cuts (PYTHIAHYDJET gen-jet loop omits them).
inline bool passesRecoMuonCuts(eventMap* em, int m, int* matchFlagR, bool checkKinematic = true) {
    if(matchFlagR[m] == 1) return false;
    if(checkKinematic){
        if(em->muPt->at(m) < muPtCut || em->muPt->at(m) > muPtMaxCut || fabs(em->muEta->at(m)) > 2.0) return false;
    }
    if(fillMu12){
        if(!isQualityMuon_tight(em->muChi2NDF->at(m),
                                em->muInnerD0->at(m),
                                em->muInnerDz->at(m),
                                em->muMuonHits->at(m),
                                em->muPixelHits->at(m),
                                em->muIsGlobal->at(m),
                                em->muIsPF->at(m),
                                em->muStations->at(m),
                                em->muTrkLayers->at(m))) return false;
    } else if(fillMu5 || fillMu7){
        if(!isQualityMuon_hybridSoft(em->muChi2NDF->at(m),
                                     em->muInnerD0->at(m),
                                     em->muInnerDz->at(m),
                                     em->muPixelHits->at(m),
                                     em->muIsTracker->at(m),
                                     em->muIsGlobal->at(m),
                                     em->muTrkLayers->at(m))) return false;
    }
    return true;
}

// Returns true if any gen-level muon (|pdgID|==13) within DR < epsilon of (muEta, muPhi)
// passes the given kinematic bounds. Defaults mean "no cut" (pass all gen muons).
// guardedBySkip=false bypasses the skipGenParticles gate — use for loops that lack it.
inline bool findGenMuonMatch(
    eventMap* em,
    double muEta, double muPhi,
    double genMuPtMin = 0.0,
    double genMuPtMax = 9999.0,
    double genMuEtaMax = 9999.0,
    bool guardedBySkip = true
) {
    if(guardedBySkip && skipGenParticles) return false;
    for(int j = 0; j < (int)em->gpptp->size(); j++){
        if(TMath::Abs(em->gppdgIDp->at(j)) != 13) continue;
        if(em->gpptp->at(j) < genMuPtMin) continue;
        if(em->gpptp->at(j) > genMuPtMax) continue;
        if(fabs(em->gpetap->at(j)) > genMuEtaMax) continue;
        if(getDr(muEta, muPhi, em->gpetap->at(j), em->gpphip->at(j)) < epsilon) return true;
    }
    return false;
}
