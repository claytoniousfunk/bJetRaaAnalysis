#pragma once
// Finds the first unmatched reco-muon within deltaR < epsilon_mm of a given jet.
// Requires config globals: muPtCut, muPtMaxCut, fillMu12, fillMu5, fillMu7,
//   doWDecayFilter, epsilon_mm.
// Updates matchFlagR[m] = 1 to prevent double-counting.
// Sets muPtRel, muPt, muEta, muPhi, muJetDr via out-params on match.
// Returns true if a muon-tagged jet was found.

inline bool findRecoMuonTag(
    eventMap* em,
    double x, double y, double z,
    int* matchFlagR,
    double& muPtRel, double& muPt, double& muEta, double& muPhi, double& muJetDr
) {
    for(int m = 0; m < em->nMu; m++){
        double muPt_m = em->muPt->at(m);
        double muEta_m = em->muEta->at(m);
        double muPhi_m = em->muPhi->at(m);
        double muJetDr_m = getDr(muEta_m, muPhi_m, y, z);
        if(matchFlagR[m] == 1) continue;
        if(muPt_m < muPtCut || muPt_m > muPtMaxCut || fabs(muEta_m) > 2.0) continue;
        if(fillMu12){
            if(!isQualityMuon_tight(em->muChi2NDF->at(m),
                                    em->muInnerD0->at(m),
                                    em->muInnerDz->at(m),
                                    em->muMuonHits->at(m),
                                    em->muPixelHits->at(m),
                                    em->muIsGlobal->at(m),
                                    em->muIsPF->at(m),
                                    em->muStations->at(m),
                                    em->muTrkLayers->at(m))) continue;
        } else if(fillMu5 || fillMu7){
            if(!isQualityMuon_hybridSoft(em->muChi2NDF->at(m),
                                         em->muInnerD0->at(m),
                                         em->muInnerDz->at(m),
                                         em->muPixelHits->at(m),
                                         em->muIsTracker->at(m),
                                         em->muIsGlobal->at(m),
                                         em->muTrkLayers->at(m))) continue;
        }
        if(doWDecayFilter){
            if(isWDecayMuon(muPt_m, x)) continue;
        }
        if(muJetDr_m < epsilon_mm){
            matchFlagR[m] = 1;
            muPtRel = getPtRel(muPt_m, muEta_m, muPhi_m, x, y, z);
            muPt = muPt_m;
            muEta = muEta_m;
            muPhi = muPhi_m;
            muJetDr = muJetDr_m;
            return true;
        }
    }
    return false;
}
