#pragma once
// Jet flavor index enum and helpers for flavor-decomposed template histogram arrays

enum JetFlavorIdx {
  kAllJets  = 0,
  kDJets    = 1,
  kUJets    = 2,
  kSJets    = 3,
  kCJets    = 4,
  kBJets    = 5,
  kBGSJets  = 6,
  kGJets    = 7,
  kXJets    = 8,
  kNJetFlavors = 9
};

static const char* kFlavorNames[kNJetFlavors] = {
  "allJets","dJets","uJets","sJets","cJets","bJets","bGSJets","gJets","xJets"
};

// Returns flavor index 1-8 for a known specific flavor, -1 for unknown.
// Never returns kAllJets (0) — allJets is always filled separately.
static int getFlavorIdx(int jetFlavorInt) {
  if(fabs(jetFlavorInt) == 1) return kDJets;
  if(fabs(jetFlavorInt) == 2) return kUJets;
  if(fabs(jetFlavorInt) == 3) return kSJets;
  if(fabs(jetFlavorInt) == 4) return kCJets;
  if(fabs(jetFlavorInt) == 5) return kBJets;
  if(jetFlavorInt == 17)      return kBGSJets;
  if(jetFlavorInt == 21)      return kGJets;
  if(jetFlavorInt == 0)       return kXJets;
  return -1;
}
