#pragma once
// Centrality bin scheme — change CENT_SCHEME to switch all scan files at once.

#define CENT_NOMINAL   1
#define CENT_ULTRAFINE 2
#define CENT_PERIPH90  3

#define CENT_SCHEME CENT_NOMINAL

#if CENT_SCHEME == CENT_NOMINAL
#include "../../headers/AnalysisSetupV2p3.h"
#include "../../headers/functions/getCentBin.h"
#define CENT_SCHEME_SUFFIX ""
#elif CENT_SCHEME == CENT_ULTRAFINE
#include "../../headers/AnalysisSetupV2p4.h"
#include "../../headers/functions/getCentBin_V2p4.h"
#define CENT_SCHEME_SUFFIX "_ultraFineCentBins"
#elif CENT_SCHEME == CENT_PERIPH90
#include "../../headers/AnalysisSetupV2p5.h"
#include "../../headers/functions/getCentBin.h"
#define CENT_SCHEME_SUFFIX "_periph90CentBins"
#endif
