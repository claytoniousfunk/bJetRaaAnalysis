# bJetRaaAnalysis — working notes

PbPb b-jet R_AA via muon-tagged jets. CMS, 5.02 TeV. The open question driving
most current work: PbPb 10-30% looks identical to 0-10% in the jets-per-Z /
PbPb-over-pp ratio, which should not happen.

## Layout

| path | what |
|---|---|
| `src/scanning/PbPb/PbPb_pfCandAnalyzer.C` | main PbPb scan (PF candidates, FastJet, mixed events) |
| `src/scanning/PbPb/PbPb_scan.C`, `src/scanning/pp/pp_scan.C` | the non-pfCand scans |
| `headers/AnalysisSetup/` | physics constants (`common.h`), `pseudoJets.h` for mixing/FastJet knobs |
| `headers/config/config_PbPb.h` | run-time flags (filters, triggers, jet-collection overrides) |
| `headers/plotting/` | shared plot helpers — **use these, do not re-paste** |
| `src/plots/` | plotting macros, one per figure family |
| `rootFiles/scanningOuput/` | scan outputs — **note the typo, "Ouput" not "Output"** |
| `figures/` | output figures (gitignored) |

**Before writing a macro for a measurement, read the reuse rules and the
measurement → macro index:** @MACROS.md

Presentations live outside the repo: `~/Documents/nuclear/GroupMeeting/2026/`
with figures in `../figures/<date>/`.

## Conventions

- **Figures are PDF.** PNG only as a throwaway for inline viewing.
- **No triangle markers.** Centre-symmetric only: circle, square, diamond, cross.
- **Colour-blind safe.** Okabe-Ito for series, kViridis/kCividis for 2D maps.
  Both are in `headers/plotting/plotStyle.h` — include it rather than pasting hex.
- **Linear axes** unless there is a reason; a log axis hides negative bins,
  which are usually the thing worth seeing.
- **Commits carry no Claude authorship trailer.**
- `bJetMuonTaggingAnalysis` (the sibling repo) is **read-only**. Fix things here.

## Shared headers — prefer these

```cpp
#include "../../../headers/plotting/coarseCent.h"   // class mapping, slice sums
#include "../../../headers/plotting/plotStyle.h"    // palette, markers, rebinning
#include "../../../headers/plotting/ratioPanel.h"   // ratio + pad splitting
#include "../../../headers/functions/divideByBinwidth.h"
```

Before these existed the coarse-class mapping was duplicated in 31 macros and
the Okabe-Ito hex codes in 52. The same ratio-panel bug was then written twice
independently. Older macros still carry local copies; leave them unless you are
already editing one.

## Traps that have cost real time

**Trigger gating and the matching denominator.** `h_vz` is filled *before*
`evtTriggerDecision` is computed, so it counts every event. A numerator gated on
`evtTriggerDecision` must be divided by `h_vz_triggerOn`, not `h_vz`. The
trigger fraction is 0.120 in 0-10% and 0.025 in 50-80%, so the error does not
cancel across centrality. Currently gated: `h_mixedMuonPtRel_recoJetPt`,
`h_realMuonPtRel_mixedFastJetPt`, `h_muonDR_inclusiveClosestJet_triggerOn`, the
data `..._triggerOn` histograms, and anything requiring `hasFastJetRecoMuonTag`.
Ungated: `h_muonDR_inclusiveClosestJet`, `h_fastJetMuonDR_inclusiveClosestFastJet`,
the injection-study histograms.

**ptRel does not depend on the jet pT scale.** `getPtRel` divides by `|jet|²`,
so the magnitude cancels; only the jet *direction* and the muon momentum matter.
A raw-vs-JEC change moves which jets land in a pT window, not the ptRel values.

**`TH1::Rebin` with an edge array returns a NEW histogram** and leaves the
original alone. Release it or it leaks per call. Every edge must land on an
existing bin boundary of the input.

**Variable bin widths need `Scale(1,"width")` before drawing**, or a wide bin
reads as a tall one — and do it *after* any integrals, which are counts.

**`FindBin(hi)` returns the bin that starts at `hi`.** Use `FindBin(hi - 1e-6)`
for a window upper edge, or every window is silently one bin too wide.

**Scan generation matters.** Several 2026-09-08 commits changed what existing
histograms contain (jet pT axis raw→JEC, tagging muon last→leading). Files from
different generations must not be compared.

Scans from 2026-09-09 onward carry a `provenance` TNamed recording the git hash,
whether the tree was dirty, and every config flag that changes the output:

```
root -l file.root
provenance->GetTitle()
```

For older files, the marker for the post-2026-09-08 generation is the presence
of `h_mixedMuonPtRel_recoJetPt`, which is booked unconditionally and is simply
empty in a same-event scan. `plotFastJetMuonPtRel_sameVsMixedEvent_coarseCent.C`
uses that marker to refuse mismatched pairs.

**`scanningOuput`** is misspelled in the repo. Paths will fail silently if you
"correct" it.

## Z normalization

Since 2026-11-09 the nominal N_Z is the **published** CMS Z yield converted with
the luminosity, not our own Z->mumu count: PbPb from PRL 127 (2021) 102002,
pp from the 5.02 TeV total cross section times the |y_Z| < 2.1 fraction
(0.763, Pythia 8). Everything is in `headers/publishedData/CMS_ZYield_5p02TeV.h`.
`BJET_NZ_SOURCE=own` restores our count (outputs `_ownNZ`). N_Z cancels in the
double ratio R_Z(b)/R_Z(incl). `calculateJetsPerZ.cc` and the legacy calo chains
still use our own count.

## ptRel fit binning

Since 2026-11-12 the CSVv2-cut ptRel fits (`fitPtRelCSVv2_2D_pp_PF.C` and everything
that includes it) use **0.2 GeV ptRel bins** by default; 0.5 GeV before. Untagged
outputs written before that date are 0.5 GeV. `BJET_PTREL_BIN=0.5` reproduces the
old fits under a `_bin0p5` suffix.

## Nominal b purity

Since 2026-11-12 the nominal b purity is the **2D muon ptRel x CSVv2 fit** (4 CSVv2
categories, HistFactory, MC tagger shapes), with its own purity systematics from
`src/systematics/systematics_CSVv2_2D_PF.C` and **no working-point and no fit-range source** (fitRangeHigh
zeroed 2026-11-16; the AN uses this nominal since then). Tag
`_csvV2_2D`, chain variant `BJET_PURITY_VARIANT=csvV2__2D` (double underscore). The
CSVv2-loose cut fit (`_csvV2loose`, nominal 2026-10-28 to 11-12) is untouched; several
results macros still default to its tag, so pass `_csvV2_2D` explicitly. The 2D model
is rejected by the data in pp and central 80-120 GeV. Recipe in MACROS.md.

## Nominal b-jet binning

Since 2026-11-16 (user: "make 7 bins nominal") the b-jet results are reported in **7 jet-pT bins**,
80 90 100 110 120 135 150 200 GeV, with the re-tuned N (pp 2; PbPb 1, 2, 3, 3 from 50-80% to 0-10%).
Results tag `_csvV2_2D7bins_nIterOpt`; the AN uses it. The **code default is still 4 bins**: every
7-bin chain step needs `root -l -b -q -e '#define BJET_PT_7BINS 1' ...` with
`BJET_PURITY_VARIANT=csvV2__2D7bins BJET_NITER=opt` (the chain refuses a tag without "7bins"),
and the purity-study macros follow the same define (`ptBinTag()` gives their outputs `_7bins`).
The 4-bin `_csvV2_2D` files remain as the previous nominal. Recipe in MACROS.md.

## Key constants (`headers/AnalysisSetup/common.h`)

`etaMax = 1.6` · `jetPtCut = 20` · `muPtCut = 15` · `muPtMaxCut = 999` ·
`epsilon_mm = 0.4` (muon-jet ΔR match) · ultraFine centrality =
17 indices (C0 inclusive + 16 slices of 5%)

Coarse classes: `0-10% = slices 1-2`, `10-30% = 3-6`, `30-50% = 7-10`,
`50-80% = 11-16`. C0 is the scan's own inclusive bin and is **not** the sum of
1-16.

## Scan workflow

Scans run on lxplus and take hours, and anything that adds or changes a
histogram needs a full rescan before it can be used. The histogram set is
largely settled as of 2026-09, so this is not usually a bottleneck; if a stretch
of work starts needing several new histograms, batch them into one scan rather
than paying a round trip each.

**Type-check before submitting.** The `#ifdef DO_FASTJET` block is the
highest-risk code here to edit, and it is now checkable locally:

```bash
tools/syntaxCheck.sh                          # PbPb_pfCandAnalyzer.C
tools/syntaxCheck.sh src/scanning/pp/pp_scan.C
```

Exit 0 is clean. It is `-fsyntax-only`, so it catches typos, wrong argument
counts, out-of-scope symbols and include-order mistakes — not logic errors and
nothing that only appears at link time. It caught a real include-order bug on
its first run (`writeProvenance.h` was included before the config headers whose
globals it reads).

Needs FastJet headers; `configure` alone suffices, no build required:

```bash
cd ~/Programs && curl -LO https://fastjet.fr/repo/fastjet-3.4.3.tar.gz
tar xzf fastjet-3.4.3.tar.gz && cd fastjet-3.4.3
./configure --prefix=$HOME/Programs/fastjet
```

The pfCand scan also needs fjcontrib (ConstituentSubtractor, for the forest
akCs4PF emulation in `headers/functions/csJetEmulation.h`). `syntaxCheck.sh`
uses `~/Programs/fastjet330` (FastJet 3.3.0 + fjcontrib 1.033, the forest's
versions) when it exists; the scan's Makefile and ACLiC drivers link
`-lfastjetcontribfragile` from `FJCONTRIB` / `FJCONTRIB_HOME`, default the
FastJet prefix.

The script generates a driver adding `using namespace std;`, because ROOT's
interpreter has it in effect and `eventMap.h` relies on that (bare
`vector<Float_t>`). Without it g++ reports dozens of spurious errors.

Jet collection is chosen at load time and swaps the *tree*, not the histograms:

```
useCaloJetsOverride  -> akPu4CaloJetAnalyzer/t     (adds _caloJets to the filename)
useFlowJetsOverride  -> akFlowPuCs4PFJetAnalyzer/t
default              -> akCs4PFJetAnalyzer/t
```

So calo and PF live in **separate files** under the same histogram names; there
are no calo-specific histograms. `useCaloJetsOverride` needs its own scan.

The PF-jet calo-match ID works the same way: `requireCaloJetMatch` (PbPb,
PYTHIAHYDJET, pp and PYTHIA configs, `src/scanning/scan_calo_jet_match.h`) drops
PF jets with no calo jet (akPu4Calo in PbPb, ak4Calo in pp) within dR < 0.2
from every reco-jet histogram and adds `_caloJetMatched` to the filename. Apply
it on both sides of an R_AA or the ID's genuine-jet loss reads as suppression. Same names, different contents — compare
matched and unmatched only across files, never assume from the key list.

## Current state of the ptRel decomposition

`S = D − T2 − T3`, where `D` is the measured μ+jet ptRel and

- `T2` = `h_injMuonPtRel_donorJetPt_inject` — mixed-event μ **reclustered into**
  the next centrality-matched event's PF candidates. Ungated, so consumers must
  scale it by the trigger fraction when everything else is gated.
- `T3` = `h_fastJetMuonPtRel_fastJetPt_PF_bkgSub_RC` — mixed μ + mixed fastJet

Known open problems, all documented in
`src/plots/muonPtRel/plotPtRelTemplateDecomposition.C`:

1. **Geometry.** T2 pairs a muon against a jet axis built without it, so its
   ⟨ptRel⟩ is ~5.1 GeV against the data's ~1.6 and it barely subtracts. The
   muon-injection study measures this; it confirms the axis effect (ΔR −65%)
   but ptRel moves little in central collisions.
2. **Muon population.** Mixed-event muons are PF candidates with no tight ID and
   no W-decay filter, so they are a looser population than the data's → the
   templates over-subtract.
3. **Double subtraction.** T2's denominator is all reco jets, genuine + fake, so
   its fake-jet part overlaps T3.
4. **50-80% b-purity moves the wrong way** under subtraction (purity rises).
5. The b-purity **fit range moves the unsubtracted baseline** (0-10%, 80-90 GeV:
   0.767 / 0.821 / 0.843 at 0-3 / 0-4 / 0-5 GeV). Current range is 0-4.

The subtraction's effect on purity is **small**: T2+T3 together are ~5% of D
in 0-10% at 80-90 GeV, and purity drops 3-8% in 0-10%, 2-4% in 10-30%, <1% in
30-50%. An earlier "~20-35% drop" (and subtracted purities of ~0.62) came from a
bin-offset bug in `plotBPurity_fakeMuonSubtractionEffect.C`, fixed 2026-10-04:
the 2026-2-12 data file books jet pT on [20,500] and the template files on
[0,500], so reusing the data's bin indices read the templates 20 GeV low.
**Always project a template's jet-pT window from its own axis.**
