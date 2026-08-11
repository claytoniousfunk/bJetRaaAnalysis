# Running PbPb_pfCandAnalyzer on condor

Run everything from `src/scanning/PbPb/`. The analyzer resolves its input list
relative to this directory.

---

## 0. Environment

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
```

Needed for `root-config` and `fastjet-config`. If you use a different view,
change `env_setup` in `condor_PbPb_pfCandAnalyzer.py` to match — the job scripts
source that string, and a mismatch means the jobs build against a different
stack than you tested with.

---

## 1. Build

```bash
make -f Makefile.pfCandAnalyzer
```

Produces `pfCandAnalyzer_PbPb`. Takes a few minutes: one translation unit, all
of ROOT + FastJet + the analysis headers, at `-O2`. No object caching, so every
build is from scratch.

**Rebuild after editing `PbPb_pfCandAnalyzer.C`, any header it includes, or
after a `git pull`.** The Makefile tracks header dependencies automatically
(`-MMD`), so `make` notices config changes such as
`doHiBinReweightToHardProbesJet80` in `headers/config/config_PbPb.h`. If you
ever suspect a stale binary, `make -f Makefile.pfCandAnalyzer clean` first —
the output directory name printed at startup encodes the active config flags
and is the quickest way to confirm which settings the binary actually has.

For a faster edit/test loop, build unoptimized:

```bash
make -f Makefile.pfCandAnalyzer CXXFLAGS="-O0 -g -DDO_FASTJET $(root-config --cflags) $(fastjet-config --cxxflags) -Wno-sign-compare -Wno-empty-body -Wno-ignored-qualifiers"
```

Rebuild at `-O2` before submitting — the mixing loop is CPU-heavy enough that it
matters for wall time.

---

## 2. Verify the FastJet linkage

```bash
ldd pfCandAnalyzer_PbPb | grep -Ei "fastjet|delphes"
```

Expect libfastjet and libfastjettools from
`/cvmfs/sft.cern.ch/lcg/releases/fastjet/...`, and **no Delphes**.

If Delphes appears, stop. The LCG views ship Delphes, which embeds its own
FastJet build and claims the `fastjet::` namespace in its rootmap;
`ClusterSequence`'s constructor is a header template, so the weak symbol binds
to Delphes' instantiation, whose `PseudoJet` layout differs. It segfaults on the
first clustering call. Avoiding that is the reason this runs as a standalone
binary instead of an ACLiC macro — see the comment block in
`main_pfCandAnalyzer.cc`.

---

## 3. Run one file interactively

```bash
time ./pfCandAnalyzer_PbPb 1
```

The argument is the 1-based index into the input list. Reference point: file 1
holds 6596 events and took **111 min** wall / 75 min CPU. Note your own timing —
it sets the JobFlavour and tells you whether to optimize before committing farm
time.

Confirm the clustering actually filled, not just that it exited cleanly:

```bash
root -l -b -q -e 'auto f=TFile::Open("<output path printed by the run>"); ((TH1D*)f->Get("h_fastJetPt_PF_bkgSub_RC_C1"))->Print();'
```

A clean exit with an empty FastJet histogram means the build lost `-DDO_FASTJET`.

---

## 4. Canary submission

Edit `condor_PbPb_pfCandAnalyzer.py`:

```python
njobs_max   = 10       # or 1 for a first look
auto_submit = False
```

```bash
python3 condor_PbPb_pfCandAnalyzer.py
condor_submit PbPb_pfCandAnalyzer_mixedEvent_fastJet/condor_submit.cfg
```

The script builds the binary first and refuses to generate jobs if the build
fails, so a compile error surfaces here rather than on every worker.

The canary is what actually tests the worker environment: AFS visibility, the
EOS write with your proxy, and the output-directory mkdir race (many jobs
racing to create the same directory — the analyzer treats "someone else made
it" as success).

---

## 5. Monitor

```bash
condor_q -nobatch                 # I = idle, R = running, H = held
condor_tail <cluster>.<proc>      # live stdout: the `evt frac:` progress lines
condor_ssh_to_job <cluster>.<proc>   # then `top` — CPU-bound vs I/O wait
```

**`.out` and `.err` files do not exist until a job completes.** CERN HTCondor
buffers stdout in the job's sandbox and transfers it back at the end. Only the
`.log` file is written live by the schedd. Use `condor_tail` for a running job;
an empty `log/` directory mid-run is normal, not a failure.

Useful checks in the `.log`:

```bash
grep -c "Job executing" PbPb_pfCandAnalyzer_mixedEvent_fastJet/log/job_*.log
```

A count above 1 means the job was evicted and restarted. There is no
checkpointing, so a restart begins from event zero — a job preempted more often
than its runtime will never finish while appearing to run for hours.

When jobs finish:

```bash
grep "Number of events" PbPb_pfCandAnalyzer_mixedEvent_fastJet/log/job_*.out
```

Input files vary in size, so per-job runtimes vary with them.

---

## 6. Full submission

Once the canary jobs produce output files:

```python
njobs_max   = None
auto_submit = True
```

```bash
python3 condor_PbPb_pfCandAnalyzer.py
```

1993 input files, one job each.

---

## 7. Merge and analyze

```bash
hadd merged.root /eos/cms/store/group/phys_heavyions/cbennett/scanningOutput/<output dataset>/PbPb_pfCandAnalyzer_output_*.root
```

Then point `fmixed_path` in `src/plots/jetPt/pseudoJets/makeFakeJetFile.C` at
the merged file. `h_vz_C*` are plain counts and add correctly, so the per-event
normalization stays valid without other changes.

---

## Settings reference

`condor_PbPb_pfCandAnalyzer.py`, top of file:

| setting | notes |
|---|---|
| `dblist` | Must match the list the analyzer picks via its `doMinBiasSample` / `doHardProbesSample` / `doSingleMuonSample` flags. The job index is matched against `ifile` inside the analyzer, so a mismatch silently processes the wrong files. |
| `env_setup` | Sourced at the top of every job script. Must match the view you built with. |
| `use_standalone_binary` | Leave `True`. `False` is the ACLiC path, which segfaults on any view shipping Delphes. |
| `force_fastjet_preload` | Only relevant on the ACLiC path. Adds `LD_PRELOAD` for libfastjet. |
| `time_flavour` | espresso 20m, microcentury 1h, longlunch 2h, workday 8h, tomorrow 1d. Longer flavours match fewer slots and queue longer — pick the shortest that comfortably covers your measured runtime. |
| `request_memory` | MB. Observed image size is ~730 MB. Larger requests match fewer slots. |
| `nsplit` | Input files per job. |
| `njobs_max` | Cap for canary runs. `None` = all files. |
| `auto_submit` | `False` prints the `condor_submit` command instead of running it. |

---

## Known issues

**Mixed-event running is I/O heavy.** The pool is rebuilt for every event:
~100 full `getEvent` calls per analyzed event, each reading the whole event
including all PF candidates. That is roughly a 100× I/O amplification and
dominates the runtime. Caching the pool per centrality class and refreshing it
every N events would cut it substantially. Not yet implemented.

**The mixed-event scan is hiBin-reweighted to HardProbes Jet80.** Within a
coarse centrality class this shifts composition toward the central edge, which
inflates the per-event jet rate by 1.05 (0–10%) up to 1.71 (50–80%) relative to
an unweighted scan. Subtracting fakes normalized this way over-subtracts by that
factor. It cancels if the fake *fraction* is formed per ultra-fine slice
(numerator and denominator carry the same weight there) and the slices are then
recombined with flat MinBias weights. Verified against the unweighted
`sameEventPFClustering` scan; the two forests themselves agree to within a few
percent, so the reweighting is the whole effect.
