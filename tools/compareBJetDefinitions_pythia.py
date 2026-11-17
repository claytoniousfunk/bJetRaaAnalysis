"""Truth-level b-jet cross section in pp at 5.02 TeV under OUR definition and under the ATLAS one of
arXiv:2204.13530 (headers/publishedData/ATLAS_bJetXsecPP_2204_13530.h), to size how much of
ours / ATLAS is definition (2026-11-16).

Pythia 8 (PyPI pythia8mc, venv ~/Programs/pythia8venv), HardQCD:all, Monash, pTHat > 40 GeV with
bias2Selection (pTHat^4), anti-kT R = 0.4 (SlowJet). b jet = b hadron (last copy, pT > 5 GeV)
within DeltaR < 0.3 of the jet axis, in every definition, so only the jet pT and the sample differ:

  ours      visible particles (muons in, neutrinos out), all b jets, |eta| < 1.6, per unit eta
  visY21    the same jets, |y| < 2.1, per unit y                 (acceptance alone)
  allY21    every final particle (neutrinos in), all b jets, |y| < 2.1   (neutrinos alone)
  atlas     every final particle, b jets holding a muon from a b-hadron decay chain (b -> mu,
            b -> c -> mu, ...; any muon pT) within DeltaR < 0.4, |y| < 2.1, divided by B = 0.206:
            what ATLAS's sigma = N_unfolded / (L B dpT dy) measures
  atlasVis  the atlas sample at the visible-jet pT (B method alone)

Output: JSON with sum w and sum w^2 per 10 GeV bin (60-300 GeV) and definition, in mb.
usage: ~/Programs/pythia8venv/bin/python -I compareBJetDefinitions_pythia.py <nEvents> <seed> <out.json>
"""
import sys, math, json
import pythia8mc as p8

nev, seed, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
py = p8.Pythia("", False)
for s in ["Beams:eCM = 5020.", "HardQCD:all = on", "PhaseSpace:pTHatMin = 40.",
          "PhaseSpace:bias2Selection = on", "PhaseSpace:bias2SelectionPow = 4.", "PhaseSpace:bias2SelectionRef = 100.",
          "Random:setSeed = on", f"Random:seed = {seed}", "Next:numberCount = 0"]:
    py.readString(s)
py.init()

PTMIN, LO, BW, NB = 50., 60., 10., 24
sjVis = p8.SlowJet(-1, 0.4, PTMIN, 3.0, 2, 2)   # visible: no neutrinos
sjAll = p8.SlowJet(-1, 0.4, PTMIN, 3.0, 1, 2)   # every final particle
defs = ["ours", "visY21", "allY21", "atlas", "atlasVis", "inclOurs", "inclY21"]
H  = {d: [0.]*NB for d in defs}
H2 = {d: [0.]*NB for d in defs}
def fill(d, pt, w):
    b = int((pt - LO)//BW)
    if 0 <= b < NB: H[d][b] += w; H2[d][b] += w*w

def hasB(i):
    a = abs(i)
    return (a//10) % 10 == 5 or (a//100) % 10 == 5 or (a//1000) % 10 == 5
def dR(y1, p1, y2, p2):
    dp = abs(p1 - p2)
    if dp > math.pi: dp = 2.*math.pi - dp
    return math.hypot(y1 - y2, dp)

sumw = 0.
for _ in range(nev):
    if not py.next(): continue
    ev = py.event; w = py.infoPython().weight(); sumw += w
    sjVis.analyze(ev); sjAll.analyze(ev)
    if sjVis.sizeJet() == 0 and sjAll.sizeJet() == 0: continue
    bh, mu = [], []
    for i in range(ev.size()):
        p = ev[i]
        if p.isHadron() and hasB(p.id()):
            last = True
            for j in range(p.daughter1(), p.daughter2() + 1):
                if j > 0 and ev[j].isHadron() and hasB(ev[j].id()): last = False; break
            if last and p.pT() > 5.: bh.append((p.y(), p.phi()))
        elif p.isFinal() and abs(p.id()) == 13:
            k, fromB = i, False
            while True:
                m = ev[k].mother1()
                if m <= 0: break
                q = ev[m]
                if q.isHadron():
                    if hasB(q.id()): fromB = True; break
                elif abs(q.id()) != 15: break
                k = m
            if fromB: mu.append((p.y(), p.phi()))
    # visible jets
    visB = []
    for j in range(sjVis.sizeJet()):
        pj = sjVis.p(j); pt, y, eta, phi = pj.pT(), pj.rap(), pj.eta(), pj.phi()
        if abs(eta) < 1.6: fill("inclOurs", pt, w)
        if abs(y) < 2.1: fill("inclY21", pt, w)
        if not any(dR(y, phi, b[0], b[1]) < 0.3 for b in bh): continue
        if abs(eta) < 1.6: fill("ours", pt, w)
        if abs(y) < 2.1: fill("visY21", pt, w)
        visB.append((y, phi, pt))
    for j in range(sjAll.sizeJet()):
        pj = sjAll.p(j); pt, y, phi = pj.pT(), pj.rap(), pj.phi()
        if abs(y) > 2.1: continue
        if not any(dR(y, phi, b[0], b[1]) < 0.3 for b in bh): continue
        fill("allY21", pt, w)
        if any(dR(y, phi, m[0], m[1]) < 0.4 for m in mu):
            fill("atlas", pt, w/0.206)
            best = min(visB, key=lambda v: dR(y, phi, v[0], v[1]), default=None)
            if best and dR(y, phi, best[0], best[1]) < 0.2: fill("atlasVis", best[2], w/0.206)

py.stat()
sig = py.infoPython().sigmaGen()   # mb, for sum of weights sumw
json.dump({"sigmaGen": sig, "sumw": sumw, "lo": LO, "bw": BW, "H": H, "H2": H2}, open(out, "w"))
