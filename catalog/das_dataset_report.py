#!/usr/bin/env python3
"""Summarize a DAS dataset and optionally cross-check its runs against a
CMS golden JSON (run-level only, not lumi-precise -- use brilcalc/filterJSON
for an exact lumi mask).

Requires dasgoclient on PATH and a valid CMS grid proxy
(voms-proxy-init --voms cms); run on lxplus/cmsconnect, not locally.
"""
import argparse
import json
import subprocess
import sys


def das_query(query: str) -> list[str]:
    try:
        result = subprocess.run(
            ["dasgoclient", "--query", query],
            capture_output=True,
            text=True,
            timeout=120,
        )
    except FileNotFoundError:
        print(
            "dasgoclient not found. Run this on lxplus/cmsconnect with a CMS "
            "environment sourced and voms-proxy-init --voms cms already done.",
            file=sys.stderr,
        )
        sys.exit(1)
    if result.returncode != 0:
        print(f"dasgoclient failed for: {query}\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    return [line for line in result.stdout.splitlines() if line.strip()]


def dataset_summary(dataset: str):
    lines = das_query(f"dataset={dataset} | grep dataset.nevents, dataset.nfiles, dataset.size")
    if not lines:
        return None
    nevents, nfiles, size = lines[0].split()
    return int(nevents), int(nfiles), int(size)


def dataset_runs(dataset: str):
    return sorted(int(r) for r in das_query(f"run dataset={dataset}"))


def load_golden_runs(golden_json_path: str):
    with open(golden_json_path) as f:
        golden = json.load(f)
    return {int(run) for run in golden}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", help="DAS dataset path, e.g. /HIHardProbes/HIRun2018A-04Apr2019-v1/AOD")
    parser.add_argument("--golden-json", help="path to a golden JSON, e.g. certificationFiles/Cert_*.txt")
    args = parser.parse_args()

    summary = dataset_summary(args.dataset)
    if summary is None:
        print(f"No dataset found: {args.dataset}")
        sys.exit(1)
    nevents, nfiles, size = summary
    print(args.dataset)
    print(f"  events: {nevents:,}   files: {nfiles:,}   size: {size / 1e9:.1f} GB")

    runs = dataset_runs(args.dataset)
    if runs:
        print(f"  runs: {len(runs)} spanning {runs[0]}-{runs[-1]}")

    if args.golden_json:
        golden_runs = load_golden_runs(args.golden_json)
        certified = [r for r in runs if r in golden_runs]
        print(f"  certified (run-level) per {args.golden_json}: {len(certified)}/{len(runs)}")
        missing = sorted(set(runs) - golden_runs)
        if missing:
            shown = missing[:20]
            suffix = " ..." if len(missing) > 20 else ""
            print(f"  runs NOT in golden JSON: {shown}{suffix}")


if __name__ == "__main__":
    main()
