#!/usr/bin/env python3
"""Build/refresh a SQLite catalog of rootFiles/scanningOuput/**.

Filenames in that tree encode dataset, cuts, and a production date as
underscore-delimited tokens (e.g.
PbPb_MinBias_Part1_..._pseudoJetCandPtMin-0.0_..._2026-6-30.root).
This walks the tree, parses each name into a "config key" (everything but
the date) plus individual tags, and records a couple of QA numbers pulled
out of each file's histograms so the result is queryable with SQL instead
of by eye.
"""
import argparse
import re
import sqlite3
import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCAN_ROOT = REPO_ROOT / "rootFiles" / "scanningOuput"
MACRO_PATH = Path(__file__).resolve().parent / "inspect_scan_file.C"
DB_PATH = Path(__file__).resolve().parent / "scan_catalog.db"

DATE_RE = re.compile(r"^\d{4}-\d{1,2}-\d{1,2}$")

SCHEMA = """
CREATE TABLE IF NOT EXISTS scan_files (
    id INTEGER PRIMARY KEY,
    path TEXT UNIQUE NOT NULL,
    dataset_group TEXT NOT NULL,
    filename TEXT NOT NULL,
    config_key TEXT NOT NULL,
    produced_on TEXT,
    size_bytes INTEGER,
    mtime TEXT,
    n_keys INTEGER,
    events_before_selection REAL,
    events_after_selection REAL
);
CREATE TABLE IF NOT EXISTS scan_file_tags (
    file_id INTEGER NOT NULL REFERENCES scan_files(id),
    tag TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_tags_tag ON scan_file_tags(tag);
CREATE INDEX IF NOT EXISTS idx_tags_file ON scan_file_tags(file_id);
"""


def parse_filename(name: str):
    stem = name[: -len(".root")]
    tokens = stem.split("_")
    dataset_group = tokens[0]
    date_idx = next((i for i, t in enumerate(tokens) if DATE_RE.match(t)), None)
    if date_idx is not None:
        y, m, d = tokens[date_idx].split("-")
        produced_on = f"{int(y):04d}-{int(m):02d}-{int(d):02d}"
        tags = tokens[:date_idx] + tokens[date_idx + 1 :]
    else:
        produced_on = None
        tags = tokens[:]
    config_key = "_".join(tags)
    return dataset_group, produced_on, tags, config_key


def inspect_root_file(path: Path):
    """Run the ROOT macro to count histogram keys and pull selection yields."""
    try:
        result = subprocess.run(
            ["root", "-l", "-b", "-q", f'{MACRO_PATH}("{path}")'],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None, None, None
    n_keys = before = after = None
    for line in result.stdout.splitlines():
        parts = line.split("|")
        if parts[0] == "NKEYS":
            n_keys = int(parts[1])
        elif parts[0] == "HIST" and parts[1] == "h_eventsBeforeSelection":
            before = float(parts[2])
        elif parts[0] == "HIST" and parts[1] == "h_eventsAfterSelection":
            after = float(parts[2])
    return n_keys, before, after


def build(with_root: bool):
    files = sorted(SCAN_ROOT.rglob("*.root"))
    conn = sqlite3.connect(DB_PATH)
    conn.executescript(SCHEMA)

    for path in files:
        rel = str(path.relative_to(REPO_ROOT))
        dataset_group, produced_on, tags, config_key = parse_filename(path.name)
        stat = path.stat()
        n_keys = before = after = None
        if with_root:
            n_keys, before, after = inspect_root_file(path)

        row = conn.execute(
            """
            INSERT INTO scan_files
                (path, dataset_group, filename, config_key, produced_on,
                 size_bytes, mtime, n_keys, events_before_selection, events_after_selection)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(path) DO UPDATE SET
                config_key = excluded.config_key,
                produced_on = excluded.produced_on,
                size_bytes = excluded.size_bytes,
                mtime = excluded.mtime,
                n_keys = COALESCE(excluded.n_keys, scan_files.n_keys),
                events_before_selection = COALESCE(excluded.events_before_selection, scan_files.events_before_selection),
                events_after_selection = COALESCE(excluded.events_after_selection, scan_files.events_after_selection)
            RETURNING id
            """,
            (
                rel,
                dataset_group,
                path.name,
                config_key,
                produced_on,
                stat.st_size,
                datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                n_keys,
                before,
                after,
            ),
        )
        file_id = row.fetchone()[0]
        conn.execute("DELETE FROM scan_file_tags WHERE file_id = ?", (file_id,))
        conn.executemany(
            "INSERT INTO scan_file_tags (file_id, tag) VALUES (?, ?)",
            [(file_id, t) for t in tags],
        )

    conn.commit()
    n_total = conn.execute("SELECT COUNT(*) FROM scan_files").fetchone()[0]
    conn.close()
    print(f"Catalogued {len(files)} files this run ({n_total} total) into {DB_PATH}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-root",
        action="store_true",
        help="Skip opening ROOT files for histogram QA metadata (faster; path/date/tags only).",
    )
    args = parser.parse_args()
    build(with_root=not args.no_root)
