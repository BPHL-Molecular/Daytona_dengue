#!/usr/bin/env python3
"""
summary_report.py — Aggregate all per-sample pipeline outputs into a single
summary TXT (tab-separated).

Usage:
    summary_report.py --qc-dir <dir> --coverage-dir <dir> --consensus-dir <dir>
                      --nextclade-dir <dir> --vadr-dir <dir> --kraken2-dir <dir>
                      --serotype-dir <dir> --screen-cov-dir <dir>
                      --trimstat-dir <dir> --phix-log-dir <dir>
                      --output <summary_report.txt>

All directories can be '.' when all files are staged flat in the work directory.

Output columns:
    sample_id, serotype, nextclade_clade, nextclade_qc_overall, kraken2_percent,
    reference, start, end,
    num_raw_reads, num_clean_reads, num_mapped_reads, percent_mapped_clean_reads,
    cov_bases_mapped, percent_genome_cov_map, mean_depth, mean_base_qual, mean_map_qual,
    assembly_length, numN, percent_ref_genome_cov,
    VADR_flag, QC_flag
"""

import argparse
import csv
import glob
import os
import re
import sys


# ---------------------------------------------------------------------------
# Loaders — each returns a dict keyed by sample_id
# ---------------------------------------------------------------------------

def load_serotype(serotype_dir):
    """Read *_serotype.txt → {sample_id: serotype_str}."""
    records = {}
    for path in glob.glob(os.path.join(serotype_dir, "*_serotype.txt")):
        sid = os.path.basename(path).replace("_serotype.txt", "")
        with open(path) as fh:
            value = fh.read().strip()
        if sid and value:
            records[sid] = value
    return records


def load_coverage(coverage_dir):
    """
    Read *.coverage.txt files (from samtools_coverage, post-ivar).
    Skips *_DENV?.coverage.txt files (those are from samtools_screen).
    Returns {sample_id: dict_of_stats}.
    """
    records = {}
    for path in glob.glob(os.path.join(coverage_dir, "*.coverage.txt")):
        fname = os.path.basename(path)
        if re.search(r'_DENV\d\.coverage\.txt$', fname):
            continue
        sid = fname.replace(".coverage.txt", "")
        with open(path) as fh:
            _header = fh.readline()
            line = fh.readline().rstrip()
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        records[sid] = {
            "reference":              cols[0],
            "start":                  cols[1],
            "end":                    cols[2],
            "num_mapped_reads":       cols[3],
            "cov_bases_mapped":       cols[4],
            "percent_genome_cov_map": cols[5],
            "mean_depth":             cols[6],
            "mean_base_qual":         cols[7],
            "mean_map_qual":          cols[8],
        }
    return records


def load_screen_coverage(screen_cov_dir):
    """
    Read *_DENV?.coverage.txt files (from samtools_screen).
    For each sample, pick the DENV reference with the highest mean depth.
    Returns {sample_id: dict_of_stats} — used as fallback for unclassified samples.
    """
    from collections import defaultdict
    sample_files = defaultdict(dict)
    for path in glob.glob(os.path.join(screen_cov_dir, "*_DENV?.coverage.txt")):
        fname = os.path.basename(path)
        m = re.match(r'(.+)_(DENV\d)\.coverage\.txt$', fname)
        if not m:
            continue
        sid, denv = m.group(1), m.group(2)
        sample_files[sid][denv] = path

    records = {}
    for sid, denv_files in sample_files.items():
        best_cols = None
        best_depth = -1.0
        for _denv, path in sorted(denv_files.items()):
            with open(path) as fh:
                fh.readline()
                line = fh.readline().rstrip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            try:
                depth = float(cols[6])
            except ValueError:
                depth = 0.0
            if depth > best_depth:
                best_depth = depth
                best_cols = cols
        if best_cols:
            records[sid] = {
                "reference":              best_cols[0],
                "start":                  best_cols[1],
                "end":                    best_cols[2],
                "num_mapped_reads":       best_cols[3],
                "cov_bases_mapped":       best_cols[4],
                "percent_genome_cov_map": best_cols[5],
                "mean_depth":             best_cols[6],
                "mean_base_qual":         best_cols[7],
                "mean_map_qual":          best_cols[8],
            }
    return records


def load_consensus(consensus_dir):
    """
    Read *.consensus.fa files.
    Returns {sample_id: {assembly_length, numN, _seq_called}}.
    """
    records = {}
    for path in glob.glob(os.path.join(consensus_dir, "*.consensus.fa")):
        sid = os.path.basename(path).replace(".consensus.fa", "")
        seq_parts = []
        with open(path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    seq_parts.append(line.rstrip())
        seq = "".join(seq_parts).upper()
        records[sid] = {
            "assembly_length": str(len(seq)),
            "numN":            str(seq.count("N")),
            "_seq_called":     len(seq) - seq.count("N"),
        }
    return records


def load_nextclade(nextclade_dir):
    """Read *_nextclade.tsv files → {seqName: row}."""
    records = {}
    for path in glob.glob(os.path.join(nextclade_dir, "*_nextclade.tsv")):
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = row.get("seqName", "").strip()
                if sid:
                    records[sid] = row
    return records


def load_vadr(vadr_dir):
    """Walk VADR result dirs → {sample_id: 'PASS'|'REVIEW'}."""
    records = {}
    for pass_list in glob.glob(os.path.join(vadr_dir, "**", "*.vadr.pass.list"), recursive=True):
        with open(pass_list) as fh:
            for line in fh:
                sid = line.strip()
                if sid:
                    records[sid] = "PASS"
    for fail_list in glob.glob(os.path.join(vadr_dir, "**", "*.vadr.fail.list"), recursive=True):
        with open(fail_list) as fh:
            for line in fh:
                sid = line.strip()
                if sid and sid not in records:
                    records[sid] = "REVIEW"
    return records


def load_kraken2(kraken2_dir):
    """
    Parse *_kraken2_report.txt files, sum % reads across dengue serotype taxa.
    Returns {sample_id: '42.35'}.
    """
    DENGUE_KEYWORDS = (
        "dengue virus 1", "dengue virus 2", "dengue virus 3", "dengue virus 4",
        "dengue virus type 1", "dengue virus type 2",
        "dengue virus type 3", "dengue virus type 4",
    )
    records = {}
    for path in glob.glob(os.path.join(kraken2_dir, "*_kraken2_report.txt")):
        sid = os.path.basename(path).replace("_kraken2_report.txt", "")
        total_pct = 0.0
        try:
            with open(path) as fh:
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 6:
                        continue
                    try:
                        pct = float(cols[0].strip())
                    except ValueError:
                        continue
                    name = cols[5].strip().lower()
                    for kw in DENGUE_KEYWORDS:
                        if kw in name:
                            total_pct += pct
                            break
        except OSError:
            pass
        records[sid] = f"{round(total_pct, 2):.2f}"
    return records


def load_trimstats(trimstat_dir):
    """
    Parse *_trimstats.txt (trimmomatic PE stdout/stderr).
    Extracts 'Input Read Pairs: N' → num_raw_reads = N * 2.
    Returns {sample_id: str(n)}.
    """
    records = {}
    for path in glob.glob(os.path.join(trimstat_dir, "*_trimstats.txt")):
        sid = os.path.basename(path).replace("_trimstats.txt", "")
        try:
            with open(path) as fh:
                content = fh.read()
            m = re.search(r'Input Read Pairs:\s+(\d+)', content)
            if m:
                records[sid] = str(int(m.group(1)) * 2)
        except OSError:
            pass
    return records


def load_phix_log(phix_log_dir):
    """
    Parse *_phix_log.txt (bbduk stderr).
    Extracts 'Result:  N reads' → num_clean_reads.
    Returns {sample_id: str(n)}.
    """
    records = {}
    for path in glob.glob(os.path.join(phix_log_dir, "*_phix_log.txt")):
        sid = os.path.basename(path).replace("_phix_log.txt", "")
        try:
            with open(path) as fh:
                content = fh.read()
            m = re.search(r'Result:\s+(\d+)\s+reads', content)
            if m:
                records[sid] = str(int(m.group(1)))
        except OSError:
            pass
    return records


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Aggregate pipeline results into summary TXT")
    parser.add_argument("--qc-dir",         required=True)
    parser.add_argument("--coverage-dir",   required=True)
    parser.add_argument("--consensus-dir",  required=True)
    parser.add_argument("--nextclade-dir",  required=True)
    parser.add_argument("--vadr-dir",       required=True)
    parser.add_argument("--kraken2-dir",    required=True)
    parser.add_argument("--serotype-dir",   required=True)
    parser.add_argument("--screen-cov-dir", required=True)
    parser.add_argument("--trimstat-dir",   required=True)
    parser.add_argument("--phix-log-dir",   required=True)
    parser.add_argument("--output",         required=True)
    args = parser.parse_args()

    serotype       = load_serotype(args.serotype_dir)
    qc             = load_qc(args.qc_dir)
    coverage       = load_coverage(args.coverage_dir)
    screen_cov     = load_screen_coverage(args.screen_cov_dir)
    consensus      = load_consensus(args.consensus_dir)
    nextclade      = load_nextclade(args.nextclade_dir)
    vadr           = load_vadr(args.vadr_dir)
    kraken2        = load_kraken2(args.kraken2_dir)
    trimstats      = load_trimstats(args.trimstat_dir)
    phix_log       = load_phix_log(args.phix_log_dir)

    all_samples = sorted(serotype.keys())

    if not all_samples:
        print("WARNING: no samples found in serotype directory — check input", file=sys.stderr)
        sys.exit(1)

    header = [
        "sample_id", "serotype",
        "nextclade_clade", "nextclade_qc_overall",
        "kraken2_percent",
        "reference", "start", "end",
        "num_raw_reads", "num_clean_reads", "num_mapped_reads", "percent_mapped_clean_reads",
        "cov_bases_mapped", "percent_genome_cov_map",
        "mean_depth", "mean_base_qual", "mean_map_qual",
        "assembly_length", "numN", "percent_ref_genome_cov",
        "VADR_flag", "QC_flag",
    ]

    rows = []
    for sid in all_samples:
        sero = serotype.get(sid, "NA")
        unclassified = (sero == "unclassified")

        cov = coverage.get(sid, screen_cov.get(sid, {})) if not unclassified else screen_cov.get(sid, {})
        con = {} if unclassified else consensus.get(sid, {})
        nc  = {} if unclassified else nextclade.get(sid, {})
        vf  = "NA" if unclassified else vadr.get(sid, "NA")
        k2  = kraken2.get(sid, "NA")
        qf  = "FAIL: Unclassified" if unclassified else qc.get(sid, "NA")

        raw_reads   = trimstats.get(sid, "NA")
        clean_reads = phix_log.get(sid, "NA")

        mapped = cov.get("num_mapped_reads", "NA")
        if mapped != "NA" and clean_reads != "NA":
            try:
                pct_mapped = f"{float(mapped) / float(clean_reads) * 100:.2f}"
            except (ValueError, ZeroDivisionError):
                pct_mapped = "NA"
        else:
            pct_mapped = "NA"

        ref_len = int(cov.get("end", 0) or 0)
        called  = con.get("_seq_called", 0)
        pct_ref = f"{(called / ref_len * 100):.4f}" if ref_len > 0 and not unclassified else "NA"

        row = {
            "sample_id":                     sid,
            "serotype":                      sero,
            "nextclade_clade":               nc.get("clade", "NA"),
            "nextclade_qc_overall":          nc.get("qc.overallStatus", "NA"),
            "kraken2_percent":               k2,
            "reference":                     cov.get("reference", "NA"),
            "start":                         cov.get("start", "NA"),
            "end":                           cov.get("end", "NA"),
            "num_raw_reads":                 raw_reads,
            "num_clean_reads":               clean_reads,
            "num_mapped_reads":              mapped,
            "percent_mapped_clean_reads":    pct_mapped,
            "cov_bases_mapped":              cov.get("cov_bases_mapped", "NA"),
            "percent_genome_cov_map":        cov.get("percent_genome_cov_map", "NA"),
            "mean_depth":                    cov.get("mean_depth", "NA"),
            "mean_base_qual":                cov.get("mean_base_qual", "NA"),
            "mean_map_qual":                 cov.get("mean_map_qual", "NA"),
            "assembly_length":               con.get("assembly_length", "NA"),
            "numN":                          con.get("numN", "NA"),
            "percent_ref_genome_cov":        pct_ref,
            "VADR_flag":                     vf,
            "QC_flag":                       qf,
        }
        rows.append(row)

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Summary report written: {args.output} ({len(rows)} samples)", file=sys.stderr)


def load_qc(qc_dir):
    """Read *_qc.tsv files → {sample_id: qc_flag}."""
    records = {}
    for path in glob.glob(os.path.join(qc_dir, "*_qc.tsv")):
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = row.get("sample_id", "").strip()
                if sid:
                    records[sid] = row.get("qc_flag", "NA").strip()
    return records


if __name__ == "__main__":
    main()
