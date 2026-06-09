#!/usr/bin/env python3
"""
summary_report.py — Aggregate all per-sample pipeline outputs into a single
summary TSV.

Usage:
    summary_report.py --qc-dir <dir> --coverage-dir <dir> --consensus-dir <dir>
                      --nextclade-dir <dir> --vadr-dir <dir> --kraken2-dir <dir>
                      --output <summary_report.tsv>

All directories can be '.' when all files are staged flat in the work directory.

Output columns:
    sample_id, reference, start, end,
    num_mapped_reads, cov_bases_mapped, percent_genome_cov_map,
    mean_depth, mean_base_qual, mean_map_qual,
    assembly_length, numN, percent_ref_genome_cov,
    QC_flag, VADR_flag,
    nextclade_clade, nextclade_qc_overall,
    nextclade_total_substitutions, nextclade_total_deletions,
    nextclade_total_insertions, nextclade_total_ns,
    kraken2_dengue_pct
"""

import argparse
import csv
import glob
import os
import sys


# ---------------------------------------------------------------------------
# Loaders — each returns a dict keyed by sample_id
# ---------------------------------------------------------------------------

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


def load_coverage(coverage_dir):
    """
    Read *.coverage.txt files (samtools coverage output).
    samtools coverage columns:
        #rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
    Returns {sample_id: dict_of_stats}.
    """
    records = {}
    for path in glob.glob(os.path.join(coverage_dir, "*.coverage.txt")):
        sid = os.path.basename(path).replace(".coverage.txt", "")
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


def load_consensus(consensus_dir):
    """
    Read *.consensus.fa files.
    Returns {sample_id: {assembly_length, numN, percent_ref_genome_cov}}.
    percent_ref_genome_cov requires reference length — passed in from coverage.
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
            "_seq_called":     len(seq) - seq.count("N"),  # used to compute pct
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Aggregate pipeline results into summary TSV")
    parser.add_argument("--qc-dir",        required=True)
    parser.add_argument("--coverage-dir",  required=True)
    parser.add_argument("--consensus-dir", required=True)
    parser.add_argument("--nextclade-dir", required=True)
    parser.add_argument("--vadr-dir",      required=True)
    parser.add_argument("--kraken2-dir",   required=True)
    parser.add_argument("--output",        required=True)
    args = parser.parse_args()

    qc        = load_qc(args.qc_dir)
    coverage  = load_coverage(args.coverage_dir)
    consensus = load_consensus(args.consensus_dir)
    nextclade = load_nextclade(args.nextclade_dir)
    vadr      = load_vadr(args.vadr_dir)
    kraken2   = load_kraken2(args.kraken2_dir)

    all_samples = sorted(set(list(qc.keys()) + list(coverage.keys()) + list(nextclade.keys())))

    if not all_samples:
        print("WARNING: no samples found — check input directories", file=sys.stderr)
        sys.exit(1)

    header = [
        "sample_id", "reference", "start", "end",
        "num_mapped_reads", "cov_bases_mapped", "percent_genome_cov_map",
        "mean_depth", "mean_base_qual", "mean_map_qual",
        "assembly_length", "numN", "percent_ref_genome_cov",
        "QC_flag", "VADR_flag",
        "nextclade_clade", "nextclade_qc_overall",
        "nextclade_total_substitutions", "nextclade_total_deletions",
        "nextclade_total_insertions", "nextclade_total_ns",
        "kraken2_dengue_pct",
    ]

    rows = []
    for sid in all_samples:
        cov = coverage.get(sid, {})
        con = consensus.get(sid, {})
        nc  = nextclade.get(sid, {})
        vf  = vadr.get(sid, "NA")
        k2  = kraken2.get(sid, "NA")
        qf  = qc.get(sid, "NA")

        ref_len = int(cov.get("end", 0) or 0)
        called  = con.get("_seq_called", 0)
        pct_ref = f"{(called / ref_len * 100):.4f}" if ref_len > 0 else "NA"

        row = {
            "sample_id":                     sid,
            "reference":                     cov.get("reference", "NA"),
            "start":                         cov.get("start", "NA"),
            "end":                           cov.get("end", "NA"),
            "num_mapped_reads":              cov.get("num_mapped_reads", "NA"),
            "cov_bases_mapped":              cov.get("cov_bases_mapped", "NA"),
            "percent_genome_cov_map":        cov.get("percent_genome_cov_map", "NA"),
            "mean_depth":                    cov.get("mean_depth", "NA"),
            "mean_base_qual":                cov.get("mean_base_qual", "NA"),
            "mean_map_qual":                 cov.get("mean_map_qual", "NA"),
            "assembly_length":               con.get("assembly_length", "NA"),
            "numN":                          con.get("numN", "NA"),
            "percent_ref_genome_cov":        pct_ref,
            "QC_flag":                       qf,
            "VADR_flag":                     vf,
            "nextclade_clade":               nc.get("clade", "NA"),
            "nextclade_qc_overall":          nc.get("qcOverallStatus", "NA"),
            "nextclade_total_substitutions": nc.get("totalSubstitutions", "NA"),
            "nextclade_total_deletions":     nc.get("totalDeletions", "NA"),
            "nextclade_total_insertions":    nc.get("totalInsertions", "NA"),
            "nextclade_total_ns":            nc.get("totalNs", "NA"),
            "kraken2_dengue_pct":            k2,
        }
        rows.append(row)

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Summary report written: {args.output} ({len(rows)} samples)", file=sys.stderr)


if __name__ == "__main__":
    main()
