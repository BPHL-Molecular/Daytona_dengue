#!/usr/bin/env python3
"""
serotype_detect.py — Determine dengue serotype by comparing BWA alignment
coverage across all four serotype references.

Usage:
    serotype_detect.py --denv1 DENV1.coverage.txt --denv2 ... --denv3 ... --denv4 ...
                       --sample-id <id> --output <serotype.txt>
                       [--detail <detail.tsv>] [--min-coverage 10]

Picks the reference with the highest coverage (% genome bases covered).
Returns 'unclassified' if the best coverage is below --min-coverage.
"""

import argparse
import sys


def parse_samtools_coverage(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 7:
                continue
            try:
                return {
                    'rname':     cols[0],
                    'numreads':  int(cols[3]),
                    'covbases':  int(cols[4]),
                    'coverage':  float(cols[5]),
                    'meandepth': float(cols[6]),
                }
            except (ValueError, IndexError):
                pass
    return {'rname': 'NA', 'numreads': 0, 'covbases': 0, 'coverage': 0.0, 'meandepth': 0.0}


def main():
    parser = argparse.ArgumentParser(
        description='Determine DENV serotype from BWA alignment coverage across all 4 references')
    parser.add_argument('--denv1',        required=True)
    parser.add_argument('--denv2',        required=True)
    parser.add_argument('--denv3',        required=True)
    parser.add_argument('--denv4',        required=True)
    parser.add_argument('--sample-id',    required=True)
    parser.add_argument('--output',       required=True)
    parser.add_argument('--detail',       required=False)
    parser.add_argument('--min-coverage', type=float, default=10.0,
                        help='Minimum %% genome coverage to call a serotype (default: 10.0)')
    args = parser.parse_args()

    candidates = {
        'DENV1': parse_samtools_coverage(args.denv1),
        'DENV2': parse_samtools_coverage(args.denv2),
        'DENV3': parse_samtools_coverage(args.denv3),
        'DENV4': parse_samtools_coverage(args.denv4),
    }

    best_serotype = max(candidates, key=lambda k: candidates[k]['coverage'])
    best_cov      = candidates[best_serotype]['coverage']
    result        = best_serotype if best_cov >= args.min_coverage else 'unclassified'

    with open(args.output, 'w') as fh:
        fh.write(result + '\n')

    if args.detail:
        with open(args.detail, 'w') as fh:
            fh.write('sample_id\tserotype\treference\tnumreads\tcovbases\tcoverage_pct\tmeandepth\tselected\n')
            for sero, stats in candidates.items():
                selected = 'YES' if (sero == best_serotype and result != 'unclassified') else ''
                fh.write(
                    f"{args.sample_id}\t{sero}\t{stats['rname']}\t"
                    f"{stats['numreads']}\t{stats['covbases']}\t"
                    f"{stats['coverage']:.4f}\t{stats['meandepth']:.2f}\t{selected}\n"
                )

    print(
        f"Serotype for {args.sample_id}: {result} "
        f"(best: {best_serotype} at {best_cov:.2f}% genome coverage)",
        file=sys.stderr
    )
    if result == 'unclassified':
        print(
            f"WARNING: {args.sample_id} did not reach minimum coverage threshold "
            f"({args.min_coverage}%) on any reference -- excluded from pipeline",
            file=sys.stderr
        )


if __name__ == '__main__':
    main()
