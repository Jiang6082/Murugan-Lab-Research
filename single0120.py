#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from collections import Counter
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats
from scipy.optimize import curve_fit
from scipy.special import factorial

# Fixed index pair requested
FIXED_FWD_INDEX = "GGACTGTG"
FIXED_REV_INDEX = "CATTCCAG"


def parse_arguments():
    """
    Keeps the same core positional args as your original script:
      overall_title wt_sequence file1 file2

    Also accepts (and ignores) any extra trailing args so you can keep your old
    command structure without breaking.
    """
    parser = argparse.ArgumentParser(
        description="Single index-pair mutation stats + log histogram (exclude 0-mutation from stats/plot)."
    )
    parser.add_argument("run_title", help="Short title for this run (used in plot titles and output files)")
    parser.add_argument("wt_sequence", help="Wild-type reference sequence (already stripped if that’s how you use it)")
    parser.add_argument("file1", help="Forward reads FASTQ/TXT")
    parser.add_argument("file2", help="Reverse reads FASTQ/TXT")

    # Accept extra args (conditions etc.) but ignore them
    parser.add_argument("condition_data", nargs="*", help="Extra args (ignored). Kept for compatibility.")

    parser.add_argument("--max-hamming", type=int, default=1,
                        help="Maximum hamming distance for index matching (default: 1)")
    parser.add_argument("--strip-length", type=int, default=8,
                        help="Number of bases to strip from each end after matching (default: 8)")

    parser.add_argument("--max-mutations-plot", type=int, default=30,
                        help="Max mutation count shown on x-axis (default: 30)")
    parser.add_argument("--align-threshold", type=int, default=10,
                        help="If initial mutation calls exceed this, do pairwise alignment (default: 10)")

    return parser.parse_args()


def reverse_complement(seq: str) -> str:
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(base, base) for base in seq[::-1])


def hamming_distance(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        return 10**9
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def read_sequences(filename: str) -> list[str]:
    sequences = []
    try:
        with open(filename, "r") as f:
            lines = f.readlines()
        # FASTQ
        if lines and lines[0].startswith("@"):
            for i in range(1, len(lines), 4):
                if i < len(lines):
                    sequences.append(lines[i].strip())
        else:
            # TXT: one sequence per line
            sequences = [line.strip() for line in lines if line.strip()]
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        sys.exit(1)
    return sequences


def check_sequence_match(seq: str, seq1: str, seq2_rc: str, max_hamming: int) -> tuple[bool, bool]:
    """Matches seq1 at start, seq2_rc at end"""
    seq1_match = False
    seq2_match = False

    if len(seq) >= len(seq1):
        start_seq = seq[:len(seq1)]
        if hamming_distance(start_seq, seq1) <= max_hamming:
            seq1_match = True

    if len(seq) >= len(seq2_rc):
        end_seq = seq[-len(seq2_rc):]
        if hamming_distance(end_seq, seq2_rc) <= max_hamming:
            seq2_match = True

    return seq1_match, seq2_match


def identify_mutations_simple(seq: str, wt_seq: str) -> list[str]:
    muts = []
    min_len = min(len(seq), len(wt_seq))
    for i in range(min_len):
        if seq[i] != wt_seq[i]:
            muts.append(f"{wt_seq[i]}{i+1}{seq[i]}")
    if len(seq) != len(wt_seq):
        muts.append(f"LengthDiff_{len(seq)}vs{len(wt_seq)}")
    return muts


def pairwise_align(seq1: str, seq2: str, match_score=1, mismatch_score=-1, gap_score=-1) -> tuple[str, str]:
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score = [[0] * cols for _ in range(rows)]

    for i in range(1, rows):
        score[i][0] = gap_score * i
    for j in range(1, cols):
        score[0][j] = gap_score * j

    for i in range(1, rows):
        for j in range(1, cols):
            diag = score[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            up = score[i - 1][j] + gap_score
            left = score[i][j - 1] + gap_score
            score[i][j] = max(diag, up, left)

    aligned1 = []
    aligned2 = []
    i, j = rows - 1, cols - 1

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            cur = score[i][j]
            diag = score[i - 1][j - 1]
            up = score[i - 1][j]
            left = score[i][j - 1]
            if cur == diag + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
                aligned1.append(seq1[i - 1])
                aligned2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif cur == up + gap_score:
                aligned1.append(seq1[i - 1])
                aligned2.append("-")
                i -= 1
            else:
                aligned1.append("-")
                aligned2.append(seq2[j - 1])
                j -= 1
        elif i > 0:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()
    return "".join(aligned1), "".join(aligned2)


def identify_mutations_with_alignment(seq: str, wt_seq: str, align_threshold: int) -> tuple[list[str], bool]:
    initial = identify_mutations_simple(seq, wt_seq)
    if len(initial) <= align_threshold:
        return initial, False

    aligned_seq, aligned_wt = pairwise_align(seq, wt_seq)
    muts = []
    wt_pos = 0

    for s, w in zip(aligned_seq, aligned_wt):
        if w != "-":
            wt_pos += 1
            if s == "-":
                muts.append(f"del{wt_pos}{w}")
            elif s != w:
                muts.append(f"{w}{wt_pos}{s}")
        elif s != "-":
            muts.append(f"ins{wt_pos}{s}")

    muts.append("ALIGNED")
    return muts, True


def filter_actual_mutations(mutations: list[str]) -> list[str]:
    """Remove flags/markers so mutation *counts* reflect real base indels/subs only."""
    return [m for m in mutations if m != "ALIGNED" and not m.startswith("LengthDiff")]


def poisson_function(x, mu):
    return np.exp(-mu) * (mu ** x) / factorial(x)


def compute_stats(vals: list[int]) -> dict:
    if not vals:
        return {}
    arr = np.array(vals, dtype=float)
    mode_res = scipy_stats.mode(arr, keepdims=True)
    mode_val = float(mode_res.mode[0]) if mode_res.count.size > 0 else float("nan")

    return {
        "n": int(len(arr)),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "std": float(np.std(arr, ddof=0)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "q1": float(np.percentile(arr, 25)),
        "q3": float(np.percentile(arr, 75)),
        "mode": mode_val,
    }


def main():
    args = parse_arguments()

    print("Reading sequence files...")
    seqs1 = read_sequences(args.file1)
    seqs2 = read_sequences(args.file2)

    if len(seqs1) != len(seqs2):
        print("Error: Files have different number of sequences")
        sys.exit(1)

    print(f"Found {len(seqs1)} sequence pairs")

    fwd = FIXED_FWD_INDEX
    rev = FIXED_REV_INDEX
    rev_rc = reverse_complement(rev)

    print("\nAnalyzing fixed index pair:")
    print(f"  Forward index: {fwd}")
    print(f"  Reverse index: {rev} (RC used for end match: {rev_rc})")
    print(f"  Max hamming distance: {args.max_hamming}")
    print(f"  Strip length: {args.strip_length} bases from each end")
    print(f"  Alignment threshold: {args.align_threshold}\n")

    wt = args.wt_sequence

    matched_total = 0
    zero_mut = 0
    nonzero_mut_counts = []

    aligned_used = 0

    for r1, r2 in zip(seqs1, seqs2):
        # Stitching rule copied from your pipeline
        r2_rc_full = reverse_complement(r2)
        if len(r2_rc_full) >= 69:
            combined = r1 + r2_rc_full[68:]
        else:
            combined = r1 + r2_rc_full

        m1, m2 = check_sequence_match(combined, fwd, rev_rc, args.max_hamming)
        if not (m1 and m2):
            continue

        matched_total += 1

        if len(combined) <= 2 * args.strip_length:
            # Too short after stripping; skip
            continue

        stripped = combined[args.strip_length:-args.strip_length]

        muts, used_align = identify_mutations_with_alignment(stripped, wt, args.align_threshold)
        if used_align:
            aligned_used += 1

        actual_muts = filter_actual_mutations(muts)
        mcount = len(actual_muts)

        if mcount == 0:
            zero_mut += 1
        else:
            nonzero_mut_counts.append(mcount)

    if matched_total == 0:
        print("No reads matched this index pair. Exiting.")
        sys.exit(0)

    zero_pct = 100.0 * zero_mut / matched_total
    stats_dict = compute_stats(nonzero_mut_counts)

    print("=" * 70)
    print("RESULTS (Fixed index pair)")
    print("=" * 70)
    print(f"Matched reads (this index pair): {matched_total}")
    print(f"Zero-mutation reads: {zero_mut} ({zero_pct:.2f}%)")

    print("\nMutation-count statistics (EXCLUDING zero-mutation reads):")
    if not stats_dict:
        print("  No nonzero-mutation reads found after filtering.")
    else:
        print(f"  N (nonzero): {stats_dict['n']}")
        print(f"  Mean:   {stats_dict['mean']:.4f}")
        print(f"  Median: {stats_dict['median']:.4f}")
        print(f"  Std:    {stats_dict['std']:.4f}")
        print(f"  Min:    {stats_dict['min']:.0f}")
        print(f"  Q1:     {stats_dict['q1']:.4f}")
        print(f"  Q3:     {stats_dict['q3']:.4f}")
        print(f"  Max:    {stats_dict['max']:.0f}")
        print(f"  Mode:   {stats_dict['mode']:.0f}")

    if aligned_used > 0:
        print(f"\nAlignment was used for {aligned_used} matched reads (high mutation counts).")

    # -------------------------
    # Plot: log-scale histogram
    # -------------------------
    if nonzero_mut_counts:
        max_x = int(args.max_mutations_plot)
        # Clip counts > max_x into max_x bin (optional; keeps x-range stable)
        clipped = [min(c, max_x) for c in nonzero_mut_counts]

        plt.figure(figsize=(7, 6))
        bins = np.arange(1, max_x + 2)  # bins covering 1..max_x inclusive

        plt.hist(clipped, bins=bins, edgecolor="black", linewidth=0.6, alpha=0.6)
        plt.yscale("log")
        plt.ylim(bottom=0.8)
        plt.xlim(1, max_x)
        
        plt.title(f"{args.run_title}\nMutation Distribution (Fixed indices, 0 excluded)")
        plt.xlabel("Number of mutations")
        plt.ylabel("Number of sequences (log scale)")
        plt.grid(True, alpha=0.3)

        # Optional Poisson overlay (like your example)
        mu = float(np.mean(nonzero_mut_counts))
        x_fit = np.arange(1, max_x + 1)
        y_fit = poisson_function(x_fit, mu) * len(nonzero_mut_counts)
        plt.plot(x_fit, y_fit, "--", linewidth=2, alpha=0.8)
        out_png = f"{args.run_title.replace(' ', '_')}_fixed_indices_mutation_hist_log.png"
        plt.tight_layout()
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"\nSaved histogram: {out_png}")
    else:
        print("\nNo nonzero-mutation reads to plot (0-mutation reads are excluded by design).")

    # Also save a small text summary next to the png for convenience
    out_txt = f"{args.run_title.replace(' ', '_')}_fixed_indices_summary.txt"
    with open(out_txt, "w") as f:
        f.write("Fixed Index Pair Mutation Summary\n")
        f.write(f"Run title: {args.run_title}\n")
        f.write(f"R1 file: {args.file1}\n")
        f.write(f"R2 file: {args.file2}\n")
        f.write(f"Forward index: {fwd}\n")
        f.write(f"Reverse index: {rev} (RC end-match: {rev_rc})\n")
        f.write(f"Max hamming: {args.max_hamming}\n")
        f.write(f"Strip length: {args.strip_length}\n")
        f.write(f"Matched reads: {matched_total}\n")
        f.write(f"Zero-mutation reads: {zero_mut} ({zero_pct:.2f}%)\n\n")
        f.write("Stats for mutation counts (excluding zeros):\n")
        if not stats_dict:
            f.write("  None (no nonzero-mutation reads)\n")
        else:
            for k in ["n", "mean", "median", "std", "min", "q1", "q3", "max", "mode"]:
                f.write(f"  {k}: {stats_dict[k]}\n")
        f.write(f"\nAlignment used for: {aligned_used} reads\n")

    print(f"Saved summary: {out_txt}")
    print("\nDone.")


if __name__ == "__main__":
    main()
