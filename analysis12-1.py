import argparse
import sys
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Analyze FASTQ/TXT files for multiple conditions with mutation rate analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example usage:
  python script.py "Overall Title" "WT_SEQUENCE" file1.fastq file2.fastq \\
    "Time_0" seq1a seq1b seq2a seq2b seq3a seq3b \\
    "Time_1" seq4a seq4b seq5a seq5b seq6a seq6b \\
    "Time_2" seq7a seq7b seq8a seq8b seq9a seq9b \\
    "Time_3" seq10a seq10b seq11a seq11b seq12a seq12b \\
    "Time_Final" seq13a seq13b seq14a seq14b seq15a seq15b
    
  Each condition can have multiple index pairs (seq_a, seq_b).
        '''
    )
    
    parser.add_argument('overall_title', help='Overall title for the combined figure')
    parser.add_argument('wt_sequence', help='Wild-type reference sequence')
    parser.add_argument('file1', help='First FASTQ/TXT file (forward reads)')
    parser.add_argument('file2', help='Second FASTQ/TXT file (reverse reads)')
    
    # Parse remaining arguments as condition groups
    parser.add_argument('condition_data', nargs='+', 
                        help='Condition name followed by index sequence pairs')
    
    parser.add_argument('--max-hamming', type=int, default=1, 
                        help='Maximum hamming distance for index matching (default: 1)')
    parser.add_argument('--strip-length', type=int, default=8, 
                        help='Number of bases to strip from each end (default: 8)')
    parser.add_argument('--log-scale', action='store_true',
                        help='Use log scale for y-axis in histograms')
    
    args = parser.parse_args()
    conditions = []
    i = 0
    while i < len(args.condition_data):
        # First item should be condition name
        condition_name = args.condition_data[i]
        i += 1
        index_sequences = []
        while i < len(args.condition_data):
            current_item = args.condition_data[i]
            if any(c not in 'ACGTN' for c in current_item.upper()) or '_' in current_item:
                # This is probably the next condition name
                break
            index_sequences.append(current_item)
            i += 1
        
        if len(index_sequences) % 2 != 0:
            parser.error(f"Condition '{condition_name}' must have pairs of index sequences (got {len(index_sequences)} sequences)")
        
        if len(index_sequences) == 0:
            parser.error(f"Condition '{condition_name}' has no index sequences")
        index_pairs = []
        for j in range(0, len(index_sequences), 2):
            index_pairs.append((index_sequences[j], index_sequences[j+1]))
        conditions.append({
            'name': condition_name,
            'index_pairs': index_pairs
        })
    
    args.conditions = conditions
    return args

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in seq[::-1])

def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def read_sequences(filename):
    sequences = []
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        if lines and lines[0].startswith('@'):
            for i in range(1, len(lines), 4):
                if i < len(lines):
                    sequences.append(lines[i].strip())
        else:
            sequences = [line.strip() for line in lines if line.strip()]
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        sys.exit(1)
    
    return sequences

def check_sequence_match(seq, seq1, seq2_rc, max_hamming=1):
    """Check if sequence matches seq1 at start and seq2_rc at end"""
    seq1_match = False
    seq2_match = False
    
    # Check beginning
    if len(seq) >= len(seq1):
        start_seq = seq[:len(seq1)]
        if hamming_distance(start_seq, seq1) <= max_hamming:
            seq1_match = True
    
    # Check end
    if len(seq) >= len(seq2_rc):
        end_seq = seq[-len(seq2_rc):]
        if hamming_distance(end_seq, seq2_rc) <= max_hamming:
            seq2_match = True
    
    return seq1_match, seq2_match

def identify_mutations(seq, wt_seq):
    """Identify mutations between sequence and wild-type"""
    mutations = []
    min_len = min(len(seq), len(wt_seq))
    
    for i in range(min_len):
        if seq[i] != wt_seq[i]:
            mutations.append(f"{wt_seq[i]}{i+1}{seq[i]}")
    
    # If sequences have different lengths, note it
    if len(seq) != len(wt_seq):
        mutations.append(f"LengthDiff_{len(seq)}vs{len(wt_seq)}")
    
    return mutations

def pairwise_align(seq1, seq2, match_score=1, mismatch_score=-1, gap_score=-1):
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    for i in range(1, rows):
        score_matrix[i][0] = gap_score * i
    for j in range(1, cols):
        score_matrix[0][j] = gap_score * j
    
    for i in range(1, rows):
        for j in range(1, cols):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = score_matrix[i-1][j] + gap_score
            insert = score_matrix[i][j-1] + gap_score
            score_matrix[i][j] = max(match, delete, insert)
    
    # Traceback to get alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = rows - 1, cols - 1
    
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score_current = score_matrix[i][j]
            score_diag = score_matrix[i-1][j-1]
            score_up = score_matrix[i-1][j]
            score_left = score_matrix[i][j-1]
            
            if score_current == score_diag + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif score_current == score_up + gap_score:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
        elif i > 0:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
    
    # Reverse the alignments
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    
    return ''.join(aligned_seq1), ''.join(aligned_seq2)

def identify_mutations_with_alignment(seq, wt_seq, align_threshold=10):
    """Identify mutations, using alignment for sequences with many differences"""
    initial_mutations = identify_mutations(seq, wt_seq)
    if len(initial_mutations) > align_threshold:
        aligned_seq, aligned_wt = pairwise_align(seq, wt_seq)
        mutations = []
        wt_pos = 0
        
        for i, (s, w) in enumerate(zip(aligned_seq, aligned_wt)):
            if w != '-':
                wt_pos += 1
                if s == '-':
                    mutations.append(f"del{wt_pos}{w}")
                elif s != w:
                    mutations.append(f"{w}{wt_pos}{s}")
            elif s != '-':
                # Insertion
                mutations.append(f"ins{wt_pos}{s}")
        
        # Add alignment flag
        mutations.append("ALIGNED")
        return mutations, True  # True indicates alignment was used
    
    return initial_mutations, False  # False indicates no alignment was used

def group_identical_sequences(sequences_with_data):
    """Group sequences that are exactly identical"""
    sequence_groups = defaultdict(list)
    
    # Group identical sequences together
    for seq, data in sequences_with_data:
        sequence_groups[seq].append(data)
    
    # Convert to list format similar to original clustering output
    unique_sequences = []
    sequence_members = []
    
    for seq, data_list in sequence_groups.items():
        unique_sequences.append(seq)
        # Reconstruct the original format with (seq, data) tuples
        members = [(seq, data) for data in data_list]
        sequence_members.append(members)
    
    return unique_sequences, sequence_members

def calculate_mutation_statistics(mutation_counts):
    """Calculate various statistics for mutation counts"""
    if not mutation_counts:
        return {}
    
    from scipy import stats as scipy_stats
    
    stats_dict = {
        'mean': np.mean(mutation_counts),
        'median': np.median(mutation_counts),
        'std': np.std(mutation_counts),
        'min': np.min(mutation_counts),
        'max': np.max(mutation_counts),
        'q1': np.percentile(mutation_counts, 25),
        'q3': np.percentile(mutation_counts, 75),
        'mode': scipy_stats.mode(mutation_counts, keepdims=True)[0][0]
    }
    
    return stats_dict

def poisson_function(x, mu):
    """Poisson distribution function for fitting"""
    from scipy.special import factorial
    return np.exp(-mu) * (mu ** x) / factorial(x)

def fit_poisson_and_calculate_rate(mutation_counts, seq_length, exclude_zeros=True):
    """Fit Poisson distribution and calculate mutation rate"""
    if exclude_zeros:
        # Exclude sequences with 0 mutations (wild-type)
        mutation_counts = [m for m in mutation_counts if m > 0]
    
    if not mutation_counts:
        return None, None, None
    
    # Calculate mean (lambda parameter for Poisson)
    mean_mutations = np.mean(mutation_counts)
    
    # Calculate per-base mutation rate
    mutation_rate_per_base = mean_mutations / seq_length
    
    # Fit Poisson distribution
    try:
        # Create histogram data for fitting
        unique, counts = np.unique(mutation_counts, return_counts=True)
        # Normalize counts to get probabilities
        total_counts = np.sum(counts)
        probabilities = counts / total_counts
        
        # Initial guess for Poisson parameter
        initial_guess = [mean_mutations]
        
        # Fit using curve_fit (though for Poisson, mean is the MLE)
        popt, _ = curve_fit(lambda x, mu: poisson_function(x, mu) * total_counts, 
                           unique, counts, p0=initial_guess, maxfev=5000)
        fitted_mu = popt[0]
    except:
        fitted_mu = mean_mutations
    
    return mean_mutations, mutation_rate_per_base, fitted_mu

def process_with_indices(sequences1, sequences2, seq1, seq2, wt1, max_hamming, strip_length):
    """Process sequences with specific index sequences"""
    if len(sequences1) != len(sequences2):
        print(f"Error: Files have different number of sequences")
        return None
    
    # Get reverse complement of seq2 for matching
    seq2_rc = reverse_complement(seq2)
    
    # Process sequences and check matches
    combined_sequences = []
    both_matched = 0
    matched_stripped_sequences = []
    
    for seq1_read, seq2_read in zip(sequences1, sequences2):
        # Create combined sequence
        seq2_rc_full = reverse_complement(seq2_read)
        if len(seq2_rc_full) >= 69:
            combined_seq = seq1_read + seq2_rc_full[68:]
        else:
            combined_seq = seq1_read + seq2_rc_full
        
        combined_sequences.append(combined_seq)
        
        # Check matches
        seq1_match, seq2_match = check_sequence_match(combined_seq, seq1, seq2_rc, max_hamming)
        
        if seq1_match and seq2_match:
            both_matched += 1
            # Strip specified number of characters from each end
            if len(combined_seq) > 2 * strip_length:
                stripped_seq = combined_seq[strip_length:-strip_length]
                matched_stripped_sequences.append(stripped_seq)
    
    if not matched_stripped_sequences:
        print(f"  No sequences with both indices matched for seq1={seq1}, seq2={seq2}")
        return None
    
    # Compare to wild-type and identify mutations
    wt_stripped = wt1
    
    sequences_with_mutations = []
    mutation_counts = []
    all_mutations_list = []
    aligned_count = 0
    
    for seq in matched_stripped_sequences:
        mutations, was_aligned = identify_mutations_with_alignment(seq, wt_stripped)
        if was_aligned:
            aligned_count += 1
            # Remove the ALIGNED flag from mutations for counting
            mutations = [m for m in mutations if m != "ALIGNED"]
        sequences_with_mutations.append((seq, mutations))
        mutation_counts.append(len(mutations))
        all_mutations_list.extend(mutations)
    
    if aligned_count > 0:
        print(f"  Used alignment for {aligned_count} sequences with high mutation counts")
    
    # Calculate mutation statistics
    mut_stats = calculate_mutation_statistics(mutation_counts)
    
    # Calculate mutation rate (excluding wild-type sequences)
    mean_mutations, mutation_rate_per_base, fitted_mu = fit_poisson_and_calculate_rate(
        mutation_counts, len(wt_stripped), exclude_zeros=True
    )
    
    # Group identical sequences
    unique_sequences, sequence_groups = group_identical_sequences(sequences_with_mutations)
    
    # Prepare sequence group data with read counts
    group_data = []
    sequence_read_counts = Counter()  # For rank abundance plot
    
    for i, (unique_seq, members) in enumerate(zip(unique_sequences, sequence_groups)):
        group_name = f"seq{chr(65 + i)}"
        consensus_seq = unique_seq
        read_count = len(members)
        
        # Track read counts for rank abundance
        sequence_read_counts[unique_seq] = read_count
        
        # Calculate mutation statistics for this group
        group_mutation_counts = [len(muts) for _, muts in members]
        group_mut_stats = calculate_mutation_statistics(group_mutation_counts) if group_mutation_counts else {}
        
        # Check if this is the WT sequence - both by sequence comparison and mutation count
        all_mutations = [mut for _, muts in members for mut in muts]
        is_wild_type = False
        
        # First check: exact sequence match
        if consensus_seq == wt_stripped:
            is_wild_type = True
            group_name = "seqWT"
        # Second check: no mutations or only length differences
        elif not all_mutations or all(mut.startswith("LengthDiff") for mut in all_mutations):
            if len(group_mutation_counts) > 0 and all(count == 0 for count in group_mutation_counts):
                is_wild_type = True
                group_name = "seqWT"
        
        group_data.append({
            'name': group_name,
            'consensus': consensus_seq,
            'count': read_count,
            'mutations': Counter([tuple(muts) for _, muts in members]),
            'mut_stats': group_mut_stats,
            'is_wild_type': is_wild_type
        })
    
    # Sort groups by read count
    group_data.sort(key=lambda x: x['count'], reverse=True)
    
    return {
        'total_sequences': len(combined_sequences),
        'both_matched': both_matched,
        'matched_stripped_sequences': matched_stripped_sequences,
        'mutation_counts': mutation_counts,
        'mut_stats': mut_stats,
        'mean_mutations_no_wt': mean_mutations,
        'mutation_rate_per_base': mutation_rate_per_base,
        'poisson_mu': fitted_mu,
        'group_data': group_data,
        'all_mutations_list': all_mutations_list,
        'sequence_read_counts': sequence_read_counts,
        'seq1': seq1,
        'seq2': seq2,
        'aligned_count': aligned_count
    }

def create_mutation_rate_histograms(all_conditions_results, wt_seq, overall_title):
    """Create histograms with Poisson fitting, excluding wild-type sequences"""
    n_conditions = len(all_conditions_results)
    
    # Create figure with subplots
    fig, axes = plt.subplots(1, n_conditions, figsize=(5 * n_conditions, 5))
    if n_conditions == 1:
        axes = [axes]
    
    # Store mutation rates for summary
    mutation_rates_summary = []
    
    for idx, (condition_name, results) in enumerate(all_conditions_results):
        ax = axes[idx]
        
        # Determine concentration labels based on index pairs
        n_samples = len(results)
        if n_samples > 1:
            # Assume different concentrations for different index pairs
            concentrations = [0, 10, 100]  # mM concentrations
            labels = [f'{conc} mM' for conc in concentrations[:n_samples]]
            colors = plt.cm.viridis(np.linspace(0.2, 0.9, n_samples))
        else:
            labels = ['']
            colors = ['blue']
        
        # Plot histograms for each concentration
        for i, (result, color, label) in enumerate(zip(results, colors, labels)):
            if result and result['mutation_counts']:
                # Exclude wild-type (0 mutations) sequences
                mutation_counts_no_wt = [m for m in result['mutation_counts'] if m > 0]
                
                if mutation_counts_no_wt:
                    # Create histogram
                    counts, bins, _ = ax.hist(mutation_counts_no_wt, bins=range(1, 31), 
                                             alpha=0.6, color=color, edgecolor='black',
                                             label=f"{label} (n={len(mutation_counts_no_wt)})", 
                                             density=False)
                    
                    # Fit Poisson distribution
                    mean_mutations = np.mean(mutation_counts_no_wt)
                    mutation_rate = mean_mutations / len(wt_seq)
                    
                    # Plot Poisson fit
                    x_fit = np.arange(1, 30)
                    y_fit = poisson_function(x_fit, mean_mutations) * len(mutation_counts_no_wt)
                    ax.plot(x_fit, y_fit, '--', color=color, alpha=0.8, linewidth=2)
                    
                    # Store results
                    mutation_rates_summary.append({
                        'condition': condition_name,
                        'concentration': label,
                        'mean_mutations': mean_mutations,
                        'mutation_rate_per_base': mutation_rate,
                        'n_sequences': len(mutation_counts_no_wt)
                    })
        
        ax.set_xlabel('Number of Mutations', fontsize=11)
        ax.set_ylabel('Number of Sequences', fontsize=11)
        ax.set_title(f'{condition_name}\n(Wild-type excluded)', fontsize=12)
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.3)
        
        # Add text with mutation rate info
        info_text = []
        for rate_info in mutation_rates_summary:
            if rate_info['condition'] == condition_name:
                info_text.append(f"{rate_info['concentration']}: μ={rate_info['mean_mutations']:.2f}, "
                               f"rate={rate_info['mutation_rate_per_base']:.4f}/bp")
        
        if info_text:
            ax.text(0.95, 0.85, '\n'.join(info_text), transform=ax.transAxes,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                   verticalalignment='top', horizontalalignment='right', fontsize=8)
    
    plt.suptitle(f'{overall_title} - Mutation Rate Analysis (WT Excluded)', fontsize=14)
    plt.tight_layout()
    
    # Save figure
    filename = f'{overall_title.replace(" ", "_")}_mutation_rates.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    return filename, mutation_rates_summary

def create_rank_abundance_plots(all_conditions_results, overall_title):
    """Create rank abundance plots for different conditions"""
    # Group conditions by type (time-based or concentration-based)
    time_conditions = {}
    concentration_conditions = defaultdict(list)
    
    for condition_name, results in all_conditions_results:
        if 'time' in condition_name.lower() or 'stage' in condition_name.lower():
            # This is a time point
            time_conditions[condition_name] = results
        else:
            # Group by concentration
            # Extract concentration value if present
            if 'mm' in condition_name.lower():
                concentration_conditions[condition_name].append((condition_name, results))
    
    # Create figure with subplots
    n_plots = len(time_conditions) + len(concentration_conditions)
    if n_plots == 0:
        n_plots = 1
        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        axes = [axes]
    else:
        fig, axes = plt.subplots(1, n_plots, figsize=(6 * n_plots, 5))
        if n_plots == 1:
            axes = [axes]
    
    plot_idx = 0
    
    # Plot time-based conditions
    for time_name, time_results in time_conditions.items():
        ax = axes[plot_idx]
        
        # Aggregate read counts across all index pairs
        read_count_distribution = Counter()
        
        for result in time_results:
            if result:
                # Count how many sequences have each read count
                for seq, count in result['sequence_read_counts'].items():
                    read_count_distribution[count] += 1
        
        if read_count_distribution:
            # Sort by read count
            sorted_counts = sorted(read_count_distribution.items())
            x_values = [x for x, y in sorted_counts]
            y_values = [y for x, y in sorted_counts]
            
            # Plot
            ax.loglog(x_values, y_values, 'o-', linewidth=2, markersize=6)
            ax.set_xlabel('Number of Reads', fontsize=11)
            ax.set_ylabel('Number of Sequences', fontsize=11)
            ax.set_title(f'Rank Abundance - {time_name}', fontsize=12)
            ax.grid(True, alpha=0.3, which='both')
        
        plot_idx += 1
    
    # Plot concentration-based conditions
    for conc_group_name, conc_conditions in concentration_conditions.items():
        if plot_idx < len(axes):
            ax = axes[plot_idx]
            
            colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(conc_conditions)))
            
            for (cond_name, cond_results), color in zip(conc_conditions, colors):
                # Aggregate read counts
                read_count_distribution = Counter()
                
                for result in cond_results:
                    if result:
                        for seq, count in result['sequence_read_counts'].items():
                            read_count_distribution[count] += 1
                
                if read_count_distribution:
                    sorted_counts = sorted(read_count_distribution.items())
                    x_values = [x for x, y in sorted_counts]
                    y_values = [y for x, y in sorted_counts]
                    
                    ax.loglog(x_values, y_values, 'o-', linewidth=2, markersize=6,
                             label=cond_name, color=color)
            
            ax.set_xlabel('Number of Reads', fontsize=11)
            ax.set_ylabel('Number of Sequences', fontsize=11)
            ax.set_title(f'Rank Abundance - {conc_group_name}', fontsize=12)
            ax.grid(True, alpha=0.3, which='both')
            ax.legend(loc='best', fontsize=9)
            
            plot_idx += 1
    
    plt.suptitle(f'{overall_title} - Rank Abundance Analysis', fontsize=14)
    plt.tight_layout()
    
    # Save figure
    filename = f'{overall_title.replace(" ", "_")}_rank_abundance.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    return filename

def export_mutated_sequences_summary(top_sequences_data, overall_title, wt_sequence):
    """Export detailed summary of top 20 mutated sequences from Mueller plots"""
    filename = f"{overall_title.replace(' ', '_')}_mutated_sequences_summary.txt"
    
    with open(filename, 'w') as f:
        f.write(f"Mutated Sequences Summary\n")
        f.write(f"Overall Title: {overall_title}\n")
        f.write(f"Wild-type Sequence: {wt_sequence}\n")
        f.write("=" * 100 + "\n\n")
        
        for i, seq_data in enumerate(top_sequences_data, 1):
            f.write(f"Seq {i}\n")
            f.write("-" * 80 + "\n")
            
            # Write full sequence
            f.write(f"Sequence: {seq_data['sequence']}\n")
            
            # Write mutations from WT
            if seq_data['mutations']:
                # Filter out ALIGNED and LengthDiff markers
                actual_mutations = [m for m in seq_data['mutations'] 
                                  if not m.startswith("LengthDiff") and m != "ALIGNED"]
                if actual_mutations:
                    f.write(f"Changes from WT: {', '.join(actual_mutations)}\n")
                else:
                    f.write("Changes from WT: None (Wild-type)\n")
            else:
                f.write("Changes from WT: None (Wild-type)\n")
            
            # Write counts at each time point
            for time_label, count in seq_data['time_counts'].items():
                f.write(f"Count at {time_label}: {count}\n")
            
            # Add relative abundance info
            f.write(f"Maximum relative abundance: {seq_data['max_abundance']:.4f}\n")
            
            f.write("\n")
    
    print(f"Exported mutated sequences summary to: {filename}")
    return filename

def create_mueller_plots_with_summary(all_conditions_results, overall_title, wt_sequence):
    """Modified version with THREE Mueller plots: all, no WT, and WT + top 20 mutants"""
    
    # First, identify which sequences are wild-type (0 mutations)
    wt_sequences = set()
    all_sequence_mutation_info = {}  # sequence -> mutation list
    
    for condition_name, results in all_conditions_results:
        for result in results:
            if result:
                for group in result['group_data']:
                    seq = group['consensus']
                    
                    # Get the actual mutations for this sequence
                    if group['mutations']:
                        most_common = group['mutations'].most_common(1)[0]
                        mutations = list(most_common[0])
                    else:
                        mutations = []
                    
                    # Store mutations for this sequence
                    all_sequence_mutation_info[seq] = mutations
                    
                    # Check if it's wild-type
                    if 'is_wild_type' in group and group['is_wild_type']:
                        wt_sequences.add(seq)
                    elif seq == wt_sequence:
                        wt_sequences.add(seq)
                    elif len(mutations) == 0:
                        wt_sequences.add(seq)
                    elif len(mutations) == 1 and mutations[0].startswith("LengthDiff"):
                        wt_sequences.add(seq)
                    else:
                        actual_mutations = [m for m in mutations if not m.startswith("LengthDiff") and m != "ALIGNED"]
                        if len(actual_mutations) == 0:
                            wt_sequences.add(seq)
    
    print(f"Identified {len(wt_sequences)} wild-type sequences")
    
    # Collect all sequence abundance data first
    time_points_dict = {}
    sequence_time_counts = defaultdict(lambda: defaultdict(int))
    sequence_abundances = defaultdict(list)
    sequence_abundances_no_wt = defaultdict(list)  # NEW: for exclude_wt plot
    
    for condition_name, results in all_conditions_results:
        if 'time' in condition_name.lower() or any(x in condition_name.lower() for x in ['0', '1', '2', '3', 'final']):
            # Extract time point number
            if 'final' in condition_name.lower():
                time_point = 4
            elif '0' in condition_name:
                time_point = 0
            elif '1' in condition_name:
                time_point = 1
            elif '2' in condition_name:
                time_point = 2
            elif '3' in condition_name:
                time_point = 3
            else:
                continue
            
            time_points_dict[time_point] = condition_name
            
            # Calculate total reads and individual sequence reads
            total_reads = 0
            total_reads_no_wt = 0  # NEW: total excluding WT
            seq_reads = {}
            
            for result in results:
                if result:
                    for group in result['group_data']:
                        seq = group['consensus']
                        seq_reads[seq] = seq_reads.get(seq, 0) + group['count']
                        total_reads += group['count']
                        sequence_time_counts[seq][condition_name] += group['count']
                        
                        # NEW: Track non-WT reads separately
                        if seq not in wt_sequences:
                            total_reads_no_wt += group['count']
            
            # Calculate relative abundances (including WT)
            for seq in seq_reads:
                relative_abundance = seq_reads[seq] / total_reads if total_reads > 0 else 0
                sequence_abundances[seq].append((time_point, relative_abundance))
                
                # NEW: Calculate relative abundances excluding WT (normalized to non-WT total)
                if seq not in wt_sequences:
                    relative_abundance_no_wt = seq_reads[seq] / total_reads_no_wt if total_reads_no_wt > 0 else 0
                    sequence_abundances_no_wt[seq].append((time_point, relative_abundance_no_wt))
    
    if not time_points_dict:
        print(f"No time-based conditions found for Mueller plot")
        return [], None
    
    # Sort time points
    sorted_times = sorted(time_points_dict.keys())
    
    # Get top mutated sequences (excluding WT) for the third plot
    mutant_abundances = {seq: max(abd for _, abd in abds) 
                        for seq, abds in sequence_abundances_no_wt.items()}  # Use normalized abundances
    top_mutants = sorted(mutant_abundances.items(), key=lambda x: x[1], reverse=True)[:20]
    
    # Create THREE plots
    filenames = []
    top_sequences_for_summary = None
    
    plot_configs = [
        (False, False, "All Sequences"),  # Include all
        (True, False, "Wild-type Excluded"),  # Exclude WT - MODIFIED BEHAVIOR
        (False, True, "Wild-type + Top 20 Mutants"),  # WT + top 20 mutants only
    ]
    
    for exclude_wt, wt_plus_top20_only, plot_subtitle in plot_configs:
        # Filter sequences based on plot type
        if wt_plus_top20_only:
            # Only include WT and top 20 mutants
            included_sequences = set([seq for seq, _ in top_mutants]) | wt_sequences
            filtered_abundances = {seq: abds for seq, abds in sequence_abundances.items()
                                  if seq in included_sequences}
        elif exclude_wt:
            # Use the normalized abundances for exclude_wt plot
            filtered_abundances = sequence_abundances_no_wt  # CHANGED: use normalized abundances
        else:
            # Include all sequences
            filtered_abundances = sequence_abundances
        
        # Get top sequences for this plot
        max_abundances = {seq: max(abd for _, abd in abds) for seq, abds in filtered_abundances.items()}
        
        if wt_plus_top20_only:
            # For WT + top 20 plot, sort WT first, then by abundance
            wt_seqs_list = [(seq, abd) for seq, abd in max_abundances.items() if seq in wt_sequences]
            mutant_seqs_list = [(seq, abd) for seq, abd in max_abundances.items() if seq not in wt_sequences]
            
            wt_seqs_list.sort(key=lambda x: x[1], reverse=True)
            mutant_seqs_list.sort(key=lambda x: x[1], reverse=True)
            
            top_sequences = wt_seqs_list + mutant_seqs_list[:20]
        else:
            top_sequences = sorted(max_abundances.items(), key=lambda x: x[1], reverse=True)[:20]
        
        # Calculate "Other" category - aggregate all non-top-20 sequences
        top_sequences_set = set(seq for seq, _ in top_sequences)
        other_abundances = []
        
        for time_point in sorted_times:
            other_abundance = 0
            for seq, time_abds in filtered_abundances.items():
                if seq not in top_sequences_set:
                    for t, abd in time_abds:
                        if t == time_point:
                            other_abundance += abd
                            break
            other_abundances.append(other_abundance)
        
        # Collect data for summary (only from the plot that excludes WT)
        if exclude_wt and top_sequences_for_summary is None:
            top_sequences_for_summary = []
            for seq, max_abd in top_sequences:
                seq_data = {
                    'sequence': seq,
                    'mutations': all_sequence_mutation_info.get(seq, []),
                    'time_counts': {},
                    'max_abundance': max_abd
                }
                
                for time_point in sorted_times:
                    time_label = time_points_dict[time_point]
                    seq_data['time_counts'][time_label] = sequence_time_counts[seq].get(time_label, 0)
                
                top_sequences_for_summary.append(seq_data)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 7))
        
        # Prepare data for stacked area chart
        sequence_data = {}
        sequence_labels = []
        
        for i, (seq, _) in enumerate(top_sequences):
            seq_abundances = dict(filtered_abundances[seq])
            sequence_data[seq] = [seq_abundances.get(t, 0) for t in sorted_times]
            
            # Create informative labels
            if seq in wt_sequences:
                label = f'WT (Seq {i+1})'
            else:
                mut_list = all_sequence_mutation_info.get(seq, [])
                actual_mutations = [m for m in mut_list if not m.startswith("LengthDiff") and m != "ALIGNED"]
                mut_count = len(actual_mutations)
                label = f'Seq {i+1} ({mut_count} mut)'
            sequence_labels.append(label)
        
        # Add "Other" category with all remaining sequences aggregated
        if sum(other_abundances) > 0:
            sequence_data['__other__'] = other_abundances
            other_seq_count = len(filtered_abundances) - len(top_sequences)
            sequence_labels.append(f'Other ({other_seq_count} seqs)')
        
        # Create stacked area chart
        x = sorted_times
        y_stack = []
        
        # Add data for top sequences
        for seq, _ in top_sequences:
            y_stack.append(sequence_data[seq])
        
        # Add "Other" category if it exists
        if sum(other_abundances) > 0:
            y_stack.append(other_abundances)
        
        y_stack = np.vstack(y_stack)
        
        # Use different colormaps for different plots
        if wt_plus_top20_only:
            # Special coloring for WT + top 20 plot
            colors = []
            for seq, _ in top_sequences:
                if seq in wt_sequences:
                    colors.append('lightgreen')
                else:
                    color_idx = len([c for c in colors if c != 'lightgreen'])
                    colors.append(plt.cm.tab20(color_idx % 20))
            if sum(other_abundances) > 0:
                colors.append('lightgray')
        elif exclude_wt:
            # For exclude WT plot - use distinct colors for top 20
            colors = list(plt.cm.tab20(np.linspace(0, 1, len(top_sequences))))
            if sum(other_abundances) > 0:
                colors.append('lightgray')
        else:
            # For all sequences plot
            colors = []
            for seq, _ in top_sequences:
                if seq in wt_sequences:
                    colors.append('lightgreen')
                else:
                    color_idx = len([c for c in colors if c != 'lightgreen'])
                    colors.append(plt.cm.tab20(color_idx % 20))
            if sum(other_abundances) > 0:
                colors.append('lightgray')
        
        # Plot stacked areas
        ax.stackplot(x, y_stack, labels=sequence_labels, colors=colors, alpha=0.8)
        
        ax.set_xlabel('Time Point', fontsize=12)
        ax.set_ylabel('Relative Abundance', fontsize=12)
        ax.set_title(f'{overall_title} - Mueller Plot ({plot_subtitle})', fontsize=14)
        
        ax.set_xticks(sorted_times)
        ax.set_xticklabels([f'T{t}' for t in sorted_times])
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add legend
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=1)
        
        # Add annotation
        if wt_plus_top20_only:
            wt_count = sum(1 for seq, _ in top_sequences if seq in wt_sequences)
            mut_count = len(top_sequences) - wt_count
            other_pct = np.mean(other_abundances) * 100 if other_abundances else 0
            ax.text(0.02, 0.98, f'Wild-type sequences: {wt_count}\n'
                               f'Top mutant sequences: {mut_count}\n'
                               f'Other sequences: {other_pct:.1f}% avg',
                   transform=ax.transAxes, fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                   verticalalignment='top', horizontalalignment='left')
        elif exclude_wt:
            other_pct = np.mean(other_abundances) * 100 if other_abundances else 0
            other_count = len(filtered_abundances) - len(top_sequences)
            ax.text(0.02, 0.98, f'Wild-type sequences excluded\n'
                               f'Top 20 mutant sequences shown\n'
                               f'Other sequences ({other_count}): {other_pct:.1f}% avg',
                   transform=ax.transAxes, fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                   verticalalignment='top', horizontalalignment='left')
        else:
            wt_count = sum(1 for seq, _ in top_sequences if seq in wt_sequences)
            other_pct = np.mean(other_abundances) * 100 if other_abundances else 0
            ax.text(0.02, 0.98, f'Wild-type sequences: {wt_count}/{len(top_sequences)}\n'
                               f'Other sequences: {other_pct:.1f}% avg',
                   transform=ax.transAxes, fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                   verticalalignment='top', horizontalalignment='left')
        
        plt.tight_layout()
        
        # Save figure
        if wt_plus_top20_only:
            filename = f'{overall_title.replace(" ", "_")}_mueller_plot_wt_plus_top20.png'
        elif exclude_wt:
            filename = f'{overall_title.replace(" ", "_")}_mueller_plot_no_wt.png'
        else:
            filename = f'{overall_title.replace(" ", "_")}_mueller_plot_all.png'
        
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        filenames.append(filename)
        print(f"  Created: {filename}")
        
        # Report statistics
        if sum(other_abundances) > 0:
            other_seq_count = len(filtered_abundances) - len(top_sequences)
            total_seqs = len(top_sequences) + other_seq_count
            avg_other = np.mean(other_abundances)*100
            print(f"    Total sequences in plot: {total_seqs} (Top {len(top_sequences)} + {other_seq_count} others)")
            print(f"    Average 'Other' category abundance: {avg_other:.1f}%")
        else:
            print(f"    Total sequences in plot: {len(top_sequences)}")
    
    return filenames, top_sequences_for_summary     

def create_combined_figure(all_conditions_results, wt_seq, overall_title, use_log_scale=False):
    """Create a combined figure with histograms, enhanced box plots, and mutation heatmaps"""
    n_conditions = len(all_conditions_results)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(6 * n_conditions, 20))
    
    # Use GridSpec for better control
    gs = gridspec.GridSpec(3, n_conditions, figure=fig, 
                          height_ratios=[1, 0.5, 1],
                          hspace=0.4, wspace=0.2)
    
    # Create subplots
    ax_histograms = [fig.add_subplot(gs[0, i]) for i in range(n_conditions)]
    ax_boxplots = [fig.add_subplot(gs[1, i]) for i in range(n_conditions)]
    ax_heatmaps = [fig.add_subplot(gs[2, i]) for i in range(n_conditions)]
    
    # Store all sequences for export
    all_sequence_data = []
    
    # Process each condition
    for idx, (condition_name, results) in enumerate(all_conditions_results):
        # Create histogram for this condition
        ax_hist = ax_histograms[idx]
        
        # Determine if this is a time-based or concentration-based condition
        is_time_based = any(keyword in condition_name.lower() for keyword in 
                           ['time', 'day', 'hour', 'stage', 'final', 'initial'])
        is_concentration_based = any(keyword in condition_name.lower() for keyword in 
                                    ['mm', 'um', 'mg', 'ug', 'conc', 'dose'])
        
        # Choose color palette
        if is_time_based:
            cmap = cm.get_cmap('viridis')
            n_samples = len(results)
            colors = [cmap(i/(n_samples-1) if n_samples > 1 else 0.5) for i in range(n_samples)]
        elif is_concentration_based:
            cmap = cm.get_cmap('Blues')
            n_samples = len(results)
            colors = [cmap(0.3 + 0.6*i/(n_samples-1) if n_samples > 1 else 0.6) for i in range(n_samples)]
        else:
            cmap = cm.get_cmap('RdYlBu_r')
            n_samples = len(results)
            colors = [cmap(i/(n_samples-1) if n_samples > 1 else 0.5) for i in range(n_samples)]
        
        # Create labels
        if n_samples > 1:
            labels = [f'{40*(i**2)-30*i} mM' for i in range(n_samples)]
        else:
            labels = ['']
        
        # Plot histograms
        max_count = 0
        all_mutation_data = []
        
        for i, (result, color, label) in enumerate(zip(results, colors, labels)):
            if result and result['mutation_counts']:
                mutation_counts = result['mutation_counts']
                
                # Exclude wild-type for mutation rate calculation
                mutation_counts_no_wt = [m for m in mutation_counts if m > 0]
                
                # Plot histogram with ALL sequences
                n, bins, patches = ax_hist.hist(mutation_counts, bins=30, range=(0, 30),
                                               alpha=0.6, edgecolor='black', linewidth=0.5,
                                               color=color, label=f"{label} (n={result['both_matched']})")
                
                max_count = max(max_count, max(n))
                all_mutation_data.append(mutation_counts)
                
                # Add Poisson fit line (excluding WT)
                if mutation_counts_no_wt:
                    mean_mut = np.mean(mutation_counts_no_wt)
                    x_fit = np.arange(1, 30)
                    y_fit = poisson_function(x_fit, mean_mut) * len(mutation_counts_no_wt)
                    ax_hist.plot(x_fit, y_fit, '--', color=color, alpha=0.8, linewidth=1)
        
        # Set log scale if requested
        if use_log_scale:
            ax_hist.set_yscale('log')
            ax_hist.set_ylim(bottom=0.8)
        
        ax_hist.set_xlabel('Number of mutations', fontsize=11)
        ax_hist.set_ylabel('Number of sequences' + (' (log scale)' if use_log_scale else ''), fontsize=11)
        ax_hist.set_title(f'{condition_name}\nMutation Distribution', fontsize=12)
        ax_hist.grid(True, alpha=0.3)
        
        if n_samples > 1:
            ax_hist.legend(loc='upper right', fontsize=8)
        
        # Create ENHANCED box plots with more detail
        ax_box = ax_boxplots[idx]
        if all_mutation_data:
            # Create box plot with additional features
            bp = ax_box.boxplot(all_mutation_data, labels=labels if n_samples > 1 else [''],
                               patch_artist=True, widths=0.6, 
                               showfliers=True,  # Show outliers
                               showmeans=True,   # Show mean as well as median
                               meanprops=dict(marker='D', markeredgecolor='red', markerfacecolor='red', markersize=6),
                               flierprops=dict(marker='o', markerfacecolor='gray', markersize=3, alpha=0.5))
            
            # Color the boxes
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            # Add detailed statistics annotations
            for i, data in enumerate(all_mutation_data):
                if data:
                    mean_val = np.mean(data)
                    median_val = np.median(data)
                    std_val = np.std(data)
                    
                    # Add text with statistics
                    y_pos = ax_box.get_ylim()[1] * 0.95
                    ax_box.text(i + 1, y_pos, 
                               f'μ={mean_val:.2f}\nσ={std_val:.2f}',
                               horizontalalignment='center',
                               fontsize=7,
                               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
            
            ax_box.set_ylabel('Number of Mutations', fontsize=10)
            ax_box.grid(True, alpha=0.3, axis='y')
            ax_box.set_title('Distribution Summary (red diamond = mean)', fontsize=10)
            
            # Set y-axis to show decimals if needed
            from matplotlib.ticker import MaxNLocator
            ax_box.yaxis.set_major_locator(MaxNLocator(integer=False, nbins=10))
        
        # Create mutation heatmap
        ax_heat = ax_heatmaps[idx]
        
        # Aggregate all sequences from all index pairs
        all_sequences = {}
        
        for result in results:
            if result:
                for group in result['group_data']:
                    seq_key = group['consensus']
                    if seq_key not in all_sequences:
                        all_sequences[seq_key] = {
                            'name': group['name'],
                            'total_count': 0,
                            'mutations': group['mutations']
                        }
                    all_sequences[seq_key]['total_count'] += group['count']
        
        # Sort and get top 20
        sorted_sequences = sorted(all_sequences.items(), 
                                key=lambda x: x[1]['total_count'], reverse=True)
        top20_sequences = sorted_sequences[:20]
        
        if top20_sequences:
            # Prepare heatmap data
            sequence_names = [f"{data['name']} ({data['total_count']})" 
                            for _, data in top20_sequences]
            seq_length = len(wt_seq)
            
            # Initialize mutation frequency matrix
            mutation_matrix = np.zeros((len(sequence_names), seq_length))
            
            # Store sequence data for export
            for i, (seq, data) in enumerate(top20_sequences):
                # Count mutations at each position
                position_counts = defaultdict(int)
                
                # Get the most common mutation pattern
                if data['mutations']:
                    most_common_mutations = data['mutations'].most_common(1)[0][0]
                else:
                    most_common_mutations = []
                
                # Store sequence info
                all_sequence_data.append({
                    'condition': condition_name,
                    'name': data['name'],
                    'sequence': seq,
                    'count': data['total_count'],
                    'mutations': list(most_common_mutations)
                })
                
                for mutations, count in data['mutations'].items():
                    for mut in mutations:
                        if not mut.startswith("LengthDiff") and mut != "ALIGNED":
                            if mut.startswith("del"):
                                pos = int(mut[3:-1]) - 1
                            elif mut.startswith("ins"):
                                pos = int(mut[3:-1]) - 1
                            else:
                                pos = int(mut[1:-1]) - 1
                            
                            if 0 <= pos < seq_length:
                                position_counts[pos] += count
                
                # Convert to binary presence
                for pos, count in position_counts.items():
                    mutation_matrix[i, pos] = 1
            
            # Plot heatmap
            im = ax_heat.imshow(mutation_matrix, cmap='YlOrRd', aspect='auto', 
                               interpolation='nearest')
            
            # X-axis labels
            x_tick_spacing = 50
            ax_heat.set_xticks(np.arange(0, seq_length, x_tick_spacing))
            ax_heat.set_xticklabels(np.arange(1, seq_length+1, x_tick_spacing))
            ax_heat.set_yticks(range(len(sequence_names)))
            ax_heat.set_yticklabels(sequence_names, fontsize=9)
            
            ax_heat.set_xlabel('Position in sequence', fontsize=11)
            if idx == 0:
                ax_heat.set_ylabel('Top 20 Sequences', fontsize=11)
            ax_heat.set_title(f'{condition_name}\nMutation Positions', fontsize=12)
            
            # Add grid
            ax_heat.set_xticks(np.arange(-0.5, seq_length, 1), minor=True)
            ax_heat.set_yticks(np.arange(-0.5, len(sequence_names), 1), minor=True)
            ax_heat.grid(which='minor', color='gray', linestyle='-', linewidth=0.1, alpha=0.3)
    
    # Add overall title
    fig.suptitle(overall_title, fontsize=16, y=0.98)
    
    # Save figure
    filename = f'{overall_title.replace(" ", "_")}_combined_analysis{"_logscale" if use_log_scale else ""}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    return filename

def export_mutation_rates_summary(mutation_rates_summary, overall_title):
    """Export mutation rates summary to text file"""
    filename = f"{overall_title.replace(' ', '_')}_mutation_rates_summary.txt"
    
    with open(filename, 'w') as f:
        f.write(f"Mutation Rate Analysis Summary\n")
        f.write(f"Overall Title: {overall_title}\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("Mutation Rates (excluding wild-type sequences):\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'Condition':<20} {'Concentration':<15} {'Mean Mutations':<15} {'Rate per Base':<15} {'N Sequences':<15}\n")
        f.write("-" * 60 + "\n")
        
        for rate_info in mutation_rates_summary:
            f.write(f"{rate_info['condition']:<20} {rate_info['concentration']:<15} "
                   f"{rate_info['mean_mutations']:<15.2f} {rate_info['mutation_rate_per_base']:<15.6f} "
                   f"{rate_info['n_sequences']:<15}\n")
    
    print(f"Exported mutation rates summary to: {filename}")

def main():
    args = parse_arguments()
    
    # Read input files once
    print("Reading sequence files...")
    sequences1 = read_sequences(args.file1)
    sequences2 = read_sequences(args.file2)
    
    if len(sequences1) != len(sequences2):
        print("Error: Files have different number of sequences")
        sys.exit(1)
    
    print(f"Found {len(sequences1)} sequence pairs")
    print(f"Processing {len(args.conditions)} conditions")
    print(f"Max hamming distance: {args.max_hamming}")
    print(f"Strip length: {args.strip_length} bases from each end")
    print(f"Log scale for histograms: {'Yes' if args.log_scale else 'No'}\n")
    
    # Process each condition
    all_conditions_results = []
    
    for condition in args.conditions:
        condition_name = condition['name']
        index_pairs = condition['index_pairs']
        
        print(f"\nProcessing condition: {condition_name}")
        print(f"  Number of index pairs: {len(index_pairs)}")
        
        # Process each index pair for this condition
        condition_results = []
        
        for i, (seq1, seq2) in enumerate(index_pairs):
            print(f"\n  Index pair {i+1}:")
            print(f"    seq1: {seq1}")
            print(f"    seq2: {seq2} (RC: {reverse_complement(seq2)})")
            
            result = process_with_indices(
                sequences1, sequences2,
                seq1, seq2, args.wt_sequence,
                args.max_hamming, args.strip_length
            )
            
            if result:
                print(f"    Matched: {result['both_matched']}")
                print(f"    Unique sequences: {len(result['group_data'])}")
                if result['mut_stats']:
                    print(f"    Mean mutations (all): {result['mut_stats']['mean']:.2f}")
                if result['mean_mutations_no_wt'] is not None:
                    print(f"    Mean mutations (no WT): {result['mean_mutations_no_wt']:.2f}")
                    print(f"    Mutation rate per base: {result['mutation_rate_per_base']:.6f}")
            
            condition_results.append(result)
        
        all_conditions_results.append((condition_name, condition_results))
    
    # Create all figures
    print("\nCreating analysis figures...")
    
    # 1. Mutation rate histograms with Poisson fitting
    print("  - Creating mutation rate histograms...")
    mut_rate_file, mutation_rates_summary = create_mutation_rate_histograms(
        all_conditions_results, args.wt_sequence, args.overall_title
    )
    
    # 2. Combined figure with enhanced box plots
    print("  - Creating combined analysis figure with enhanced box plots...")
    combined_file = create_combined_figure(
        all_conditions_results, args.wt_sequence, args.overall_title, args.log_scale
    )
    
    # 3. Rank abundance plots
    print("  - Creating rank abundance plots...")
    rank_abundance_file = create_rank_abundance_plots(
        all_conditions_results, args.overall_title
    )
    
    # 4. Mueller plots (FOUR versions now)
    print("  - Creating Mueller plots (4 versions)...")
    mueller_files, top_sequences_data = create_mueller_plots_with_summary(
        all_conditions_results, args.overall_title, args.wt_sequence
    )
    
    # 5. Export mutation rates summary
    export_mutation_rates_summary(mutation_rates_summary, args.overall_title)
    
    # 6. Export mutated sequences summary
    print("  - Exporting mutated sequences summary...")
    if top_sequences_data:
        mutated_seq_summary_file = export_mutated_sequences_summary(
            top_sequences_data, args.overall_title, args.wt_sequence
        )
    else:
        mutated_seq_summary_file = None
        print("    No sequence data available for summary")
    print("\nAnalysis complete!")
    print(f"Output files created:")
    print(f"  - {mut_rate_file}")
    print(f"  - {combined_file}")
    print(f"  - {rank_abundance_file}")
    for mueller_file in mueller_files:
        print(f"  - {mueller_file}")
    print(f"  - {args.overall_title.replace(' ', '_')}_mutation_rates_summary.txt")
    if mutated_seq_summary_file:
        print(f"  - {mutated_seq_summary_file}")
    
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()