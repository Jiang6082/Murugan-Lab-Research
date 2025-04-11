#!/usr/bin/env python3

import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import seaborn as sns
from concurrent.futures import ProcessPoolExecutor
import glob
import gzip
import re
from Bio import Align

class FastqMutationAnalyzer:
    def __init__(self, reference_sequence=None):
        """Initialize the analyzer with an optional reference sequence."""
        self.reference_sequence = reference_sequence
        self.mutations_per_read = []
        self.mutation_types = []
        self.mutation_positions = []
        self.read_lengths = []
        self.reverse_complemented = []  # Track which reads were reverse complemented
        self.revcomp_count = 0  # Count how many reads were reverse complemented
        self.filtered_count = 0  # Count how many reads were filtered out
        
    def set_reference(self, reference_sequence):
        """Set the reference sequence for alignment."""
        self.reference_sequence = reference_sequence
    
    def detect_file_type(self, file_path):
        """Detect if file is a FASTQ or txt file containing FASTQ data."""
        is_gzipped = file_path.endswith('.gz')
        
        # Check file extension
        if file_path.lower().endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
            return 'fastq', is_gzipped
        elif file_path.lower().endswith('.txt'):
            return 'txt', False
        else:
            # Try to peek at the file to determine format
            opener = gzip.open if is_gzipped else open
            mode = "rt" if is_gzipped else "r"
            
            try:
                with opener(file_path, mode) as f:
                    # Read first few lines to check format
                    lines = [next(f) for _ in range(4)]
                    
                    # Check if it looks like FASTQ format
                    if lines[0].startswith('@') and lines[2].startswith('+'):
                        return 'fastq', is_gzipped
                    else:
                        # See if it's txt with FASTQ data
                        if any(line.startswith('@') for line in lines):
                            return 'txt', False
            except:
                pass
                
            # Default to txt if can't determine
            return 'txt', is_gzipped
    
    def parse_txt_as_fastq(self, file_path, is_gzipped=False):
        """Parse a txt file that contains FASTQ-like data."""
        opener = gzip.open if is_gzipped else open
        records = []
        
        with opener(file_path, "rt") as handle:
            lines = handle.readlines()
            
        i = 0
        while i < len(lines):
            # Look for a FASTQ record pattern
            if i+3 < len(lines) and lines[i].startswith('@') and lines[i+2].startswith('+'):
                # This looks like a FASTQ record
                header = lines[i].strip()
                seq = lines[i+1].strip()
                qual_header = lines[i+2].strip()
                qual = lines[i+3].strip()
                
                # Create SeqRecord object
                record = SeqRecord(
                    Seq(seq),
                    id=header[1:].split()[0],  # Remove @ and take first part as ID
                    name=header[1:].split()[0],
                    description=header[1:],
                    letter_annotations={"phred_quality": [ord(c)-33 for c in qual]}
                )
                records.append(record)
                i += 4
            else:
                # Try to find the next record
                i += 1
                
        return records
    
    def extract_reference_from_file(self, file_path, method="consensus"):
        """Extract a reference sequence from the file, which could be FASTQ or txt."""
        file_type, is_gzipped = self.detect_file_type(file_path)
        reads = []
        
        # Get reads based on file type
        if file_type == 'fastq':
            opener = gzip.open if is_gzipped else open
            with opener(file_path, "rt") as handle:
                if method == "first":
                    # Use the first read as reference
                    for record in SeqIO.parse(handle, "fastq"):
                        self.reference_sequence = str(record.seq)
                        print(f"Reference sequence set from first read ({len(self.reference_sequence)} bp)")
                        return self.reference_sequence
                
                elif method == "longest" or method == "consensus":
                    # Need to collect reads
                    count = 0
                    for record in SeqIO.parse(handle, "fastq"):
                        reads.append(str(record.seq))
                        count += 1
                        if method == "consensus" and count >= 100:
                            break
        else:  # txt file
            records = self.parse_txt_as_fastq(file_path, is_gzipped)
            if method == "first" and records:
                self.reference_sequence = str(records[0].seq)
                print(f"Reference sequence set from first read ({len(self.reference_sequence)} bp)")
                return self.reference_sequence
            elif method == "longest" or method == "consensus":
                reads = [str(record.seq) for record in records]
                if method == "consensus":
                    reads = reads[:1000]  # Limit to first 1000 for consensus
        
        # Process reads based on method
        if method == "longest" and reads:
            self.reference_sequence = max(reads, key=len)
            print(f"Reference sequence set from longest read ({len(self.reference_sequence)} bp)")
            return self.reference_sequence
            
        elif method == "consensus" and reads:
            # Find the most common length
            lengths = Counter([len(read) for read in reads])
            most_common_length = lengths.most_common(1)[0][0]
            
            # Keep only reads of the most common length
            filtered_reads = [read for read in reads if len(read) == most_common_length]
            
            # Build consensus
            consensus = ""
            for i in range(most_common_length):
                bases = Counter([read[i] for read in filtered_reads if i < len(read)])
                if bases:
                    consensus += bases.most_common(1)[0][0]
                else:
                    consensus += "N"  # Default if no data
            
            self.reference_sequence = consensus
            print(f"Reference sequence set using consensus ({len(self.reference_sequence)} bp)")
            # print(f"The reference sequence is: {self.reference_sequence}")
            return self.reference_sequence
        
        raise ValueError(f"Could not extract reference sequence from {file_path} using method {method}")
    
    def get_reverse_complement(self, seq):
        """Generate reverse complement of a sequence."""
        complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                           'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                           'N': 'N', 'n': 'n'}
        
        # Generate complement and reverse
        complement = ''.join([complement_dict.get(base, base) for base in seq])
        reverse_complement = complement[::-1]
        
        return reverse_complement
    
    def analyze_file(self, file_path, max_reads=None, try_revcomp=False, mutation_threshold=100):
        """Analyze mutations in a file (FASTQ or txt containing FASTQ data)."""
        if not self.reference_sequence:
            self.extract_reference_from_file(file_path)
            
        read_count = 0
        processed_count = 0
        file_type, is_gzipped = self.detect_file_type(file_path)
        
        print(f"Analyzing {file_path} as {file_type} file...")
        
        if file_type == 'fastq':
            opener = gzip.open if is_gzipped else open
            with opener(file_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    # Analyze the read with appropriate mode
                    processed = self._analyze_read(record, try_revcomp, mutation_threshold)
                    
                    read_count += 1
                    if processed:
                        processed_count += 1
                    
                    if read_count % 1000 == 0:
                        print(f"Processed {read_count} reads, kept {processed_count}")
                    
                    if max_reads and read_count >= max_reads:
                        break
        else:  # txt file
            records = self.parse_txt_as_fastq(file_path, is_gzipped)
            for record in records:
                # Analyze the read with appropriate mode
                processed = self._analyze_read(record, try_revcomp, mutation_threshold)
                
                read_count += 1
                if processed:
                    processed_count += 1
                
                if read_count % 1000 == 0:
                    print(f"Processed {read_count} reads, kept {processed_count}")
                
                if max_reads and read_count >= max_reads:
                    break
        
        filtered_count = read_count - processed_count        
        print(f"Completed analysis of {read_count} reads from {file_path}")
        print(f"Kept {processed_count} reads, filtered out {filtered_count} reads")
        
        if try_revcomp:
            print(f"Used reverse complement for {self.revcomp_count} reads that exceeded mutation threshold of {mutation_threshold}")
        
    def analyze_multiple_files(self, file_paths, max_reads_per_file=None, try_revcomp=False, mutation_threshold=100):
        """Analyze mutations in multiple files."""
        for file_path in file_paths:
            if os.path.exists(file_path):
                self.analyze_file(file_path, max_reads=max_reads_per_file, 
                                 try_revcomp=try_revcomp, mutation_threshold=mutation_threshold)
            else:
                print(f"Warning: File {file_path} not found!")
    
    def _perform_alignment(self, read_seq):
        """Perform alignment and return mutation count and details"""
        try:
            # Create a pairwise aligner with similar scoring to the original
            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2.0
            aligner.mismatch_score = -1.0
            aligner.open_gap_score = -2.0
            aligner.extend_gap_score = -0.5
            
            # Perform alignment
            alignments = aligner.align(self.reference_sequence, read_seq)
            
            # Use the best alignment
            if len(alignments) > 0:
                best_alignment = alignments[0]
                
                # Extract alignment operations and build gapped sequences
                ref_aligned = ""
                read_aligned = ""
                
                # Get alignment coordinates
                ref_start, ref_end = best_alignment.aligned[0][0][0], best_alignment.aligned[0][-1][1]
                read_start, read_end = best_alignment.aligned[1][0][0], best_alignment.aligned[1][-1][1]
                
                # Add any leading gaps
                if ref_start > 0:
                    read_aligned += read_seq[:ref_start]
                    ref_aligned += "-" * ref_start
                    
                if read_start > 0:
                    ref_aligned += self.reference_sequence[:read_start]
                    read_aligned += "-" * read_start
                    
                # Process the aligned blocks
                ref_pos = max(0, read_start)
                read_pos = max(0, ref_start)
                
                # Build aligned sequences based on the cigar-like operations
                for i, (ref_block, read_block) in enumerate(zip(best_alignment.aligned[0], best_alignment.aligned[1])):
                    # Handle any gaps between blocks
                    if ref_pos < ref_block[0]:
                        ref_aligned += self.reference_sequence[ref_pos:ref_block[0]]
                        read_aligned += "-" * (ref_block[0] - ref_pos)
                        ref_pos = ref_block[0]
                    
                    if read_pos < read_block[0]:
                        read_aligned += read_seq[read_pos:read_block[0]]
                        ref_aligned += "-" * (read_block[0] - read_pos)
                        read_pos = read_block[0]
                    
                    # Add aligned segment
                    ref_aligned += self.reference_sequence[ref_block[0]:ref_block[1]]
                    read_aligned += read_seq[read_block[0]:read_block[1]]
                    ref_pos = ref_block[1]
                    read_pos = read_block[1]
                
                # Add any trailing sequence
                if ref_pos < len(self.reference_sequence):
                    ref_aligned += self.reference_sequence[ref_pos:]
                    read_aligned += "-" * (len(self.reference_sequence) - ref_pos)
                
                if read_pos < len(read_seq):
                    read_aligned += read_seq[read_pos:]
                    ref_aligned += "-" * (len(read_seq) - read_pos)
                
                # Verify equal lengths before analyzing mutations
                if len(ref_aligned) != len(read_aligned):
                    # print(f"Warning: Aligned sequences have different lengths: ref={len(ref_aligned)}, read={len(read_aligned)}")
                    # Truncate to shorter length if needed
                    min_len = min(len(ref_aligned), len(read_aligned))
                    ref_aligned = ref_aligned[:min_len]
                    read_aligned = read_aligned[:min_len]
                
                # Collect mutation data
                mutations = 0
                mutation_types = []
                mutation_positions = []
                
                for i in range(len(ref_aligned)):
                    if ref_aligned[i] != read_aligned[i] and ref_aligned[i] != '-' and read_aligned[i] != '-':
                        mutations += 1
                        # Record mutation type (substitution)
                        mutation_types.append(f"{ref_aligned[i]}>{read_aligned[i]}")
                        # Estimate position in reference
                        ref_pos = len(ref_aligned[:i].replace('-', ''))
                        mutation_positions.append(ref_pos)
                    
                    elif ref_aligned[i] == '-':  # Insertion in read
                        mutations += 1
                        mutation_types.append(f"ins_{read_aligned[i]}")
                        # Position is the previous reference base
                        ref_pos = len(ref_aligned[:i].replace('-', ''))
                        mutation_positions.append(ref_pos)
                    
                    elif read_aligned[i] == '-':  # Deletion in read
                        mutations += 1
                        mutation_types.append(f"del_{ref_aligned[i]}")
                        # Position is this reference base
                        ref_pos = len(ref_aligned[:i].replace('-', ''))
                        mutation_positions.append(ref_pos)
                
                return mutations, mutation_types, mutation_positions
            
            # If no alignment found
            return 0, [], []
        
        except Exception as e:
            print(f"Error during alignment: {e}")
            return 0, [], []
    
    def _analyze_read(self, record, try_revcomp=False, mutation_threshold=100):
        """
        Analyze a single read for mutations using PairwiseAligner.
        
        Returns:
            bool: True if the read was processed and kept, False if filtered
        """
        read_seq = str(record.seq)
        
        # Skip if reference is not set or read is empty
        if not self.reference_sequence or not read_seq:
            return False
        
        # First, try the original sequence
        mutations, mutation_types, mutation_positions = self._perform_alignment(read_seq)
        
        # Handle based on mode and mutation count
        if try_revcomp and mutations > mutation_threshold:
            # Try reverse complement mode
            revcomp_seq = self.get_reverse_complement(read_seq)
            revcomp_mutations, revcomp_types, revcomp_positions = self._perform_alignment(revcomp_seq)
            
            # If reverse complement has fewer mutations, use it instead
            if revcomp_mutations < mutations:
                mutations = revcomp_mutations
                mutation_types = revcomp_types
                mutation_positions = revcomp_positions
                used_revcomp = True
                self.revcomp_count += 1
            else:
                used_revcomp = False
                
            # Store results
            self.read_lengths.append(len(read_seq))
            self.mutations_per_read.append(mutations)
            self.mutation_types.extend(mutation_types)
            self.mutation_positions.extend(mutation_positions)
            self.reverse_complemented.append(used_revcomp)
            return True
            
        elif not try_revcomp and mutations > mutation_threshold:
            # Filter mode - discard this read
            self.filtered_count += 1
            return False
        
        else:
            # Normal mode or mutations below threshold
            self.read_lengths.append(len(read_seq))
            self.mutations_per_read.append(mutations)
            self.mutation_types.extend(mutation_types)
            self.mutation_positions.extend(mutation_positions)
            self.reverse_complemented.append(False)
            return True
        
    def plot_results(self, output_dir='.'):
        """Generate plots of the analysis results."""
        os.makedirs(output_dir, exist_ok=True)
        
        # Create a DataFrame for easier analysis
        data = {
            'mutations_per_read': self.mutations_per_read,
            'read_length': self.read_lengths,
            'reverse_complemented': self.reverse_complemented
        }
        df = pd.DataFrame(data)
        
        # Plot 1: Histogram of mutations per read
        plt.figure(figsize=(10, 6))
        bin_size = 1  # Change to 5 if 1 is too granular
        
        # Calculate the bins based on the range of data
        max_mutations = max(df['mutations_per_read']) if len(df['mutations_per_read']) > 0 else 10
        min_mutations = min(df['mutations_per_read']) if len(df['mutations_per_read']) > 0 else 0
        bins = range(int(min_mutations), 100, bin_size)
        
        # Create the histogram with specific bins
        sns.histplot(df['mutations_per_read'], bins=bins, kde=True)
        plt.title('Distribution of Mutations per Read')
        plt.xlim(0, 25)
        plt.xlabel('Number of Mutations')
        plt.ylabel('Frequency')
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(output_dir, 'mutations_per_read_histogram.png'))
        plt.close()
        
        # Plot 2: Mutation types
        if self.mutation_types:
            plt.figure(figsize=(12, 8))
            mutation_counts = Counter(self.mutation_types)
            
            # Sort by frequency
            labels, values = zip(*sorted(mutation_counts.items(), key=lambda x: x[1], reverse=True))
            
            # If there are too many categories, limit to top N
            if len(labels) > 20:
                labels = labels[:20]
                values = values[:20]
                plt.title('Top 20 Mutation Types')
            else:
                plt.title('Mutation Types')
                
            plt.bar(labels, values)
            plt.xlabel('Mutation Type')
            plt.ylabel('Frequency')
            plt.xticks(rotation=90)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'mutation_types_barplot.png'))
            plt.close()
        
        # Plot 3: Mutation positions along reference
        if self.mutation_positions:
            plt.figure(figsize=(12, 6))
            sns.histplot(self.mutation_positions, bins=100, kde=True)
            plt.title('Mutation Positions along Reference Sequence')
            plt.xlabel('Position in Reference')
            plt.ylabel('Frequency')
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(output_dir, 'mutation_positions_histogram.png'))
            plt.close()
        
        # Plot 4: Mutations vs read length
        # plt.figure(figsize=(10, 6))
        # sns.scatterplot(data=df, x='read_length', y='mutations_per_read', 
        #                hue='reverse_complemented' if any(df['reverse_complemented']) else None,
        #                alpha=0.5)
        # plt.title('Mutations vs Read Length')
        # plt.xlabel('Read Length')
        # plt.ylabel('Number of Mutations')
        # plt.grid(True, alpha=0.3)
        # if any(df['reverse_complemented']):
        #     plt.legend(title='Reverse Complemented')
        # plt.savefig(os.path.join(output_dir, 'mutations_vs_length_scatter.png'))
        # plt.close()
        
        # Plot 5: Mutation rate per read
        # Avoid division by zero
        # df['mutation_rate'] = df.apply(
        #     lambda row: row['mutations_per_read'] / row['read_length'] if row['read_length'] > 0 else 0, 
        #     axis=1
        # )
        
        # plt.figure(figsize=(10, 6))
        # if any(df['reverse_complemented']):
        #     # Split by whether reverse complemented or not
        #     sns.histplot(data=df, x='mutation_rate', hue='reverse_complemented', kde=True, common_norm=False)
        #     plt.title('Mutation Rate per Read (by Orientation)')
        # else:
        #     sns.histplot(df['mutation_rate'], kde=True)
        #     plt.title('Mutation Rate per Read')
        # plt.xlabel('Mutation Rate (mutations/bp)')
        # plt.ylabel('Frequency')
        # plt.grid(True, alpha=0.3)
        # plt.savefig(os.path.join(output_dir, 'mutation_rate_histogram.png'))
        # plt.close()

        # Plot 6: Mutation matrix (substitutions as 4x4 matrix + indels)
        if self.mutation_types:
            plt.figure(figsize=(10, 8))
            
            # Initialize the substitution matrix
            bases = ['A', 'T', 'G', 'C']
            substitution_matrix = np.zeros((4, 4), dtype=int)  # Use int dtype here
            insertions = {"A": 0, "T": 0, "G": 0, "C": 0, "total": 0}
            deletions = {"A": 0, "T": 0, "G": 0, "C": 0, "total": 0}
            
            # Parse mutation types to fill the matrix
            for mutation in self.mutation_types:
                if '>' in mutation:  # Substitution
                    try:
                        ref, alt = mutation.split('>')
                        if ref in bases and alt in bases:
                            i = bases.index(ref)
                            j = bases.index(alt)
                            substitution_matrix[i, j] += 1
                    except:
                        continue
                elif mutation.startswith('ins_'):  # Insertion
                    base = mutation[4:]
                    if base in bases:
                        insertions[base] += 1
                    insertions["total"] += 1
                elif mutation.startswith('del_'):  # Deletion
                    base = mutation[4:]
                    if base in bases:
                        deletions[base] += 1
                    deletions["total"] += 1
            
            # Create the heatmap for substitutions
            ax = plt.subplot(1, 2, 1)
            sns.heatmap(substitution_matrix, annot=True, fmt="d", cmap="YlGnBu",
                        xticklabels=bases, yticklabels=bases, cbar_kws={'label': 'Count'})
            plt.title('Nucleotide Substitutions')
            plt.xlabel('To')
            plt.ylabel('From')
            
            # Create the barplot for indels
            ax2 = plt.subplot(1, 2, 2)
            indel_data = {
                'A ins': insertions.get('A', 0),
                'T ins': insertions.get('T', 0),
                'G ins': insertions.get('G', 0),
                'C ins': insertions.get('C', 0),
                'A del': deletions.get('A', 0),
                'T del': deletions.get('T', 0),
                'G del': deletions.get('G', 0),
                'C del': deletions.get('C', 0)
            }
            
            # Create bar colors - blue for insertions, red for deletions
            colors = ['blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red']
            
            # Create the bar plot
            bars = plt.bar(indel_data.keys(), indel_data.values(), color=colors)
            plt.title('Insertions and Deletions')
            plt.xticks(rotation=45)
            plt.ylabel('Count')
            
            # Add totals as text
            plt.figtext(0.5, 0.02, f"Total insertions: {insertions['total']} | Total deletions: {deletions['total']}", 
                        ha='center', fontsize=10, bbox={"facecolor":"white", "alpha":0.5, "pad":5})
            
            plt.tight_layout(rect=[0, 0.05, 1, 0.95])  # Adjust layout to make room for the text
            plt.savefig(os.path.join(output_dir, 'mutation_matrix.png'))
            plt.close()       
        
        # Additional plot for reverse complement if used
        if any(df['reverse_complemented']):
            # Comparison of mutations before/after revcomp
            plt.figure(figsize=(10, 6))
            revcomp_stats = {
                'Original': len(df[~df['reverse_complemented']]),
                'Reverse Complemented': len(df[df['reverse_complemented']])
            }
            plt.bar(revcomp_stats.keys(), revcomp_stats.values())
            plt.title('Number of Reads Using Original vs Reverse Complement')
            plt.ylabel('Count')
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(output_dir, 'reverse_complement_usage.png'))
            plt.close()
        
        # Save summary statistics
        stats = {
            'Total reads analyzed': len(self.mutations_per_read),
            'Reads filtered out': self.filtered_count,
            'Average mutations per read': np.mean(self.mutations_per_read) if self.mutations_per_read else 0,
            'Median mutations per read': np.median(self.mutations_per_read) if self.mutations_per_read else 0,
            'Max mutations in a read': max(self.mutations_per_read) if self.mutations_per_read else 0,
            'Average read length': np.mean(self.read_lengths) if self.read_lengths else 0,
            'Total mutations identified': len(self.mutation_types),
            'Unique mutation types': len(set(self.mutation_types)),
            'Reads using reverse complement': sum(self.reverse_complemented)
        }
        
        with open(os.path.join(output_dir, 'mutation_stats_summary.txt'), 'w') as f:
            for key, value in stats.items():
                f.write(f"{key}: {value}\n")
        
        # Save detailed results to CSV
        # df.to_csv(os.path.join(output_dir, 'mutation_data.csv'), index=False)
        
        # Save mutation types counts
        # if self.mutation_types:
        #     pd.DataFrame(Counter(self.mutation_types).most_common(), 
        #                 columns=['mutation_type', 'count']).to_csv(
        #                     os.path.join(output_dir, 'mutation_types.csv'), index=False)
        
        print(f"Results saved to {output_dir}")
        return stats

def main():
    parser = argparse.ArgumentParser(description='Analyze mutations in FASTQ files or text files containing FASTQ data')
    parser.add_argument('input_files', nargs='+', help='Input file(s) to analyze. Can be FASTQ or text files. Supports glob patterns.')
    parser.add_argument('--reference', '-r', help='Reference sequence (file in FASTA format or string)')
    parser.add_argument('--output', '-o', default='mutation_analysis_results', help='Output directory')
    parser.add_argument('--max-reads', '-m', type=int, default=None, help='Maximum number of reads to analyze per file')
    parser.add_argument('--extract-ref', choices=['first', 'consensus', 'longest'], default='consensus',
                        help='Method to extract reference if not provided')
    parser.add_argument('--try-revcomp', action='store_true', 
                        help='Try reverse complement for highly mutated sequences')
    parser.add_argument('--mutation-threshold', type=int, default=100, 
                        help='Mutation threshold: with --try-revcomp, sequences with more mutations than this will be tried as reverse complement; without --try-revcomp, sequences with more mutations will be filtered out')
    
    args = parser.parse_args()
    
    # Expand glob patterns in input_files
    expanded_files = []
    for pattern in args.input_files:
        expanded = glob.glob(pattern)
        if expanded:
            expanded_files.extend(expanded)
        else:
            print(f"Warning: No files match pattern '{pattern}'")
    
    if not expanded_files:
        print("Error: No valid input files provided")
        return
    
    # Create the analyzer
    analyzer = FastqMutationAnalyzer()
    
    # Set the reference sequence
    if args.reference:
        if os.path.exists(args.reference):
            # Check if reference file is gzipped
            is_gzipped = args.reference.endswith('.gz')
            opener = gzip.open if is_gzipped else open
            # Read from FASTA file
            with opener(args.reference, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    analyzer.set_reference(str(record.seq))
                    print(f"Reference sequence set from file ({len(analyzer.reference_sequence)} bp)")
                    break
        else:
            # Use the provided string as reference
            analyzer.set_reference(args.reference)
            print(f"Reference sequence set from input ({len(analyzer.reference_sequence)} bp)")
    else:
        # Will extract reference from first file
        print(f"No reference provided. Will extract from first file using '{args.extract_ref}' method.")
    
    # Print mode information
    if args.try_revcomp:
        print(f"Mode: Try reverse complement for sequences with >{args.mutation_threshold} mutations")
    else:
        print(f"Mode: Filter out sequences with >{args.mutation_threshold} mutations")
    
    # Analyze the files
    analyzer.analyze_multiple_files(expanded_files, args.max_reads, 
                                   try_revcomp=args.try_revcomp, 
                                   mutation_threshold=args.mutation_threshold)
    
    # Plot and save results
    stats = analyzer.plot_results(args.output)
    
    # Print summary
    print("\nAnalysis Summary:")
    for key, value in stats.items():
        print(f"{key}: {value}")

if __name__ == "__main__":
    main()