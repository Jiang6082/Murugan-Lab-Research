# MutationAnalyzer-041025

This is for analyzing sequences in any given fastq files/txt snippets, and outputs graphs such as mutation positions, types of mutations, and mutation counts. It also outputs a stats summary txt file for a more direct interpretation. 

In v2, the revComp argument is updated alongside a new argument calledm mutation_threshold. Instead of the original modes, there are now three options to choose from: 
- try: considers the reverse complements for sequences that havea # of mutations above the threshold
- record: only records down data from the sequences with above-threshold mutations
- filter: filter out sequences with above-threshold mutations from the output.
