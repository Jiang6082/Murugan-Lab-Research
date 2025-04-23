# MutationAnalyzer-041025

This is for analyzing sequences in any given fastq files/txt snippets, and outputs graphs such as mutation positions, types of mutations, and mutation counts. It also outputs a stats summary txt file for a more direct interpretation. 

Update in v2:
the revComp argument is updated alongside a new argument calledm mutation_threshold. Instead of the original modes, there are now three options to choose from: 
- try: considers the reverse complements for sequences that havea # of mutations above the threshold
- record: only records down data from the sequences with above-threshold mutations
- filter: filter out sequences with above-threshold mutations from the output.

Updates in v3:
-  the try option has been updated, such that it only considers the sequeces' reverse complements that makes sense(have a mutation ct. below mutation_threshold). It if doesn't make sense whatsoever, it is filtered out.
-  percent_identity is introduced in _perform_alignment to reduce irrelevant alignments. Now the function only count mutations if alignment quality is decent(custom).
-  Poisson plot is now included as part of the plot_results function
-  The mutation positions graph's y-axis has been limited(custom, currently at 95th percentile).
-  The file can now also take in .seq files(just pure sequences, no fastq format data).
