# ViterbiGeneFinder
Build:
g++ main.cc viterbi.cc dna_sequence.cc -std=c++0x

Run:
./a.out <sequence_file>

Performs Viterbi Training to find the most probable state path based on the passed in emission
probability matrix for each state. By looking for areas in the sequence where we emit a disproporitionate
number of Cs and Gs to find "CpG islands" we can predict the start regions of genes.
