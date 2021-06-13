/*
 * This file describes the Viterbi class which finds the most probable path
 * of states based on a given emission and transition matrix.
 *
 * Author: danielin@uw.edu
 *
 */

#ifndef VITERBI_H_
#define VITERBI_H_

#include <deque>
#include <vector>

#include "dna_sequence.h"

namespace viterbi {

class Viterbi {
 public:
  // S x S matrix of transition probabilities between each state.
  typedef std::vector<std::vector<double>> TransitionMatrix;

  // S x A matrix of the probability of emitting each character in an alphabet
  // in each state.
  typedef std::vector<std::vector<double>> EmissionMatrix;

  // Constructor. Builds the  Viterbi class using the provided
  // 'transition_matrix' and 'emission_matrix'.
  Viterbi(const TransitionMatrix& transition_matrix,
          const EmissionMatrix& emission_matrix);

  // Perform 'num_iterations' of viterbi training using the given
  // 'emission_sequence'.
  std::deque<int> ViterbiTrain(const std::vector<char>& emission_sequence,
                               int num_iterations);

  // Print the island "hits" for State 2 in 'probable_path'. Print the location
  // at length of the first 'max_hits' islands.
  void PrintHits(const std::deque<int>& probable_path,
                 const int max_hits,
                 const std::vector<char>& emission_sequence);

  // Fill 'log_prob_matrix' based on frquency of the counts in 'count_matrix'
  // for each row.
  void FillLogProbabilityMatrix(
    const std::vector<std::vector<int>>& count_matrix,
    std::vector<std::vector<double>> *log_prob_matrix);

  // Return the sequence of states corresponding to the most probable
  // corresponding states to emit 'emission_sequence'.
  std::deque<int> MostProbableStatePath(
    const std::vector<char>& emission_sequence);

  // Given the 'current_state' and the index we looking at in the sequence use
  // 'state_emission_prob_table' to get the most probable previous state.
  // The pair returned is of form <State Number, log(probability)>.
  inline std::pair<int, double> MaxTransitionStateProbability(
    const std::vector<std::vector<double>>& state_emission_prob_table,
    int current_state,
    int sequence_index);


 private:
  // S x S State transition log probability matrix.
  TransitionMatrix transition_matrix_;

  // S x A Emission log probability per state matrix.
  EmissionMatrix emission_matrix_;
};

} // namespace

#endif // VITERBI_H_
