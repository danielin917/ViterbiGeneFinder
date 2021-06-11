/*
 * This file contains the implementation of the Viterbi class.
 *
 * danielin@uw.edu
 *
 */

#include <deque>
#include <limits>
#include <iostream>
#include <string>
#include <vector>

#include "viterbi.h"

using namespace std;

static double kNegativeMaxDouble = -1 * numeric_limits<double>::max();

namespace viterbi {

//-----------------------------------------------------------------------------

Viterbi::Viterbi(const TransitionMatrix& transition_matrix,
                 const EmissionMatrix& emission_matrix) :
  transition_matrix_(transition_matrix),
  emission_matrix_(emission_matrix) {

}

//-----------------------------------------------------------------------------

deque<int> Viterbi::ViterbiTrain(const vector<char>& emission_sequence,
                                 const int num_iterations) {

  deque<int> probable_path;
  for (int ii = 0; ii < num_iterations - 1; ++ii) {
    probable_path = MostProbableStatePath(emission_sequence);

    PrintHits(probable_path, 5 /* max_hits */, emission_sequence);
    // We should just have one more entry for the initial state.
    assert(probable_path.size() == (emission_sequence.size() + 1));
    assert(probable_path[0] == 0);

    // We are not training the the begin state transition so kick it out.
    probable_path.pop_front();

    // S x S matrix that will count how many times we expect that we transfered
    // between 2 states.
    vector<vector<int>> transition_count_matrix(
      transition_matrix_.size(),
      vector<int>(transition_matrix_[0].size(), 0));

    // S x A matrix that will count how many times we emitted each value in the
    // alphabet in each particular state.
    vector<vector<int>> emission_count_matrix(
      emission_matrix_.size(), vector<int>(emission_matrix_[0].size(), 0));

    ++emission_count_matrix[probable_path[0]][emission_sequence[0]];
    for (int ii = 1; ii < emission_sequence.size(); ++ii) {
      ++emission_count_matrix[probable_path[ii]][emission_sequence[ii]];
      ++transition_count_matrix[probable_path[ii - 1]][probable_path[ii]];
    }
    FillLogProbabilityMatrix(emission_count_matrix, &emission_matrix_);
    FillLogProbabilityMatrix(transition_count_matrix, &transition_matrix_);
  }

  // Do final print with all hits.
  probable_path = MostProbableStatePath(emission_sequence);
  PrintHits(probable_path, numeric_limits<int>::max(), emission_sequence);
  return probable_path;
}

//-----------------------------------------------------------------------------

void Viterbi::PrintHits(const deque<int>& probable_path,
                        const int max_hits,
                        const vector<char>& emission_sequence) {

  bool in_state_2 = false;
  int num_hits = 0;
  // We start after the beginning sequence.
  int start = -1;
  for (int ii = 1; ii < probable_path.size(); ++ii) {
    int next_state = probable_path[ii];
    if (next_state == 2 && !in_state_2) {
      if (num_hits < max_hits) {
        start = ii;
      }
      ++num_hits;
      in_state_2 = true;
    } else if (next_state != 2 && in_state_2) {
      if (start != -1) {
        cout << "Start: " << start << " End: " << ii
             << " Length: " << ii - start << endl;
         start = -1;
      }
      in_state_2 = false;
    }
  }

  // If we end in the island we must print this if needed.
  if (start != -1) {
    cout << "Start: " << start << " End: " << probable_path.size()
         << " Length: " << probable_path.size() - start << endl;
     start = -1;
  }
  cout << "Number of hits: " << num_hits << endl;
}

//-----------------------------------------------------------------------------

void Viterbi::FillLogProbabilityMatrix(
  const vector<vector<int>>& count_matrix,
  vector<vector<double>> *const log_prob_matrix) {

  // We start at 1 since we are not touching the begin state row.
  for (int ii = 1; ii < count_matrix.size(); ++ii) {
    double row_sum = 0;
    for (int kk = 0; kk < count_matrix[ii].size(); ++kk) {
      row_sum += count_matrix[ii][kk];
    }

    for (int kk = 0; kk < count_matrix[ii].size(); ++kk) {
      (*log_prob_matrix)[ii][kk] = log((1.0 * count_matrix[ii][kk]) / row_sum);
    }
  }
}

//-----------------------------------------------------------------------------

deque<int> Viterbi::MostProbableStatePath(
  const vector<char>& emission_sequence) {

  // Print the current transition_matrix_ and emission_matrix_;
  cout << "Transition Matrix" << endl;
  for (int ii = 0; ii < transition_matrix_.size(); ++ii) {
    for (int kk = 0; kk < transition_matrix_[ii].size(); ++kk) {
      string prob_string = to_string(exp(transition_matrix_[ii][kk]));
      if (prob_string.size() > 10) {
        prob_string = prob_string.substr(0, 10);
      } else if( prob_string.size() < 10) {
        prob_string = string(10 - prob_string.size(), ' ') + prob_string;
      }
      cout << prob_string << " ";
    }
    cout << endl;
  }

  cout << "Emission Matrix" << endl;
  for (int ii = 0; ii < emission_matrix_.size(); ++ii) {
    for (int kk = 0; kk < emission_matrix_[ii].size(); ++kk) {
      string prob_string = to_string(exp(emission_matrix_[ii][kk]));
      if (prob_string.size() > 10) {
        prob_string = prob_string.substr(0, 10);
      } else if(prob_string.size() < 10) {
        prob_string = string(10 - prob_string.size(), ' ') + prob_string;
      }
      cout << prob_string << " ";
    }
    cout << endl;
  }

  const int num_states = transition_matrix_.size();

  // Sequence length x State table of probability of each state at a given
  // emission index. One extra column is added to front to represent initial
  // state.
  vector<vector<double>> state_emission_prob_table(
    num_states,
    vector<double>(emission_sequence.size() + 1, kNegativeMaxDouble));

  // Corresponding backtrace table.
  vector<vector<int>> backtrace_table(
    num_states,
    vector<int>(emission_sequence.size(), -1));

  // Initial state probability is 1 at beginning so log(1) = 0;
  state_emission_prob_table[0][0] = 0;

  for (int ii = 0; ii < emission_sequence.size(); ++ii) {
    // We can never transition into begin state so don't fill begin state row.
    for (int kk = 1; kk < num_states; ++kk) {
      int prev_state = -1;
      double max_log_prev_probability = -1;
      tie(prev_state, max_log_prev_probability) =
        MaxTransitionStateProbability(state_emission_prob_table,
                                      kk /* current_state */,
                                      ii /* sequence_index */);

      assert(emission_sequence[ii] < emission_matrix_[kk].size());

      // Add the emission log probability to the probability log probability
      // for this state at the current index.
      state_emission_prob_table[kk][ii + 1] =
        max_log_prev_probability + emission_matrix_[kk][emission_sequence[ii]];
      backtrace_table[kk][ii]= prev_state;
    }
  }

  // Now find the highest probability end state.
  int end_state = -1;
  double best_log_prob = kNegativeMaxDouble;
  for (int ii = 0; ii < num_states; ++ii) {
    if (state_emission_prob_table[ii].back() >
          best_log_prob) {
      end_state = ii;
      best_log_prob = state_emission_prob_table[ii].back();
    }
  }

  // Now do backtrace.
  deque<int> probable_path;
  probable_path.push_front(end_state);
  int next_state = end_state;
  double log_prob_sum = 0;
  for (int ii = backtrace_table[0].size() - 1; ii >= 0; --ii) {
    log_prob_sum += state_emission_prob_table[next_state][ii + 1];
    next_state = backtrace_table[next_state][ii];
    probable_path.push_front(next_state);
  }
  // Assert that the first state in backtrace is the begin state.
  assert(next_state == 0);
  cout << "Log probability sum: " << log_prob_sum << endl;
  return probable_path;
}

//-----------------------------------------------------------------------------

pair<int, double> Viterbi::MaxTransitionStateProbability(
  const vector<vector<double>>& state_emission_prob_table,
  const int current_state,
  const int sequence_index) {
  assert(current_state >= 0);
  assert(sequence_index >= 0);
  assert(sequence_index < state_emission_prob_table[0].size());
  assert(current_state < transition_matrix_[0].size());

  const int num_states = transition_matrix_.size();
  assert(!transition_matrix_.empty());

  int prev_state = -1;

  // Must use negative max as min() will not be a negative number.
  double max_log_probability = kNegativeMaxDouble;
  for (int ii = 0; ii < num_states; ++ii) {
    if (ii == 0 && sequence_index > 0) {
      // No need to check transition from begin state after first iteration.
      continue;
    }
    // The sequence index is one ahead of the probability table because of the
    // initial state. So below we are checking the probability of the previous
    // state not the current state.
    const double log_probability =
      state_emission_prob_table[ii][sequence_index] +
        transition_matrix_[ii][current_state];
    if (log_probability > max_log_probability) {
      max_log_probability = log_probability;
      prev_state = ii;
    }
  }
  assert(prev_state > -1);
  return pair<int, double>(prev_state, max_log_probability);
}

//-----------------------------------------------------------------------------

} // namespace.
