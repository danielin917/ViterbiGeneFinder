/*
 * This file contains the driver for the viterbi path finder.
 *
 * Author: Dan Lin, danielin@uw.edu
 */

#include <fstream>
#include <iostream>
#include <vector>

#include "dna_sequence.h"
#include "viterbi.h"

using namespace viterbi;
using namespace std;

bool ExtractSequencesFromFile(const string& filename,
                              vector<DNASequence> *dna_sequences);

bool ExtractNumberSequencesFromFile(const string& filename,
                                    vector<char> *const dna_sequences);

// Below used for loaded Dice testing./////////////
Viterbi::TransitionMatrix dice_transition_matrix =
  /*            Beg                                   1           2    */
  /* Beg */ { { -1 * numeric_limits<double>::max(), log(.5),  log(0.5) },
  /*  1  */   { -1 * numeric_limits<double>::max(), log(.95), log(.05) },
  /*  2  */   { -1 * numeric_limits<double>::max(), log(.1),   log(.9) } };

Viterbi::EmissionMatrix dice_emission_matrix =
  /*              1                     2             3               4              5              6  */
  /* Beg */ { {   0,                    0,            0,              0,             0,             0,  },
  /*  1  */   { log(.166666), log(.166666), log(.166666), log(.16666666), log(.1666666), log(.1666666)  },
  /*  2  */   { log(.1),      log(.1),      log(.1),      log(.1),              log(.1),       log(.5)  } };
////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Not enough args provided." << endl
         << "Format: ./a.out <sequence_file>" << endl;
    return -1;
  }
  string test_filename(argv[1]);

  cout << "Reading input" << endl;
  vector<DNASequence> dna_sequences;
  ExtractSequencesFromFile(test_filename, &dna_sequences);
  cout << "Extracted " << dna_sequences.size() << " sequences" << endl;

  Viterbi::TransitionMatrix transition_matrix =
    /*            Beg                                   1              2    */
    /* Beg */ { { -1 * numeric_limits<double>::max(), log(.9999), log(.0001) },
    /*  1  */   { -1 * numeric_limits<double>::max(), log(.9999), log(.0001) },
    /*  2  */   { -1 * numeric_limits<double>::max(), log(.01),   log(.99) } };

 Viterbi::EmissionMatrix emission_matrix =
    /*              A            C        G         T  */
    /* Beg */ { {   0,           0,       0,        0  },
    /*  1  */   { log(.25), log(.25), log(.25), log(.25)  },
    /*  2  */   { log(.20), log(.30), log(.30), log(.20)  } };

  vector<char> emission_sequence;
  for (int ii = 0; ii < dna_sequences[0].sequence().size(); ++ii) {
    emission_sequence.push_back(
      static_cast<int>(CharToNucleotide(dna_sequences[0].sequence()[ii])));
  }
  vector<DNASequence>().swap(dna_sequences);

  //ExtractNumberSequencesFromFile(test_filename, &emission_sequence);

  Viterbi v(transition_matrix, emission_matrix);
  deque<int> viterbi_path =
    v.ViterbiTrain(emission_sequence, 10 /* num_iterations */);

// Print state path.
/*
  for (int ii = 0; ii < viterbi_path.size(); ++ii) {
    cout << viterbi_path[ii] << " ";
  }
  cout << endl;
*/
  return 0;
}

bool ExtractNumberSequencesFromFile(const string& filename,
                                    vector<char> *const emission_sequence) {

  ifstream ifs;
  ifs.open(filename.c_str(), std::ifstream::in);
  if (!ifs.good()) {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }

  string sequence;
  string identifier;
  while (true) {
    char peek_char = ifs.peek();
    if (peek_char == EOF) {
      // All input has been read.
      if (!sequence.empty()) {
        for (int ii = 0; ii < sequence.size(); ++ii) {
          emission_sequence->push_back(sequence[ii] - '0' - 1);
        }
      }
      break;
    }

    string line_str;
    getline(ifs, line_str);
    if (peek_char == '>') {
      if (!sequence.empty()) {
        for (int ii = 0; ii < sequence.size(); ++ii) {
          emission_sequence->push_back(sequence[ii] - '0' - 1);
        }
        sequence = "";
        identifier = "";
      }
      // Get the identifier set after '>' and prior to ':'.
      int end = line_str.find(":");
      identifier = line_str.substr(1, end - 1);
      continue;
    }

    sequence += line_str;
  }
  return true;
}

bool ExtractSequencesFromFile(const string& filename,
                              vector<DNASequence> *const dna_sequences) {
  assert(dna_sequences);

  ifstream ifs;
  ifs.open(filename.c_str(), std::ifstream::in);
  if (!ifs.good()) {
    cerr << "Unable to open file: " << filename << endl;
    return false;
  }

  string sequence;
  string identifier;
  while (true) {
    char peek_char = ifs.peek();
    if (peek_char == EOF) {
      // All input has been read.
      if (!sequence.empty()) {
        dna_sequences->push_back(DNASequence(sequence, identifier));
      }
      break;
    }

    string line_str;
    getline(ifs, line_str);
    if (peek_char == '>') {
      if (!sequence.empty()) {
        dna_sequences->push_back(DNASequence(sequence, identifier));
        sequence = "";
        identifier = "";
      }
      // Get the identifier set after '>' and prior to ':'.
      int end = line_str.find(":");
      identifier = line_str.substr(1, end - 1);
      continue;
    }

    line_str = DNASequence::MakeValidDNASequence(line_str);
    sequence += line_str;
  }
  return true;
}
