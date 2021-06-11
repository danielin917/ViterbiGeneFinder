/*
 *
 */

#include <iostream>

#include "dna_sequence.h"

using namespace std;

//-----------------------------------------------------------------------------

Nucleotide CharToNucleotide(const char c) {
  const char upperc = toupper(c);
  switch (upperc) {
    case 'A':
      return Nucleotide::kA;
    case 'C':
      return Nucleotide::kC;
    case 'G':
      return Nucleotide::kG;
    case 'T':
      return Nucleotide::kT;
    default:
      cerr << "Invalid nucleotide char " << c << endl;
      assert(false);
      break;
  }
  __builtin_unreachable();
}

//-----------------------------------------------------------------------------

char NucleotideToChar(Nucleotide n) {
  switch (n) {
    case Nucleotide::kA:
      return 'A';
    case Nucleotide::kC:
      return 'C';
    case Nucleotide::kG:
      return 'G';
    case Nucleotide::kT:
      return 'T';
    default:
      cerr << "Invalid nucleotide " << static_cast<int>(n) << endl;
      assert(false);
      break;
  }
  __builtin_unreachable();
}

//-----------------------------------------------------------------------------

DNASequence::DNASequence() {

}

//-----------------------------------------------------------------------------

DNASequence::DNASequence(const string& sequence, const string& identifier) :
  sequence_(sequence),
  identifier_(identifier) {
}

//-----------------------------------------------------------------------------

string DNASequence::MakeValidDNASequence(const string& line) {
  string sequence = line;
  for (int ii = 0; ii < sequence.size(); ++ii) {
    if (!IsValidNucleotide(sequence[ii])) {
      // Hack just make all invalid chars a T.
      sequence[ii] = 'T';
    }
  }
  return sequence;
}

//-----------------------------------------------------------------------------

bool DNASequence::IsValidDNASequence(const string& sequence) {
  for (int ii = 0; ii < sequence.size(); ++ii) {
    if (!IsValidNucleotide(sequence[ii])) {
      return false;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------

bool DNASequence::IsValidNucleotide(const char c) {
  const char upper_c = toupper(c);
  if (upper_c != 'A' &&
      upper_c != 'C' &&
      upper_c != 'G' &&
      upper_c != 'T') {
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

void DNASequence::PrintWithColumns() const {
  cout << identifier_ << endl;
  for (int ii = 0; ii < sequence_.size(); ++ii) {
    if (ii % 60 == 0 && ii != 0) {
      cout << endl;
    }
    cout << sequence_[ii];
  }
  cout << endl;
}

//-----------------------------------------------------------------------------
