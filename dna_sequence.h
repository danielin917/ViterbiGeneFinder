/*
 *
 *
 */

#ifndef DNA_SEQUENCE_H_
#define DNA_SEQUENCE_H_

#include <string>

//-----------------------------------------------------------------------------

enum class Nucleotide {
  kA = 0 /* Adenine */,
  kC = 1 /* Cytosine */,
  kG = 2 /* Guanine */,
  kT = 3 /* Thymine */,
  kNumNucleotides = 4
};

//-----------------------------------------------------------------------------

Nucleotide CharToNucleotide(const char c);
char NucleotideToChar(Nucleotide n);

//-----------------------------------------------------------------------------

class DNASequence {
 public:
  // Constructor.
  DNASequence();

  // Construct sequence based on given 'sequence' and 'identifier'.
  DNASequence(const std::string& sequence, const std::string& identifier);

  // Accesors.
  const std::string& sequence() const { return sequence_; }
  const std::string& identifier() const { return identifier_; }

  // Prints the DNA sequence out in columns.
  void PrintWithColumns() const;

  // Substitutes all unknown characters in 'line' with a T.
  static std::string MakeValidDNASequence(const std::string& line);

  // Returns true if 'sequence' is a a valid DNA sequence.
  static bool IsValidDNASequence(const std::string& sequence);

  // Whether 'c' translates to a valid nucleotide {A, C, G, T}.
  static bool IsValidNucleotide(char c);

 private:
  // Name or identifier for this sequence.
  std::string identifier_;

  // Vector containing a sequence of nucleotides.
  std::string sequence_;
};

//-----------------------------------------------------------------------------

#endif
