// *-*-C++-*-*
#ifndef ROBJECTWRITER_H
#define ROBJECTWRITER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

///class to write an R object
class RObjectWriter :public std::ofstream{

public:

  ///no-arg c'tor
  RObjectWriter();
  ///c'tor with filename, equivalent to default c'tor + open
  RObjectWriter(const char*);
  ///c'tor with filename, equivalent to default c'tor + open
  RObjectWriter(const std::string&);
  ///open file and start writing R object
  void open(const char*);
  ///open file and start writing R object
  void open(const std::string&);
  ///finish writing R object by writing dimensions and dimnames
  ///supply a vector of dimensions and a vector of labels for at least the first dimension
  void close(const std::vector<int>& dims, const std::vector< std::vector<std::string> >& dimnames);
  //write a comment, prefixed by a comment character
  RObjectWriter& comment(const char* c);
  /**
     manipulator to start new line.
     Call **before a new line** not after the old one. 
     First call starts a new line. Every subsequent call writes a comma then starts a new line.
     If you just want a new line, use endl instead.
  */
  friend void newline(RObjectWriter& R);

  friend void comma(RObjectWriter& R);

  ///stream insertion operator to accept newline manipulator
  ///this allows "R << newline;"
  friend RObjectWriter& operator<<( RObjectWriter & R, void (*man)(RObjectWriter&));

private:
  bool needNewLineComma;///< used to determine if a comma is needed before a newline
  bool needComma;///<used to determine if a comma is needed before a new element
  std::string comment_cache;

  ///overloaded close function, private to prevent closing without writing dimensions
  void close();
  ///write dimensions of R object
  void WriteDimensions(const std::vector<int>&, const std::vector< std::vector<std::string> >&);
};

//declarations of friend functions

void newline(RObjectWriter& R);
void comma(RObjectWriter& R);

std::ostream& newline(std::ostream& os);
std::ostream& comma(std::ostream& os);

RObjectWriter& operator<<( RObjectWriter & R, void (*man)(RObjectWriter&));

//define a comma followed by a space, for convenience. 
#define COMMA ", "

#endif
