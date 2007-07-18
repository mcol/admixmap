#include "bclib/StringSplitter.h"
#include <cctype>
#include <stdexcept>

BEGIN_BCLIB_NAMESPACE

//definitions of struct members of StringSplitter
struct StringSplitter::StringSplitterState
{
  virtual ~StringSplitterState(){};
  virtual StringSplitterState *onWhite(StringSplitter *_this, char c) = 0;
  virtual StringSplitterState *onGraph(StringSplitter *_this, char c) = 0;
  virtual StringSplitterState *onQuote(StringSplitter *_this) = 0;
  virtual StringSplitterState *onEndLn(StringSplitter *_this) = 0;
};

struct StringSplitter::StateWhite : public StringSplitter::StringSplitterState
{
  virtual StringSplitterState *onWhite(StringSplitter *_this, char c);
  virtual StringSplitterState *onGraph(StringSplitter *_this, char c);
  virtual StringSplitterState *onQuote(StringSplitter *_this);
  virtual StringSplitterState *onEndLn(StringSplitter *_this);
};

struct StringSplitter::StateInWord : public StringSplitter::StringSplitterState
{
  virtual StringSplitterState *onWhite(StringSplitter *_this, char c);
  virtual StringSplitterState *onGraph(StringSplitter *_this, char c);
  virtual StringSplitterState *onQuote(StringSplitter *_this);
  virtual StringSplitterState *onEndLn(StringSplitter *_this);
};

struct StringSplitter::StateInString : public StringSplitter::StringSplitterState
{
  virtual StringSplitterState *onWhite(StringSplitter *_this, char c);
  virtual StringSplitterState *onGraph(StringSplitter *_this, char c);
  virtual StringSplitterState *onQuote(StringSplitter *_this);
  virtual StringSplitterState *onEndLn(StringSplitter *_this);
};

//static member declarations
StringSplitter::StateWhite StringSplitter::STATE_WHITE;
StringSplitter::StateInWord StringSplitter::STATE_IN_WORD;
StringSplitter::StateInString StringSplitter::STATE_IN_STRING;

//
//  StateWhite functions.
//
StringSplitter::StringSplitterState *
StringSplitter::StateWhite::onWhite(StringSplitter *, char)
{
  return &STATE_WHITE;
}

StringSplitter::StringSplitterState *
StringSplitter::StateWhite::onQuote(StringSplitter *_this)
{
  _this->addChar('"');
  return &STATE_IN_STRING;
}

StringSplitter::StringSplitterState *
StringSplitter::StateWhite::onGraph(StringSplitter *_this, char c)
{
  _this->addChar(c);
  return &STATE_IN_WORD;
}

StringSplitter::StringSplitterState *
StringSplitter::StateWhite::onEndLn(StringSplitter *)
{
  return 0;
}


//
//  StateInWord functions.
//
StringSplitter::StringSplitterState *
StringSplitter::StateInWord::onWhite(StringSplitter *_this, char)
{
  _this->wordComplete();
  return &STATE_WHITE;
}

StringSplitter::StringSplitterState *
StringSplitter::StateInWord::onQuote(StringSplitter *_this)
{
  _this->addChar('"');
  return &STATE_IN_WORD;
}

StringSplitter::StringSplitterState *
StringSplitter::StateInWord::onGraph(StringSplitter *_this, char c)
{
  _this->addChar(c);
  return &STATE_IN_WORD;
}

StringSplitter::StringSplitterState *
StringSplitter::StateInWord::onEndLn(StringSplitter *_this)
{
  _this->wordComplete();
  return 0;
}

//
//  StateInString functions.
//
StringSplitter::StringSplitterState *
StringSplitter::StateInString::onWhite(StringSplitter *_this, char c)
{
  _this->addChar(c);
  return &STATE_IN_STRING;
}

StringSplitter::StringSplitterState *
StringSplitter::StateInString::onQuote(StringSplitter *_this)
{
  _this->addChar('"');
  _this->wordComplete();
  return &STATE_WHITE;
}

StringSplitter::StringSplitterState *
StringSplitter::StateInString::onGraph(StringSplitter *_this, char c)
{
  _this->addChar(c);
  return &STATE_IN_STRING;    
}

StringSplitter::StringSplitterState *
StringSplitter::StateInString::onEndLn(StringSplitter *)
{
  throw std::runtime_error("Unexpected end of line while looking for closing quote.");
}


//
// StringSplitter functions
//
StringSplitter::StringSplitter()
{
}

StringSplitter::~StringSplitter()
{
}

const Vector_s& StringSplitter::split(const char *p, char delim)
{
  result_.clear();
  current_.clear();

  for (StringSplitterState *state = &STATE_WHITE; state; ++p)
    {
      if (isspace(*p) || *p == delim) {
	state = state->onWhite(this, *p);
      } else if (*p == '"') {
	state = state->onQuote(this);
      } else if (*p == 0) {
	state = state->onEndLn(this);
      } else {
	state = state->onGraph(this, *p);        
      }
    }

  return result_;
}

const Vector_s& StringSplitter::split(const std::string& str, char delim)
{
  return split(str.c_str(), delim);
}

void StringSplitter::addChar(char c)
{
  current_ += c;
}

void StringSplitter::wordComplete()
{
  result_.push_back(current_);
  current_.clear();
}
// can use this version if empty strings are not allowed, and use " as delimiter
void StringSplitter::Tokenize(const std::string& str,
			      std::vector<std::string>& tokens,
			      const std::string& delimiters = " ")
{
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos)
    {
      if(pos==lastPos)tokens.push_back("");
      else
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}

/**
   Splits a string into tokens.
   delimiters is a string of characters on which to split.

   Example: "1, 2, 3, 4" becomes ["1" "2" "3" "4"] if delimiters is ", "
*/
// void StringSplitter::Tokenize(const std::string& str,
// 			      std::vector<std::string>& tokens,
// 			      const std::string& delimiters = " ")
// {
//   std::string::size_type lastPos = 0;
//   std::string::size_type pos     = 0;
//   std::string::size_type paren = 0;
//   std::string::size_type closeparen = 0;
//   do    {
//     // Skip delimiters.  Note the "not_of"
//     lastPos = str.find_first_not_of(delimiters, pos);
    
//     // Find next "non-delimiter"
//     pos = str.find_first_of(delimiters, lastPos);
    
//     paren = str.find_first_of("\"", lastPos);
//     if(std::string::npos != paren ){
//       closeparen = str.find_first_of("\"", paren+1);
//       if(std::string::npos == closeparen) 
// 	throw std::string("missing closing quotes");
//       if(paren < pos && pos < closeparen){//skip delimiters within quotes
// 	//++lastPos;
// 	pos =closeparen+1;
//       }
//     }
//     if(std::string::npos != pos || std::string::npos != lastPos)
//       // Found a token, add it to the vector.
//       tokens.push_back(str.substr(lastPos, pos - lastPos));
//   }
  
//   while (std::string::npos != pos || std::string::npos != lastPos);
// }

//possible alternative

void StringSplitter::QuickTokenize(const std::string& text, std::vector<std::string>& tokens, const std::string& delim) {
  tokens.clear();
  std::string::size_type b(text.find_first_not_of(delim));
  while (b != std::string::npos) {
    std::string::size_type e(text.find_first_of(delim, b));
    tokens.push_back(text.substr(b, e - b));
    b = text.find_first_not_of(delim, (e<text.size() ? e: text.size()) );
  }
} 
END_BCLIB_NAMESPACE
