#include "bclib/StringSplitter.h"
#include <cctype>
#include <stdexcept>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

BEGIN_BCLIB_NAMESPACE
using std::string;
using std::vector;

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

const Vector_s& StringSplitter::split(const string& str, char delim)
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

#ifdef HAVE_BOOST_H
// can use this version if multiple separators are not to be merged
// if multiple separators are to be merged, use the default 
// TokenizerFunction model which is char_separator<char>
void StringSplitter::Tokenize(const string& str,
			      vector<string>& tokens,
			      const string separators, bool, bool)
{
  // default separators are space or tab
  // multiple separators are not merged
  // text delimiters may be single quote, double quote or missing
  using namespace std;
  using namespace boost;
  escaped_list_separator<char> sep("\\", separators, "\"\'");
  tokenizer<escaped_list_separator<char> > tok(str, sep);
  for (tokenizer<escaped_list_separator<char> >::iterator tok_iter = tok.begin();
       tok_iter != tok.end(); ++tok_iter) 
      tokens.push_back(*tok_iter);
}

#else

/**
   Tokenizes a string. ie splits a string into 'tokens', determined by separators, sep.
   Substrings within quotes are preserved as tokens. If merge=true, consecutive separators 
   are merged (count as a single separator). Otherwise, two consecutive separators result in a token
   consisting of an empty string.
   //TODO: option to ignore quotes - can use QuickTokenize instead for now
*/
void StringSplitter::Tokenize(const string& str, vector<string>& tokens, const string sep, bool merge, bool dequote){

  string::size_type pos = 0, pos2 = 0;
  const string::size_type SIZE = str.size();
  if(SIZE ==0) return;

  do{
    if (sep.find(str[pos],0)!=string::npos){
      // is separator
      if(merge){//move to next non-separator
	pos = pos2 = str.find_first_not_of(sep, pos);
      }else{//add empty string and move to next character
	tokens.push_back("");
	++pos;
	++pos2;
      }
    }else if(str[pos]=='\"'){
      //is quote - find closing quote
      pos2 = str.find_first_of("\"", pos+1);
      if(pos2 == string::npos)
	throw string("ERROR: mismatched quotes in StringSplitter::Tokenize\n");
      //add contents of quotes
      if(dequote)
	tokens.push_back(str.substr(pos+1, pos2-pos-1));
      else
	tokens.push_back(str.substr(pos, pos2+1-pos));
      if(sep.find(str[pos2+1],0)!=string::npos){//ignore sep character immediately after closing quote
	if(merge)//move to next non-separator
	  pos=str.find_first_not_of(sep, pos2+1);
	else//move to next character after separator
	  pos=pos2+2;
      }
      else//otherwise move to next character after closing quote
	pos=pos2+1;

    }else{
      //is other character
      pos2 = str.find_first_of(sep+"\"", pos);
      tokens.push_back(str.substr(pos, pos2-pos));
      if(sep.find(str[pos2],0)!=string::npos){//ignore sep character immediately after closing quote
	if(merge)
	  pos=str.find_first_not_of(sep, pos2);
	else
	  pos=pos2+1;
      }
      else
	pos=pos2;
    }
  }while(pos < SIZE);

}

#endif

void StringSplitter::QuickTokenize(const string& text, vector<string>& tokens, const string delim) {
  tokens.clear();
  string::size_type b(text.find_first_not_of(delim));
  while (b != string::npos) {
    string::size_type e(text.find_first_of(delim, b));
    tokens.push_back(text.substr(b, e - b));
    b = text.find_first_not_of(delim, (e<text.size() ? e: text.size()) );
  }
} 

END_BCLIB_NAMESPACE
