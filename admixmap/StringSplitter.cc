#include "StringSplitter.h"
#include <cctype>

struct StringSplitter::StringSplitterState
{
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


StringSplitter::StateWhite StringSplitter::STATE_WHITE;
StringSplitter::StateInWord StringSplitter::STATE_IN_WORD;
StringSplitter::StateInString StringSplitter::STATE_IN_STRING;

//
//  StateWhite methods.
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
//  StateInWord methods.
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
//  StateInString methods.
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
    throw std::runtime_error("Unexpected end of line looking up for closing quote.");
}

StringSplitter::StringSplitter()
{
}

StringSplitter::~StringSplitter()
{
}

const Vector_s& StringSplitter::split(const char *p)
{
    result_.clear();
    current_.clear();

    for (StringSplitterState *state = &STATE_WHITE; state; ++p)
    {
        if (isspace(*p) || *p == ',' || *p == ' ') {
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

const Vector_s& StringSplitter::split(const std::string& str)
{
    return split(str.c_str());
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
