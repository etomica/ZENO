// Generated by Bisonc++ V6.01.00 on Mon, 20 Jul 2020 17:42:23 -0400

// base/comment

// $insert class.ih
#include "XyzParser.ih"

// The FIRST element of SR arrays shown below uses `d_type', defining the
// state's type, and `d_lastIdx' containing the last element's index. If
// d_lastIdx contains the REQ_TOKEN bitflag (see below) then the state needs
// a token: if in this state d_token is Reserved__::UNDETERMINED__, nextToken() will be
// called

// The LAST element of SR arrays uses `d_token' containing the last retrieved
// token to speed up the (linear) seach.  Except for the first element of SR
// arrays, the field `d_action' is used to determine what to do next. If
// positive, it represents the next state (used with SHIFT); if zero, it
// indicates `ACCEPT', if negative, -d_action represents the number of the
// rule to reduce to.

// `lookup()' tries to find d_token in the current SR array. If it fails, and
// there is no default reduction UNEXPECTED_TOKEN__ is thrown, which is then
// caught by the error-recovery function.

// The error-recovery function will pop elements off the stack until a state
// having bit flag ERR_ITEM is found. This state has a transition on errTok__
// which is applied. In this errTok__ state, while the current token is not a
// proper continuation, new tokens are obtained by nextToken(). If such a
// token is found, error recovery is successful and the token is
// handled according to the error state's SR table and parsing continues.
// During error recovery semantic actions are ignored.

// A state flagged with the DEF_RED flag will perform a default
// reduction if no other continuations are available for the current token.

// The ACCEPT STATE never shows a default reduction: when it is reached the
// parser returns ACCEPT(). During the grammar
// analysis phase a default reduction may have been defined, but it is
// removed during the state-definition phase.

// So:
//      s_x[] = 
//      {
//                  [_field_1_]         [_field_2_]
//
// First element:   {state-type,        idx of last element},
// Other elements:  {required token,    action to perform},
//                                      ( < 0: reduce, 
//                                          0: ACCEPT,
//                                        > 0: next state)
//      }

// base/declarations

namespace // anonymous
{
    char const author[] = "Frank B. Brokken (f.b.brokken@rug.nl)";

    enum Reserved__
    {
        UNDETERMINED__   = -2,
        EOF__            = -1,
        errTok__         = 256
    };
    enum StateType       // modify statetype/data.cc when this enum changes
    {
        NORMAL,
        ERR_ITEM,
        REQ_TOKEN,
        ERR_REQ,    // ERR_ITEM | REQ_TOKEN
        DEF_RED,    // state having default reduction
        ERR_DEF,    // ERR_ITEM | DEF_RED
        REQ_DEF,    // REQ_TOKEN | DEF_RED
        ERR_REQ_DEF // ERR_ITEM | REQ_TOKEN | DEF_RED
    };    
    enum StateTransition
    {
        ACCEPT__   = 0,     // `ACCEPT' TRANSITION
    };

    struct PI__     // Production Info
    {
        size_t d_nonTerm; // identification number of this production's
                            // non-terminal 
        size_t d_size;    // number of elements in this production 
    };

    struct SR__     // Shift Reduce info, see its description above
    {
        union
        {
            int _field_1_;      // initializer, allowing initializations 
                                // of the SR s_[] arrays
            StateType d_type;
            int       d_token;
        };
        union
        {
            int _field_2_;

            int d_lastIdx;          // if negative, the state uses SHIFT
            int d_action;           // may be negative (reduce), 
                                    // postive (shift), or 0 (accept)
        };
    };

    // $insert staticdata
    
    enum                        // size to expand the state-stack with when
    {                           // full
        STACK_EXPANSION__ = 10
    };

// Productions Info Records:
PI__ const s_productionInfo[] = 
{
     {0, 0}, // not used: reduction values are negative
     {265, 0}, // 1: snapshots ->  <empty>
     {265, 2}, // 2: snapshots ->  snapshots snapshot
     {266, 6}, // 3: snapshot (NEWLINE) ->  atom_count_integer NEWLINE #0001 comment NEWLINE atoms
     {268, 0}, // 4: #0001 ->  <empty>
     {269, 0}, // 5: atoms ->  <empty>
     {269, 3}, // 6: atoms (NEWLINE) ->  atoms atom NEWLINE
     {270, 4}, // 7: atom ->  word real real real
     {271, 1}, // 8: atom_count_integer (ATOM_COUNT_INT) ->  ATOM_COUNT_INT
     {272, 1}, // 9: real (FLOAT) ->  FLOAT
     {272, 1}, // 10: real (INT) ->  INT
     {273, 1}, // 11: word (WORD) ->  WORD
     {273, 1}, // 12: word (FLOAT) ->  FLOAT
     {273, 1}, // 13: word (INT) ->  INT
     {267, 0}, // 14: comment ->  <empty>
     {267, 2}, // 15: comment ->  comment word
     {274, 1}, // 16: snapshots_$ ->  snapshots
};

// State info and SR__ transitions for each state.


SR__ s_0[] =
{
    { { DEF_RED}, {  2} },             
    { {     265}, {  1} }, // snapshots
    { {       0}, { -1} },             
};

SR__ s_1[] =
{
    { { REQ_TOKEN}, {        5} },                      
    { {       266}, {        2} }, // snapshot          
    { {       271}, {        3} }, // atom_count_integer
    { {       261}, {        4} }, // ATOM_COUNT_INT    
    { {     EOF__}, { ACCEPT__} },                      
    { {         0}, {        0} },                      
};

SR__ s_2[] =
{
    { { DEF_RED}, {  1} }, 
    { {       0}, { -2} }, 
};

SR__ s_3[] =
{
    { { REQ_TOKEN}, { 2} },           
    { {       260}, { 5} }, // NEWLINE
    { {         0}, { 0} },           
};

SR__ s_4[] =
{
    { { DEF_RED}, {  1} }, 
    { {       0}, { -8} }, 
};

SR__ s_5[] =
{
    { { DEF_RED}, {  2} },         
    { {     268}, {  6} }, // #0001
    { {       0}, { -4} },         
};

SR__ s_6[] =
{
    { { DEF_RED}, {   2} },           
    { {     267}, {   7} }, // comment
    { {       0}, { -14} },           
};

SR__ s_7[] =
{
    { { REQ_TOKEN}, {  6} },           
    { {       260}, {  8} }, // NEWLINE
    { {       273}, {  9} }, // word   
    { {       259}, { 10} }, // WORD   
    { {       258}, { 11} }, // FLOAT  
    { {       257}, { 12} }, // INT    
    { {         0}, {  0} },           
};

SR__ s_8[] =
{
    { { DEF_RED}, {  2} },         
    { {     269}, { 13} }, // atoms
    { {       0}, { -5} },         
};

SR__ s_9[] =
{
    { { DEF_RED}, {   1} }, 
    { {       0}, { -15} }, 
};

SR__ s_10[] =
{
    { { DEF_RED}, {   1} }, 
    { {       0}, { -11} }, 
};

SR__ s_11[] =
{
    { { DEF_RED}, {   1} }, 
    { {       0}, { -12} }, 
};

SR__ s_12[] =
{
    { { DEF_RED}, {   1} }, 
    { {       0}, { -13} }, 
};

SR__ s_13[] =
{
    { { REQ_DEF}, {  6} },         
    { {     270}, { 14} }, // atom 
    { {     273}, { 15} }, // word 
    { {     259}, { 10} }, // WORD 
    { {     258}, { 11} }, // FLOAT
    { {     257}, { 12} }, // INT  
    { {       0}, { -3} },         
};

SR__ s_14[] =
{
    { { REQ_TOKEN}, {  2} },           
    { {       260}, { 16} }, // NEWLINE
    { {         0}, {  0} },           
};

SR__ s_15[] =
{
    { { REQ_TOKEN}, {  4} },         
    { {       272}, { 17} }, // real 
    { {       258}, { 18} }, // FLOAT
    { {       257}, { 19} }, // INT  
    { {         0}, {  0} },         
};

SR__ s_16[] =
{
    { { DEF_RED}, {  1} }, 
    { {       0}, { -6} }, 
};

SR__ s_17[] =
{
    { { REQ_TOKEN}, {  4} },         
    { {       272}, { 20} }, // real 
    { {       258}, { 18} }, // FLOAT
    { {       257}, { 19} }, // INT  
    { {         0}, {  0} },         
};

SR__ s_18[] =
{
    { { DEF_RED}, {  1} }, 
    { {       0}, { -9} }, 
};

SR__ s_19[] =
{
    { { DEF_RED}, {   1} }, 
    { {       0}, { -10} }, 
};

SR__ s_20[] =
{
    { { REQ_TOKEN}, {  4} },         
    { {       272}, { 21} }, // real 
    { {       258}, { 18} }, // FLOAT
    { {       257}, { 19} }, // INT  
    { {         0}, {  0} },         
};

SR__ s_21[] =
{
    { { DEF_RED}, {  1} }, 
    { {       0}, { -7} }, 
};


// State array:
SR__ *s_state[] =
{
  s_0,  s_1,  s_2,  s_3,  s_4,  s_5,  s_6,  s_7,  s_8,  s_9,
  s_10,  s_11,  s_12,  s_13,  s_14,  s_15,  s_16,  s_17,  s_18,  s_19,
  s_20,  s_21,
};

} // anonymous namespace ends


// $insert namespace-open
namespace xyz_parser
{

// $insert polymorphicCode
namespace Meta__
{

size_t const *t_nErrors;
// $insert idoftag
char const *idOfTag__[] = {
    "STRING_TYPE",
    "INT_TYPE",
    "FLOAT_TYPE",
    "<undefined>"
};

size_t const *s_nErrors__;

Base::~Base()
{}

}   // namespace Meta__

// If the parsing function call (i.e., parse()' needs arguments, then provide
// an overloaded function.  The code below doesn't rely on parameters, so no
// arguments are required.  Furthermore, parse uses a function try block to
// allow us to do ACCEPT and ABORT from anywhere, even from within members
// called by actions, simply throwing the appropriate exceptions.


// base/base1
XyzParserBase::XyzParserBase()
:
    d_token(Reserved__::UNDETERMINED__),
    // $insert baseclasscode
    d_requiredTokens__(0)
{
    Meta__::t_nErrors = &d_nErrors__;
}

// base/clearin
void XyzParserBase::clearin__()
{
    d_nErrors__ = 0;
    d_stackIdx = -1;
    d_stateStack.clear();
    d_token = Reserved__::UNDETERMINED__;
    d_next = TokenPair{ Reserved__::UNDETERMINED__, STYPE__{} };
    d_recovery = false;
    d_acceptedTokens__ = d_requiredTokens__;
    d_val__ = STYPE__{};

    push__(0);
}

// base/debugfunctions

void XyzParserBase::setDebug(bool mode)
{
    d_actionCases__ = false;
    d_debug__ = mode;
}

void XyzParserBase::setDebug(DebugMode__ mode)
{
    d_actionCases__ = mode & ACTIONCASES;
    d_debug__ =       mode & ON;
}

// base/lex
void XyzParserBase::lex__(int token)
{
    d_token = token;

    if (d_token <= 0)
        d_token = Reserved__::EOF__;

    d_terminalToken = true;
}

// base/lookup
int XyzParserBase::lookup__() const
{
    // if the final transition is negative, then we should reduce by the rule
    // given by its positive value.

    SR__ const *sr = s_state[d_state];
    SR__ const *last = sr + sr->d_lastIdx;

    for ( ; ++sr != last; )           // visit all but the last SR entries
    {
        if (sr->d_token == d_token)
            return sr->d_action;
    }

    if (sr == last)   // reached the last element
    {
        if (sr->d_action < 0)   // default reduction
        {
            return sr->d_action;                
        }

        // No default reduction, so token not found, so error.
        throw UNEXPECTED_TOKEN__;
    }

    // not at the last element: inspect the nature of the action
    // (< 0: reduce, 0: ACCEPT, > 0: shift)

    int action = sr->d_action;


    return action;
}

// base/pop
void XyzParserBase::pop__(size_t count)
{
    if (d_stackIdx < static_cast<int>(count))
    {
        ABORT();
    }

    d_stackIdx -= count;
    d_state = d_stateStack[d_stackIdx].first;
    d_vsp = &d_stateStack[d_stackIdx];

}

// base/poptoken
void XyzParserBase::popToken__()
{
    d_token = d_next.first;
    d_val__ = std::move(d_next.second);

    d_next.first = Reserved__::UNDETERMINED__;
}

// base/push
void XyzParserBase::push__(size_t state)
{
    size_t currentSize = d_stateStack.size();
    if (stackSize__() == currentSize)
    {
        size_t newSize = currentSize + STACK_EXPANSION__;
        d_stateStack.resize(newSize);
    }

    ++d_stackIdx;
    d_stateStack[d_stackIdx] = 
                    StatePair{ d_state = state, std::move(d_val__) };

    d_vsp = &d_stateStack[d_stackIdx];

    if (d_stackIdx == 0)
    {
    }
    else
    {
    }
}

// base/pushtoken
void XyzParserBase::pushToken__(int token)
{
    d_next = TokenPair{ d_token, std::move(d_val__) };
    d_token = token;
}

// base/redotoken
void XyzParserBase::redoToken__()
{
    if (d_token != Reserved__::UNDETERMINED__)
        pushToken__(d_token);
}

// base/reduce
void XyzParserBase::reduce__(int rule)
{
    PI__ const &pi = s_productionInfo[rule];

    d_token = pi.d_nonTerm;
    pop__(pi.d_size);

    d_terminalToken = false;
}

// base/shift
void XyzParserBase::shift__(int action)
{
    push__(action);
    popToken__();               // token processed

    if (d_recovery and d_terminalToken)
    {
        d_recovery = false;
        d_acceptedTokens__ = 0;
    }
}

// base/startrecovery
void XyzParserBase::startRecovery__()
{
    int lastToken = d_token;                // give the unexpected token a
                                            // chance to be processed
                                            // again.

    pushToken__(Reserved__::errTok__);      // specify errTok__ as next token
    push__(lookup__());                     // push the error state

    d_token = lastToken;                    // reactivate the unexpected
                                            // token (we're now in an
                                            // ERROR state).

    d_recovery = true;
}

// base/top
inline size_t XyzParserBase::top__() const
{
    return d_stateStack[d_stackIdx].first;
}

// derived/errorrecovery
void XyzParser::errorRecovery__()
{
    // When an error has occurred, pop elements off the stack until the top
    // state has an error-item. If none is found, the default recovery
    // mode (which is to abort) is activated. 
    //
    // If EOF is encountered without being appropriate for the current state,
    // then the error recovery will fall back to the default recovery mode.
    // (i.e., parsing terminates)



    if (d_acceptedTokens__ >= d_requiredTokens__)// only generate an error-
    {                                           // message if enough tokens 
        ++d_nErrors__;                          // were accepted. Otherwise
        error();                                // simply skip input
    }

    // get the error state
    while (not (s_state[top__()][0].d_type & ERR_ITEM))
    {
        pop__();
    }

    // In the error state, looking up a token allows us to proceed.
    // Continuation may be require multiple reductions, but eventually a
    // terminal-token shift is used. See nextCycle__ for details.

    startRecovery__();
}

// derived/executeaction
void XyzParser::executeAction__(int production)
try
{
    if (token__() != Reserved__::UNDETERMINED__)
        pushToken__(token__());     // save an already available token
    switch (production)
    {
        // $insert actioncases
        
        case 2:
#line 19 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
            d_val__ = std::move(vs__(-1));
        }
        break;

        case 3:
#line 28 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
            d_val__ = std::move(vs__(-5));
        }
        break;

        case 4:
#line 23 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
         addSnapshot(vs__(-1).get<Tag__::INT_TYPE>());
         }
        break;

        case 6:
#line 35 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
            d_val__ = std::move(vs__(-2));
        }
        break;

        case 7:
#line 40 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
         addAtom(vs__(-3).get<Tag__::STRING_TYPE>(), vs__(-2).get<Tag__::FLOAT_TYPE>(), vs__(-1).get<Tag__::FLOAT_TYPE>(), vs__(0).get<Tag__::FLOAT_TYPE>());
         }
        break;

        case 8:
#line 47 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
         d_val__ = stoi(d_scanner.matched());
         }
        break;

        case 9:
#line 54 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
         d_val__ = stod(d_scanner.matched());
         }
        break;

        case 10:
#line 59 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
         d_val__ = stod(d_scanner.matched());
         }
        break;

        case 11:
#line 66 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
         d_val__ = d_scanner.matched();
         }
        break;

        case 12:
#line 71 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
         d_val__ = d_scanner.matched();
         }
        break;

        case 13:
#line 76 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
         d_val__ = d_scanner.matched();
         }
        break;

        case 15:
#line 85 "/home/dcj/Projects/zeno/cpp/zeno-cli/src/XyzParser/xyz_grammar"
        {
         addCommentWord(vs__(0).get<Tag__::STRING_TYPE>());
         }
        break;

    }
}
catch (std::exception const &exc)
{
    exceptionHandler(exc);
}

// derived/nextcycle
void XyzParser::nextCycle__()
try
{
    if (s_state[state__()]->d_type & REQ_TOKEN)
        nextToken__();              // obtain next token


    int action = lookup__();        // lookup d_token in d_state

    if (action > 0)                 // SHIFT: push a new state
    {
        shift__(action);
        return;
    }

    if (action < 0)            // REDUCE: execute and pop.
    {

        if (recovery__())
            redoToken__();
        else
            executeAction__(-action);
                                            // next token is the rule's LHS
        reduce__(-action); 
        return;
    }

    if (recovery__())
        ABORT();
    else 
        ACCEPT();
}
catch (ErrorRecovery__)
{
    if (not recovery__())
        errorRecovery__();
    else
    {
        if (token__() == Reserved__::EOF__)
            ABORT();
        popToken__();               // skip the failing token
    }
}


// derived/nexttoken
void XyzParser::nextToken__()
{ 
    // If d_token is Reserved__::UNDETERMINED__ then if savedToken__() is
    // Reserved__::UNDETERMINED__ another token is obtained from lex(). Then
    // savedToken__() is assigned to d_token.

                                    // no need for a token: got one already
    if (token__() != Reserved__::UNDETERMINED__) 
    {
        return;                             
    }

    if (savedToken__() != Reserved__::UNDETERMINED__)
    {
        popToken__();               // consume pending token
    }
    else
    {
        ++d_acceptedTokens__;       // accept another token (see
                                    // errorRecover())
        lex__(lex());
        print__();
    }
    print();
}

// derived/print
void XyzParser::print__()
{
// $insert print
}

// derived/parse
int XyzParser::parse()
try 
{
    // The parsing algorithm:
    // Initially, state 0 is pushed on the stack, and all relevant variables
    // are initialized by Base::clearin__.
    //
    // Then, in an eternal loop:
    //
    //  1. If a state is a REQ_TOKEN type, then the next token is obtained
    //     from nextToken().  This may very well be the currently available
    //     token. When retrieving a terminal token d_terminal is set to true.
    //
    //  2. lookup() is called, d_token is looked up in the current state's
    //     SR_ array.
    //
    //  4. Depending on the result of the lookup() function the next state is
    //     shifted on the parser's stack, a reduction by some rule is applied,
    //     or the parsing function returns ACCEPT(). When a reduction is
    //     called for, any action that may have been defined for that
    //     reduction is executed.
    //
    //  5. An error occurs if d_token is not found, and the state has no
    //     default reduction.

    clearin__();                            // initialize, push(0)

    while (true)
    {
// $insert prompt
        nextCycle__();
    }
}
catch (Return__ retValue)
{
    return retValue or d_nErrors__;
}


// derived/tail

// $insert namespace-close
}
