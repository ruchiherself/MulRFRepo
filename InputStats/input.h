/*
 * Copyright (C) 2007 Andre Wehe, Mukul Bansal, Oliver Eulenstein
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#ifndef INPUT_H
#define INPUT_H

// legal characters
// std::string legalChars4Name = "abcdefghijklmnopqrtsuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
#define legalChar4Name(c) (\
        ((c>='a') && (c<='z')) || \
        ((c>='A') && (c<='Z')) || \
        ((c>='0') && (c<='9')) | \
        (c=='+') | \
        (c=='-') | \
        (c=='.') | \
        (c=='_') \
    )
#define legalChar4Number(c) (\
        ((c>='0') && (c<='9')) || \
        (c=='.') || \
        (c=='-') || \
        (c=='+') || \
        (c=='e') \
    )

#include <stack>
#include <string>
#include <sstream>

namespace NS_input {

using namespace std;

inline bool legalstring(const std::string &str) {
    for (unsigned int i=0,iEE=str.length(); i<iEE; ++i) {
        if (!legalChar4Name(str[i])) return false;
    }
    return true;
}

inline std::string getlegalstring(const std::string &str) {
    if (legalstring(str)) return str;
    std::ostringstream os;
    os << '"';
    for (unsigned int i=0,iEE=str.length(); i<iEE; ++i) {
        if (str[i] != '"') os << str[i];
    }
    os << '"';
    return os.str();
}

// interface for callbacks for comments
class InputCommentCallbackBase {
public:
    virtual void commentFound(string &str) = 0;
    virtual ~InputCommentCallbackBase() {};
};

// handles intup from a stream
class Input {
public:
    std::string legalChars4Name;
    istream *in;
    int prevline, line;
    int prevcolumn, column;
    stack<char> st;
    bool skipComments;

    Input(istream *in = NULL, int line = 1, int column = 1) : in(in) {
        legalChars4Name = "abcdefghijklmnopqrtsuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
        this->line = line;
        prevline = line;
        this->column = column;
        prevcolumn = column;
        skipComments = true;
        callbackflag = false;
    }

    // set the input stream
    inline void setInputStream(istream &is) {
        in = &is;
    }

    // checks if c is a whitespace
    inline bool whitespace(char &c) {
        return ((c==' ') || (c=='\t') || (c=='\n') || (c=='\r') || (c=='\f'));
    }

    // gets the next character
    inline bool nextChar(char &c) {
        if (!st.empty()) {
            c = st.top();
            st.pop();
            return true;
        }
        if (!in->get(c)) return false;
        prevcolumn = column++;
        prevline = line;
        if (c == '\n') {
            line++;
            column = 1;
        }
        return true;
    }

    // gets the next character that is not a whitespace
    inline bool nextAnyChar(char &c) {
        while (nextChar(c)) {
            if ((skipComments) && (c=='[')) {
                ostringstream comment;
                comment << c;
                for (;;) {
                    if (!nextChar(c)) return false;
                    comment << c;
                    if (c==']') {
                        if (!nextChar(c)) return false;
                        break;
                    }
                }
                string str = comment.str();
                if (callbackflag) {
                    commentcallback->commentFound(str);
                }
            }
            if (!whitespace(c)) return true;
        }
        return false;
    }

    // skip following white spaces
    inline bool skipWhiteSpaces() {
        char c;
        if (!nextAnyChar(c)) return false;
        pushBack(c);
        return true;
    }

    // puts 1 character back onto the stream
    inline void pushBack(const char c) {
        st.push(c);
    }

    // reads a string
    inline string getString() {
        ostringstream name;
        char c;
        while (nextChar(c)) {
            if (whitespace(c)) {
                pushBack(c);
                break;
            }
            name << c;
        };
        return name.str();
    }

    // read a comment
    inline string getComment() {
        setSkipComments(false);
        ostringstream comment;
        char c;
        if (nextAnyChar(c)) {
            if (c != '[') pushBack(c);
            else {
                comment << '[';
                while (nextChar(c) && (c != ']')) {
                    comment << c;
                }
                comment << ']';
            }
        }
        setSkipComments(true);
        return comment.str();
    }

    // read a name
    inline string getName() {
        ostringstream name;
        char c;
        if (!nextAnyChar(c)) return "";
        if ((c == '\'') || (c == '"')) {
            while (nextChar(c) && (c != '\'') && (c != '"')) {
                name << c;
            }
        } else {
            do {
                if (legalChar4Name(c)) {
                    name << c;
                } else {
                    pushBack(c);
                    break;
                }
            } while (nextChar(c));
        }
        return name.str();
    }

    // read a numeric value
    inline double readNumber() {
        double value;
        if (!readNumber(value)) ERROR_exit("cannot read numeric value at " << getPos());
        return value;
    }

    // read a numeric value
    template<typename T>
    inline bool readNumber(T &value) {
        ostringstream name;
        char c;
        if (!nextAnyChar(c)) return 0;
        if ((c == '\'') || (c == '"')) {
            while (nextChar(c) && (c != '\'') && (c != '"')) {
                name << c;
            }
        } else {
            do {
                if (legalChar4Number(c)) {
                    name << c;
                } else {
                    pushBack(c);
                    break;
                }
            } while (nextChar(c));
        }
        // trim leading whitespace
        string str2 = name.str();
        char const delims[] = " \t\r\n\f";
        string::size_type notwhite = str2.find_first_not_of(delims);
        str2.erase(0, notwhite);

        // trim trailing whitespace
        notwhite = str2.find_last_not_of(delims);
        str2.erase(notwhite+1);

        // convert to float
        istringstream ss(str2);
        ss >> value;
        if (ss.fail()) return false;

        // look for junk after the value
        string junk;
        ss >> junk;
        if (junk.size() != 0) return false;

        return true;
    }

    // return the current postion as a string
    inline string getPos() {
        ostringstream os;
        os << "line " << line << " column " << column;
        return os.str();
    }

    // return the previous postion as a string
    inline string getLastPos() {
        ostringstream os;
        os << "line " << prevline << " column " << prevcolumn;
        return os.str();
    }

    // ignore comments when reading input stream or not
    inline void setSkipComments(bool flag) {
        skipComments = flag;
    }

    // set a callback that gets executed when a comment is read
    InputCommentCallbackBase *commentcallback;
    bool callbackflag;
    inline void setCommentCallback(InputCommentCallbackBase *obj) {
        callbackflag = (obj != NULL);
        commentcallback = obj;
    }
};

} // end namespace

#endif
