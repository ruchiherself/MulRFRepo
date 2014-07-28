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

#ifndef ARGUMENT_H

#include <vector>
#include <map>
#include <string.h>
#include "common.h"
//using namespace std;

class Argument {
    public: class arg_type {
        public: std::string key;
        public: std::string value;
        public: bool visited;

        //Constructors 
        public: arg_type() {
            key = "";
            value = "";
            visited = false;
        }
        public: arg_type(std::string key) {
            this->key = key;
            value = "";
            visited = false;
        }
        public: arg_type(std::string key, std::string value) {
            this->key = key;
            this->value = value;
            visited = false;
        }
        
        // specialization for std::string
        inline void convert(std::string &var) const {
            var = value;
        }
        
        // convert a character std::string into another datatype
        template<class T>
        void convert(T &var) const {
            std::istringstream ist(value);
            if (!(ist >> var))
                ERROR_exit("invalid value of argument " << key << ' ' << value);
        }
        
        // true if argument has a value
        public: inline bool hasValue() const {
            return (!value.empty());
        }
        public: inline bool operator < (const arg_type& ref) const {
            return key < ref.key;
        }
        public: inline bool operator > (const arg_type& ref) const {
            return key > ref.key;
        }
        public: inline bool operator == (const arg_type& ref) const {
            return key == ref.key;
        }
    };

    //Data Member of Argument class
    private: std::map<std::string, arg_type> args;

    // adds arguments from the command line
    public: inline void add(const int argc, char *argsv[]) {
        if (argc <= 1)
            return;
        std::string key = "";
        std::string value = "";
        for (int i=1; i<argc; i++) {
            std::string str(argsv[i]);
            if (str.length() < 2)
                ERROR_exit("argument too short");

            if (str[0] == '-') { // -t
                const size_t pos = str.find("=");
                if (pos != std::string::npos) {
                    key = str.substr(0,pos);
                    value = str.substr(pos+1);
                } else {
                    key = str;
                }
                for (;;) {
                    ++i;
                    if (i >= argc) break;
                    std::string str(argsv[i]);
                    if (str[0] == '-') {
                        --i;
                        break;
                    }
                        if (!value.empty()) value += " ";
                        value = value + str;
                }
                add_insert(key,value);
                key = "";
                value = "";
            }
        }
    }
    
    private:
    template<class T>
    inline void add_insert(std::string key, T value) {
        if (args.find(key) != args.end())
            ERROR_exit("more than one argument " << key);
        args[key] = arg_type(key, value);
        key = "";
    }

    public:
	// searches for a particular argument
	inline arg_type *find(const std::string &key) {
		std::map<std::string, arg_type>::iterator itr = args.find(key);
		if (itr == args.end()) return NULL;
		arg_type &arg = itr->second;
		arg.visited = true;
		return &arg;
	}

	// searches for any argument
	inline arg_type *findAny(std::vector<std::string> &list) {
		for (std::vector<std::string>::iterator itr = list.begin(); itr != list.end(); itr++) {
			arg_type *arg = find(*itr);
			if (arg != NULL) return arg;
		}
		return NULL;
	}

// 	// convert a character std::string into another datatype
// 	template<class T> void convert(T &var) const;

	// return unused arguments
	inline std::vector<arg_type*> unusedArgs() {
		std::vector<arg_type*> unused;
		for (std::map<std::string, arg_type>::iterator itr = args.begin(); itr != args.end(); itr++) {
			if (!itr->second.visited) unused.push_back(&itr->second);
		}
		return unused;
	}

    inline void unusedArgsError() {
        std::vector<arg_type*> unused = unusedArgs();
        if (unused.empty()) return;
        MSG("unknown argument(s):");
        for (unsigned int i=0,iEE=unused.size(); i<iEE; ++i) {
            MSG(unused[i]->key);
        }
        ERROR_EXIT_COMMAND;
    }

    inline bool existArg(const std::string &k) {
        const Argument::arg_type *arg = find(k);
        if (arg == NULL) return false;
        if (arg->hasValue())
            MSG_exit(k << " does not needs a value");
        return true;
    }

    inline bool existArg2(const std::string &k1, const std::string &k2) {
        const Argument::arg_type *arg1 = find(k1);
        const Argument::arg_type *arg2 = find(k2);
        const Argument::arg_type *arg = arg1 != NULL ? arg1 : arg2;
        if (arg == NULL)
            return false;
        if (arg->hasValue())
            MSG_exit(k1 << " [" << k2 << "] does not needs a value");
        return true;
    }

    template<class T>
    inline bool getArgVal(const std::string &k, T &value) {
        const Argument::arg_type *arg = find(k);
        if (arg == NULL) 
            MSG_exit(k << " is not set");
        if (!arg->hasValue())
            MSG_exit(k << " needs a value");
        arg->convert(value);
        return true;
    }

    template<class T>
    inline bool getArgVal2(const std::string &k1, const std::string &k2, T &value) {
        const Argument::arg_type *arg1 = find(k1);
        const Argument::arg_type *arg2 = find(k2);
        if ((arg1 != NULL) && (arg2 != NULL))
            MSG_exit(k1 << " and " << k2 << " are exclusive");
        const Argument::arg_type *arg = arg1 != NULL ? arg1 : arg2;
        if (arg == NULL)
            MSG_exit(k1 << " [" << k2 << "] is not set");
        if (!arg->hasValue())
            MSG_exit(k1 << " [" << k2 << "] needs a value");
        arg->convert(value);
        return true;
    }

    template<class T>
    inline bool existArgVal(const std::string &k, T &value) {
        const Argument::arg_type *arg = find(k);
        if (arg == NULL) return false;
        if (!arg->hasValue()) MSG_exit(k << " needs a value");
        arg->convert(value);
        return true;
    }

    template<class T>
    inline bool existArgVal2(const std::string &k1, const std::string &k2, T &value) {
        const Argument::arg_type *arg1 = find(k1);
        const Argument::arg_type *arg2 = find(k2);
        if ((arg1 != NULL) && (arg2 != NULL))
            MSG_exit(k1 << " and " << k2 << " are exclusive");
        const Argument::arg_type *arg = arg1 != NULL ? arg1 : arg2;
        if (arg == NULL)
            return false;
        if (!arg->hasValue())
            MSG_exit(k1 << " [" << k2 << "] needs a value");
        arg->convert(value);
        return true;
    }
};

#endif
