#ifndef UTIL_HH
#define UTIL_HH

#include <sstream>
#include <string>
#include <set>
#include <queue>
#include <boost/dynamic_bitset.hpp>
//#include "rmq.h"

namespace util {
    // return the minium of 2 values
    template <class X>
    inline const X& min(const X &a, const X &b) {
        return (b<a) ? b : a;
    }

    // return the minium of 2 values
    template <class X>
    inline const X& max(const X &a, const X &b) {
        return (b>a) ? b : a;
    }

    // swap the values or 2 variables
    template <class X>
    inline void swap(X &a, X &b) {
        X temp = a;
        a = b;
        b = temp;
    }
    // convert a printable data type into a string
    template<typename T>
    inline std::string stringify (const T& x) {
        std::ostringstream out;
        if (!(out << x))
        return "";
        return out.str();
    }

    // triplet container class
    template<typename T1, typename T2, typename T3> class triplet {
    public:
        typedef T1 first_type;
        typedef T2 second_type;
        typedef T3 third_type;

        triplet (const T1 &a, const T2 &b, const T3 &c) : first (a), second (b), third (c) {}

        triplet () : first (T1 ()), second (T2 ()), third (T3 ()) {}

        template<typename U1, typename U2, typename U3> triplet (const triplet<U1, U2, U3> &t)
        : first (t.first), second (t.second), third (t.third) {}

    public:
        T1 first;
        T2 second;
        T3 third;
    };

    // dummy class
    class empty {};

    // vector class
    template<class T>
    class vector { // lightweight vector similar to std::vector
        public: typedef T value_type;
        private: value_type *data;
        private: unsigned int data_size;
        private: unsigned int index_end;
        public: vector() { data = NULL; data_size = 0; index_end = 0; }
        public: vector(const vector<T> &r) { // copy constructor
            copy(r);
        }
        public: vector& operator=(const vector<T>& r) { // assign operator
            if (this != &r) {
                this->free();
                copy(r);
            }
            return *this;
        }
        protected: inline void copy(const vector<T> &r) {
            if (r.data_size != 0) {
                index_end = r.index_end;
                data_size = r.data_size;
                data = new value_type[data_size];
                memcpy(data, r.data, data_size * sizeof(value_type));
            }
        }
        public: ~vector() { free(); }
        private: inline void free() { if (data != NULL) delete [] data; }
        public: inline void set_min_size(const unsigned int s) {
            if (data_size < s) {
                free();
                data_size = s;
                data = new value_type[data_size];
            }
            index_end = 0;
        }
        public: inline void push_back(const value_type &item) { data[index_end++] = item; }
        public: inline value_type &operator[](const unsigned int i) { return data[i]; }
        public: inline unsigned int size() { return index_end; }
        public: inline bool empty() { return index_end == 0; }
        public: inline void clear() { index_end = 0; }
    };

    // split FIFO queue - allows to insert at the end or at a fixed point in the middle
    template<typename T>
    class splitqueue {
        protected: std::queue<T> a,b;
        protected: size_t a_max;
        public: splitqueue() { a_max = 0; }
        public: splitqueue(size_t &p) : a_max(p) {}
        public: inline size_t size() {
            return a.size() + b.size();
        }
        public: inline void set_split(std::size_t p) {
            a_max = p;
        }
        public: inline T& front() {
            if (!a.empty()) return a.front();
            return b.front();
        }
        public: inline void push(const T v) {
            if (a.size() < a_max) {
                a.push(v);
            } else {
                b.push(v);
            }
        }
        public: inline void push_middle(const T v) {
            a.push(v);
        }
        public: inline void pop() {
            if (a.empty()) {
                if (b.empty()) return;
                b.pop();
            } else {
                a.pop();
                if (a.size() < a_max) {
                    if (b.empty()) return;
                    a.push(b.front());
                    b.pop();
                }
            }
        }
    };

    class lookup {
        protected: boost::dynamic_bitset<> tab;
        protected: std::vector<unsigned int> vals;
        public: inline void set_min_size(std::size_t s) {
            clear();
            if (tab.size() < s) {
                tab.resize(s);
                vals.reserve(s);
            }
        }
        public: inline void clear() {
            for (unsigned int i=0,iEE=vals.size(); i<iEE; ++i) {
                tab[vals[i]] = 0;
            }
            vals.clear();
        }
        public: inline void add(const unsigned int v) {
            if (tab[v] == 0) {
                tab[v] = 1;
                vals.push_back(v);
            }
        }
        public: inline bool exist(const unsigned int v) {
            return (tab[v] != 0);
        }
        public: inline std::size_t item_size() {
            return vals.size();
        }
    };

	    // trim str and store it into value
    template<class T>
    bool convert(std::string str, T &value) {
        // trim leading whitespace
        std::string str2 = str;
        char const delims[] = " \t\r\n\f";
        std::string::size_type notwhite = str2.find_first_not_of(delims);
        str2.erase(0, notwhite);

        // trim trailing whitespace
        notwhite = str2.find_last_not_of(delims);
        str2.erase(notwhite+1);

        // convert to value
        std::istringstream ss(str2);
        ss >> value;
        if (ss.fail()) return false;

        // look for junk after the value
        std::string junk;
        ss >> junk;
        if (junk.size() != 0) return false;
        return true;
    }

    //convert seconds into days, hours, minutes, and seconds format: Added by Ruchi
    inline void convertTime(int input_seconds,int &d, int &h, int &m, int &s){
        const int HOURS_IN_DAY = 24;
        const int MINS_IN_HOUR = 60;
        const int SECS_IN_MIN = 60;        

        s = input_seconds % SECS_IN_MIN;
        m = (input_seconds/SECS_IN_MIN) % MINS_IN_HOUR;
        h =  ((input_seconds/SECS_IN_MIN)/MINS_IN_HOUR) % HOURS_IN_DAY;
        d = ((input_seconds/SECS_IN_MIN)/MINS_IN_HOUR)/HOURS_IN_DAY;
        
        return;
    }

    // extract RGB from a string of the format "#rrggbb"
    inline bool extractRGB(const std::string &str, unsigned int &r, unsigned int &g, unsigned int &b) {
        // trim leading whitespace
        std::string str2 = str;
        char const delims[] = " \t\r\n\f";
        std::string::size_type notwhite = str2.find_first_not_of(delims);
        str2.erase(0, notwhite);

        // trim trailing whitespace
        notwhite = str2.find_last_not_of(delims);
        str2.erase(notwhite+1);

        {
            std::string str = std::string("0x") + str2.substr(1,2);
            std::stringstream ss;
            ss << std::hex << str;
            ss >> r;
            if (ss.fail()) return false;
        }
        {
            std::string str = std::string("0x") + str2.substr(3,2);
            std::stringstream ss;
            ss << std::hex << str;
            ss >> g;
            if (ss.fail()) return false;
        }
        {
            std::string str = std::string("0x") + str2.substr(5,2);
            std::stringstream ss;
            ss << std::hex << str;
            ss >> b;
            if (ss.fail()) return false;
        }
        return true;
    }
    // concatenate 2 strings and places a delimiter if needed
    inline std::string concatenateWithDelimiter(const std::string &l_str, const std::string &delim, const std::string &r_str) {
        if (l_str.empty()) return r_str;
        if (r_str.empty()) return l_str;
        return l_str + delim + r_str;
    }

}

#endif //#ifndef UTIL_HH
