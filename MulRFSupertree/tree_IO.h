/*
 * Copyright (C) 2009 Andre Wehe
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

#ifndef TREE_IO_H
#define TREE_IO_H

#include "common.h"
#include "input.h"
#include "tree_traversal.h"
#include "util.h"
#include <iostream>
#include <iomanip>
#include <stack>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#ifdef NOHASH
#include <map>
#else
#include <boost/unordered_map.hpp>
#endif


namespace aw {

using namespace std;

#ifdef NOHASH
typedef std::map<unsigned int,std::string> idx2name;
#else
typedef boost::unordered_map<unsigned int,std::string> idx2name;
#endif

// stores one or more node/edge weights
// example
// weights<double> a=10.3;
// weight<int> b; b[0]=23; b[1]=21;
template<class VALUE>
class weights_type {
    public: typedef VALUE value_type;
    private: std::vector<VALUE> data;
    public: inline VALUE &operator[](const unsigned int i) {
        if (i >= data.size()) data.resize(i+1);
        return data[i];
    }
    public: inline weights_type& operator=(const VALUE r) {
        (*this)[0] = r;
        return *this;
    }
    public: inline size_t size() {
        return data.size();
    }
    // required function for tree_IO
    public: inline void push_back(const VALUE v) {
        data.push_back(v);
    }
    // required function for tree_IO
    public: inline void display_in_newick(std::ostream &os) {
        BOOST_FOREACH(const VALUE &v, data) {
            os << ':' << v;            
        }
    }
};
// map from node_id to weight(s)
// VALUE is the data type precision of the weight
template<class VALUE>
class idx2weight_type
#ifdef NOHASH
: public std::map<unsigned int,weights_type<VALUE> >
#else
: public boost::unordered_map<unsigned int,weights_type<VALUE> >
#endif
{
    public: typedef weights_type<VALUE> data_type;
};

// statndard weight data type is double
typedef idx2weight_type<double> idx2weight; // default is double;
typedef idx2weight_type<double> idx2weight_double;
typedef idx2weight_type<float> idx2weight_float;
typedef idx2weight_type<int> idx2weight_int;
typedef idx2weight_type<std::string> idx2weight_string;

// read the tree from a string stream (newick formatted e.g. ((name1,name2),name3);)
template<class TREE, class WEIGHTS>
bool stream2tree(std::istream &is, TREE &tree, idx2name &names, WEIGHTS &weights, float &t_w) {
    const unsigned int default_root = tree.node_size();
    char c;
    NS_input::Input input = &is;
    std::string rooting = input.getComment();
    std::string weighting;
    if(!rooting.empty())
        weighting = input.getComment();
    if (!input.nextAnyChar(c)) return false;
    std::stack<unsigned int> parents;
    unsigned int lin = UINT_MAX; // Last Internal Node
    for (;;) {
        if (c == ';') {
            lin = UINT_MAX;
            if (!parents.empty()) ERROR_return/*ERROR_exit*/("problem in the tree expression " << input.getLastPos());
            break;
        } else
        if (c == '(') { // new subtree
            const unsigned int subtree_root = tree.new_node();
            if (!parents.empty()) { // add an edge to the parent unless it is the root node
                tree.add_edge(parents.top(), subtree_root);
            }
            parents.push(subtree_root);
        } else
        if (c == ',') { // sibling
            lin = UINT_MAX;
        } else
        if (c == ')') { // subtree completed
            lin = parents.top();
            parents.pop();
        } else
        if (c == ':') { // sibling
            std::string w_str = input.getName();
            if (w_str.empty()) ERROR_return("cannot read weight/branch length value " << input.getLastPos());
            typename WEIGHTS::data_type::value_type w;
            if (util::convert(w_str,w)) weights[lin].push_back(w);
            else ERROR_return("cannot read weight/branch length value " << input.getLastPos());
//             double weight;
//             if (input.readNumber(weight)) weights[lin].push_back(weight);
//             else ERROR_return("cannot read weight/branch length value " << input.getLastPos());
        } else { // name
            if (lin == UINT_MAX) { // name for a leaf or internal node?
                lin = tree.new_node();
                if (!parents.empty()) { // parents can be empty when the tree is only 1 node "i.e. [&R] species;"
                    tree.add_edge(parents.top(), lin);
                }
            }
            input.pushBack(c);
            std::string name = input.getName();
            if (name.empty()) ERROR_return("problem in the tree expression " << input.getLastPos());
            
            names.insert(idx2name::value_type(lin, name));
        }
        if (!input.nextAnyChar(c)) return false;
    }
    tree.root = default_root;
    //Setting tree weight values
    if (rooting.find("WEIGHT")!=string::npos) {
        int posE = rooting.find("=");
        t_w = boost::lexical_cast<float>(rooting.substr(posE+1,rooting.length()-10));
    } else if (weighting.find("WEIGHT")!=string::npos) {
        int pos = weighting.find("=");
        t_w = boost::lexical_cast<float>(weighting.substr(pos+1,weighting.length()-10));
    }
    else {
        t_w = 1.0f;
    }
    if(t_w<0.0f || t_w > 1.0f)
        ERROR_exit("Negative tree weight!");
    return true;
}

template<class TREE>
inline bool stream2tree(std::istream &is, TREE &tree, idx2name &names, float &t_w) {
    idx2weight weights;
    return stream2tree(is, tree, names, weights, t_w);
}

inline bool stream2constr(std::istream &is, std::vector<std::string> &list) {    
    char c;
    NS_input::Input input = &is;
    if (!input.nextAnyChar(c)) return false;
    unsigned int lin = UINT_MAX; // Last   
    for (;;) {        
        if (c == ';') {          
            break;
        } else
        if (c == ',') {
            //do nothing
        } else { // name                   
            input.pushBack(c);
            std::string name = input.getName();            
            if (name.empty()) ERROR_return("problem in the constraint file " << input.getLastPos());
            list.push_back(name);
        } 
        if (!input.nextAnyChar(c)) return false;        
    }    
    
    return true;
}

template<class TREE>
inline bool stream2tree(std::istream &is, TREE &tree) {
    idx2name names;
    return stream2tree(is, tree, names);
}

enum tree_format{ROOTED,UNROOTED};


// output a tree in newick format
template<class TREE, class WEIGHTS>
bool tree2newick(std::ostream &os, TREE &tree, idx2name &names, WEIGHTS &weights, unsigned int root, tree_format flag) {
    if (tree.empty()) return false; // tree is empty
    if (flag == UNROOTED) os << "[&U]";
    const unsigned int save_root = tree.root;
    tree.root = root;    
    TREE_DFS2(v,tree) {        
        if (!tree.is_leaf(v.idx)) {
            if (v.direction == PREORDER) os << '(';
            if (v.direction == INORDER) os << ',';
            if (v.direction == POSTORDER) os << ')';
        }
        
        if (v.direction == POSTORDER) {            
            { // print name
                idx2name::iterator itr = names.find(v.idx);
                if (itr != names.end()) os << NS_input::getlegalstring(itr->second);
            }
            { // print weight(s)
                typename WEIGHTS::iterator itr = weights.find(v.idx);
                if (itr != weights.end()) itr->second.display_in_newick(os);
            }            
        }
    }
    os << ';';
    tree.root = save_root;
    return true;
}

// write a tree to a string stream
// example:
template<class TREE, class WEIGHTS>
inline bool tree2newick(std::ostream &os, TREE &tree, idx2name &names, WEIGHTS &weights) {
    if (tree.is_unrooted()) { // unrooted tree
        return tree2newick(os, tree, names, weights, 0, UNROOTED);
    } else { // rooted tree
        return tree2newick(os, tree, names, weights, tree.root, ROOTED);
    }
}
template<class TREE>
inline bool tree2newick(std::ostream &os, TREE &tree, idx2name &names) {
    idx2weight weights;
    return tree2newick(os, tree, names, weights);
}
template<class TREE>
inline bool tree2newick(std::ostream &os, TREE &tree) {
    idx2name names;
    return tree2newick(os, tree, names);
}

// read multiple elements from a string stream
// example: bird,frog,plant;
template<class T>
bool stream2names(std::istream &is, T &names) {
    char c;
    NS_input::Input input = &is;
    for (;;) {
        std::string name = input.getName();
        if (!name.empty()) names.push_back(name);
        if (!input.nextAnyChar(c)) return false;
        if (c == ';') break;
        else if (c == ',') {}
        else input.pushBack(c);
    }
    return true;
}

// write multiple elements to a string stream
// example: bird,frog,plant;
template<class T>
bool names2stream(std::ostream &os, T &names) {
    bool first = true;
    BOOST_FOREACH(const std::string &s, names) {
        if (!first) os << ','; first = false;
        os << s;
    } os << ';';
    return true;
}

// read multiple pairwise elements from a string stream
// example 1:
//   bird  ACGTS
//   frog  AATGG
//   plant TTGCG
//   ;
//
// example 2:
//   bird  ACGTS,frog AATGG,plant TTGCG;
inline bool stream2table(std::istream &is, std::vector<std::pair<std::string,std::string> > &tab) {
    char c;
    NS_input::Input input = &is;
    for (;;) {
        std::string key = input.getName();
        if (key.empty()) return false;
        std::string val = input.getName();
        if (val.empty()) return false;
        tab.push_back(std::pair<std::string,std::string>(key,val));
        if (!input.nextAnyChar(c)) return false;
        if (c == ';') break;
        else if (c == ',') {}
        else input.pushBack(c);
    }
    return true;
}

// write multiple pairwise elements to a string stream
// example 1:
//   bird  ACGTS
//   frog  AATGG
//   plant TTGCG
//   ;
//
// example 2: with delimiter1 = " " and delimiter2 = ","
//   bird ACGTS,frog AATGG,plant TTGCG;
template<class T>
bool table2stream(std::ostream &os, T &tab, std::string delimiter1 = " ", std::string delimiter2 = "") {
    BOOST_FOREACH(typename T::value_type &w, tab) {
        os << w.first << delimiter1 << w.second << delimiter2 << std::endl;
    } os << ';';
    return true;
}


}

#endif
