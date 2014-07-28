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

#ifndef TREE_SUBTREE_INFO_H
#define TREE_SUBTREE_INFO_H

#include "common.h"
#include "tree.h"
#include "tree_traversal.h"
#include <vector>
#include <boost/foreach.hpp>

namespace aw {

using namespace std;

// provides constant time access to subtree relational information
// accepts unrooted trees and rooted trees
// - subtree containment
// - subtree size nodes/leaves/internal nodes
// O(n) precomputation
// O(1) lookup
class SubtreeInfo {
    protected: class NodeVals {
        public: unsigned int begin, end;
//         public: NodeVals() { begin=UINT_MAX; end=UINT_MAX; }
    };
    protected: vector<NodeVals> ranges;
    protected: vector<unsigned int> leaves;
    protected: vector<unsigned int> internal;
    protected: unsigned int node_num;
    protected: unsigned int leaf_num;
    protected: unsigned int internal_num;
    protected: inline void free() {
        ranges.clear();
        leaves.clear();
        internal.clear();
    }
    public: template<class TREE> bool create(TREE &t) {
        free();
        if (t.node_size() == 0) return true;
        ranges.resize(t.node_size());
        leaves.resize(t.node_size());
        internal.resize(t.node_size());
        unsigned int oldroot = t.root;
        if (!t.is_rooted()) t.root = 0;
        {
            node_num = 0;
            leaf_num = 0;
            internal_num = 0;
            TREE_DFS2(v,t) {
                switch (v.direction) {
                    case PREORDER: {
                        ranges[v.idx].begin = node_num;
                        leaves[v.idx] = leaf_num;
                        internal[v.idx] = internal_num;
                        if (t.is_leaf(v.idx)) leaf_num++;
                        else internal_num++;
                    }break;
                    case INORDER: {
                        node_num++;
                    } break;
                    case POSTORDER: {
                        ranges[v.idx].end = node_num;
                        leaves[v.idx] = leaf_num - leaves[v.idx];
                        internal[v.idx] = internal_num - internal[v.idx];
                    } break;
                    default: break;
                }
            }
        }
        t.root = oldroot;
        return true;
    }

    // number of leafnodes of a subtree
    // u = subtree root node; p = parent of u
    public: inline unsigned int leaf_size(const unsigned int u, const unsigned int p) {
        if (!in_range(p,u)) return leaves[u];
        else return leaf_num - leaves[p];
    }

    // number of internal nodes of a subtree
    // u = subtree root node; p = parent of u
    public: inline unsigned int internal_node_size(const unsigned int u, const unsigned int p) {
        if (!in_range(p,u)) return internal[u];
        else return internal_num - internal[p];
    }

    // number of nodes of a subtree
    // u = subtree root node; p = parent of u
    public: inline unsigned int subtree_size(const unsigned int u, const unsigned int p) {
        return (leaf_size(u,p) + internal_node_size(u,p));
    }

    // true if node u is contained in subtree defined by v,pv
    // u,v = subtree root nodes; pv = parent of v;
    public: inline bool is_contained(const unsigned int u, const unsigned int v, const unsigned int pv) {
        if (in_range(pv,v)) {
            return !in_range(u,pv);
        }
        return in_range(u,v);
    }

    // true if the range of u is enclosed in the range of v
    protected: inline bool in_range(const unsigned int u, const unsigned int v) {
//         P(u);
//         P(ranges[u].begin);
//         P(ranges[u].end);
//         P(v);
//         P(ranges[v].begin);
//         P(ranges[v].end);
//         P(((ranges[u].begin >= ranges[v].begin) && (ranges[u].end <= ranges[v].end)));
        if (v == NONODE) return true;
        if (u == NONODE) return false;
        return ((ranges[u].begin >= ranges[v].begin) && (ranges[u].end <= ranges[v].end));
    }
};

// provides constant time access to subtree relational information for rooted trees
// - subtree containment
// - parent, sibling
// - subtree size
// - number of nodes
// O(n) precomputation
// O(1) lookup
template<class TREE>
class SubtreeInfoRooted {
    protected: unsigned int node_size;
    protected: unsigned int *parents;
    protected: struct NodeVals {
        unsigned int begin, end;
    };
    protected: NodeVals *ranges;
    protected: TREE *tree_ptr;
    public: SubtreeInfoRooted() {
        init();
    }
    public: ~SubtreeInfoRooted() {
        free();
    }
    public: SubtreeInfoRooted(const SubtreeInfoRooted &r) { // copy constructor
        init();
        copy(r);
    }
    public: SubtreeInfoRooted& operator=(const SubtreeInfoRooted& r) { // assign operator
        if (this != &r) {
            free();
            copy(r);
        }
        return *this;
    }
    protected: inline void init() {
        parents = NULL;
        ranges = NULL;
        tree_ptr = NULL;
        node_size = 0;
    }
    protected: inline void copy(const SubtreeInfoRooted &r) {
        tree_ptr = r.tree_ptr;
        if (r.tree_ptr != NULL) {
            node_size = r.node_size;
            parents = new unsigned int[node_size]; memcpy(parents, r.parents, node_size * sizeof(unsigned int));
            ranges = new NodeVals[node_size]; memcpy(ranges, r.ranges, node_size * sizeof(NodeVals));
        }
    }
    protected: inline void free() {
        if (parents != NULL) delete [] parents;
        if (ranges != NULL) delete [] ranges;
        init();
    }
    public: inline void clear() {
        free();
    }
    public: bool create(TREE &t) {
        free();
        tree_ptr = &t;
        if (t.is_unrooted()) {
            WARNING("rooted tree expected");
            return false;
        }
        const unsigned int size = t.node_size();
        if (node_size < size) {
            free();
            tree_ptr = &t;
            node_size = size;
            parents = new unsigned int[node_size];
            ranges = new NodeVals[node_size];
        }
        unsigned int node_num = 0;
        TREE_DFS2(v,t) {
            switch (v.direction) {
                case PREORDER: {
                    parents[v.idx] = v.parent;
                    ranges[v.idx].begin = node_num;
                } break;
                case INORDER: {
                    node_num++;
                } break;
                case POSTORDER: {
                    ranges[v.idx].end = node_num;
                } break;
                default: break;
            }
        }
        return true;
    }

    // true if node u is contained in subtree with root v
    public: inline bool is_contained(const unsigned int u, const unsigned int v) {
        return in_range(u,v);
    }

    // the parent node
    public: inline unsigned int parent(const unsigned int u) {
        if (u == NONODE) return NONODE;
        return parents[u];
    }

    // sibling node
    // in case of multipe siblings the first one is picked
    public: inline unsigned int sibling_binary(const unsigned int u) {
        const unsigned int p = parent(u);
        const unsigned int pp = parent(p);
        cout<< "parent of u=p:"<<p<< "parent of p=pp"<<pp;

        BOOST_FOREACH(const unsigned int &adj,tree_ptr->adjacent(p)) {
            if ((adj != pp) && (adj != u)) {
                 cout<< "sibling of"<<u<<adj;
                return adj;
            }cout<< "parent of u=p:"<<p<< "parent of p=pp"<<pp;

        }
        return NONODE;
    }

    // subtree node size (including root of the subtree)
    public: inline unsigned int subtree_size(const unsigned int u) {
        return ranges[u].end - ranges[u].begin;
    }

//     // return the nodes along the path from node u to v
//     public: inline void get_path(const unsigned int u, const unsigned int v, std::vector<unsigned int> &path) {
//         unsigned int t = u;
//         for (t = u; !is_contained(v,t); t = parent(t)) path.push_back(t);
//         if (t != u) path.push_back(t);
//         std::vector<unsigned int> rpath;
//         for (t = v; !is_contained(u,t); t = parent(t)) rpath.push_back(t);
//         for (unsigned int i=rpath.size(),iEE=0; i>iEE; --i) path.push_back(rpath[i-1]);
//     }

    // true if the range of u is enclosed in the range of v
    protected: inline bool in_range(const unsigned int u, const unsigned int v) {
        if (v == NONODE) return true;
        if (u == NONODE) return false;
        return ((ranges[v].begin <= ranges[u].begin) && (ranges[u].end <= ranges[v].end));
    }
};

template<class TREE>
class SubtreeParent {
    protected: unsigned int parents_size;
    protected: unsigned int *parents;
    protected: typedef TREE tree_type;
    protected: tree_type *tree_ptr;
    public: SubtreeParent() {
        init();
    }
    public: SubtreeParent(const SubtreeParent &r) { // copy constructor
        init();
        copy(r);
    }
    public: SubtreeParent& operator=(const SubtreeParent& r) { // assign operator
        if (this != &r) {
            this->free();
            copy(r);
        }
        return *this;
    }
    protected: inline void init() {
        tree_ptr = NULL;
        parents = NULL;
        parents_size = 0;
    }
    protected: inline void copy(const SubtreeParent &r) {
        tree_ptr = r.tree_ptr;
        if (r.tree_ptr != NULL) {
            parents_size = r.parents_size;
            parents = new unsigned int[parents_size];
            memcpy(parents, r.parents, parents_size * sizeof(unsigned int));
        }
    }
    public: ~SubtreeParent() {
        free();
    }
    protected: inline void free() {
        if (parents != NULL) delete [] parents;
        init();
    }
    public: inline void create(tree_type &t) {
        free();
        tree_ptr = &t;
        update();
    }

    public: inline void tPtrUpdate(tree_type &t) {        
        tree_ptr = &t;        
    }

    public: inline tree_type tPtrReturn() {   //Added by Ruchi
        return *tree_ptr;
    }

    public: inline void updateTree(const unsigned int x, const unsigned int x1, const unsigned int y, const unsigned int y1) {   //Added by Ruchi
        tree_ptr->remove_edge(x,x1);
        tree_ptr->remove_edge(y,y1);
        tree_ptr->add_edge(x,y1);
        tree_ptr->add_edge(y,x1);
    }


    public: inline void update() {
        TREE &t = *tree_ptr;
        if (t.is_unrooted()) {
            ERROR_exit("rooted tree expected");
        }
        const unsigned int size = t.node_size();
        if (parents_size < size) {
            tree_type *ptr = tree_ptr;
            free();
            tree_ptr = ptr;
            parents_size = size;
            parents = new unsigned int[parents_size];
        }
        TREE_PREORDER2(v,t) {
            parents[v.idx] = v.parent;
        }
    }
    public: inline void update(const unsigned int v, const unsigned int vp) {
        parents[v] = vp;
    }

    // the parent node
    public: inline unsigned int parent(const unsigned int u) {
        if (u == NONODE) return NONODE;
        return parents[u];
    }
    // sibling node
    // in case of multipe siblings the first one is picked
    public: inline unsigned int sibling_binary(const unsigned int u) {
        const unsigned int p = parent(u);        
        const unsigned int pp = parent(p);   
       
        //BOOST_FOREACH(const unsigned int &adj,tree_ptr->adjacent_vector(p))
            //MSG("Adjacent of "<<p<<":"<<adj<<" ");
        
        BOOST_FOREACH(const unsigned int &adj,tree_ptr->adjacent(p)) {            
            if ((adj != pp) && (adj != u)) {
                //MSG("Returning sibling of"<<u<<" "<<adj<<" ");
                return adj;
            }            
        }
        return NONODE;
    }
    
    public: inline bool operator==(const SubtreeParent<TREE> &r) {
        if (parents_size != r.parents_size) {
            return false;
        }
        if (tree_ptr != r.tree_ptr) {
            return false;
        }
        for (unsigned int i=0; i<parents_size; ++i) {
            if (parents[i] != r.parents[i]) {
                P(i << ' ' << parents[i] << "!=" << r.parents[i]);
                return false;
            }
        }
        return true;
    }
    public: inline bool operator!=(const SubtreeParent<TREE> &r) {
        return !(*this == r);
    }
};

// provides constant time access to subtree relational information
// accepts unrooted trees and rooted trees
// - subtree node size
// O(n) precomputation
// O(1) lookup
class SubtreeSizes {
    protected: class NodeVals {
        public: unsigned int begin, end;
    };
    protected: vector<NodeVals> ranges;
    protected: unsigned int node_num;
    protected: inline void free() {
        ranges.clear();
    }
    public: inline void clear() {
        free();
    }
    public: template<class TREE> bool create(TREE &t) {
        free();
        if (t.node_size() == 0) return true;
        ranges.resize(t.node_size());
        unsigned int oldroot = t.root;
        if (!t.is_rooted()) t.root = 0;
        node_num = 0;
        TREE_DFS2(v,t) {
            switch (v.direction) {
                case PREORDER: {
                    ranges[v.idx].begin = node_num;
                } break;
                case INORDER: {
                    node_num++;
                } break;
                case POSTORDER: {
                    ranges[v.idx].end = node_num;
                } break;
                default: break;
            }
        }
        t.root = oldroot;
        return true;
    }

    // number of nodes of a subtree
    // u = subtree root node; p = parent of u
    public: inline unsigned int subtree_size(const unsigned int u, const unsigned int p) {
        return (!in_range(p,u)) ? ranges[u].end - ranges[u].begin : node_num - ranges[u].end + ranges[u].begin;
    }

    // true if node u is contained in subtree defined by v,pv
    // u,v = subtree root nodes; pv = parent of v;
    public: inline bool is_contained(const unsigned int u, const unsigned int v, const unsigned int pv) {
        if (in_range(pv,v)) {
            return !in_range(u,pv);
        }
        return in_range(u,v);
    }

    // true if the range of u is enclosed in the range of v
    protected: inline bool in_range(const unsigned int u, const unsigned int v) {
        if (v == NONODE) return true;
        if (u == NONODE) return false;
        return ((ranges[u].begin >= ranges[v].begin) && (ranges[u].end <= ranges[v].end));
    }
};

// provides constant time access to subtree sizes
// accepts rooted trees
// - subtree node size
// O(n) precomputation
// O(1) lookup
class SubtreeSizesRooted {
    protected: class NodeVals {
        public: unsigned int begin, end;
    };
    protected: vector<NodeVals> ranges;
    protected: vector<unsigned int> leaves;
    protected: vector<unsigned int> internal;
    protected: unsigned int node_num;
    protected: unsigned int leaf_num;
    protected: unsigned int internal_num;
    protected: inline void free() {
        ranges.clear();
        leaves.clear();
        internal.clear();
    }
    public: inline void clear() {
        free();
    }
    public: template<class TREE> bool create(TREE &t) {
        free();
        if (t.node_size() == 0) return true;
        if (t.is_unrooted()) {
            WARNING("rooted tree expected");
            return false;
        }
        ranges.resize(t.node_size());
        leaves.resize(t.node_size());
        internal.resize(t.node_size());
        unsigned int oldroot = t.root;
        if (!t.is_rooted()) t.root = 0;
        {
            node_num = 0;
            leaf_num = 0;
            internal_num = 0;
            TREE_DFS2(v,t) {
                switch (v.direction) {
                    case PREORDER: {
                        ranges[v.idx].begin = node_num;
                        leaves[v.idx] = leaf_num;
                        internal[v.idx] = internal_num;
                        if (t.is_leaf(v.idx)) leaf_num++;
                        else internal_num++;
                    }break;
                    case INORDER: {
                        node_num++;
                    } break;
                    case POSTORDER: {
                        ranges[v.idx].end = node_num;
                        leaves[v.idx] = leaf_num - leaves[v.idx];
                        internal[v.idx] = internal_num - internal[v.idx];
                    } break;
                    default: break;
                }
            }
        }
        t.root = oldroot;
        return true;
    }

    // number of leafnodes of a subtree
    // u = subtree root node
    public: inline unsigned int leaf_size(const unsigned int u) {
        return leaves[u];
    }

    // number of internal nodes of a subtree
    // u = subtree root node
    public: inline unsigned int internal_node_size(const unsigned int u) {
        return internal[u];
    }

    // number of nodes of a subtree
    // u = subtree root node
    public: inline unsigned int subtree_size(const unsigned int u) {
        return (leaf_size(u) + internal_node_size(u));
    }

    // true if node u is contained in subtree rooted at v
    // u,v = subtree root nodes
    public: inline bool is_contained(const unsigned int u, const unsigned int v) {
        return in_range(u,v);
    }

    // true if the range of u is enclosed in the range of v
    protected: inline bool in_range(const unsigned int u, const unsigned int v) {
        if (v == NONODE) return true;
        if (u == NONODE) return false;
        return ((ranges[u].begin >= ranges[v].begin) && (ranges[u].end <= ranges[v].end));
    }
};

// provides constant time access to cluster sizes
// accepts rooted trees only
// - subtree leafset size (duplicated labels are counted only once per leafset)
// O(n) precomputation
// O(1) lookup
class ClusterSizesRooted {
    protected: class NodeVals {
        public: unsigned int begin, end;
    };
    protected: vector<NodeVals> ranges;
    protected: unsigned int node_num;
    protected: inline void free() {
        ranges.clear();
    }
    public: template<class TREE> bool create(TREE &t, aw::idx2name &taxa) {
        free();
        if (t.node_size() == 0) return true;
        if (t.is_unrooted()) {
            WARNING("rooted tree expected");
            return false;
        }
        ranges.resize(t.node_size());
        const unsigned int oldroot = t.root;
        if (!t.is_rooted()) t.root = 0;
        // find left and right nodes with identical labels
        std::vector<unsigned int> next_left(t.node_size(),NONODE);
        boost::unordered_map<std::string, unsigned int> last_occurrence;
        TREE_INORDER(v,t) if (t.is_leaf(v)) {
            const std::string &name = taxa[v];
            boost::unordered_map<std::string, unsigned int>::iterator itr = last_occurrence.find(name);
            if (itr != last_occurrence.end()) {
                const unsigned int &last = itr->second;
                next_left[v] = last;
            }
            last_occurrence[name] = v;
        }
        // create LCA look-up table
        aw::LCA t_lca; t_lca.clear(); t_lca.create(t);
        // determine the clusters where multi-labels become relevant
        std::vector<unsigned int> dec_size(t.node_size(),0);
        TREE_FOREACHNODE(v,t) {
            const unsigned int &u = next_left[v];
            if (u != NONODE) {
                const unsigned int &ca = t_lca.lca(v,u);
                ++dec_size[ca];
            }
        }
        TREE_POSTORDER2(v,t) {
            BOOST_FOREACH(const unsigned int &c,t.children(v.idx,v.parent)) dec_size[v.idx] += dec_size[c];
        }
        // compute cluster sizes
        node_num = 0;
        TREE_DFS2(v,t) {
            switch (v.direction) {
                case PREORDER: {
                    ranges[v.idx].begin = node_num;
                } break;
                case INORDER: {
                    if (t.is_leaf(v.idx)) node_num++;
                } break;
                case POSTORDER: {
                    ranges[v.idx].end = node_num - dec_size[v.idx];
                } break;
                default: break;
            }
        }
        t.root = oldroot;
        return true;
    }

    // number of nodes of a subtree
    // u = subtree root node
    public: inline unsigned int subtree_size(const unsigned int u) {
        return ranges[u].end - ranges[u].begin;
    }

    // true if node u is contained in subtree rooted at v
    // u,v = subtree root nodes
    public: inline bool is_contained(const unsigned int u, const unsigned int v) {
        return in_range(u,v);
    }

    // true if the range of u is enclosed in the range of v
    protected: inline bool in_range(const unsigned int u, const unsigned int v) {
        if (v == NONODE) return true;
        if (u == NONODE) return false;
        return ((ranges[u].begin >= ranges[v].begin) && (ranges[u].end <= ranges[v].end));
    }
};

} // end namespace

#endif
