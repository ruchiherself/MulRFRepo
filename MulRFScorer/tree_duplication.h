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

#ifndef TREE_DUPLICATIONS_H
#define TREE_DUPLICATIONS_H

#include "common.h"
#include "tree.h"
#include "tree_traversal.h"
#include "tree_LCA.h"
#include "tree_name_map.h"
#include "tree_LCA_mapping.h"
#include "tree_subtree_info.h"
#include <limits.h>
#include <boost/random.hpp>

namespace aw {

using namespace std;

template<class TREE>
class TreeClusters {
    protected: unsigned int node_size;
    protected: unsigned int *clusters;
    protected: typedef TREE tree_type;
    protected: tree_type *s_tree_ptr;
    protected: tree_type *g_tree_ptr;
    protected: aw::TreetaxaMap *s_nmap;
    protected: aw::TreetaxaMap *g_nmap;
    public: TreeClusters() {
        init();
    }
    public: TreeClusters(const TreeClusters &r) { // copy constructor
        init();
        copy(r);
    }
    public: TreeClusters& operator=(const TreeClusters& r) { // assign operator
        if (this != &r) {
            this->free();
            copy(r);
        }
        return *this;
    }
    protected: inline void init() {
        s_tree_ptr = NULL;
        g_tree_ptr = NULL;
        g_nmap = NULL;
        s_nmap = NULL;
        clusters = NULL;
        node_size = 0;
    }
    protected: inline void copy(const TreeClusters &r) {
        s_tree_ptr = r.s_tree_ptr;
        g_tree_ptr = r.g_tree_ptr;
        g_nmap = r.g_nmap;
        s_nmap = r.s_nmap;
        if (r.s_tree_ptr != NULL) {
            node_size = r.node_size;
            clusters = new unsigned int[node_size];
            memcpy(clusters, r.clusters, node_size * sizeof(unsigned int));
        }
    }
    public: ~TreeClusters() {
        free();
    }
    protected: inline void free() {
        if (clusters != NULL) delete [] clusters;
        init();
    }

    public: inline void create(tree_type &st, tree_type &gt, aw::TreetaxaMap &gmap, aw::TreetaxaMap &smap, aw::LCAmapping &map) {
        free();
        g_tree_ptr = &gt;
        s_tree_ptr = &st;
        s_nmap = &smap;
        g_nmap = &gmap;
        
        unsigned int size = st.node_size();
        node_size = size;
        clusters = new unsigned int[size];

        unsigned int count;
        TREE_POSTORDER2(v,st){
            if (!st.is_leaf(v.idx)) {
                count = 0;
                BOOST_FOREACH(const unsigned int &c,st.children(v.idx,v.parent)) {                    
                    count = count + cluster(c); }                
                clusters[v.idx] = count;
                
            }
            else {
                if(map.mapping(v.idx)!=NONODE) clusters[v.idx] = 1;
                else clusters[v.idx] = 0;  }
            //std::cout<<"clst of "<<v.idx<<" "<<clusters[v.idx];
        }
    }

    public: inline void stPtrUpdate(tree_type &st) {
        s_tree_ptr = &st;
    }

    public: inline void gtPtrUpdate(tree_type &gt) {
        g_tree_ptr = &gt;
    }

    public: inline void smapUpdate(aw::TreetaxaMap &smap) {
        s_nmap = &smap;
    }

    public: inline void gmapUpdate(aw::TreetaxaMap &gmap) {
        g_nmap = &gmap;
    }

//    public: inline void update() {
//        TREE &st = *s_tree_ptr;
//        TREE &gt = *g_tree_ptr;
//        aw::TreetaxaMap &smap = *s_nmap;
//        aw::TreetaxaMap &gmap = *g_nmap;
//
//
//        unsigned int size = st.node_size;
//
//        if (node_size < size) {
//            tree_type *sptr = s_tree_ptr;
//            tree_type *gptr = g_tree_ptr;
//            free();
//            s_tree_ptr = sptr;
//            g_tree_ptr = gptr;
//            node_size = size;
//            clusters = new unsigned int[node_size];
//        }
//
//        unsigned int count;
//        TREE_POSTORDER2(v,st){
//            if (!st.is_leaf(v.idx)) {
//                count = 0;
//                BOOST_FOREACH(const unsigned int &c,st.children(v.idx,v.parent))
//                    count = count + st.return_clstSz(c);
//                clusters[v.idx] = count; }
//            else {
//                const unsigned int v_gid = smap.gid(v.idx);
//                if(gmap.exists(v_gid)) clusters[v.idx] = 1;
//                else clusters[v.idx] = 0;  }
//        }
//    }

    public: inline void update(const unsigned int v, const unsigned int clst) {
        clusters[v] = clst;
    }

    // the parent node
    public: inline unsigned int cluster(const unsigned int u) {
        if (u == NONODE) return NONODE;
        return clusters[u];
    }

};


// true if it is a gene duplication
inline bool compute_duplication(const unsigned int s_map, const unsigned int c_map_0, const unsigned int c_map_1) {
    if ((c_map_0 == NONODE) || (c_map_1 == NONODE)) return false;
    return ((c_map_0 == s_map) || (c_map_1 == s_map));
}
inline bool compute_duplication(LCAmapping &g_map, const unsigned int s, const unsigned int c0, const unsigned int c1) {
    const unsigned int &c_map_0 = g_map.mapping(c0);
    const unsigned int &c_map_1 = g_map.mapping(c1);
    if ((c_map_0 == NONODE) || (c_map_1 == NONODE)) return false;
    const unsigned int &s_map = g_map.mapping(s);

    //cout<<"mapping=="<<c_map_0<<c_map_1<<s_map;
    return ((c_map_0 == s_map) || (c_map_1 == s_map));
}

inline bool compute_duplication(Tree t, LCAmapping &g_map, const unsigned int s, const unsigned int ps) {
    //find children of s
    unsigned int ch[2];
    t.children(s, ps, ch);

    //cout<<"\n Node "<<s<<" parent "<<ps;
    //cout<<"Children of "<<s<<" are "<<ch[0]<<" "<<ch[1];

    bool res = compute_duplication(g_map, s, ch[0], ch[1]);
    return res;
}

// compute the gene duplications induced by a gene tree
template<class TREE>
inline unsigned int compute_duplications(TREE &g_tree, LCAmapping &g_map) {
    unsigned int dups = 0;
    TREE_POSTORDER2(v,g_tree) {
        if (!g_tree.is_leaf(v.idx)) {
            unsigned int ch[2];
            g_tree.children(v.idx,v.parent,ch);
            if (compute_duplication(g_map,v.idx,ch[0],ch[1])) ++dups;
            // const unsigned int &c_map_0 = g_map.mapping(ch[0]);
            // const unsigned int &c_map_1 = g_map.mapping(ch[1]);
            // if ((c_map_0 == NONODE) || (c_map_1 == NONODE)) continue;
            // const unsigned int &s_map = g_map.mapping(v.idx);
            // if ((c_map_0 == s_map) || (c_map_1 == s_map)) {
            //     ++dups;
            // }
        }
    }
    return dups;
}

// compute the RF score for input tree and supertree :By ruchi
template<class TREE>
inline float compute_rf_score(TREE &s_tree, TREE &g_tree, LCAmapping &s_map, std::pair<unsigned int,unsigned int> &node_count,unsigned int s_int, float wt) {
    unsigned int score = 0;
    unsigned int gmap;

    //check later
    //MSG("Setting initial score");
    TREE_FOREACHNODE(v,g_tree) g_tree.init_score(v);    

    TREE_POSTORDER2(v,s_tree)
        if (!s_tree.is_leaf(v.idx) && s_tree.root!=v.idx) {
            gmap = s_map.mapping(v.idx);            
            if(gmap!=NONODE && (s_tree.return_clstSz(v.idx) == g_tree.return_clstSz(gmap)))
                g_tree.incr_score(gmap,1);                  
        }
    TREE_FOREACHNODE(v,g_tree) {        
        if (!g_tree.is_leaf(v) && g_tree.root!=v)
            if(g_tree.return_score(v) == 0)
                score = score + 2;
    }
    
    score =  score + s_int - node_count.first;
     
    return score*wt;
}

// compute the RF score for input tree and supertree WITH TREECLUSTER :By ruchi
template<class TREE>
inline unsigned int compute_rf_score(TREE &s_tree, TREE &g_tree, LCAmapping &s_map, TreeClusters<aw::Tree> &sclst, unsigned int g_inodes, unsigned int s_inodes) {
    unsigned int score = 0;
    unsigned int gmap;

    TREE_POSTORDER2(v,s_tree) {        
        if(!s_tree.is_leaf(v.idx) && s_tree.root!=v.idx) {           
            std::vector<unsigned int> ch;
            s_tree.children(v.idx,v.parent,ch);
            int ct = 0;
            for(unsigned int j=0; j<ch.size(); ++j)
                if(s_map.mapping(ch[j])!=NONODE) ct +=1;           
            if(ct<2) continue;
            gmap = s_map.mapping(v.idx);            
            if(sclst.cluster(v.idx) != g_tree.return_clstSz(gmap))
                score = score + 2;
        }
    }    
    score =  score - g_inodes + s_inodes;
    //MSG("score is:"<<score);
    return score;
}

// compute the RF score for input tree and supertree WITH TREECLUSTER :By ruchi
template<class TREE>
inline unsigned int compute_rf_score(TREE &s_tree, TREE &g_tree, LCAmapping &s_map, unsigned int v, unsigned int pv, TreeClusters<aw::Tree> &sclst) {
    unsigned int gmap;

    //check later
    //MSG("Setting initial score");
    //TREE_FOREACHNODE(v,s_tree) s_tree.init_score(v);

    if (!s_tree.is_leaf(v) && s_tree.root!=v) {
        std::vector<unsigned int> ch; s_tree.children(v,pv,ch);
        if(s_map.mapping(ch[0])==NONODE || s_map.mapping(ch[1])==NONODE) {
            //std::cout<<" c & map"<<ch[0]<<ch[1]<<" "<<s_map.mapping(ch[0])<<s_map.mapping(ch[1]);
            return 0;   }
        gmap = s_map.mapping(v);
        //std::cout<<" gmap "<<gmap;
        //std::cout<<"sclst.cluster(v) "<<sclst.cluster(v);
        //std::cout<<"g_tree.return_clstSz(gmap)"<<g_tree.return_clstSz(gmap);
        if(sclst.cluster(v) != g_tree.return_clstSz(gmap))
            return 2;
    }  

    return 0;
}

// compute the gene duplications induced by multiple gene trees
template<class TREE>
inline unsigned int compute_duplications(std::vector<TREE> &g_trees, std::vector<LCAmapping> &g_maps) {
    unsigned int dups = 0;
    for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
        dups += compute_duplications(g_trees[i],g_maps[i]);
    }
    return dups;
}

// compute the gene duplications induced by a gene tree
template<class STREE,class GTREE>
inline unsigned int compute_duplications(STREE &s_tree, idx2name &s_names, GTREE &g_tree, idx2name &g_names) {
    TaxaMap taxamap;
    taxamap.insert(s_names);
    taxamap.insert(g_names);
    TreetaxaMap s_nmap; s_nmap.create(s_names,taxamap);
    TreetaxaMap g_nmap; g_nmap.create(g_names,taxamap);
    LCA lca; lca.create(s_tree);
    LCAmapping lca_map; lca_map.create(lca,s_nmap,g_tree,g_nmap);
    const unsigned int dups = compute_duplications(g_tree,lca_map);
    return dups;
}

// compute the gene duplications induced by multiple gene trees
template<class STREE,class GTREE>
inline unsigned int compute_duplications(STREE &s_tree, TreetaxaMap &s_nmap, std::vector<GTREE> &g_trees, std::vector<TreetaxaMap> &g_nmaps) {
    LCA s_lca; s_lca.create(s_tree);
    unsigned int dups = 0;
    for (unsigned int i=0,iEE=g_nmaps.size(); i<iEE; ++i) {
        GTREE &g_tree = g_trees[i];
        TreetaxaMap &g_nmap = g_nmaps[i];
        LCAmapping lca_map; lca_map.create(s_lca,s_nmap,g_nmap,g_tree);
        dups += compute_duplications(g_tree,lca_map);
    }
    return dups;
}

// compute the gene duplications induced by multiple gene trees
template<class STREE,class GTREE>
inline unsigned int compute_duplications(STREE &s_tree, idx2name &s_names, std::vector<GTREE> &g_trees, std::vector<idx2name> &g_names) {
    TaxaMap taxamap;
    taxamap.insert(s_names);
    TreetaxaMap s_nmap; s_nmap.create(s_names,taxamap);
    std::vector<TreetaxaMap> g_nmaps(g_names.size());
    for (unsigned int i=0,iEE=g_names.size(); i<iEE; ++i) {
        idx2name &n = g_names[i];
        taxamap.insert(n);
        g_nmaps[i].create(n,taxamap);
    }
    const unsigned int dups = compute_duplications(s_tree,s_nmap,g_trees,g_nmaps);
    return dups;
}

#ifndef AW_RANDOMGEN
boost::mt19937 rng;
#endif

// find the best SPR move on the subtree - there can be multiple equal ones, then only one of them is returned
// return
//   location = edge(u,v) for the location with lowest duplications
//   duplications = lowest duplications
// not thread-safe
template<class STREE,class GTREE>
bool bestSPRlocation(
    const unsigned int subtree, const unsigned int subtree_parent, STREE &s_tree,
    std::vector<GTREE> &g_trees, std::vector<aw::LCAmapping> &g_lmaps,
    std::pair<unsigned int, unsigned int> &location, unsigned int &duplications
) {
    const unsigned int &subtree_left = subtree;
    const unsigned int current_root = s_tree.root;
    if (subtree_left == current_root) return false; // someone wants me to place the tree inside the tree -> can't do that
    const bool move = !(subtree_parent == current_root); // unless the subtree is already placed at the root ...
    std::vector<unsigned int> subtree_parent_adj;
    s_tree.children(subtree_parent,subtree_left,subtree_parent_adj);
    if (move) { // move the subtree to the root
        // if (subtree_parent_adj.size() != 2) ERROR_exit("something is wrong here");
        s_tree.spr_to_root(subtree_left,subtree_parent);
    }
    bool location_change;
    { // Bansal, Eulenstein, Wehe algorithm to determine best SPR for subtree_left
        // update current LCA mapping
        aw::LCA s_lca; s_lca.create(s_tree); // compute LCAs for the species tree
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) g_lmaps[i].update_LCA_internals(s_lca,g_trees[i]); // update the LCA mapping for internal nodes of the gene trees
        static aw::SubtreeInfoRooted<aw::Tree> s_info; s_info.create(s_tree);
        // determine locations and gene duplication changes
        const unsigned int &s_root = subtree_parent;
        const unsigned int subtree_right = *s_tree.children(subtree_parent,subtree_left).begin();
        std::vector<unsigned int> dups_inc(s_tree.node_size(),0);
        std::vector<unsigned int> dups_dec(s_tree.node_size(),0);
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            aw::Tree &g_tree = g_trees[i];
            static std::vector<unsigned int> g_lmap2; // store the secondary LCA mappings of a gene tree (reuseable)
            if (g_tree.node_size() > g_lmap2.size()) g_lmap2.resize(g_tree.node_size());
            aw::LCAmapping &g_lmap = g_lmaps[i];
            TREE_DFS2(v,g_tree) {
                switch (v.direction) {
                    case aw::PREORDER: {
                        const unsigned int &v_map = g_lmap.mapping(v.idx);
                        if (v_map != s_root) { // no need to traverse further into the subtree
                            v.skip();
                        }
                    } break;
                    case aw::POSTORDER: {
                        const unsigned int &v_map = g_lmap.mapping(v.idx);
                        if (v_map == s_root) {
                            unsigned int ch[2]; g_tree.children(v.idx,v.parent,ch);
                            const unsigned int ch_map[2] = { g_lmap.mapping(ch[0]),g_lmap.mapping(ch[1]) };
                            bool case0;
                            // missing lca mapping, usually caused by missing taxa in the species tree
                            // no gene duplications will ever be lost or gained
                            if ((case0 = ch_map[0] == NONODE) || (ch_map[1] == NONODE)) {
                                g_lmap2[v.idx] = g_lmap2[ch[case0 ? 1 : 0]];
                                break;
                            }
                            const bool to_root_0 = (ch_map[0] == s_root);
                            const bool to_root_1 = (ch_map[1] == s_root);
                            // parent and both children map to the root
                            // no gene duplications will every be lost or gained
                            if (to_root_0 && to_root_1) {
//                                 const unsigned int ch_map2[2] = { g_lmap2[ch[0]],g_lmap2[ch[1]] };
//                                 const unsigned int v_map2 = s_lca.lca(ch_map2[0],ch_map2[1]);
                                const unsigned int &ch_map2_0 = g_lmap2[ch[0]];
                                const unsigned int &ch_map2_1 = g_lmap2[ch[1]];
                                const unsigned int v_map2 = s_lca.lca(ch_map2_0,ch_map2_1);
                                g_lmap2[v.idx] = v_map2;
                                break;
                            }
                            const bool in_r_0 = s_info.is_contained(ch_map[0],subtree_right);
                            const bool in_r_1 = s_info.is_contained(ch_map[1],subtree_right);
                            // parent and one child map to the root, one child maps into the right subtree
                            if ((case0 = to_root_0 && in_r_1) || (to_root_1 && in_r_0)) {
//                                 const unsigned int ch_map2[2] = {
//                                     case0 ? g_lmap2[ch[0]] : ch_map[0],
//                                     case0 ? ch_map[1] : g_lmap2[ch[1]]
//                                 };
//                                 const unsigned int &s_l = ch_map2[case0 ? 0 : 1]; // virtually maps to the right subtree
//                                 const unsigned int &s_r = ch_map2[case0 ? 1 : 0]; // maps to the right subtree
//                                 const unsigned int &s_p = g_lmap2[v.idx] = s_lca.lca(ch_map2[0],ch_map2[1]); // virtually maps to the lca of both child mappings
                                const unsigned int &ch_map2_0 = case0 ? g_lmap2[ch[0]] : ch_map[0];
                                const unsigned int &ch_map2_1 = case0 ? ch_map[1] : g_lmap2[ch[1]];
                                const unsigned int &s_l = case0 ? ch_map2_0 : ch_map2_1; // virtually maps to the right subtree
                                const unsigned int &s_r = case0 ? ch_map2_1 : ch_map2_0; // maps to the right subtree
                                const unsigned int &s_p = g_lmap2[v.idx] = s_lca.lca(ch_map2_0,ch_map2_1); // virtually maps to the lca of both child mappings
                                const unsigned int &s_pp = s_info.parent(s_p); // parent of s_p
                                unsigned int s_pch[2]; s_tree.children(s_p, s_pp, s_pch); // children of s_p
                                if (s_r == s_p) { // duplication remains because of the right child mapping
                                } else
                                if (s_info.is_contained(s_l,s_pch[0])) {
                                    const unsigned int &s_ll = s_pch[0]; // left child of s_p
                                    ++dups_dec[s_ll];
                                } else
                                if (s_info.is_contained(s_l,s_pch[1])) {
                                    const unsigned int &s_ll = s_pch[1]; // left child of s_p
                                    ++dups_dec[s_ll];
                                }
                                break;
                            }
                            const bool in_l_0 = s_info.is_contained(ch_map[0],subtree_left);
                            const bool in_l_1 = s_info.is_contained(ch_map[1],subtree_left);
                            // parent maps to the root, one child maps into the left subtree, one child maps into the right subtree
                            // then gene duplication occurs when moving below the mapped node in the right subtree
                            if ((case0 = in_l_0 && in_r_1) || (in_l_1 && in_r_0)) {
                                const unsigned int &v_map2 = ch_map[case0 ? 1 : 0];
                                g_lmap2[v.idx] = v_map2;
                                ++dups_inc[v_map2];
                                break;
                            }
                            // parent and one child map to the root, one child maps into the left subtree
                            // then no gene duplications will every be lost
                            if ((case0 = to_root_0 && in_l_1) || (to_root_1 && in_l_0)) {
                                const unsigned int &v_map2 = g_lmap2[ch[case0 ? 0 : 1]];
                                g_lmap2[v.idx] = v_map2;
                                break;
                            }
                            // something is wrong
                            ERROR_exit("something is wrong");
                        }
                    } break;
                    default: break;
                }
            }
        }
        { // accumulate gene duplication changes
            unsigned int d = duplications = compute_duplications(g_trees,g_lmaps); // compute initial gene duplications
            // static std::vector<std::pair<unsigned int,unsigned int> > candidates;
            typedef std::pair<unsigned int,unsigned int> candidates_item;
            static util::vector<candidates_item> candidates; candidates.set_min_size(s_tree.node_size());
            for (aw::Tree::iterator_dfs v=s_tree.begin_dfs(subtree_right,s_root),vEE=s_tree.end_dfs(); v!=vEE; ++v) {
                switch (v.direction) {
                    case aw::PREORDER: {
                        d -= dups_dec[v.idx];
                        if (d < duplications) {
                            duplications = d;
                            candidates.clear();
                        }
                        if (d == duplications) {
                            candidates.push_back(std::pair<unsigned int,unsigned int>(v.idx,v.parent));
                        }
                        // { // brute force verification of the algorithm
                        //     aw::Tree s_tree2 = s_tree;
                        //     s_tree2.spr_from_root(subtree_left,s_info.sibling_binary(subtree_left),v.idx,v.parent);
                        //     aw::LCA s_lca2; s_lca2.create(s_tree2);
                        //     std::vector<aw::LCAmapping> g_lmaps2 = g_lmaps;
                        //     for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) g_lmaps2[i].update_LCA_internals(s_lca2,g_trees[i]); // update the LCA mapping for internal nodes of the gene trees
                        //     unsigned int dd=compute_duplications(g_trees,g_lmaps2);
                        //     if (d != dd) std::cout << v.idx << ':' << d << '(' << dd << ") ";
                        // }
                        d += dups_inc[v.idx];
                    } break;
                    case aw::POSTORDER: {
                        d -= dups_inc[v.idx];
                        d += dups_dec[v.idx];
                    } break;
                    default: break;
                }
            }
            { // any better location found?
                location_change = true;
                const unsigned int u = subtree_parent_adj[0];
                const unsigned int v = (subtree_parent_adj.size() == 1) ? s_root : subtree_parent_adj[1];
                for (unsigned int i=0,iEE=candidates.size(); i<iEE; ++i) {
                    candidates_item &w = candidates[i];
                    if (
                        ((u == w.first) && (v == w.second)) ||
                        ((v == w.first) && (u == w.second))
                    ) location_change = false;
                }
            }
            boost::uniform_int<> range(0,candidates.size()-1);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(aw::rng, range);
            location = candidates[die()];
            candidates.clear();
        }
    }
    if (move) { // move the subtree back to its original location
        s_tree.spr_from_root(subtree_left,current_root,subtree_parent_adj[0],subtree_parent_adj[1]);
    }
    return location_change;
}

// move a subtree to the root; SPR move to the root; tree is assumed to be fully binary
template<class TREE>
inline bool move2root_binary(TREE &t, const unsigned int u, const unsigned int pu) {
    // unusual
    if ((u == NONODE) || (pu == NONODE)) {
        ERROR_exit("Something is wrong here");
    }
    // necessary?
    if (pu == t.root) return false;
    // prune subtree with parent
    unsigned int ch[2];
    t.children(pu,u,ch);
    
    t.remove_edge(pu,ch[0]);
    //cout<<"|Remove edge "<<pu<<" "<<ch[0];
    t.remove_edge(pu,ch[1]);
    //cout<<"|Remove edge "<<pu<<" "<<ch[1];
    t.add_edge(ch[0],ch[1]);
    //cout<<"|Add edge "<<ch[0]<<" "<<ch[1];
    // regraft subtree to root
    t.add_edge(pu,t.root);
    //cout<<"|Add edge "<<pu<<" "<<t.root<<"|";
    // define the new root
    t.root = pu;
    return true;
}

// move a subtree from the root into an edge; SPR move from the root; tree is assumed to be fully binary
template<class TREE>
inline bool move2edge_binary(TREE &t, const unsigned int subtree, const unsigned int parent, const unsigned int v, const unsigned int pv) {
    // unusual
    if ((subtree == NONODE) || (parent == NONODE) || (v == NONODE) || (pv == NONODE)) {
        ERROR_exit("Something is wrong here");
    }
    // necessary?
    if (pv == NONODE) return false;
    if (parent == pv) return false;
    if (parent == v) return false;
    // prune subtree with parent/root

    std::vector<unsigned int> adj_par;  t.children(parent,subtree,adj_par);
    //const unsigned int su = *t.children(parent,subtree).begin();
    t.remove_edge(parent,adj_par[0]);    
    //cout<<"|Removed edge- "<<parent<<" "<<adj_par[0];

    // define the new root
    if(t.root == parent) {
        t.root = pv;
        //cout<<"|Root is "<<t.root;
     }        
    else{
        t.remove_edge(parent,adj_par[1]);
        //cout<<"|Removed edge- "<<parent<<" "<<adj_par[1];
        t.add_edge(adj_par[0],adj_par[1]);
        //cout<<"|Add edge- "<<adj_par[0]<<" "<<adj_par[1];
    }       

     t.add_edge(parent,pv);
     //cout<<"|Add edge- "<<parent<<" "<<pv<<"|";
    // regraft subtree into edge
     t.remove_edge(v,pv);
     //cout<<"|Removed edge- "<<v<<" "<<pv;
    
     t.add_edge(parent,v);
     //cout<<"|Add edge- "<<parent<<" "<<v<<"|";
    
    return true;
}

// move a subtree back to the root form the edge; Reverse of the above opeartion; tree is assumed to be fully binary :added by ruchi
template<class TREE>
inline bool REVmove2edge_binary(TREE &t, const unsigned int subtree, const unsigned int parent, const unsigned int v, const unsigned int pv, const unsigned int ppv) {
    // unusual
    if ((subtree == NONODE) || (parent == NONODE) || (v == NONODE) || (pv == NONODE)) {
        ERROR_exit("Something is wrong here");
    }
    // necessary?
    if (pv == NONODE) return false;
    if (parent == pv) return false;
    if (parent == v) return false;

    t.remove_edge(parent,pv);
    //cout<<"|Remove edge- "<<parent<<" "<<pv;
    t.remove_edge(parent,v);
    //cout<<"|Remove edge- "<<parent<<" "<<v;

    if(t.root == pv){
        t.root = parent;
        //cout<<"|New root- "<<t.root;
    }
    else {
        //was added by ruchi
        t.add_edge(ppv,parent);
        //cout<<"|Add edge- "<<ppv<<" "<<parent;
        t.remove_edge(ppv,pv);
        //cout<<"|Remove edge- "<<ppv<<" "<<pv;
    }

    // prune subtree with parent/root
    //const unsigned int su = *t.children(pu,u).begin();
    t.add_edge(parent,pv);
    //cout<<"|Add edge- "<<parent<<" "<<pv;
    // define the new root
    //t.root = su;
    // regraft subtree into edge
    t.add_edge(v,pv);
    //cout<<"|Add edge- "<<pv<<" "<<v<<"|";
    
    return true;
}






// store gene tree traversal - for better memory access
class g_tree_stride {
    public: class u_node {
        public: unsigned int node;
        public: union {
            struct {
                unsigned int v2;
                std::size_t skip;
            } preorder;
            struct {
                unsigned int ch[2];
                // unsigned int parent;
            } postorder;
        };
        public: inline aw::traversal_states direction() { // encode PREORDER/POSTORDER into children
            return (preorder.v2 == node) ? aw::PREORDER : aw::POSTORDER;
        };
    };
    private: std::vector<u_node> data;
    public: std::size_t size() { return data.size(); }
    private: std::size_t node_num;
    public: std::size_t node_size() { return node_num; }
    public: inline u_node &operator[](std::size_t i) { return data[i]; }
    public: template<class TREE> inline void create(TREE &t) {
        data.clear();
        node_num = t.node_size();
        data.reserve(node_num * 2);
        aw::SubtreeSizes t_info; t_info.create(t);
        TREE_DFS2(v,t) {
            switch (v.direction) {
                case aw::PREORDER: {
                    data.resize(data.size()+1);
                    u_node &n = data.back();
                    n.node = n.preorder.v2 = v.idx;
                    n.preorder.skip = 2 * t_info.subtree_size(v.idx,v.parent);
                } break;
                case aw::POSTORDER: {
                    data.resize(data.size()+1);
                    u_node &n = data.back();
                    n.node = v.idx;
                    // n.postorder.parent = v.parent;
                    for (unsigned int i=t.children(v.idx,v.parent,n.postorder.ch); i<2; ++i) {
                        n.postorder.ch[i] = aw::NONODE;
                    }
                } break;
                default: break;
            }
        }
    }
};

// compute the gene duplcation changes for an SPR-down operation
// input:
//   s_tree - species tree where SPR-subtree is at root
//   s_parent - parents of species nodes
//   s_lca - LCAs of species nodes
//   subtree_left - prune subtree_left
//   subtree_right - regraft into subtree
//   parent - parent of subtree_left and subtree_right
//   g_stride - preordered gene tree nodes
//   g_lmap - gene tree LCA mapping
// output:
//   dups_inc - duplication increasment @ node; initialized to 0
//   dups_dec - duplication deccreasment @ node; initialized to 0
template<class STREE>
void compute_dupchanges4SPR(STREE &s_tree, SubtreeParent<aw::Tree> &s_parent, LCA &s_lca,
                            const unsigned int parent, const unsigned int subtree_left, const unsigned int subtree_right,
                            g_tree_stride &g_stride, LCAmapping &g_lmap,
                            util::vector<unsigned int> &dups_inc, util::vector<unsigned int> &dups_dec
) {
    // Bansal, Eulenstein, Wehe algorithm to determine best SPR
    static util::vector<unsigned int> g_lmap2;
    g_lmap2.set_min_size(g_stride.node_size());
    for (unsigned int j=0,jEE=g_stride.size(); j<jEE; ++j) {
        aw::g_tree_stride::u_node &n = g_stride[j];
        switch (n.direction()) {
            case aw::PREORDER: {
                const unsigned int &v_map = g_lmap.mapping(n.node);
                if (v_map == NONODE) break;
                if (v_map != parent) { // no need to traverse further into the subtree
                    j += n.preorder.skip - 1;
                }
            } break;
            case aw::POSTORDER: {
                const unsigned int &v_map = g_lmap.mapping(n.node);
                if (v_map == NONODE) break;
                if (v_map == parent) {
                    const unsigned int * const ch = n.postorder.ch;
                    const unsigned int &ch_map_0 = g_lmap.mapping(ch[0]);
                    const unsigned int &ch_map_1 = g_lmap.mapping(ch[1]);
                    bool case0;
                    // missing lca mapping, usually caused by missing taxa in the species tree
                    // no gene duplications will ever be lost or gained
                    if ((case0 = ch_map_0 == aw::NONODE) || (ch_map_1 == aw::NONODE)) {
                        g_lmap2[n.node] = g_lmap2[ch[case0 ? 1 : 0]];
                        break;
                    }
                    const bool to_root_0 = (ch_map_0 == parent);
                    const bool to_root_1 = (ch_map_1 == parent);
                    // parent and both children map to the root
                    // no gene duplications will every be lost or gained
                    if (to_root_0 && to_root_1) {
                        const unsigned int &ch_map2_0 = g_lmap2[ch[0]];
                        const unsigned int &ch_map2_1 = g_lmap2[ch[1]];
                        const unsigned int v_map2 = s_lca.lca(ch_map2_0,ch_map2_1);
                        g_lmap2[n.node] = v_map2;
                        break;
                    }
                    const bool in_r_0 = s_lca.lca(ch_map_0,subtree_right) == subtree_right;
                    const bool in_r_1 = s_lca.lca(ch_map_1,subtree_right) == subtree_right;
                    // parent and one child map to the root, one child maps into the right subtree
                    if ((case0 = to_root_0 && in_r_1) || (to_root_1 && in_r_0)) {
                        const unsigned int &ch_map2_0 = case0 ? g_lmap2[ch[0]] : ch_map_0;
                        const unsigned int &ch_map2_1 = case0 ? ch_map_1 : g_lmap2[ch[1]];
                        const unsigned int &s_l = case0 ? ch_map2_0 : ch_map2_1; // virtually maps to the right subtree
                        const unsigned int &s_r = case0 ? ch_map2_1 : ch_map2_0; // maps to the right subtree
                        const unsigned int &s_p = g_lmap2[n.node] = s_lca.lca(ch_map2_0,ch_map2_1); // virtually maps to the lca of both child mappings
                        const unsigned int &s_pp = s_parent.parent(s_p); // parent of s_p
                        unsigned int s_pch[2]; s_tree.children(s_p, s_pp, s_pch); // children of s_p
                        if (s_r == s_p) { // duplication remains because of the right child mapping
                        } else
                        if (s_lca.lca(s_l,s_pch[0]) == s_pch[0]) {
                            const unsigned int &s_ll = s_pch[0]; // left child of s_p
                            ++dups_dec[s_ll];
                        } else
                        if (s_lca.lca(s_l,s_pch[1]) == s_pch[1]) {
                            const unsigned int &s_ll = s_pch[1]; // left child of s_p
                            ++dups_dec[s_ll];
                        }
                        break;
                    }
                    const bool in_l_0 = s_lca.lca(ch_map_0,subtree_left) == subtree_left;
                    const bool in_l_1 = s_lca.lca(ch_map_1,subtree_left) == subtree_left;
                    // parent maps to the root, one child maps into the left subtree, one child maps into the right subtree
                    // then gene duplication occurs when moving below the mapped node in the right subtree
                    if ((case0 = in_l_0 && in_r_1) || (in_l_1 && in_r_0)) {
                        const unsigned int &v_map2 = case0 ? ch_map_1 : ch_map_0;
                        g_lmap2[n.node] = v_map2;
                        ++dups_inc[v_map2];
                        break;
                    }
                    // parent and one child map to the root, one child maps into the left subtree
                    // then no gene duplications will every be lost
                    if ((case0 = to_root_0 && in_l_1) || (to_root_1 && in_l_0)) {
                        const unsigned int &v_map2 = g_lmap2[ch[case0 ? 0 : 1]];
                        g_lmap2[n.node] = v_map2;
                        break;
                    }
                    // something is wrong
                    ERROR_exit("something is wrong");
                }
            } break;
            default: break;
        }
    }
}

// accumulate gene duplication changes and compute the gene duplication score for the SPR moves
// input:
//   s_tree - species tree where SPR-subtree is at root
//   subtree,parent - (right)-subtree with possible regrafting edges
//   dups - current gene duplication score when SPR-subtree is at root
//   dups_inc - duplication increasment @ node
//   dups_dec - duplication decreasment @ node
// output:
//   new_location - SPR-edge for lowest gene duplications
//   dups - gene duplications when SPR is moved to new_location
//   ambiguity - number of locations with same lowest gene duplications
//   dups_inc dups_dec - will be reset to 0 for used entries
template<class STREE>
inline void accumulate_dup_changes(STREE &s_tree, const unsigned int subtree, const unsigned int parent,
                                   util::vector<unsigned int> &dups_inc, util::vector<unsigned int> &dups_dec,
                                   std::pair<unsigned int, unsigned int> &new_location, unsigned int &dups, unsigned int &ambiguity
) {
    unsigned int d = dups;
    typedef std::pair<unsigned int,unsigned int> candidates_item;
    static util::vector<candidates_item> candidates;
    candidates.set_min_size(s_tree.node_size());
    for (typename STREE::iterator_dfs v=s_tree.begin_dfs(subtree,parent),vEE=s_tree.end_dfs(); v!=vEE; ++v) {
        switch (v.direction) {
            case aw::PREORDER: {
                d -= dups_dec[v.idx];
                if (d < dups) {
                    dups = d;
                    candidates.clear();
                }
                if (d == dups) {
                    candidates.push_back(std::pair<unsigned int,unsigned int>(v.idx,v.parent));
                }
                d += dups_inc[v.idx];
            } break;
            case aw::POSTORDER: {
                d -= dups_inc[v.idx];
                d += dups_dec[v.idx];
                dups_inc[v.idx] = 0;
                dups_dec[v.idx] = 0;
            } break;
            default: break;
        }
    }
    ambiguity = candidates.size();
// P(candidates.size()-1);
    boost::uniform_int<> range(0,candidates.size()-1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(aw::rng, range);
    new_location = candidates[die()];
    candidates.clear();
}

} // namespace end

// compute the gene duplcation changes for an SPR-down operation
// input:
//   s_tree - species tree where SPR-subtree is at root
//   s_parent - parents of species nodes
//   s_lca - LCAs of species nodes
//   subtree_left - prune subtree_left
//   subtree_right - regraft into subtree
//   parent - parent of subtree_left and subtree_right
//   g_stride - preordered gene tree nodes
//   g_lmap - gene tree LCA mapping
// output:
//   dups_inc - duplication increasment @ node; initialized to 0
//   dups_dec - duplication deccreasment @ node; initialized to 0
template<class STREE>
void compute_dupchanges4SPR_singletaxon(STREE &s_tree, aw::SubtreeParent<aw::Tree> &s_parent, aw::LCA &s_lca,
                            const unsigned int parent, const unsigned int subtree_left, const unsigned int subtree_right,
                            aw::g_tree_stride &g_stride, aw::LCAmapping &g_lmap,
                            util::vector<unsigned int> &dups_inc, util::vector<unsigned int> &dups_dec
) {
    // Bansal, Eulenstein, Wehe algorithm to determine best SPR
    static util::vector<unsigned int> g_lmap2; g_lmap2.set_min_size(g_stride.node_size());
    for (unsigned int j=0,jEE=g_stride.size(); j<jEE; ++j) {
        aw::g_tree_stride::u_node &n = g_stride[j];
        switch (n.direction()) {
            case aw::PREORDER: {
                const unsigned int &v_map = g_lmap.mapping(n.node);
                if (v_map == aw::NONODE) break;
                if (v_map != parent) { // no need to traverse further into the subtree
                    j += n.preorder.skip - 1;
                }
            } break;
            case aw::POSTORDER: {
                const unsigned int &v_map = g_lmap.mapping(n.node);
                if (v_map == aw::NONODE) break;
                if (v_map == parent) {
                    const unsigned int * const ch = n.postorder.ch;
                    const unsigned int &ch_map_0 = g_lmap.mapping(ch[0]);
                    const unsigned int &ch_map_1 = g_lmap.mapping(ch[1]);
                    bool case0;
                    // missing lca mapping, usually caused by missing taxa in the species tree
                    // no gene duplications will ever be lost or gained
                    if ((case0 = ch_map_0 == aw::NONODE) || (ch_map_1 == aw::NONODE)) {
                        g_lmap2[n.node] = g_lmap2[ch[case0 ? 1 : 0]];
                        break;
                    }
                    const bool to_root_0 = (ch_map_0 == parent);
                    const bool to_root_1 = (ch_map_1 == parent);
                    // parent and both children map to the root
                    // no gene duplications will every be lost or gained
                    if (to_root_0 && to_root_1) {
                        const unsigned int &ch_map2_0 = g_lmap2[ch[0]];
                        const unsigned int &ch_map2_1 = g_lmap2[ch[1]];
                        const unsigned int v_map2 = s_lca.lca(ch_map2_0,ch_map2_1);
                        g_lmap2[n.node] = v_map2;
                        break;
                    }
                    const bool in_r_0 = s_lca.lca(ch_map_0,subtree_right) == subtree_right;
                    const bool in_r_1 = s_lca.lca(ch_map_1,subtree_right) == subtree_right;
                    // parent and one child map to the root, one child maps into the right subtree
                    if ((case0 = to_root_0 && in_r_1) || (to_root_1 && in_r_0)) {
                        const unsigned int &ch_map2_0 = case0 ? g_lmap2[ch[0]] : ch_map_0;
                        const unsigned int &ch_map2_1 = case0 ? ch_map_1 : g_lmap2[ch[1]];
                        const unsigned int &s_l = case0 ? ch_map2_0 : ch_map2_1; // virtually maps to the right subtree
                        const unsigned int &s_r = case0 ? ch_map2_1 : ch_map2_0; // maps to the right subtree
                        const unsigned int &s_p = g_lmap2[n.node] = s_lca.lca(ch_map2_0,ch_map2_1); // virtually maps to the lca of both child mappings
                        const unsigned int &s_pp = s_parent.parent(s_p); // parent of s_p
                        unsigned int s_pch[2]; s_tree.children(s_p, s_pp, s_pch); // children of s_p
                        if (s_r == s_p) { // duplication remains because of the right child mapping
                        } else
                        if (s_lca.lca(s_l,s_pch[0]) == s_pch[0]) {
                            const unsigned int &s_ll = s_pch[0]; // left child of s_p
                            ++dups_dec[s_ll];
                        } else
                        if (s_lca.lca(s_l,s_pch[1]) == s_pch[1]) {
                            const unsigned int &s_ll = s_pch[1]; // left child of s_p
                            ++dups_dec[s_ll];
                        }
                        break;
                    }
                    const bool in_l_0 = s_lca.lca(ch_map_0,subtree_left) == subtree_left;
                    const bool in_l_1 = s_lca.lca(ch_map_1,subtree_left) == subtree_left;
                    // parent maps to the root, one child maps into the left subtree, one child maps into the right subtree
                    // then gene duplication occurs when moving below the mapped node in the right subtree
                    if ((case0 = in_l_0 && in_r_1) || (in_l_1 && in_r_0)) {
                        const unsigned int &v_map2 = case0 ? ch_map_1 : ch_map_0;
                        g_lmap2[n.node] = v_map2;
                        ++dups_inc[v_map2];
                        break;
                    }
                    // parent and one child map to the root, one child maps into the left subtree
                    // then no gene duplications will every be lost
                    if ((case0 = to_root_0 && in_l_1) || (to_root_1 && in_l_0)) {
                        const unsigned int &v_map2 = g_lmap2[ch[case0 ? 0 : 1]];
                        g_lmap2[n.node] = v_map2;
                        break;
                    }
                    // something is wrong
                    ERROR_exit("something is wrong");
                }
            } break;
            default: break;
        }
    }
}

#endif
