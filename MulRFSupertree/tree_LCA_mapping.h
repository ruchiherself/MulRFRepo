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

#ifndef TREE_LCAMAPPING_H
#define TREE_LCAMAPPING_H

#include "common.h"
#include "tree_traversal.h"
#include "tree_LCA.h"
#include "tree_name_map.h"

namespace aw {

using namespace std;

// compute the LCA mapping between a gene tree and a species tree
class LCAmapping {
private: std::vector<unsigned int> _map;
    public: LCAmapping() { }
    public: ~LCAmapping() { }
    // create the LCA mapping between 2 trees
    public: template<class TREE>
    inline void create(LCA &s_lca, TreetaxaMap &s_map, TreetaxaMap &g_map, TREE &g_tree) {
        this->free();
        // reserve space for data
        update_LCA_leaves(s_map,g_map,g_tree);
        // create internal mapping
        if (!g_tree.is_rooted()) ERROR_exit("rooted tree expected"); // LCA mapping for rooted gene trees only
        update_LCA_internals(s_lca,g_tree);
    }
    // update the LCA mapping of all leaf nodes between 2 trees
    public: template<class TREE> inline void update_LCA_leaves(TreetaxaMap &s_map, TreetaxaMap &g_map, TREE &g_tree) {
        // create leaf mapping
        _map.resize(g_tree.node_size());
        TREE_FOREACHLEAF(v,g_tree) {
            unsigned int gid = g_map.gid(v);            
//            if(s_map.exists(gid)) _map[v] = s_map.id(gid);
//            else _map[v] = NONODE;
        }        
    }

    // update the LCA mapping of all leaf nodes between 2 trees
    public: template<class TREE> inline void update_LCA_leaves(TreetaxaMap &g_map, TreetaxaMap &s_map, TREE &g_tree, TREE &s_tree) {
        // create leaf mapping
        _map.resize(s_tree.node_size());
        unsigned int ggid_count = s_map.unq_leaves();

        std::set<unsigned int> ggids;
        s_map.unq_gids(ggids);

        if(ggid_count!=ggids.size())  ERROR_exit("Error!");

        for(std::set<unsigned int>::iterator it = ggids.begin(); it!=ggids.end(); it++) {
            unsigned int ggid = *it;
            std::vector<unsigned int> s_ids;
            s_map.ids(ggid,s_ids);
            std::vector<unsigned int> g_ids;
            g_map.ids(ggid,g_ids);
            unsigned int j=0;
            for(unsigned int i=0; i<s_ids.size(); ++i){
                if(j<g_ids.size())
                    _map[s_ids[i]] = g_ids[j++];
                else
                    _map[s_ids[i]] = NONODE;
            }
        }
    
//        for(unsigned int ggid=0; ggid<ggid_count; ++ggid){
//            std::vector<unsigned int> s_ids;
//            s_map.ids(ggid,s_ids);
//            std::vector<unsigned int> g_ids;
//            g_map.ids(ggid,g_ids);
//            unsigned int j=0;
//            for(unsigned int i=0; i<s_ids.size(); ++i){
//
//                if(j<g_ids.size())
//                    _map[s_ids[i]] = g_ids[j++];
//                else
//                    _map[s_ids[i]] = NONODE;
//                //std::cout<<"-"<<s_ids[i]<<" "<<_map[s_ids[i]];
//            }
//        }
    }
    
    // update the LCA leaf mapping for one gene tree node
    public: inline void update_LCA_leaf(TreetaxaMap &s_map, TreetaxaMap &g_map, const unsigned int gene_id) {
        if (gene_id >= _map.size()) _map.resize(gene_id+1);
        std::vector<unsigned int> m;
        const unsigned int &v = gene_id;
        g_map.mapping(v,s_map,m);
        switch (m.size()) {
            case 0: {
                _map[v] = NONODE;
            } break;
            case 1: {
                _map[v] = m[0];
            } break;
            default: {
                ERROR_exit("invalid leaf mapping");
            } break;
        }
    }
    // update the LCA mapping of a single internal node
    public: template<class TREE> inline void update_LCA_internal(LCA &s_lca, TREE &g_tree, const unsigned int gene_id, const unsigned int parent) {
        unsigned int v_map = NONODE;
        //std::cout<<"--gd "<<gene_id<<"--";
        BOOST_FOREACH(const unsigned int &c,g_tree.children(gene_id,parent)) {
            v_map = s_lca.lca(v_map,_map[c]);            
        }
        //MSG_nonewline("<<"<<gene_id<<" "<<v_map<<">>");
        _map[gene_id] = v_map;
    }
    // update the LCA mapping of a single internal node
    public: inline void update_LCA_internal_binary(LCA &s_lca, const unsigned int gene_id, const unsigned int ch0, const unsigned int ch1) {
        _map[gene_id] = s_lca.lca(_map[ch0],_map[ch1]);
    }
    // update the LCA mapping of all internal nodes between 2 trees
    public: template<class TREE> inline void update_LCA_internals(LCA &s_lca, TREE &g_tree) {
        if (!g_tree.is_rooted()) ERROR_exit("rooted tree expected"); // LCA mapping for rooted gene trees only
        TREE_POSTORDER2(v,g_tree) {
            if (!g_tree.is_leaf(v.idx)) {
                update_LCA_internal(s_lca, g_tree, v.idx,v.parent);
                //std::cout<<" *lca "<<v.idx<<" -> "<<_map[v.idx];
            }
        }
    }

    // set the LCA mapping for one gene tree node - manually
    public: inline void set_LCA(const unsigned int gene_id, const unsigned int species_id) {
        if (gene_id >= _map.size()) _map.resize(gene_id+1);
        _map[gene_id] = species_id;
        //cout<<"\n new mapping of"<<gene_id<<" is " <<species_id;
    }

    // lca mapping
    public: inline unsigned int &mapping(const unsigned int gene_id) {
        if(gene_id == NONODE) ERROR_exit("Error in tree_LCA_mapping");  //added by ruchi
        return _map[gene_id];
    }
    public: inline unsigned int &operator[](const unsigned int gene_id) {
        return _map[gene_id];
    }

    public: inline void clear() {
        free();
    }
    protected: inline void free() {
        _map.clear();
    }
};

} // namespace end

#endif
