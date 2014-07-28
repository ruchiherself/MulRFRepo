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

#ifndef TREE_SUBTREE_PROPERTIES_H
#define TREE_SUBTREE_PROPERTIES_H

#include "common.h"
#include "tree_traversal.h"
#include "tree_LCA.h"
#include <limits.h>
#include <vector>

namespace aw {

using namespace std;

class NodeDistance {
    protected: vector<unsigned int> levels;
    protected: LCA *lca_map;
    protected: bool external_lca;
    public: NodeDistance() {
        init();
    }
    public: NodeDistance(const NodeDistance &r) { // copy constructor
        init();
        copy(r);
    }
    public: NodeDistance& operator=(const NodeDistance& r) { // assign operator
        if (this != &r) {
            this->free();
            copy(r);
        }
        return *this;
    }
    protected: inline void init() {
        lca_map = NULL;
        external_lca = false;
    }
    protected: inline void copy(const NodeDistance &r) {
        levels = r.levels;
        external_lca = r.external_lca;
        if (r.lca_map != NULL) {
            lca_map = new LCA;
            *lca_map = *(r.lca_map);
        }
    }
    public: ~NodeDistance() {
        this->free();
    }
    protected: inline void free() {
        if ((lca_map != NULL) && (external_lca == false)) delete lca_map;
        init();
    }
    public: void clear() {
        free();
    }
    public: template<class TREE> bool create(TREE &t) {
        this->free();
        lca_map = new LCA;
        lca_map->create(t);
        levels.clear();
        levels.resize(t.node_size());
        TREE_DFS2(v,t) {
            levels[v.idx] = v.lvl;
        }
        return true;
    }
    public: template<class TREE> bool create(TREE &t, LCA &lca) {
        this->free();
        external_lca = true;
        lca_map = &lca;
        levels.clear();
        levels.resize(t.node_size());
        TREE_DFS2(v,t) {
            levels[v.idx] = v.lvl;
        }
        return true;
    }
    public: inline unsigned int distance(const unsigned int u, const unsigned int v) {
        const unsigned int lca = lca_map->lca(u,v);
        const unsigned int &lvl_u = levels[u];
        const unsigned int &lvl_v = levels[v];
        const unsigned int &lvl_lca = levels[lca];
        return (lvl_u - lvl_lca) + (lvl_v - lvl_lca);
    }
};

} // end namespace

#endif
