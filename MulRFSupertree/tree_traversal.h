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

#ifndef TREE_TRAVERSAL_H
#define TREE_TRAVERSAL_H

#include "common.h"
#include "tree.h"

#define TREE_DFS2(VAL, tree) \
    for (aw::Tree::iterator_dfs VAL=(tree).begin_dfs(),itrEE=(tree).end_dfs(); VAL!=itrEE; ++(VAL))

#define TREE_EULERTOUR2(VAL, tree) \
    for (aw::Tree::iterator_eulertour VAL=(tree).begin_eulertour(),itrEE=(tree).end_eulertour(); VAL!=itrEE; ++(VAL))

#define TREE_PREORDER2(VAL, tree) \
    for (aw::Tree::iterator_preorder VAL=(tree).begin_preorder(),itrEE=(tree).end_preorder(); VAL!=itrEE; ++(VAL))

#define TREE_INORDER2(VAL, tree) \
    for (aw::Tree::iterator_inorder VAL=(tree).begin_inorder(),itrEE=(tree).end_inorder(); VAL!=itrEE; ++(VAL))

#define TREE_POSTORDER2(VAL, tree) \
    for (aw::Tree::iterator_postorder VAL=(tree).begin_postorder(),itrEE=(tree).end_postorder(); VAL!=itrEE; ++(VAL))

// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------

#define TREE_FOREACHNODE(VAL, tree) \
    for (unsigned int VAL=UINT_MAX; VAL==UINT_MAX; VAL=0) \
    for (aw::Tree::iterator itr=(tree).begin(),itrEE=(tree).end(); aw::_tt_compare23(itr,itrEE,VAL); ++itr)

#define TREE_FOREACHLEAF(VAL, tree) \
    for (unsigned int VAL=UINT_MAX; VAL==UINT_MAX; VAL=0) \
    for (aw::Tree::iterator_leafnodes itr=(tree).begin_leafnodes(),itrEE=(tree).end_leafnodes(); aw::_tt_compare23(itr,itrEE,VAL); ++itr)

#define TREE_FOREACHINTERNAL(VAL, tree) \
    for (unsigned int VAL=UINT_MAX; VAL==UINT_MAX; VAL=0) \
    for (aw::Tree::iterator_internalnodes itr=(tree).begin_internalnodes(),itrEE=(tree).end_internalnodes(); aw::_tt_compare23(itr,itrEE,VAL); ++itr)

// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------

#define TREE_DFS(VAL, tree) \
    for (unsigned int VAL=UINT_MAX; VAL==UINT_MAX; VAL=0) \
    for (aw::Tree::iterator_dfs itr=(tree).begin_dfs(),itrEE=(tree).end_dfs(); aw::_tt_compare23(itr,itrEE,VAL); ++itr)

#define TREE_EULERTOUR(VAL, tree) \
    for (unsigned int VAL=UINT_MAX; VAL==UINT_MAX; VAL=0) \
    for (aw::Tree::iterator_eulertour itr=(tree).begin_eulertour(),itrEE=(tree).end_eulertour(); aw::_tt_compare23(itr,itrEE,VAL); ++itr)

#define TREE_PREORDER(VAL, tree) \
    for (unsigned int VAL=UINT_MAX; VAL==UINT_MAX; VAL=0) \
    for (aw::Tree::iterator_preorder itr=(tree).begin_preorder(),itrEE=(tree).end_preorder(); aw::_tt_compare23(itr,itrEE,VAL); ++itr)

#define TREE_INORDER(VAL, tree) \
    for (unsigned int VAL=UINT_MAX; VAL==UINT_MAX; VAL=0) \
    for (aw::Tree::iterator_inorder itr=(tree).begin_inorder(),itrEE=(tree).end_inorder(); aw::_tt_compare23(itr,itrEE,VAL); ++itr)

#define TREE_POSTORDER(VAL, tree) \
    for (unsigned int VAL=UINT_MAX; VAL==UINT_MAX; VAL=0) \
    for (aw::Tree::iterator_postorder itr=(tree).begin_postorder(),itrEE=(tree).end_postorder(); aw::_tt_compare23(itr,itrEE,VAL); ++itr)

// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------

namespace aw {

// function for the simple traversal macros
template<class T> inline bool _tt_compare23(T &itr, T &end, unsigned int &v) {
    if (itr == end) return false;
    v = *itr;
    return true;
}

} // end of namespace

#endif
