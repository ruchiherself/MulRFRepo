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

#ifndef TREE_NAME_MAP_H
#define TREE_NAME_MAP_H

#include "common.h"
#include "tree_IO.h"
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <boost/foreach.hpp>

#ifndef NOHASH
#include <boost/unordered_map.hpp>
#endif

namespace aw {

using namespace std;

#ifdef NOHASH
typedef std::map<unsigned int, unsigned int> id2id;
typedef std::map<unsigned int, std::vector<unsigned int> > id2ids;
#else
typedef boost::unordered_map<unsigned int, unsigned int> id2id;
typedef boost::unordered_map<unsigned int, std::vector<unsigned int> > id2ids;
#endif

// map identical node names
class TaxaMap {
    private: std::vector<std::string> gid2name; // [global id]:taxon name as string
    #ifdef NOHASH
    private: std::map<std::string, unsigned int> name2gid; // [taxon name]:global id
    #else
    private: boost::unordered_map<std::string, unsigned int> name2gid; // [taxon name]:global id
    #endif
    
    // number taxa
    public: inline unsigned int size() {
        return gid2name.size();
    }
    // assign taxa global IDs
    public: inline void insert(aw::idx2name &names4tree) {
        if (gid2name.size() < gid2name.size()) gid2name.reserve(gid2name.size());
        BOOST_FOREACH(aw::idx2name::value_type &w,names4tree) {
            const std::string &taxon = w.second;
            insert(taxon);            
        }
    }
    // assign taxa global IDs
    public: template<class T> inline void insert(T &names4tree) {
        if (gid2name.size() < gid2name.size()) gid2name.reserve(gid2name.size());
        BOOST_FOREACH(const std::string &taxon,names4tree) {
            insert(taxon);
        }
    }
    // assign a taxon a global ID
    public: inline void insert(const std::string &taxon) {
        if (name2gid.find(taxon) == name2gid.end()) {
            gid2name.push_back(taxon);
            name2gid[taxon] = gid2name.size()-1;
            //std::cout<<"<"<<taxon<<" "<<name2gid[taxon]<<">";
        }        
    }
    
    // return taxon name
    public: inline const std::string &taxon(const unsigned int gid) {
        return gid2name[gid];
    }
    // return global ID
    public: inline const unsigned int &gid(const std::string name) {
        return name2gid[name];
    }
    // assign a taxon a global ID
    public: inline bool exist(const std::string &taxon) {
        return (name2gid.find(taxon) != name2gid.end());
    }
};

// map taxa of a tree to a TaxaMap
class TreetaxaMap {
    #ifdef NOHASH
    private: typedef std::map<unsigned int,unsigned int> id2gid_type;
    #else
    private: typedef boost::unordered_map<unsigned int,unsigned int> id2gid_type;
    #endif
    private: id2gid_type id2gid; // [node id]:global id
    #ifdef NOHASH
    private: typedef std::multimap<unsigned int, unsigned int> gid2id_type;
    #else
    private: typedef boost::unordered_multimap<unsigned int, unsigned int> gid2id_type;
    #endif
    private: gid2id_type gid2id; // [global id]:list of node id
    private: unsigned int uniq_leaf;   //stores the number of unique leaves int the tree  :RUCHI


    // create the mapping
    public: inline void create(aw::idx2name &names4tree, TaxaMap &m) {
        // mapping from tree nodes to global IDs
        init_unq_leaves();
        BOOST_FOREACH(aw::idx2name::value_type &w,names4tree) {
            const aw::idx2name::key_type &id = w.first;
            std::string &taxon = w.second;
            const unsigned int &gid = m.gid(taxon);
            if(!exists(gid)) ++uniq_leaf;
            id2gid[id] = gid;
            gid2id.insert(pair<unsigned int,unsigned int>(gid,id));
        }
    }    

    //return the number of unique leaves :RUCHI
    public: inline unsigned int unq_leaves() {
        return uniq_leaf;
    }

    //initialize the number of unique leaves :RUCHI
    public: inline void init_unq_leaves() {
        uniq_leaf = 0;
    }

    //return the number of ids for a gid :RUCHI
    public: inline int ids_count(const unsigned int gid) {
        int cnt = 0;
        std::pair<gid2id_type::iterator,gid2id_type::iterator> ret = gid2id.equal_range(gid);
        //cout<<"\n"<<"mapping=>"<<gid<<" ";
        for (gid2id_type::iterator it=ret.first; it!=ret.second; ++it) {
            cnt++;
        }
        return cnt;
    }

    //check if the gid exists
    public: inline bool exists(const unsigned int gid){
        if(gid2id.find(gid) == gid2id.end())
            return false;
        else
            return true;
    }


    // remove a single mapping
    public: inline void remove(const unsigned int id) {
        const unsigned int gid = id2gid[id];
        gid2id.erase(gid);
        id2gid.erase(id);
    }
    // insert a single mapping
//    public: inline void add(const unsigned int id, const unsigned int gid) {
//        id2gid[id] = gid;
//        gid2id.insert(pair<unsigned int,unsigned int>(gid,id));
//
//    }
    // insert a single mapping   0 = unique, 1 = notunique
    public: inline void add(const unsigned int id, const unsigned int gid, const unsigned int unq) {
        if(unq == 0) ++uniq_leaf;
        id2gid[id] = gid;
        gid2id.insert(pair<unsigned int,unsigned int>(gid,id));

    }
    // insert a single mapping
    public: inline void add(const unsigned int id, std::string &taxon, TaxaMap &m) {
        const unsigned int &gid = m.gid(taxon);
        id2gid[id] = gid;
        gid2id.insert(pair<unsigned int,unsigned int>(gid,id));
    }
    // return the ids with global id
    public: inline void ids(const unsigned int gid, std::vector<unsigned int> &_ids) {
        std::pair<gid2id_type::iterator,gid2id_type::iterator> ret = gid2id.equal_range(gid);
        //cout<<"\n"<<"mapping=>"<<gid<<" ";
        for (gid2id_type::iterator it=ret.first; it!=ret.second; ++it) {
            const unsigned int &id = it->second;
            _ids.push_back(id);
            //cout<<id<<" ";
        }
    }

    // return one ids with global id: Ruchi
    public: inline unsigned int one_id(const unsigned int gid) {
        std::pair<gid2id_type::iterator,gid2id_type::iterator> ret = gid2id.equal_range(gid);
        //cout<<"\n"<<"mapping=>"<<gid<<" ";
        for (gid2id_type::iterator it=ret.first; it!=ret.second; ++it) 
            return it->second;        
    }

    // return the global id corresponding to id
    public: inline unsigned int gid(const unsigned int id) {
        return id2gid[id];
    }
    // return the ids matching a list of global ids
    public: template<class T> inline void gids2ids(T &gids, std::set<unsigned int> &ids_cont) {
        BOOST_FOREACH(const unsigned int &gid,gids) {
            std::vector<unsigned int> m;
            ids(gid,m);
            BOOST_FOREACH(const unsigned int &id,m) {
                ids_cont.insert(id);
            }
        }
    }

    // return the ids matching a list of global ids
    public: template<class T> inline void gids2ids(T &gids, std::vector<unsigned int> &ids_cont) {
        BOOST_FOREACH(const unsigned int &gid,gids) {
            std::vector<unsigned int> m;
            ids(gid,m);
            BOOST_FOREACH(const unsigned int &id,m) {
                ids_cont.push_back(id);
            }
        }
    }
    // return nodes in the dest tree where node_id maps to
    public: inline void mapping(const unsigned int &node_id, TreetaxaMap &dest, std::vector<unsigned int> &mapped) {
        const unsigned int &gid = id2gid[node_id];
        dest.ids(gid,mapped);
    }
    // return the global IDs that exist in 2 trees
    public: inline void intersection(TreetaxaMap &dest, std::set<unsigned int> &m) {
        std::set<unsigned int> gid_src;
        BOOST_FOREACH(const id2gid_type::value_type &w,id2gid) {
            const unsigned int &gid = w.second;
            gid_src.insert(gid);
        }
        std::set<unsigned int> gid_dest;
        BOOST_FOREACH(const id2gid_type::value_type &w,dest.id2gid) {
            const unsigned int &gid = w.second;
            gid_dest.insert(gid);
        }
        std::set_intersection(gid_src.begin(), gid_src.end(), gid_dest.begin(), gid_dest.end(), std::inserter(m,m.begin()));
    }
    // return the global IDs that do not exist in dest_tree
    public: inline void difference(TreetaxaMap &dest, std::set<unsigned int> &m) {
        std::set<unsigned int> gid_src;
        BOOST_FOREACH(const id2gid_type::value_type &w,id2gid) {
            const unsigned int &gid = w.second;
            gid_src.insert(gid);
        }
        std::set<unsigned int> gid_dest;
        BOOST_FOREACH(const id2gid_type::value_type &w,dest.id2gid) {
            const unsigned int &gid = w.second;
            gid_dest.insert(gid);
        }
        std::set_difference(gid_src.begin(), gid_src.end(), gid_dest.begin(), gid_dest.end(), std::inserter(m,m.begin()));
    }
    // return taxa
    public: inline void taxa(TaxaMap &m, aw::idx2name &names4tree) {
        BOOST_FOREACH(const id2gid_type::value_type &w,id2gid) {
            const unsigned int &id = w.first;
            const unsigned int &gid = w.second;
            const std::string &name = m.taxon(gid);
            names4tree[id] = name;
        }
    }
    // return all global IDs
    public: inline void gids(std::vector<unsigned int> &t_gids) {
        t_gids.reserve(t_gids.size() + id2gid.size());
        BOOST_FOREACH(const id2gid_type::value_type &w,id2gid) {
            const unsigned int &gid = w.second;
            t_gids.push_back(gid);
        }
    }
    // return all global IDs /SET :RUCHI
    public: inline void unq_gids(std::set<unsigned int> &gids) {
        BOOST_FOREACH(const id2gid_type::value_type &w,id2gid) {
            const unsigned int &gid = w.second;
            gids.insert(gid);
        }
    }
    // return the ids matching a global ID
    public: inline void gid2ids(const unsigned int gid, std::vector<unsigned int> &ids_cont) {
        BOOST_FOREACH(gid2id_type::value_type &w,gid2id.equal_range(gid)) {
            ids_cont.push_back(w.second);
        }
    }
};

} // namespace end

#endif
