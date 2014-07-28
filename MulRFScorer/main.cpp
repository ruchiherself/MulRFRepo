/* 
 * File:   main.cpp
 * Author: ruchi
 *
 * Created on April 2, 2014, 3:05 PM
 */

#include <stdlib.h>
#include "argument.h"
#include "common.h"
#include "tree.h"
#include "tree_IO.h"
#include "gauge.h"
#include "tree_traversal.h"
#include "tree_LCA.h"
#include "tree_LCA_mapping.h"
#include "tree_name_map.h"
#include "tree_duplication.h"
#include <boost/foreach.hpp>
#include <boost/progress.hpp>
#include "boost/tuple/tuple.hpp"

#ifndef NOHASH
#include <boost/unordered_set.hpp>
#endif

static const unsigned int NONODE = UINT_MAX;
typedef boost::unordered_map<unsigned int,int> gid2ctype;

/*
 * 
 */
int main(int ac, char* av[]) {
        {
        std::ostringstream os; os << "command:";
        for (int i = 0; i < ac; i++) os << ' ' << av[i];
        MSG(os.str());
    }
    std::string trees_filename;    
    std::string output_filename;
    bool stree_first = true;
    
    {
        Argument a; a.add(ac, av);
        // help
        if (a.existArg2("-h","--help")) {
            MSG("options:");
            MSG("  -i [ --input ] arg      input trees (file in NEWICK format)");
            MSG("  -o [ --output ] arg     write the trees into a file");
            MSG("  -h [ --help ]           produce help message");
            MSG("");
            MSG("example:");
            MSG("  " << av[0] << " -i inputF.newick -o outputF.newick");
            exit(0);
        }

        // input trees
        if (a.existArgVal2("-i", "--input", trees_filename)) MSG("input file: " << trees_filename) else MSG("using standard input");
        
        // output file
        if (a.existArgVal2("-o", "--output", output_filename)) MSG("output file: " << output_filename);
        
        
                
        // unknown arguments?
        a.unusedArgsError();
    }
    // -----------------------------------------------------------------------------------

        //create output stream
    std::ofstream ouput_fs;
    if (!output_filename.empty()) {
        ouput_fs.open(output_filename.c_str());
        if (!ouput_fs) ERROR_exit("cannot write file '" << output_filename << "'");
    }
    std::ostream &output = output_filename.empty() ? std::cout : ouput_fs;

    // read trees
    // output:
    aw::Tree s_tree;
    aw::idx2name s_taxa;
    std::vector<aw::Tree> g_trees;
    std::vector<aw::Tree> rs_trees;
    std::vector<aw::idx2name> g_taxa;
    std::vector<aw::LCAmapping> s_lmaps;
    aw::TreetaxaMap s_nmap;
    std::vector<aw::LCA> g_lca;
    std::vector<float> g_weights;
    {
        const std::string filename = trees_filename;
        { // read trees
            std::ifstream ifs;
            std::istream &is = filename.empty() ? std::cin : ifs;
            if (!filename.empty()) {
                ifs.open(filename.c_str());
                if (!ifs) ERROR_exit("cannot read file '" << filename << "'");
            }
            if (stree_first) {
                float t_w = 1.0f;
                if (!aw::stream2tree(is, s_tree, s_taxa,t_w)) ERROR_exit("No species tree found in file '" << filename << "'");
            }
            MSG_nonewline("Reading input trees: ");
            aw::gauge_exp g; aw::gauge_init(&g);
            for (;;) {
                aw::Tree t;
                aw::idx2name t_names;
                float t_w = 1.0f;
                if (!aw::stream2tree(is, t, t_names,t_w)) break;
                g_taxa.push_back(t_names);
                g_trees.push_back(t);                
                g_weights.push_back(t_w);
                aw::gauge_inc(&g);
            }
            aw::gauge_end(&g);
        }
        MSG("Input trees: " << g_trees.size());
        if (g_trees.empty()) ERROR_exit("No input trees found in file '" << filename << "'");
    }

    if (stree_first) { // check input for binary
            TREE_FOREACHNODE(v,s_tree) {
                const unsigned int d = s_tree.degree(v);
                bool isroot = false;
                if(s_tree.is_rooted()) isroot = (s_tree.root == v);
                else isroot = (0 == v);
                if (d == 0) continue; // single isolated node
                if ((d == 1) && (!isroot)) continue; // leaf
                if ((d == 3) && (!isroot)) continue; // binary interior node but not root
                if ((d == 2) && isroot) continue; // root node with children
                ERROR_exit("Initial species tree"<<" is not binary");
            }
        MSG("Initial species tree pass binary test");
    }

    // map taxa labels
    aw::TaxaMap taxamap;  //to store all taxon and global id
    std::vector<aw::TreetaxaMap> g_nmaps; //for mapping taxamap and idx2name (of a tree)
    {
        g_nmaps.resize(g_taxa.size());
        for (unsigned int i=0,iEE=g_taxa.size(); i<iEE; ++i) {
            aw::idx2name &n = g_taxa[i];
            taxamap.insert(n);
            g_nmaps[i].create(n,taxamap);
        }
        MSG("Taxa: " << taxamap.size());
    }

    std::vector<std::pair<unsigned int,unsigned int> > g_nodes;  //pair <internal node,leaf count>
    std::vector<unsigned int> root_leaf;
    { // gene tree nodes
        unsigned int c = 0, t = 0, t1 = 0, c1 = 0;
        for (unsigned int k=0,kEE=g_trees.size(); k<kEE; ++k){
            bool first = true;
            TREE_PREORDER2(v,g_trees[k]) {
                ++c;
                if (g_trees[k].is_leaf(v.idx)) {
                    ++t;
                    if(first) { root_leaf.push_back(g_nmaps[k].gid(v.idx)); first = false; }
                }
            }
            if(g_trees[k].degree(0)>2) ++c;
            t1 +=t; c1 +=c;
            g_nodes.push_back(std::pair<unsigned int,unsigned int>(c-t-1,t));
            t = 0; c = 0;
        }
        MSG("Input tree nodes: " << c1 << " (" << t1 << " taxa)");
    }

    if (stree_first) {
        aw::idx2name &n = s_taxa;
        if(n.size()<taxamap.size())  ERROR_exit("Error: Initial species tree doesn't have all leaves!!");
        taxamap.insert(n);
        s_nmap.create(n,taxamap);
    }

    std::vector<unsigned int> leaves_trim;  //if need to trim some leaves from the species tree

    //counting maximum number of copies of a gene node in a genetree
    gid2ctype gid2cnt;    //<gid,count>
    {
        for(unsigned int i=0,iKK=taxamap.size(); i<iKK; ++i) gid2cnt[i]=0;
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            TREE_FOREACHLEAF(v,g_trees[i]) {
                unsigned int gid_v = g_nmaps[i].gid(v);
                int gid_cnt = g_nmaps[i].ids_count(gid_v);
                if(gid2cnt[gid_v] < gid_cnt) gid2cnt[gid_v] = gid_cnt;
        }   }

        for(unsigned int i=0,iKK=taxamap.size(); i<iKK; ++i)
          if(gid2cnt[i]==0) {
              WARNING("Supertree leaf "<<taxamap.taxon(i)<<" not in any of the input trees!");
              leaves_trim.push_back(i);
              gid2cnt.erase(i);
          }
    }

    //actually trimming leaves
    BOOST_FOREACH(const unsigned int &c, leaves_trim) {
        unsigned int loc = s_nmap.one_id(c);
        if(s_tree.degree(loc)!=1) {
            MSG("Singleton node!");
            continue;
        }
        if(!s_tree.trim_leaf(loc))
            ERROR_exit("Error trimming leaf...");
    }

    //Extending the supertree for RF computation...
    if (stree_first) {
        BOOST_FOREACH(const gid2ctype::value_type &w, gid2cnt) {
            unsigned int ggid = w.first;
            int l_cnt = w.second;
            if(l_cnt==0) ERROR_exit("Error");
            if(l_cnt==1) continue;
            unsigned int sid;
            sid = s_nmap.one_id(ggid); //since initial supertree is not multilabeled
            std::vector<unsigned int> new_leaves;
            s_tree.extend_leaf(sid,l_cnt,new_leaves);
            std::string taxon = taxamap.taxon(ggid);
            BOOST_FOREACH(const unsigned int &c,new_leaves)
                s_nmap.add(c,ggid,1);
        }
    }

    {   // root the trees by one leaf
        for (unsigned int k=0; k<g_trees.size(); ++k){
            unsigned int k_id = g_nmaps[k].one_id(root_leaf[k]);
            g_trees[k].rootBy(k_id);
        }
    }

    {   //Calculate cluster size for input trees
        unsigned int count;
        for (unsigned int k=0; k<g_trees.size(); ++k)
            TREE_POSTORDER2(v,g_trees[k]) {
                if (g_trees[k].is_leaf(v.idx))
                    g_trees[k].update_clst(v.idx,1);
                else {  count = 0;
                    BOOST_FOREACH(const unsigned int &c,g_trees[k].children(v.idx,v.parent))
                        count = count + g_trees[k].return_clstSz(c);
                    g_trees[k].update_clst(v.idx,count); }
            }
    }

    {   //Root Supertree + define initial LCA leaf Mapping
        s_lmaps.clear();
        s_lmaps.resize(g_trees.size());
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            rs_trees.push_back(s_tree);
            s_lmaps[i].update_LCA_leaves(g_nmaps[i],s_nmap,g_trees[i],rs_trees[i]);
            std::vector<unsigned int> ch;
            unsigned int g_rt;
            g_trees[i].adjacent(0,ch);
            if(g_trees[i].is_leaf(ch[0])) g_rt = ch[0]; else g_rt = ch[1];
            ch.clear();
            unsigned int ggid = g_nmaps[i].gid(g_rt);
            s_nmap.ids(ggid,ch);
            unsigned int s_id;
            BOOST_FOREACH(const unsigned int &c,ch)
                if(s_lmaps[i].mapping(c)==g_rt) {
                    s_id = c;   break;  }
            rs_trees.back().rootBy(s_id);
        }
    }

    std::vector<unsigned int> rs_int_nodes;  //internal nodes in each rs_tree
    rs_int_nodes.clear();
    {   //store internal nodes count for each rs_tree
        unsigned int temp;
        for (unsigned int k=0; k<rs_trees.size(); ++k) {
            temp = g_nmaps[k].unq_leaves() - 2;
            BOOST_FOREACH(const gid2ctype::value_type &w, gid2cnt) {
                unsigned int ggid = w.first;
                int g_count = g_nmaps[k].ids_count(ggid);
                if(g_count > 1) temp++;
            }
            rs_int_nodes.push_back(temp);
        }
    }

    {   //Computing cluster size for supertrees: computed based on leaf mapping
        unsigned int count;
        for (unsigned int k=0, kEE=rs_trees.size(); k<kEE; ++k){
            TREE_POSTORDER2(v,rs_trees[k]) {
                if (!rs_trees[k].is_leaf(v.idx)) {
                    count = 0;
                    BOOST_FOREACH(const unsigned int &c,rs_trees[k].children(v.idx,v.parent))
                        count = count + rs_trees[k].return_clstSz(c);
                    rs_trees[k].update_clst(v.idx,count);
                }
                else {
                    unsigned int gid = s_nmap.gid(v.idx);
                    if(s_lmaps[k].mapping(v.idx)!=NONODE)    //:FOR MUL-TREES
                        rs_trees[k].update_clst(v.idx,1);
                    else rs_trees[k].update_clst(v.idx,0);
                }
            }
        }
    }

    {   g_lca.clear();  //store lca if it is done first time
        aw::LCA lca;
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            lca.create(g_trees[i]);
            g_lca.push_back(lca); }
    }

    std::vector<float> g_scr;
    float scr = 0.0f;
    {   for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i)
            s_lmaps[i].update_LCA_internals(g_lca[i],rs_trees[i]);
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            std::pair<unsigned int,unsigned int> p = g_nodes[i];
            g_scr.push_back(aw::compute_rf_score(rs_trees[i],g_trees[i],s_lmaps[i],p,rs_int_nodes[i],g_weights[i]));
            scr = scr + g_scr[i] ;
        }
        MSG_nonewline("\nMulRF Score: "<<std::fixed<<std::setprecision(2)<<scr);
    }

    {   //preprocessing of s_tree
        BOOST_FOREACH(const gid2ctype::value_type &w, gid2cnt) {
            unsigned int ggid = w.first;
            int l_cnt = w.second;
            if(l_cnt==1) continue;
            unsigned int sid,adj;
            sid = s_nmap.one_id(ggid); //since initial supertree is not multilabeled
            std::vector<unsigned int> ch;
            s_tree.adjacent(sid,ch);
            unsigned int fakeInt = ch[0];
            BOOST_FOREACH(const unsigned int &c, s_tree.adjacent(fakeInt))
                if(!s_tree.is_leaf(c)) { adj = c;  break; }

            s_tree.disconnect_node(fakeInt);
            s_tree.add_edge(sid,adj);
        }

        output <<"[ Species Tree: Unrooted RF Score = "<<std::fixed<<std::setprecision(2)<<scr<<"]"<< std::endl;
        aw::tree2newick(output,s_tree,s_taxa); output << std::endl;
    }

    
        
    for(int mn=0, mnEE=g_trees.size(); mn<mnEE; ++mn) {
         output <<"\n[ Gene Tree "<<mn<< " MulRF Score = "<<std::fixed<<std::setprecision(2)<<g_scr[mn]<<"]"<< std::endl;
         output<<"[&WEIGHT="<<std::fixed<<std::setprecision(2)<<g_weights[mn]<<"]";
         aw::tree2newick(output,g_trees[mn],g_taxa[mn]); output << std::endl; 
    }

    return (EXIT_SUCCESS);
}

