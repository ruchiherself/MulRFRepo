/*
 * A project to create RF super tree of Multi-labaled trees
 * File:   main.cpp
 * Author: ruchi
 *
 * Created on March 1, 2012, 3:21 PM
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
#include "tree_node_distance.h"
#include "tree_duplication.h"
#include "rf_compute.h"
#include <boost/foreach.hpp>
#include <boost/progress.hpp>
#include "boost/tuple/tuple.hpp"


#ifndef NOHASH
#include <boost/unordered_set.hpp>
#endif

static const unsigned int NONODE = UINT_MAX;
long double EPSILON = 0.00001;

struct chEdge { unsigned int x, y, px; };

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
    std::string c_filename;
    bool stree_first = false;
    std::string output_filename;    
    bool alltrees = false;
    bool initialtree = false;
    bool inputrees = false;
    bool constr = false;    
    unsigned int seed = std::time(0);
    unsigned int SPR_rounds = 0; 
    {
        Argument a; a.add(ac, av);
        // help
        if (a.existArg2("-h","--help")) {
            MSG("options:");
            MSG("  -i [ --input ] arg      input trees (file in NEWICK format)");
            MSG("       --stree            first input tree is a starting species tree");
            MSG("  -c [--constraints] arg   constraints");
            MSG("  -o [ --output ] arg     write the trees into a file (file in NEWICK format)");
            MSG("       --alltrees         synonym for --initialtree and --inputrees");
            MSG("       --initialtree      output the initial species tree");
            MSG("       --inputrees        output the input trees");            
            MSG("       --seed arg         random generator seed");            
            MSG("  -h [ --help ]           produce help message");
            MSG("");
            MSG("example:");
            MSG("  " << av[0] << " -i inputF.newick -o outputF.newick");
            exit(0);
        }
        
        // input trees
        if (a.existArgVal2("-i", "--input", trees_filename)) MSG("input file: " << trees_filename) else MSG("using standard input");
        // starting species tree
        stree_first = a.existArg("--stree");
        // constraints 
        if (a.existArgVal2("-c","--constraints", c_filename)) {
            MSG("constraints file: " << c_filename);
            constr = true;
        }        
        if(constr && stree_first)  WARNING("Constraints doesn't work with starting species tree!");
      
        // output file
        if (a.existArgVal2("-o", "--output", output_filename)) MSG("output file: " << output_filename);
        // output initial species tree
        initialtree = a.existArg("--initialtree");
        // output gene trees
        inputrees = a.existArg("--inputrees");        
        // output all trees
        alltrees = a.existArg("--alltrees");
        if (alltrees) {
            inputrees = true;
            initialtree = true;
        }
        // random seed
        a.existArgVal("--seed", seed);
        MSG("seed: " << seed);
        aw::rng.seed(static_cast<unsigned int>(seed));
        // unknown arguments?
        a.unusedArgsError();
    }
    // -----------------------------------------------------------------------------------

    int t3 = clock();

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
    std::vector<std::vector<std::string> > c_taxa;
    std::vector<aw::Tree> rs_trees;
    std::vector<aw::idx2name> g_taxa;    
    std::vector<aw::LCAmapping> s_lmaps;
    aw::TreetaxaMap s_nmap;
    std::vector<aw::LCA> g_lca;
    std::vector<float> g_weights;
    {
        // read trees -------------------------------------
        {
            const std::string filename = trees_filename;
            std::ifstream ifs;
            std::istream &is = filename.empty() ? std::cin : ifs;
            if (!filename.empty()) {
                ifs.open(filename.c_str());
                if (!ifs) ERROR_exit("cannot read file '" << filename << "'");
            }
            if (stree_first) {
                float t_w = 1.0f;  //we will not use this weight
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
            MSG("Input trees: " << g_trees.size());
            if (g_trees.empty()) ERROR_exit("No input trees found in file '" << filename << "'");
         }

        // reading constriants file ----------------------------
         if (constr) {   
            std::ifstream ifs;
            std::istream &is = c_filename.empty() ? std::cin : ifs;
            if (!c_filename.empty()) {
                ifs.open(c_filename.c_str());
                if (!ifs) ERROR_exit("cannot read file '" << c_filename << "'");
            }                        
            for (;;) {
                std::vector<std::string> constr_list;
                if (!aw::stream2constr(is, constr_list)) break;
                if(constr_list.size()<=1)
                    ERROR_exit("Too small constraint!");
                c_taxa.push_back(constr_list);
            }

            MSG("Constraints: " << c_taxa.size());
            if (c_taxa.empty()) ERROR_exit("No constraints found in file '" << c_filename << "'");
        }       
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
  
    // checking constraints
    boost::unordered_map<unsigned int, unsigned int> gid2c; //<global id, order of its list>
    {       
        for (unsigned int i=0,iEE=c_taxa.size(); i<iEE; ++i) {
            std::vector<std::string> ls = c_taxa[i];
            BOOST_FOREACH(const std::string &st,ls) {
                if(!taxamap.exist(st))
                    ERROR_exit("Constraints error, unique taxa: "<<st);
                unsigned int gid = taxamap.gid(st);

                if(gid2c.find(gid) == gid2c.end())
                    gid2c.insert(std::pair<unsigned int,unsigned int>(gid,i));
                else
                    ERROR_exit("Constraints error, overlap!");                
            }                         
        }        
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
          if(gid2cnt[i]==0) ERROR_exit("Error: Some leaf of the supertree is not in any of the input trees!!");
    }

    int cst[c_taxa.size()];  //FOR CONSTRAINTS.... to store each clades's root node in the species tree
    for(unsigned int i=0; i<c_taxa.size(); ++i)
        cst[i] = NONODE;
      
    // build starting tree using leaf adding ************************************************************************************************
    if (!stree_first) {        
        int t1 = clock();
        
        MSG("Building initial species tree...");
        std::vector<unsigned int> s_inodes,g_inodes;  //internal node in s_tree, g_tree
        std::queue<unsigned int> taxa_queue;

        if(true) { // random taxa order
            std::vector<unsigned int> nodes; nodes.reserve(taxamap.size());
            for (unsigned int i=0,iEE=taxamap.size(); i<iEE; ++i) nodes.push_back(i);
            boost::uniform_int<> range(0,nodes.size()-1);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(aw::rng, range);
            for (unsigned int i=0,iEE=nodes.size(); i<iEE; ++i) {
                const unsigned int j=die();
                util::swap(nodes[i],nodes[j]);
            }
            BOOST_FOREACH(const unsigned int &i,nodes) {
                taxa_queue.push(i);
            }
        }

         if(false){   // taxa order based on taxa connectivity  NOT IN USE RIGHT NOW
            std::vector<unsigned int> linkage(taxamap.size(),0);  //store sum of num of leaves in each tree it appears
            {
                for (unsigned int k=0,kEE=g_trees.size(); k<kEE; ++k) {
                    aw::Tree &g_tree = g_trees[k];
                    aw::TreetaxaMap &g_nmap = g_nmaps[k];
                    std::pair<unsigned int,unsigned int> p = g_nodes[k];
                    const unsigned int taxa_size = p.second;
                    TREE_PREORDER2(v,g_tree) if (g_tree.is_leaf(v.idx)) {
                        const unsigned v_gid = g_nmap.gid(v.idx);
                        linkage[v_gid] += taxa_size;  }
                }
            }

            std::vector<unsigned int> &priority = linkage;  //to get highly used taxa first
            typedef std::multimap<unsigned int, unsigned int> a_type;  a_type a;
            for (unsigned int i=0,iEE=priority.size(); i<iEE; ++i) a.insert(a_type::value_type(priority[i],i));
            //taxa_queue has taxons in high linkage to low linkage order
            for (a_type::reverse_iterator i=a.rbegin(),iEE=a.rend(); i!=iEE; ++i) {taxa_queue.push(i->second);
            }
        }

        // create all species nodes
        unsigned int nodeCount = 0;
        nodeCount = 2*gid2cnt.size() - 1;
        BOOST_FOREACH(const gid2ctype::value_type &w, gid2cnt) {
            if(w.second > 1) nodeCount +=w.second; 
        }
        for (unsigned int i=0,iEE=(nodeCount); i<iEE; ++i) s_tree.new_node(); // leaf + internal nodes

        //start building the real tree now...............................
        unsigned int node_c = 0;
        { // create an initial tree by 2 unique leaves            
            unsigned int id0, id1;
            id0 = taxa_queue.front(); taxa_queue.pop();
            id1 = taxa_queue.front(); taxa_queue.pop();
            int leav0 = gid2cnt[id0];
            int leav1 = gid2cnt[id1];            
            
            unsigned int p, c0, c1;
            p = node_c++;
            s_tree.root = p;
            s_nmap.init_unq_leaves();
            if(leav0>1) {
                c0 = node_c++;  s_tree.set_fake(c0);
                bool unqFirst = true;
                for(int i=0; i<leav0; ++i) {
                    unsigned int l = node_c++;
                    s_tree.add_edge(c0,l);
                    if(unqFirst) { s_nmap.add(l,id0,0); unqFirst = false; }
                    else s_nmap.add(l,id0,1);
                    s_taxa.insert(aw::idx2name::value_type(l, taxamap.taxon(id0)));  }
            } else {
                c0 = node_c++;
                s_nmap.add(c0,id0,0);
                s_taxa.insert(aw::idx2name::value_type(c0, taxamap.taxon(id0))); }

            if(constr && gid2c.find(id0) != gid2c.end()) {  //FOR CONSTRAINTS ..................
                unsigned int list = gid2c[id0];
                s_tree.set_constr(c0,list);
                cst[list] = c0;                
            }           

            if(leav1>1) {
                c1 = node_c++;  s_tree.set_fake(c1);
                bool unqFirst = true;
                for(int i=0; i<leav1; ++i) {
                    unsigned int l = node_c++;
                    s_tree.add_edge(c1,l);
                    if(unqFirst) { s_nmap.add(l,id1,0); unqFirst = false; }
                    else s_nmap.add(l,id1,1);
                    s_taxa.insert(aw::idx2name::value_type(l, taxamap.taxon(id1)));  }
            } else {
                c1 = node_c++;
                s_nmap.add(c1,id1,0);
                s_taxa.insert(aw::idx2name::value_type(c1, taxamap.taxon(id1))); }
            if(constr && gid2c.find(id1) != gid2c.end()) {
                unsigned int list = gid2c[id1];
                s_tree.set_constr(c1,list);                
                cst[list] = c1;                
            }            
            s_tree.add_edge(p,c0); s_tree.add_edge(p,c1);

            if(constr && s_tree.constr_num(c0)!= NONODE && s_tree.constr_num(c0)==s_tree.constr_num(c1)) {  //FOR CONSTRAINTS ..................
                unsigned int list = gid2c[id1];
                s_tree.set_constr(p,s_tree.constr_num(c0));
                s_tree.set_constr(c0,NONODE);
                s_tree.set_constr(c1,NONODE);
                cst[list] = p;                
            }
        }
                
        {   //Calculate cluster size for input trees
            unsigned int count;
            for (unsigned int k=0, kEEE=g_trees.size(); k<kEEE; ++k) {
                unsigned int inodes = 0;
                TREE_POSTORDER2(v,g_trees[k]) {
                    if (g_trees[k].is_leaf(v.idx)){
                        unsigned int ggid = g_nmaps[k].gid(v.idx);
                        if(s_nmap.exists(ggid)) g_trees[k].update_clst(v.idx,1);
                        else g_trees[k].update_clst(v.idx,0);

                    }
                    else {  count = 0; unsigned int ncount = 0;
                        BOOST_FOREACH(const unsigned int &c,g_trees[k].children(v.idx,v.parent)) {
                            count = count + g_trees[k].return_clstSz(c);
                            if(g_trees[k].return_clstSz(c)>0)  ++ncount;
                        }
                        g_trees[k].update_clst(v.idx,count);
                        if(ncount>1) ++inodes;
                    }
                }
                if(inodes>0)  --inodes; //subtracting for root
                g_inodes.push_back(inodes);
            }
        }       
        
        {   //MSG("Updating LCA leaf mapping");
            s_lmaps.resize(g_trees.size());
            for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i)
                s_lmaps[i].update_LCA_leaves(g_nmaps[i],s_nmap,g_trees[i],s_tree);
        }
        
        std::vector< aw::TreeClusters<aw::Tree> > s_clst; s_clst.resize(g_trees.size());
        //updating s_tree cluster + s_inodes
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            s_clst[i].create(s_tree,g_trees[i],g_nmaps[i],s_nmap,s_lmaps[i]);
            unsigned int inodes = 0;
            TREE_POSTORDER2(v,s_tree) {
                unsigned int count = 0;
                BOOST_FOREACH(const unsigned int &c,s_tree.children(v.idx,v.parent))
                    if(s_clst[i].cluster(c)>0) ++count;
                if(count>1) ++inodes;
            }
            if(inodes>0)  --inodes;            
            s_inodes.push_back(inodes);
        }

        std::vector<unsigned int> g_scr;
        float scr = 0;
        {   aw::LCA lca;
            for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {                
                lca.create(g_trees[i]);
                g_lca.push_back(lca);    // compute LCAs for the input tree
                s_lmaps[i].update_LCA_internals(g_lca[i],s_tree); }
        }        
        
        // precompute parents
        aw::SubtreeParent<aw::Tree> s_parent; s_parent.create(s_tree);
        std::vector<aw::SubtreeParent<aw::Tree> > g_parents(g_trees.size());
        for (unsigned int k=0,kEE=g_trees.size(); k<kEE; ++k) g_parents[k].create(g_trees[k]);        

        //Adding remaning leaves-----------------------------------------------------------------------------------------
        while(taxa_queue.size()!=0) { 
            const unsigned int gid = taxa_queue.front(); taxa_queue.pop();            
            const unsigned int p = node_c++;
            unsigned int gid_list = c_taxa.size();  //this gid's list
            unsigned int c;            
            bool in_clade = false;
          
            if(gid2cnt[gid]>1) {
                c = node_c++;   s_tree.set_fake(c);
                bool unqFirst = true;
                for(int i=0; i<gid2cnt[gid]; ++i) {
                    unsigned int l = node_c++;
                    s_tree.add_edge(c,l);
                    if(unqFirst) { s_nmap.add(l,gid,0); unqFirst=false; }
                    else s_nmap.add(l,gid,1);
                    s_taxa.insert(aw::idx2name::value_type(l, taxamap.taxon(gid)));  }
            } else {
                c = node_c++;
                s_nmap.add(c,gid,0);
                s_taxa.insert(aw::idx2name::value_type(c, taxamap.taxon(gid))); }            

            s_tree.add_edge(c,p);  unsigned int old_root = s_tree.root;            
            
            if(constr && gid2c.find(gid) != gid2c.end() && cst[gid2c[gid]] != NONODE) {  //if new gid should go in an existing cade               
                gid_list = gid2c[gid];
                unsigned int clade_rt = cst[gid_list];
                unsigned int p_clade_rt = s_parent.parent(clade_rt);                
                s_tree.remove_edge(clade_rt,p_clade_rt); s_tree.add_edge(clade_rt,p); s_tree.add_edge(p_clade_rt,p);                
                s_parent.update(clade_rt,p); s_parent.update(p,p_clade_rt); s_parent.update(c,p);   //update parents
                in_clade = true;
            }
            else { //new gid should start from the root                
                s_tree.add_edge(s_tree.root,p);                
                s_tree.root = p;
                s_parent.update(old_root,p); s_parent.update(p,NONODE); s_parent.update(c,p);   //update parents 
            }            

            //updating parents for newly added leaves
            std::vector<unsigned int> sids;
            s_nmap.ids(gid,sids);            
            if(sids.size()>1)
                for(unsigned int s=0; s<sids.size(); ++s)
                    s_parent.update(sids[s],c);
            s_parent.tPtrUpdate(s_tree);            

            //update LCAs
            if(!in_clade){               
                for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {                   
                    if(g_nmaps[i].exists(gid)){
                        //updating s_inodes number...
                        std::vector<unsigned int> gids;
                        g_nmaps[i].ids(gid,gids);
                        if(gids.size()>1) {  //updating number of internal nodes in the new s copies
                            unsigned int old_s = s_inodes[i]; s_inodes.erase(s_inodes.begin()+i); s_inodes.insert(s_inodes.begin()+i,old_s+2);
                        }
                        else {
                            unsigned int old_s = s_inodes[i]; s_inodes.erase(s_inodes.begin()+i); s_inodes.insert(s_inodes.begin()+i,old_s+1);
                        }                        
                        
                        //mapping leaf nodes...
                        unsigned int j=0;
                        for(unsigned int k=0; k<sids.size(); ++k){
                            if(j<gids.size())
                                s_lmaps[i].set_LCA(sids[k],gids[j++]);
                            else
                                s_lmaps[i].set_LCA(sids[k],NONODE);
                        }
                        //mapping internal nodes...                        
                        if(sids.size()>1)
                            s_lmaps[i].update_LCA_internal(g_lca[i],s_tree,c,p);
                        s_lmaps[i].set_LCA(s_tree.root, g_lca[i].lca(s_lmaps[i].mapping(old_root),s_lmaps[i].mapping(c)));
                    }
                    else {
                        if(sids.size()>1)
                            for(unsigned int s=0; s<sids.size(); ++s)
                                s_lmaps[i].set_LCA(sids[s],NONODE);
                        s_lmaps[i].set_LCA(c,NONODE);
                        s_lmaps[i].set_LCA(s_tree.root, s_lmaps[i].mapping(old_root));
                    }
                }
            }
            else { //when the leaf was NOT added above root
                for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
                    if(!g_nmaps[i].exists(gid)){
                        if(sids.size()>1)
                            for(unsigned int s=0; s<sids.size(); ++s)
                                s_lmaps[i].set_LCA(sids[s],NONODE);
                        s_lmaps[i].set_LCA(c,NONODE);                       
                        unsigned int clade_rt = cst[gid_list];
                        s_lmaps[i].set_LCA(p, s_lmaps[i].mapping(clade_rt));
                    }
                    else {  //we will have to do all mapping etc from scratch now
                        //updating s_inodes number...
                        std::vector<unsigned int> gids;
                        g_nmaps[i].ids(gid,gids);
                        if(gids.size()>1) {  //updating number of internal nodes in the new s copies
                            unsigned int old_s = s_inodes[i];
                            s_inodes.erase(s_inodes.begin()+i);
                            s_inodes.insert(s_inodes.begin()+i,old_s+2);
                        }
                        else {
                            unsigned int old_s = s_inodes[i];
                            s_inodes.erase(s_inodes.begin()+i);
                            s_inodes.insert(s_inodes.begin()+i,old_s+1);
                        }                        
                        //mapping leaf nodes...
                        unsigned int j=0;
                        for(unsigned int k=0; k<sids.size(); ++k){ // map some leaves and leave others...
                            if(j<gids.size())
                                s_lmaps[i].set_LCA(sids[k],gids[j++]);
                            else
                                s_lmaps[i].set_LCA(sids[k],NONODE);
                        }
                        //mapping internal nodes...
                        s_lmaps[i].update_LCA_internals(g_lca[i],s_tree);
                    }
                }
            }
            

            //Updating clusters for s_tree & input trees + g_inodes
            for (unsigned int k=0, kEE=g_trees.size(); k<kEE; ++k) {
                //redoing clusters for s_tree
                s_clst[k].create(s_tree,g_trees[k],g_nmaps[k],s_nmap,s_lmaps[k]);

                //clusters for gene trees + g_inodes
                if(g_nmaps[k].exists(gid)) {
                    unsigned int inodes = 0;
                    TREE_POSTORDER2(v,g_trees[k]) {
                        if(g_trees[k].is_leaf(v.idx)) {
                            if(g_nmaps[k].gid(v.idx)==gid)
                                g_trees[k].update_clst(v.idx,1); }
                        else {
                            int count = 0; unsigned int ncount = 0;
                            BOOST_FOREACH(const unsigned int &l,g_trees[k].children(v.idx,v.parent)) {
                                count = count + g_trees[k].return_clstSz(l);
                                if(g_trees[k].return_clstSz(l)>0) ++ncount;          }
                            g_trees[k].update_clst(v.idx,count);                            
                            if(ncount>1) ++inodes;
                        }
                    }
                    g_inodes.erase(g_inodes.begin()+k);
                    if(inodes>0) inodes--;                    
                    g_inodes.insert(g_inodes.begin()+k,inodes);
                }
            }            

            {   g_scr.clear();  scr = 0;
                for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {                    
                    g_scr.push_back(aw::compute_rf_score(s_tree,g_trees[i],s_lmaps[i],s_clst[i],g_inodes[i],s_inodes[i]));
                    scr = scr + g_scr[i]*g_weights[i] ;
                }
            }
         
            unsigned int subtree = c;
            unsigned int psubtree = s_parent.parent(subtree);
            aw::Tree best_tree = s_tree, rnd_tree = s_tree;
            float best_score = scr;
            unsigned int itr_start = s_parent.sibling_binary(subtree);
            unsigned int itr_par = psubtree;
            unsigned int last_node = NONODE;

            //MOVE DOWN LOOP..............................
            for (aw::Tree::iterator_dfs m=rnd_tree.begin_dfs(itr_start,itr_par),mEE=rnd_tree.end_dfs(); m!=mEE; ++m) {                
                if(m.idx == itr_start) {
                    if(constr && !in_clade && rnd_tree.constr_num(m.idx)!=NONODE)  break;                   
                    continue;  }
                                
                if(rnd_tree.is_fake(m.parent)) continue;  //FOR MULTree

                if(constr && !in_clade && last_node==NONODE && s_tree.constr_num(m.parent) != NONODE)  {
                    last_node = m.parent;
                    continue;
                }

                if(constr && !in_clade && last_node == m.idx) last_node = NONODE;
                if(constr && !in_clade && last_node != NONODE)  continue;                
                float rf_new = 0, rf_old = 0;

                switch (m.direction) {
                    case aw::PREORDER: {
                        if(constr && !in_clade && s_tree.constr_num(m.idx) != NONODE)  last_node = m.idx;
                        
                        for (unsigned int n=0,nEE=g_trees.size(); n<nEE; ++n) {
                            unsigned int old_par_score = aw::compute_rf_score(s_tree, g_trees[n], s_lmaps[n], psubtree, s_parent.parent(psubtree), s_clst[n]);                            
                            rf_old += old_par_score*g_weights[n];
                            unsigned int old_pm_score = aw::compute_rf_score(s_tree, g_trees[n], s_lmaps[n], m.parent, psubtree, s_clst[n]);                            
                            rf_old += old_pm_score*g_weights[n];                        }
                        aw::move2edge_binary(s_tree, subtree, psubtree, m.idx, m.parent);
                        s_parent.update(m.parent,s_parent.parent(psubtree));
                        s_parent.update(psubtree,m.parent);
                        s_parent.update(m.idx,psubtree);
                        s_parent.tPtrUpdate(s_tree);

                        for (unsigned int n=0,nEE=g_trees.size(); n<nEE; ++n) {                            
                             s_lmaps[n].set_LCA(m.parent, s_lmaps[n].mapping(psubtree));                             
                             unsigned int subt_map = s_lmaps[n].mapping(subtree);
                             unsigned int midx_map = s_lmaps[n].mapping(m.idx);                             
                             s_lmaps[n].set_LCA(psubtree, g_lca[n].lca(subt_map, midx_map));                             
                             s_clst[n].update(m.parent,s_clst[n].cluster(psubtree));                             
                             s_clst[n].update(psubtree,s_clst[n].cluster(m.idx)+s_clst[n].cluster(subtree));                             
                             unsigned int new_par_score = aw::compute_rf_score(s_tree, g_trees[n], s_lmaps[n], psubtree, m.parent, s_clst[n]);
                             rf_new += new_par_score*g_weights[n];
                             unsigned int new_pm_score = aw::compute_rf_score(s_tree, g_trees[n], s_lmaps[n], m.parent, s_parent.parent(m.parent), s_clst[n]);                             
                             rf_new += new_pm_score*g_weights[n];
                        }
                        scr = scr - (rf_old - rf_new);
                        if(fabs(best_score-scr) > EPSILON){
                            best_score = scr;
                            best_tree = s_tree;                            
                        }
                    } break;
                    case aw::POSTORDER: {                       
                        
                        for (unsigned int n=0,nEE=g_trees.size(); n<nEE; ++n) {
                            unsigned int old_par_score = aw::compute_rf_score(s_tree, g_trees[n], s_lmaps[n], psubtree, s_parent.parent(psubtree), s_clst[n]);
                            rf_old += old_par_score*g_weights[n];
                            unsigned int old_pm_score = aw::compute_rf_score(s_tree, g_trees[n], s_lmaps[n], m.parent, s_parent.parent(m.parent), s_clst[n]);
                            rf_old += old_pm_score*g_weights[n];
                        }
                        aw::REVmove2edge_binary(s_tree, subtree, psubtree, m.idx, m.parent, s_parent.parent(m.parent));
                        s_parent.update(psubtree,s_parent.parent(m.parent));
                        s_parent.update(m.parent,psubtree);
                        s_parent.update(m.idx,m.parent);
                        s_parent.tPtrUpdate(s_tree);

                        for (unsigned int n=0,nEE=g_trees.size(); n<nEE; ++n) {
                             s_lmaps[n].set_LCA(psubtree, s_lmaps[n].mapping(m.parent));
                             unsigned int m_map = s_lmaps[n].mapping(m.idx);
                             unsigned int subsib_map = s_lmaps[n].mapping(s_parent.sibling_binary(m.idx));
                             s_lmaps[n].set_LCA(m.parent, g_lca[n].lca(m_map, subsib_map));
                             s_clst[n].update(psubtree,s_clst[n].cluster(m.parent));
                             s_clst[n].update(m.parent,s_clst[n].cluster(m.idx)+s_clst[n].cluster(s_parent.sibling_binary(m.idx)));
                             unsigned int new_par_score = aw::compute_rf_score(s_tree, g_trees[n], s_lmaps[n], psubtree,  s_parent.parent(psubtree), s_clst[n]);
                             rf_new += new_par_score*g_weights[n];
                             unsigned int new_pm_score = aw::compute_rf_score(s_tree, g_trees[n], s_lmaps[n], m.parent, psubtree, s_clst[n]);
                             rf_new += new_pm_score*g_weights[n];
                        }
                        scr = scr - (rf_old - rf_new);
                    } break;
                    default: break;
                }
            }
            s_tree = best_tree;
            s_parent.create(s_tree);                        

            if(constr) {
                if(in_clade && (s_parent.parent(cst[gid_list])==p)) { //update the clade info                   
                    s_tree.set_constr(cst[gid_list],NONODE);                    
                    cst[gid_list] = p;
                    s_tree.set_constr(p,gid_list);
                }
                else if(!in_clade && gid2c.find(gid) != gid2c.end()) {
                    unsigned int list = gid2c[gid];
                    cst[list] = c;
                    s_tree.set_constr(c,list);
                }               
            }
       
            s_clst.clear(); s_clst.resize(g_trees.size());
            for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i)
                s_clst[i].create(s_tree,g_trees[i],g_nmaps[i],s_nmap,s_lmaps[i]);

            {   float scr1 = 0; g_scr.clear();
                for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i)                   
                    s_lmaps[i].update_LCA_internals(g_lca[i],s_tree); 
                for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) { 
                    g_scr.push_back(aw::compute_rf_score(s_tree,g_trees[i],s_lmaps[i],s_clst[i],g_inodes[i],s_inodes[i]));
                    scr1 = scr1 + g_scr[i]*g_weights[i] ;
                }
            
                //if(fabs(scr1 - best_score) > EPSILON) ERROR_exit("Scores doesn't match!!");

            }          
        }        

        s_tree.rootInit();
        //aw::tree2newick(output,s_tree,s_taxa); output << std::endl;
        unsigned int num = NONODE;

        //setting in_clade for each node of clades...
        for (aw::Tree::iterator_dfs m=s_tree.begin_dfs(),mEE=s_tree.end_dfs(); m!=mEE; ++m) {
            if(m.idx==0) continue;
            switch (m.direction) {
                case aw::PREORDER: {
                    if(num==NONODE && s_tree.constr_num(m.idx)!=NONODE)
                        num = s_tree.constr_num(m.idx);
                    s_tree.set_in_cld(m.idx,num);
                } break;
                case aw::POSTORDER: {
                    if(num!=NONODE && s_tree.constr_num(m.idx)!=NONODE)
                        num = NONODE;
                } break;
                default: {continue;} break;
            }            
        }
      
        int t2 = clock();
        {
            long ttime = ((long)(t2-t1))/CLOCKS_PER_SEC;
            int d, h, m, s;
            util::convertTime(ttime,d,h,m,s);

            MSG_nonewline("Elapsed time in building initial tree: ");
            if(d!=0) MSG_nonewline(d<<"d ");
            if(h!=0) MSG_nonewline(h<<"h ");
            if(m!=0) MSG_nonewline(m<<"m ");
            MSG(s<<"s ");
        }
    } //Initial tree is ready ***************************************************************************************************************************************

    //Extending the supertree for RF computation...
    if (stree_first) {        
        BOOST_FOREACH(const gid2ctype::value_type &w, gid2cnt) {
            unsigned int ggid = w.first;
            int l_cnt = w.second;
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

    // Write the tree in the output file if asked for
    if (initialtree) {
        output << "[ Initial Species Tree ]" << std::endl;        

        {   //preprocessing of s_tree
            aw::Tree temp_stree = s_tree;
            BOOST_FOREACH(const gid2ctype::value_type &w, gid2cnt) {
                unsigned int ggid = w.first;
                int l_cnt = w.second;
                if(l_cnt==1) continue;
                unsigned int sid,adj;
                sid = s_nmap.one_id(ggid); //since initial supertree is not multilabeled
                std::vector<unsigned int> ch;
                temp_stree.adjacent(sid,ch);
                unsigned int fakeInt = ch[0];
                BOOST_FOREACH(const unsigned int &c, temp_stree.adjacent(fakeInt))
                    if(!temp_stree.is_leaf(c)) { adj = c;  break; }

                temp_stree.disconnect_node(fakeInt);
                temp_stree.add_edge(sid,adj);
            }            
            aw::tree2newick(output,temp_stree,s_taxa); output << std::endl;
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
            unsigned int g_rt;  //g_tree root leaf
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

    std::vector<unsigned int> g_scr;
    float scr = 0;
    {   for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i)
            s_lmaps[i].update_LCA_internals(g_lca[i],rs_trees[i]);
        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
            std::pair<unsigned int,unsigned int> p = g_nodes[i];
            g_scr.push_back(aw::compute_rf_score(rs_trees[i],g_trees[i],s_lmaps[i],p,rs_int_nodes[i]));            
            scr = scr + g_scr[i]*g_weights[i] ;
        }
        MSG_nonewline("\nCurrent RF Score: "<<std::fixed<<std::setprecision(2)<< scr);
    }

    aw::Tree bestTree = s_tree; //to store best tree in one SPR neighborhood
    float bestScore = scr;

    //***********************************************     SPR START     ***********************************************************************
    for(;;)
    {
        SPR_rounds++;
        std::vector<bool> treeEft;
        unsigned int lost_node = NONODE;
        unsigned int x, px, y;
        std::vector<struct chEdge> spr_edge;   //round robin

        //s_tree is not changed in any way -- CLADES are preserved....
        //Starting unrooted SPR....
        for (aw::Tree::iterator_dfs v=s_tree.begin_dfs(),vEE=s_tree.end_dfs(); v!=vEE; ++v) {
            //get a unique edge in unrooted s_tree            
            if(v.idx==0) continue;            
            if(s_tree.is_fake(v.parent) && s_tree.is_leaf(v.idx)) {  continue;} //escepting taking edges of extended stars :FOR MUL-TREE
            switch (v.direction) {
                case aw::PREORDER: {
                    if(v.parent==0 && lost_node==NONODE) {
                        lost_node = v.idx; continue; }
                    else if(v.parent==0 && lost_node!=NONODE) {
                        if(s_tree.is_fake(v.idx) && s_tree.is_leaf(lost_node)) continue; //escepting taking edges of extended stars :FOR MUL-TREE
                        x = lost_node; y = v.idx;  px = v.parent; }
                    else {
                        x = v.idx; y = v.parent; px = y;}
                } break;
                case aw::POSTORDER: {  continue; } break;
                default: {continue;} break;
            }
            chEdge a;            
            a.x = x; a.px = px; a.y = y;
            spr_edge.push_back(a);
        }

        //for randomization
        boost::uniform_int<> range(0,spr_edge.size()-1);
        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > die(aw::rng, range);
        for (unsigned int i=0,iEE=spr_edge.size(); i<iEE; ++i) {
            const unsigned int j=die();
            chEdge temp;
            temp = spr_edge[i];
            spr_edge[i] = spr_edge[j];
            spr_edge[j] = temp;
        }

        int spr_side = die()%2;   //for more randomization
        bool reg_x = true;
        if(spr_side==1)
            reg_x = false;

        for(int sd = 0; sd<2; ++sd) {
            for (unsigned int qi=0,qiEE=spr_edge.size(); qi<qiEE; ++qi) {
                x = spr_edge[qi].x; y = spr_edge[qi].y; px = spr_edge[qi].px;               

                aw::Tree us_tree = s_tree;
                us_tree.delRoot();

                //Storing leaves below x in s_tree with tag 'x' and others with 'y'
                //y is parent of x or both are siblings (for edge having root (0)) in s_tree
                boost::unordered_map<unsigned int, char> slid2char;
                unsigned int reg_leaf = NONODE;
                {
                    TREE_FOREACHLEAF(vl,s_tree)  slid2char[vl]='y';
                    for (aw::Tree::iterator_dfs v1=s_tree.begin_dfs(x,px),vEE1=s_tree.end_dfs(); v1!=vEE1; ++v1) {
                        if(v1.idx == px) break;
                        if(s_tree.is_leaf(v1.idx)) slid2char[v1.idx]='x';   }
                }

                //for generalization
                unsigned int prn_side, rgft_side;
                char my_char, oth_char;
                if(reg_x) {
                    prn_side=x; rgft_side=y; my_char = 'x'; oth_char = 'y';
                } else {
                    prn_side=y; rgft_side=x; my_char = 'y'; oth_char = 'x';}

                if(constr) {  //FOR CONSTRAINT ----------------------------------------------------------------
                    if(px!=y) {
                        if(s_tree.constr_num(x)!= NONODE && s_tree.constr_num(y)!=NONODE) { continue;
                        } else if(s_tree.constr_num(x)!=NONODE) {
                            if(prn_side==y) continue;
                        } else if(s_tree.constr_num(y)!=NONODE) {
                            if(prn_side==x) continue;
                        }
                    }
                    else {
                        if(s_tree.constr_num(x) != NONODE) {  //x side can be pruned but no y side since x makes a clade
                            if(s_tree.constr_num(y) != NONODE)  ERROR_exit("ERROR");
                            if(prn_side==y)  continue;
                        } else if(s_tree.constr_num(y) != NONODE) {  //y side can be pruned but no x side since y makes a clade
                            if(s_tree.constr_num(x) != NONODE)  ERROR_exit("ERROR");
                            if(prn_side==x)  continue;
                        } else if(s_tree.in_cld(x)==s_tree.in_cld(y) && s_tree.constr_num(x)==NONODE && s_tree.constr_num(y)==NONODE) {
                            //LIMITED regraft
                            if(prn_side==y) continue;
                        }
                    }
                }

                TREE_FOREACHLEAF(v2,s_tree) { // find a leaf that is not multiple
                    bool flgg = false;
                    if(slid2char[v2]== oth_char) {
                        std::vector<unsigned int> ch;
                        s_tree.adjacent(v2,ch);
                        if(ch.size()>1) ERROR_exit("Leaf has more than one adjacent nodes!");
                        if(!s_tree.is_fake(ch[0])) {
                            reg_leaf=v2;
                            flgg = true; }
                        else {
                            reg_leaf = ch[0];
                            flgg = true; }
                    }
                    if(flgg) {
                        if(constr) {  //FOR CONSTRAINT ----------------------------------------------------------------
                            if(s_tree.constr_num(prn_side)!=NONODE) {   //can go any where but in an another clade
                                if(s_tree.in_cld(reg_leaf)!=NONODE)
                                    continue;
                            } else {
                                if(s_tree.in_cld(prn_side)!=s_tree.in_cld(reg_leaf))
                                    continue;
                            }
                        }
                        break;
                    }
                }
                
                if(reg_leaf==NONODE) ERROR_exit("Uninitialised reg_leaf");
                int round = us_tree.spr_to_edge(prn_side,rgft_side,reg_leaf);   //Regraaft XX above reg_leaf in YY

                if(round == 0) {
                    treeEft.clear();  rs_trees.clear();
                    std::vector<char> reroot (g_trees.size());
                    for (unsigned int k=0,kEEE=g_trees.size(); k<kEEE; ++k) {
                        bool noX = false, noY = false;
                        unsigned int rootAt, old_root;
                        TREE_FOREACHLEAF(w,g_trees[k]) {
                            unsigned int gid = g_nmaps[k].gid(w);
                            unsigned int sid = s_nmap.one_id(gid);
                            if(!noX && slid2char[sid]==my_char) noX = true;
                            if(!noY && slid2char[sid]==oth_char) { rootAt = w; noY = true;}
                            if(noX && noY) break;   //ADDED 9th SEPT
                        }
                        rs_trees.push_back(us_tree);

                        if(!noX || !noY) {
                            treeEft.push_back(false);  continue;   } //NO Need to do for this round of this tree
                        else treeEft.push_back(true);

                        BOOST_FOREACH(const unsigned int &w, g_trees[k].adjacent(0))
                            if(g_trees[k].is_leaf(w)) old_root = w;  //Assuming input trees have more than 2 leaf3

                        //check if we really need to reroot input tree
                        unsigned int gid = g_nmaps[k].gid(old_root);
                        unsigned int sid = s_nmap.one_id(gid);
                        if(slid2char[sid]==oth_char){
                            rootAt = old_root;
                            reroot[k] = 'N';   }
                        else  reroot[k] = 'Y';

                         //rooting s_tree & g_tree by same leaf
                        if(reroot[k]=='Y') g_trees[k].rootBy(rootAt);
                        unsigned int gRootAt = g_nmaps[k].gid(rootAt);                        
                        std::vector<unsigned int> child;
                        s_nmap.ids(gRootAt,child);
                        std::vector<unsigned int> ch1;
                        rs_trees.back().adjacent(child[0],ch1);
                        if(ch1.size()>1) ERROR_exit("Leaf has more than one adjacent nodes!");
                        BOOST_FOREACH(const unsigned int &c,child){    //:FOR MUL-TREES
                            if(s_lmaps[k].mapping(c)==rootAt) 
                                rs_trees.back().addRoot(c,ch1[0]);                            
                        }                        
                    }

                    {   //Computing cluster size for supertrees: computed based on leaf mapping
                        unsigned int count;
                        for (unsigned int k=0, kEE=rs_trees.size(); k<kEE; ++k){
                            if(!treeEft[k]) continue;
                            TREE_POSTORDER2(v,rs_trees[k]) {                                
                                if (!rs_trees[k].is_leaf(v.idx)) {                                      
                                    count = 0;
                                    BOOST_FOREACH(const unsigned int &c,rs_trees[k].children(v.idx,v.parent))
                                        count = count + rs_trees[k].return_clstSz(c);
                                    rs_trees[k].update_clst(v.idx,count);  }
                                else {
                                    if(s_lmaps[k].mapping(v.idx)!=NONODE)    //:FOR MUL-TREES
                                        rs_trees[k].update_clst(v.idx,1);
                                    else rs_trees[k].update_clst(v.idx,0);  }
                            }
                        }
                    }
                    
                    {   //Calculate cluster size for input trees
                        unsigned int count;
                        for (unsigned int k=0, kEE=g_trees.size(); k<kEE; ++k){
                            if(!treeEft[k] || (reroot[k]=='N')) continue;                            
                            TREE_POSTORDER2(v, g_trees[k])
                                if (g_trees[k].is_leaf(v.idx))
                                    g_trees[k].update_clst(v.idx,1);
                                else {  count = 0;
                                    BOOST_FOREACH(const unsigned int &c,g_trees[k].children(v.idx,v.parent))
                                        count += g_trees[k].return_clstSz(c);
                                    g_trees[k].update_clst(v.idx,count); }
                        }
                    }
                    
                    std::vector<unsigned int> g_score;
                    float score = 0;
                    {
                        aw::LCA lca;
                        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
                            if(!treeEft[i] || (reroot[i]=='N')) continue;
                            lca.create(g_trees[i]);
                            g_lca.erase(g_lca.begin()+i);
                            g_lca.insert(g_lca.begin()+i,lca);
                        }
                        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
                            if(!treeEft[i]) continue;
                            s_lmaps[i].update_LCA_internals(g_lca[i],rs_trees[i]); }
                        for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
                            if(!treeEft[i]){ g_score.push_back(g_scr[i]);
                                score = score + g_scr[i]; continue; }
                            std::pair<unsigned int,unsigned int> p = g_nodes[i];
                            g_score.push_back(aw::compute_rf_score(rs_trees[i],g_trees[i],s_lmaps[i],p,rs_int_nodes[i]));
                            score = score + g_score[i]*g_weights[i] ;
                        }
                    }

                    //for easy parent-child relationship in rs_trees
                    std::vector<aw::SubtreeParent<aw::Tree> > rs_parents(rs_trees.size());
                    for (unsigned int k=0,kEE=g_trees.size(); k<kEE; ++k) {
                        if(!treeEft[k]) continue;
                        rs_parents[k].create(rs_trees[k]);
                    }

                    //rooting us_tree for iteration
                    unsigned int reg_leaf_adj;
                    BOOST_FOREACH(const unsigned int &w, us_tree.adjacent(rgft_side))
                        if(w!=reg_leaf && w!=prn_side) reg_leaf_adj = w;

                    us_tree.addRoot(reg_leaf,rgft_side);  //root it for traversal
                    if((bestScore-score) > EPSILON) {
                    bestTree = us_tree; bestScore = score; }
                    aw::SubtreeParent<aw::Tree> us_parent; us_parent.create(us_tree);

                    unsigned int last_a, last_b, last_c, a1, b1, c1;
                    std::string last_dir;
                    bool fake = false;

                    if(us_tree.is_fake(reg_leaf_adj) || us_tree.is_leaf(reg_leaf_adj)) continue;

                    //*************************     Starting MOVE-DOWN thing     **************************************************************************************
                    for (aw::Tree::iterator_dfs p=us_tree.begin_dfs(reg_leaf_adj,rgft_side),pEE=us_tree.end_dfs(); p!=pEE; ++p) {                        
                        if(p.idx == reg_leaf_adj) continue;
                        
                        if(fake && !us_tree.is_fake(p.idx)) continue;

                        if(constr) {  //FOR CONSTRAINT..............
                            if(us_tree.constr_num(prn_side)!=NONODE) { //x is root of a clade
                                if(us_tree.in_cld(p.idx)!=NONODE && us_tree.constr_num(p.idx)==NONODE)
                                    continue;
                            } else if(us_tree.constr_num(prn_side)==NONODE && us_tree.in_cld(prn_side)!=NONODE) { //x inside a clade
                                if(us_tree.in_cld(prn_side)!=us_tree.in_cld(p.idx))
                                    continue;
                            } else {  //x is no where clade
                                if(us_tree.in_cld(p.idx)!=NONODE && us_tree.constr_num(p.idx)==NONODE)
                                    continue;
                            }
                        }

                        //Moving subtree X from {a1,b1} to edge {b1,c1}
                        switch (p.direction) {
                            case aw::PREORDER: {                                
                                if(p.parent != reg_leaf_adj) {
                                    if(last_dir=="PRE") {
                                        a1 = last_b; b1 = last_c; c1 = p.idx;
                                    } else {
                                        a1 = last_c; b1 = last_b; c1 = p.idx;
                                    }
                                } else { //in the start of traversal
                                    a1 = reg_leaf; b1 = reg_leaf_adj; c1 = p.idx;
                                }
                                last_dir = "PRE";
                                if(us_tree.is_fake(p.idx))
                                    fake = !fake;                                
                            } break;

                            case aw::POSTORDER: {                                
                                if(last_dir=="PRE") {
                                    a1 = last_c; b1 = last_b; c1 = last_a;
                                } else {
                                    a1 = last_b; b1 = last_c; c1 = us_parent.parent(last_c);
                                    if(c1 == rgft_side) c1 = reg_leaf;
                                }
                                last_dir = "POST";
                                if(us_tree.is_fake(p.idx)) fake = !fake;
                            } break;
                            default: {continue;} break;
                        }
                        
                        last_a = a1; last_b = b1; last_c = c1;
                        if(us_tree.is_fake(b1) && us_tree.is_leaf(c1))
                            ERROR_exit("Wrong move-down");

                        //find the score of each tree when regrafted x-subtree at edge {b1,c1} from {a1,b1}
                        for (unsigned int i=0,iEE=rs_trees.size(); i<iEE; ++i) {
                            if(!treeEft[i]) continue;
                            if(rs_parents[i].parent(c1)==b1 && rs_parents[i].parent(b1)==rgft_side) {
                                unsigned int real_a1;
                                if(rs_parents[i].parent(rgft_side)==a1) real_a1 = a1;
                                else if(rs_parents[i].parent(rgft_side)==0) real_a1=0;
                                else ERROR_exit("Error in the tree");

                                unsigned int sib_c1 = rs_parents[i].sibling_binary(c1);
                                rs_trees[i].moveSub(real_a1,b1,c1,rgft_side);  //update tree
                                rs_parents[i].tPtrUpdate(rs_trees[i]); //update parent-child relationships
                                rs_parents[i].update(c1,rgft_side);
                                rs_parents[i].update(rgft_side,b1);
                                rs_parents[i].update(b1,real_a1);

                                //update lca and score
                                unsigned int old_b1_map = s_lmaps[i].mapping(b1);
                                unsigned int old_yy_map = s_lmaps[i].mapping(rgft_side);
                                unsigned int new_yy_map = g_lca[i].lca(s_lmaps[i].mapping(prn_side),s_lmaps[i].mapping(c1));
                                unsigned int old_b1_map_scr = g_trees[i].return_score(old_b1_map);
                                g_score[i] = g_score[i] + rc::old_map_chg(rs_trees[i],g_trees[i],b1,c1,sib_c1,old_b1_map,old_b1_map_scr);
                                g_trees[i].update_score(old_b1_map,old_b1_map_scr);
                                s_lmaps[i].set_LCA(b1,old_yy_map);
                                s_lmaps[i].set_LCA(rgft_side,new_yy_map);
                                unsigned int new_yy_map_scr = g_trees[i].return_score(new_yy_map);
                                g_score[i] = g_score[i] + rc::new_map_chg(rs_trees[i],g_trees[i],rgft_side,prn_side,c1,new_yy_map,new_yy_map_scr);
                                g_trees[i].update_score(new_yy_map,new_yy_map_scr);

                                //update clusters
                                rs_trees[i].update_clst(b1,rs_trees[i].return_clstSz(rgft_side));
                                rs_trees[i].update_clst(rgft_side,rs_trees[i].return_clstSz(prn_side)+rs_trees[i].return_clstSz(c1));
                            } else if(rs_parents[i].parent(rgft_side)==b1 && rs_parents[i].parent(c1)==b1) {
                                rs_trees[i].moveSub(a1,b1,c1,rgft_side);  //update tree
                                rs_parents[i].tPtrUpdate(rs_trees[i]); //update parent-child relationships
                                rs_parents[i].update(c1,rgft_side);
                                rs_parents[i].update(rgft_side,b1);
                                rs_parents[i].update(a1,b1);

                                //update lca and score
                                unsigned int old_yy_map = s_lmaps[i].mapping(rgft_side);
                                unsigned int new_yy_map = g_lca[i].lca(s_lmaps[i].mapping(prn_side),s_lmaps[i].mapping(c1));
                                unsigned int old_yy_map_scr = g_trees[i].return_score(old_yy_map);
                                g_score[i] = g_score[i] + rc::old_map_chg(rs_trees[i],g_trees[i],rgft_side,prn_side,a1,old_yy_map,old_yy_map_scr);
                                g_trees[i].update_score(old_yy_map,old_yy_map_scr);
                                s_lmaps[i].set_LCA(rgft_side,new_yy_map);
                                unsigned int new_yy_map_scr = g_trees[i].return_score(new_yy_map);
                                g_score[i] = g_score[i] + rc::new_map_chg(rs_trees[i],g_trees[i],rgft_side,prn_side,c1,new_yy_map,new_yy_map_scr);
                                g_trees[i].update_score(new_yy_map,new_yy_map_scr);

                                //update clusters
                                rs_trees[i].update_clst(rgft_side,rs_trees[i].return_clstSz(prn_side)+rs_trees[i].return_clstSz(c1));
                            }  else if(rs_parents[i].parent(a1)==rgft_side && rs_parents[i].parent(rgft_side)==b1) {
                                unsigned int real_c1;
                                if(rs_parents[i].parent(b1)==c1) real_c1 = c1;
                                else if(rs_parents[i].parent(b1)==0) real_c1=0;
                                else ERROR_exit("Error in the tree");

                                unsigned int sib_yy = rs_parents[i].sibling_binary(rgft_side);
                                rs_trees[i].moveSub(a1,b1,real_c1,rgft_side);  //update tree
                                rs_parents[i].tPtrUpdate(rs_trees[i]); //update parent-child relationships
                                rs_parents[i].update(b1,rgft_side);
                                rs_parents[i].update(rgft_side,real_c1);
                                rs_parents[i].update(a1,b1);

                                //update lca and score
                                unsigned int old_yy_map = s_lmaps[i].mapping(rgft_side);
                                unsigned int old_b1_map = s_lmaps[i].mapping(b1);
                                unsigned int new_b1_map = g_lca[i].lca(s_lmaps[i].mapping(sib_yy),s_lmaps[i].mapping(a1));
                                unsigned int old_yy_map_scr = g_trees[i].return_score(old_yy_map);
                                g_score[i] = g_score[i] + rc::old_map_chg(rs_trees[i],g_trees[i],rgft_side,prn_side,a1,old_yy_map,old_yy_map_scr);
                                g_trees[i].update_score(old_yy_map,old_yy_map_scr);
                                s_lmaps[i].set_LCA(rgft_side,old_b1_map);
                                s_lmaps[i].set_LCA(b1,new_b1_map);
                                unsigned int new_b1_map_scr = g_trees[i].return_score(new_b1_map);
                                g_score[i] = g_score[i] + rc::new_map_chg(rs_trees[i],g_trees[i],b1,sib_yy,a1,new_b1_map,new_b1_map_scr);
                                g_trees[i].update_score(new_b1_map,new_b1_map_scr);

                                //update clusters
                                rs_trees[i].update_clst(rgft_side,rs_trees[i].return_clstSz(b1));
                                rs_trees[i].update_clst(b1,rs_trees[i].return_clstSz(sib_yy)+rs_trees[i].return_clstSz(a1));
                            }
                            else  ERROR_exit("SOME ERROR");
                        }

                        score = 0;
                        for (unsigned int mm=0,mmEE=g_trees.size(); mm<mmEE; ++mm)
                            score = score + g_score[mm]*g_weights[mm] ;

                        if((bestScore-score) > EPSILON) {
                            for(int mn=0, mnEE=treeEft.size(); mn<mnEE; ++mn)
                                if(treeEft[mn]) {
                                    bestTree = rs_trees[mn];
                                    break; }
                            bestScore = score; }
                    }                    
                }                
            }
            reg_x = !reg_x;
        }

        MSG_nonewline('\r');
        MSG_nonewline("Current RF Score: "<<std::fixed<<std::setprecision(2)<< bestScore);

        s_tree = bestTree;       

        {   //should root s_tree at right place: not below fake internal node
            std::vector<unsigned int> ch;
            s_tree.adjacent(0,ch);
            if(s_tree.is_fake(ch[0]) && s_tree.is_leaf(ch[1]) || s_tree.is_fake(ch[1]) && s_tree.is_leaf(ch[0])) {
                unsigned int nf;
                if(s_tree.is_fake(ch[0])) nf = ch[0]; else nf = ch[1];
                ch.clear();
                BOOST_FOREACH(const unsigned int &c, s_tree.adjacent(nf))
                    if(!s_tree.is_leaf(c) && c!=0)
                         s_tree.rootBy(nf,c);
            }
        }

        if(bestScore == 0 || bestScore==scr)  break;  //exit if no improvement or score is already zero
        
        { // root the supertree copies again by one leaf
            rs_trees.clear();  unsigned int rt;
            for (unsigned int k=0, iEE=g_trees.size(); k<iEE; ++k) {               
                BOOST_FOREACH(const unsigned int &w, g_trees[k].adjacent(0))
                    if(g_trees[k].is_leaf(w)) rt = w;                
                unsigned int ggid = g_nmaps[k].gid(rt);
                std::vector<unsigned int> ch;
                s_nmap.ids(ggid,ch);
                unsigned int s_id;
                BOOST_FOREACH(const unsigned int &c,ch)
                    if(s_lmaps[k].mapping(c)==rt)
                        s_id = c;                
                rs_trees.push_back(s_tree);
                rs_trees.back().rootBy(s_id);                
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
                        rs_trees[k].update_clst(v.idx,count);      }
                    else {                                               
                        if(s_lmaps[k].mapping(v.idx)!=NONODE)    //:FOR MUL-TREES
                            rs_trees[k].update_clst(v.idx,1);
                        else rs_trees[k].update_clst(v.idx,0);   }                    
                }
            }
        }

        g_scr.clear(); scr = 0;
        {   //no need to do LCA computations again...
            for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
                s_lmaps[i].update_LCA_internals(g_lca[i],rs_trees[i]); }
            for (unsigned int i=0,iEE=g_trees.size(); i<iEE; ++i) {
                std::pair<unsigned int,unsigned int> p = g_nodes[i];
                g_scr.push_back(aw::compute_rf_score(rs_trees[i],g_trees[i],s_lmaps[i],p,rs_int_nodes[i]));
                scr = scr + g_scr[i]*g_weights[i] ;
            }            
            if(fabs(scr-bestScore)>EPSILON) ERROR_exit("SCR and bestScore doesn't match!!!");
        }        
    }

    MSG("\nSPR neighborhood searches: "<<SPR_rounds);

    {   //outputing input trees and output super tree
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

            //if constr then reroot
            if(constr) {
                if(s_tree.root==0) {
                    TREE_FOREACHLEAF(lf, s_tree) {  //Added on 29th May
                        if(s_tree.in_cld(lf)==NONODE)
                            s_tree.rootBy(lf);
                    }
                }
            }
            
            output <<"[ Species Tree: Unrooted RF Score = "<<std::fixed<<std::setprecision(2)<<bestScore<<"]"<< std::endl;
            aw::tree2newick(output,s_tree,s_taxa); output << std::endl;
        }
        
        if(inputrees) {            
            for(int mn=0, mnEE=g_trees.size(); mn<mnEE; ++mn) {
                 output <<"\n[ Gene Tree "<<mn<< " MulRF Score = "<<std::fixed<<std::setprecision(2)<<g_scr[mn]*g_weights[mn]<<"]"<< std::endl;
                 output<<"[&WEIGHT="<<std::fixed<<std::setprecision(2)<<g_weights[mn]<<"]";
                 aw::tree2newick(output,g_trees[mn],g_taxa[mn]); output << std::endl; } }
    }

    int t4 = clock();
    {   //for timing...
        long ttime = ((long)(t4-t3))/CLOCKS_PER_SEC;
        int d,h,m,s;
        util::convertTime(ttime,d,h,m,s);
        MSG_nonewline("Total elapsed time: ");
        if(d!=0) MSG_nonewline(d<<"d ");
        if(h!=0) MSG_nonewline(h<<"h ");
        if(m!=0) MSG_nonewline(m<<"m ");
        MSG_nonewline(s<<"s ");
        output <<"\n[ Time "<<d<<"d "<<h<<"h "<<m<<"m "<<s<<"s "<<"]"<< std::endl;
    }


}

