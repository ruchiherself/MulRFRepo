/* 
 * File:   rf_compute.h
 * Author: ruchi
 *
 * Created on December 3, 2010, 1:47 PM
 */

#ifndef _RF_COMPUTE_H
#define	_RF_COMPUTE_H

#include "tree.h"
#include "tree_subtree_info.h"
#include "tree_LCA_mapping.h"


#include <stdio.h>
#include <vector>

namespace rc {
    static const unsigned int NONODE = UINT_MAX;

    using namespace std;

    //rs_tree and g_tree are pass-by-reference for optimization
    inline signed int old_map_chg(aw::Tree &rs_tree, aw::Tree &g_tree, unsigned int s, unsigned int c, unsigned int d, unsigned int s_map, unsigned int &s_map_score){
        signed int score_chg = 0;
        unsigned int n_score;        

        if(s_map!=NONODE && !g_tree.is_leaf(s_map) && g_tree.root!=s_map) {
            //og2_score = g_tree.return_score(g2);
            n_score = s_map_score;
            //cout<<"c clstSz(c):"<<c<<" "<<rs_tree.return_clstSz(c)<<"d clstSz(d)"<<d<<" "<<rs_tree.return_clstSz(d)<<"g2 clstSz(g2)"<<g2<<" "<<g_tree.return_clstSz(g2);
            if((rs_tree.return_clstSz(c) + rs_tree.return_clstSz(d))== g_tree.return_clstSz(s_map))
                n_score = n_score - 1;
            //P(g2<<" "<<n_score<<" "<<og2_score);
            //cout<<"\n < new_score and ols_score"<<n_score<<og2_score<<">";
            if(n_score==0 && s_map_score>0){
                score_chg = score_chg + 2;
                //cout<<" -score_dec- ";
            }
            s_map_score = n_score;
        }        
        return score_chg;
    }

    //rs_tree and g_tree are pass-by-reference for optimization
    inline signed int new_map_chg(aw::Tree &rs_tree, aw::Tree &g_tree, unsigned int s, unsigned int c, unsigned int d, unsigned int s_map, unsigned int &s_map_score){
        signed int score_chg = 0;
        unsigned int n_score;

        //MSG_nonewline("s,g1,g2,og1_score,og2_score:"<<s<<"-"<<g1<<"-"<<g2<<"-"<<og1_score<<"-"<<og2_score);

        if(s_map!=NONODE && !g_tree.is_leaf(s_map) && g_tree.root!=s_map) {
            //og2_score = g_tree.return_score(g2);
            n_score = s_map_score;
            //cout<<"c clstSz(c):"<<c<<" "<<rs_tree.return_clstSz(c)<<"d clstSz(d)"<<d<<" "<<rs_tree.return_clstSz(d)<<"g2 clstSz(g2)"<<g2<<" "<<g_tree.return_clstSz(g2);
            if((rs_tree.return_clstSz(c) + rs_tree.return_clstSz(d))== g_tree.return_clstSz(s_map))
                n_score = n_score + 1;
            //P(g2<<" "<<n_score<<" "<<og2_score);
            //cout<<"\n < new_score and ols_score"<<n_score<<og2_score<<">";
            if(n_score>0 && s_map_score==0){
                score_chg = score_chg - 2;
                //cout<<" -score_dec- ";
            }
            s_map_score = n_score;
        }       
        return score_chg;
    }

}

#endif	/* _RF_COMPUTE_H */

