/* 
 * File:   main.cpp
 * Author: ruchi
 *
 * Created on April 2, 2014, 3:28 PM
 */

#include <stdlib.h>
#include "argument.h"
#include "common.h"
#include "tree.h"
#include "tree_IO.h"
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

    {
        Argument a; a.add(ac, av);
        // help
        if (a.existArg2("-h","--help")) {
            MSG("options:");
            MSG("  -i [ --input ] arg      input trees (file in NEWICK format)");
            MSG("  -o [ --output ] arg     output file");
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
    std::vector<aw::Tree> g_trees;
    std::vector<aw::idx2name> g_taxa;    
    
    {
        const std::string filename = trees_filename;
        { // read trees
            std::ifstream ifs;
            std::istream &is = filename.empty() ? std::cin : ifs;
            if (!filename.empty()) {
                ifs.open(filename.c_str());
                if (!ifs) ERROR_exit("cannot read file '" << filename << "'");
            }
            
            MSG_nonewline("Reading input trees: ");
            
            for (;;) {
                aw::Tree t;
                aw::idx2name t_names;
                float t_w = 1.0f;
                if (!aw::stream2tree(is, t, t_names,t_w)) break;
                g_taxa.push_back(t_names);
                g_trees.push_back(t);
                output <<"Input gene tree of taxa = "<<t_names.size()<< " and weight = "<<std::fixed<<std::setprecision(2)<<t_w<< std::endl;                
            }            
        }
        MSG("Input trees: " << g_trees.size());
        if (g_trees.empty()) ERROR_exit("No input trees found in file '" << filename << "'");
    }

    
    return (EXIT_SUCCESS);
}

