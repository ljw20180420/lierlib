#ifndef PARSER_H
#define PARSER_H

#include "tools.h"
#include <algorithm>
#include <filesystem>
#include <thread>
#include <cstring>
#include <stdlib.h>
#include <iostream>

void print_help(int argc, char **argv)
{
    for(int i=1; i<argc; ++i)
        if(!strcmp(argv[i],"-help"))
        {
            std::cout << "###Basic Usage\n"
            << "rearrangement -file input_file -ref_file reference_file -ALIGN_MAX 1 -THR_NUM 6\n"

            << "\n### Parameters\n"
            << "\n# Help Related\n"
            << "-help: Display help.\n"

            << "\n# Input Files\n"
            << "-file: Input file of NGS reads.\n"
            << "-ref_file: Reference file.\n"

            << "\n# Aligning Engine\n"
            << "-DIVCON: Use divide-and-conquer linear space method. This parameter is valid only if -ALIGN_MAX is 1. (default: false)\n"

            << "\n# Aligning Parameters\n"
            << "-s0: Mismatching score. (default: -3)\n"
            << "-s1: Matching score. (default: +1)\n"
            << "-u: Gap-extending penalty. (default: -2)\n"
            << "-v: Gap-opening penalty. (default: -5)\n"
            << "-alg_type: Method to scoring unaligned parts between aligned segments in query sequence and reference (-alg_type local|imbed|contain|local_imbed). (default: local)\n"          
            << "-per_thres: Percentage threshold of repeated reads. The same read may repeat many times in the file. Percentage of a read is calculated by (repeated # of the read in the file)/(total # of reads in the file). Reads with percentages less than -per_thres are excluded from analyses. (default: 0.0)\n"

            << "\n# Output Options\n"
            << "-mode: The method to explain the indel (-mode overlapping|nonoverlapping). (default: nonoverlapping)\n"
            << "-ALIGN_MAX: The maximally reserved number of best alignments for each read (Each read may have several best alignments to the reference). (default: 5)\n"

            << "\n# Analyze Indels\n"
            << "-indel : Post-analysis for alignments. (default : false)\n"

            << "\n# EM Configs\n"
            << "-EM : Predict the cut distribution by EM algorithm. (default: false)\n"
            << "-ini_alpha: Path to the file used to initialize the distribution of the left ligation ends. If not specified, a uniform distribution over the ligation range will be used. Note that the file must be consistent with the information in reference sequences.\n"
            << "-ini_beta: Similar as -ini_alpha but for the right ligation ends.\n"
            << "-ini_pi: Path to the file used to initialize the distribution of the middle insertions. If not specified, a uniform distribution over all possible middle insertions upper to length specified by -MID_MAX will be used.\n"
            << "-MID_MAX: The maximal length of middle insertion. (default: 0)\n"
            << "-thres: The terminal threshold for EM algorithm. (default: 0.000001)\n"

            << "\n# Parallel Options\n"
            << "-SEQ_BATCH: Set the number of reads passed in one time. (default: 300)\n"  
            << "-THR_NUM: Thread number. (default: the hardware available threads)\n";
            exit(0);
        }
    return;
}

struct Command_content
{
    // Input Files
    std::string file="";
    std::string ref_file="";
    // Aligning Engine
    bool DIVCON=false;
    // Aligning Parameters
    int s0=-3;
    int s1=1;
    int u=-2;
    int v=-5;
    std::string alg_type="local";
    double per_thres=0.0;
    // Output Options
    std::string mode="nonoverlapping";
    int ALIGN_MAX=5;
    // Analyze Indels
    bool indel=false;
    // EM Configs
    bool EM=false;
    std::string ini_alpha="";
    std::string ini_beta="";
    std::string ini_pi="";
    int MID_MAX=0;
    double thres=1e-6;
    // Parallel Options
    int SEQ_BATCH=300;
    int THR_NUM;
};

Command_content command(int argc, char **argv)
{
    Command_content cc;
    cc.THR_NUM = std::thread::hardware_concurrency();

    for(int i=1; i<argc-1; ++i)
    {
        // Input Files
        if(!strcmp(argv[i],"-file"))
            cc.file=argv[i+1];
        if(!strcmp(argv[i],"-ref_file"))
            cc.ref_file=argv[i+1];
        // Aligning Engine
        if(!strcmp(argv[i],"-DIVCON"))
        {
            if(!strcasecmp(argv[i+1],"false"))
                cc.DIVCON=false;
            else
            {
                if(!strcasecmp(argv[i+1],"true"))
                    cc.DIVCON=true;
                else
                {
                    std::cerr << "-DIVCON must be true or false\n";
                    exit(-1);
                }
            }
        }
        // Aligning Parameters
        if(!strcmp(argv[i],"-s0"))
            cc.s0=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-s1"))
            cc.s1=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-u"))
            cc.u=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-v"))
            cc.v=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-alg_type"))
            cc.alg_type=argv[i+1];
        if(!strcmp(argv[i],"-per_thres"))
            cc.per_thres=atof(argv[i+1]);
        // Output Options
        if(!strcmp(argv[i],"-mode"))
            cc.mode=argv[i+1];
        if(!strcmp(argv[i],"-ALIGN_MAX"))
            cc.ALIGN_MAX=atoi(argv[i+1]);
        // Analyze Indels
        if(!strcmp(argv[i],"-indel"))
        {
            if(!strcasecmp(argv[i+1],"false"))
                cc.indel=false;
            else
            {
                if(!strcasecmp(argv[i+1],"true"))
                    cc.indel=true;
                else
                {
                    std::cerr << "-indel must be true or false\n";
                    exit(-1);
                }
            }
        }
        // EM Configs
        if(!strcmp(argv[i],"-EM"))
        {
            if(!strcasecmp(argv[i+1],"false"))
                cc.EM=false;
            else
            {
                if(!strcasecmp(argv[i+1],"true"))
                    cc.EM=true;
                else
                {
                    std::cerr << "-EM must be true or false\n";
                    exit(-1);
                }
            }
        }
        if(!strcmp(argv[i],"-ini_alpha"))
            cc.ini_alpha=argv[i+1];
        if(!strcmp(argv[i],"-ini_beta"))
            cc.ini_beta=argv[i+1];
        if(!strcmp(argv[i],"-ini_pi"))
            cc.ini_pi=argv[i+1];
        if(!strcmp(argv[i],"-MID_MAX"))
            cc.MID_MAX=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-thres"))
            cc.thres=atof(argv[i+1]);
        // Parallel Options
        if(!strcmp(argv[i],"-SEQ_BATCH"))
            cc.SEQ_BATCH=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-THR_NUM"))
            cc.THR_NUM=atoi(argv[i+1]);
    }
    
    if (cc.file.empty())
    {
        std::cerr << "-file is not specified\n";
        exit(-1);
    }
    if (cc.ref_file.empty())
    {
        std::cerr << "-ref_file is not specified\n";
        exit(-1);
    }
    if (cc.mode!="overlapping" && cc.mode!="nonoverlapping")
    {
        std::cerr << "-mode must be overlapping or nonoverlapping\n";
        exit(-1);
    }
    
    return cc;
}

#endif