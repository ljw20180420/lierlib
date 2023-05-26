#ifndef MODE_MUTATE_H
#define MODE_MUTATE_H

#include <iostream>
#include <random>
#include <fstream>


std::string inner_mutate(std::string seg)
{
    std::string mut;
    double mut_rate;
    int row=-1;
    while(row<int(seg.size()-1))
    {
        mut_rate=uni_real_dis(generator);
        if(mut_rate<MUT_RATE)
            mut.push_back(nuc[nuc_dis(generator)]);
        else if(mut_rate<2*MUT_RATE)
            ++row;
        else
        {
            ++row;
            mut_rate=uni_real_dis(generator);
            if(mut_rate<MUT_RATE)
                mut.push_back(nuc[nuc_dis(generator)]);
            else
                mut.push_back(seg[row]);
        }
    }
    while(row<int(seg.size()))
    {
        mut_rate=uni_real_dis(generator);
        if(mut_rate<MUT_RATE)
            mut.push_back(nuc[nuc_dis(generator)]);
        else
            ++row;
    }
    return mut;
}

std::string random_seg(int seg_len)
{
    std::string seg;
    for(int i=0; i<seg_len; ++i)
        seg.push_back(nuc[nuc_dis(generator)]);
    return seg;
}

std::string random_cle(std::string seg)
{
    int r1=cle_dis(generator);
    int r2=seg.size()-cle_dis(generator);
    if(r1<r2)
        return seg.substr(r1, r2-r1);
    else
        return std::string();
}


void mode_mutate(const int SEG_NUM, const int SEG_UP, const int SEG_LOW, const int CLE_UP, const int INS_UP, const double MUT_RATE, const int READ_NUM)
{
    std::default_random_engine generator;
    std::uniform_int_distribution<int> seg_dis(SEG_LOW,SEG_UP);
    std::uniform_int_distribution<int> cle_dis(0,CLE_UP);
    std::uniform_int_distribution<int> ins_dis(0,INS_UP);
    std::uniform_int_distribution<int> nuc_dis(0,3);
    std::uniform_real_distribution<double> uni_real_dis(0.0,1.0);
    const std::string nuc("ACGT");


    std::vector<std::string> ref;
    int cum_len=0;
    std::ofstream fout_truth;
    std::ofstream fout("reference");
    for(int i=0; i<SEG_NUM; ++i)
    {
        ref.push_back(random_seg(seg_dis(generator)));
        int right_exp=cum_len+CLE_UP/2;
        cum_len+=ref.back().size();
        int left_exp=cum_len-CLE_UP/2;
        fout << right_exp << '\t' << right_exp-CLE_UP/6 << '\t' << right_exp+CLE_UP/6 << '\n' << ref.back() << '\n' << left_exp << '\t' << left_exp-CLE_UP/6 << '\t' << left_exp+CLE_UP/6 << '\n';
    }
    fout.close();
    
    fout.open("local");
    fout_truth.open("local_truth");
    for(int r=0; r<READ_NUM; ++r)
    {
        std::string read;
        std::string read_truth;
        std::string tmp;
        tmp=random_seg(ins_dis(generator));
        read.append(tmp);
        read_truth.append(tmp);
        read_truth.push_back('|');
        for(int i=0; i<ref.size(); ++i)
        {
            tmp=inner_mutate(random_cle(ref[i]));
            read.append(tmp);
            read_truth.append(tmp);
            read_truth.push_back('|');
            tmp=random_seg(ins_dis(generator));
            read.append(tmp);
            read_truth.append(tmp);
            if(i!=ref.size()-1)
                read_truth.push_back('|');
        }
        fout << read << '\n';
        fout_truth << read_truth << '\n';
    }
    fout.close();
    fout_truth.close();

    fout.open("local_imbed");
    fout_truth.open("local_imbed_truth");
    for(int r=0; r<READ_NUM; ++r)
    {
        std::string read;
        std::string read_truth;
        std::string tmp;
        for(int i=0; i<ref.size(); ++i)
        {
            tmp=inner_mutate(random_cle(ref[i]));
            read.append(tmp);
            read_truth.append(tmp);
            if(i!=ref.size()-1)
            {
                read_truth.push_back('|');
                tmp=random_seg(ins_dis(generator));
                read.append(tmp);
                read_truth.append(tmp);
                read_truth.push_back('|');
            }
        }
        fout << read << '\n';
        fout_truth << read_truth << '\n';
    }
    fout.close();
    fout_truth.close();
    
    fout.open("imbed");
    fout_truth.open("imbed_truth");
    for(int r=0; r<READ_NUM; ++r)
    {
        std::string read;
        std::string read_truth;
        std::string tmp;
        for(int i=0; i<ref.size(); ++i)
        {
            tmp=inner_mutate(random_cle(ref[i]));
            read.append(tmp);
            read_truth.append(tmp);
            if(i!=ref.size()-1)
                read_truth.push_back('|');
        }
        fout << read << '\n';
        fout_truth << read_truth << '\n';
    }
    fout.close();
    fout_truth.close();

    fout.open("contain");
    fout_truth.open("contain_truth");
    for(int r=0; r<READ_NUM; ++r)
    {
        std::string read;
        std::string read_truth;
        std::string tmp;
        tmp=random_seg(ins_dis(generator));
        read.append(tmp);
        read_truth.append(tmp);
        read_truth.push_back('|');
        for(int i=0; i<ref.size(); ++i)
        {
            tmp=inner_mutate(ref[i]);
            read.append(tmp);
            read_truth.append(tmp);
            read_truth.push_back('|');
            tmp=random_seg(ins_dis(generator));
            read.append(tmp);
            read_truth.append(tmp);
            if(i!=ref.size()-1)
                read_truth.push_back('|');
        }
        fout << read << '\n';
        fout_truth << read_truth << '\n';
    }
    fout.close();
    fout_truth.close();
}

#endif