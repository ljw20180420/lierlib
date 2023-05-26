#ifndef OUTPUT_H
#define OUTPUT_H

#include "structures.h"
#include <fstream>
#include <map>

void indel_label_fun(std::vector<Align>& aligns, std::vector<int>& left_exp, std::vector<int>& right_exp, std::string& file)
{
    std::vector<std::vector<std::array<int,3>>> LPRPILs(aligns.size());
    std::array<int,3> LPRPIL;
    for(auto & align : aligns)
    {        
        for(int i=1; i<right_exp.size(); ++i)
        {
            LPRPIL[0]=align.left[i-1]-left_exp[i-1];
            LPRPIL[1]=align.right[i]-right_exp[i];
            LPRPIL[2]=align.MID[i].size();
            LPRPILs[align.index-1].push_back(LPRPIL);
        }
    }
    std::ofstream fout(file+".indel");
    for(auto & LPRPIL_row : LPRPILs)
    {
        for(int i=0; i<LPRPIL_row.size(); ++i)
        {
            fout << LPRPIL_row[i][0] << '\t' << LPRPIL_row[i][1] << '\t' << LPRPIL_row[i][2];
            if(i!=LPRPIL_row.size()-1)
                fout << '\t';                
        }
        fout << '\n';
    }
    fout.close();
}

#endif