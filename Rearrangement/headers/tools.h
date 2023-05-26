#ifndef TOOLS_H
#define TOOLS_H

#include "structures.h"
#include <cmath>
#include <tuple>

std::vector<std::string> my_split(std::string line, std::string sc)
{
    std::vector<std::string> parts;
    std::string part;
    for(unsigned i=0; i<line.size(); ++i)
    {
        if(sc.find(line[i]) != std::string::npos)
        {
            parts.push_back(std::move(part));
            part.clear();
        }
        else
            part.push_back(line[i]);
    }
    parts.push_back(std::move(part));
    return parts;
}


void div_num_fun(std::vector<Align>::iterator first, std::vector<Align>::iterator last)
{
    if (first==last) 
        return;
    std::vector<Align>::iterator tracer = first;
    double index_num=1;
    while (++first != last)
    {
        if((first-1)->index!=first->index)
        {
            for(;tracer<first;++tracer)
            {
                tracer->num/=index_num;
            }
            index_num=1;
        }
        else
            index_num+=1;
    }
    for(;tracer<first;++tracer)
    {
        tracer->num/=index_num;
    }
}

class MY_LESS
{
    std::string& ref;
public:
    MY_LESS(std::string& x) : ref(x) {}
    bool operator()(const std::tuple<int,int,std::string>& LRM1, const std::tuple<int,int,std::string>& LRM2) const
    {
        int i=std::min(std::get<0>(LRM1),std::get<0>(LRM2));
        int up1=ref.size()+std::get<0>(LRM1)+std::get<2>(LRM1).size()-std::get<1>(LRM1);
        int up2=ref.size()+std::get<0>(LRM2)+std::get<2>(LRM2).size()-std::get<1>(LRM2);
        while(i<up1 && i<up2)
        {
            char C1, C2;
            if(i<std::get<0>(LRM1)) 
                C1=ref[i];
            else if(i<std::get<0>(LRM1)+std::get<2>(LRM1).size())
                C1=std::get<2>(LRM1)[i-std::get<0>(LRM1)];
            else
                C1=ref[i-std::get<0>(LRM1)-std::get<2>(LRM1).size()+std::get<1>(LRM1)];
            if(i<std::get<0>(LRM2)) 
                C2=ref[i];
            else if(i<std::get<0>(LRM2)+std::get<2>(LRM2).size())
                C2=std::get<2>(LRM2)[i-std::get<0>(LRM2)];
            else 
                C2=ref[i-std::get<0>(LRM2)-std::get<2>(LRM2).size()+std::get<1>(LRM2)];
            
            if(C1<C2)
                return true;
            else if(C1>C2)
                return false;
            ++i;
        }
        if(up1<up2)
            return true;
        else
            return false;
    }
};

#endif