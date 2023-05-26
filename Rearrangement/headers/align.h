#ifndef ALIGN_H
#define ALIGN_H

#include "structures.h"
#include <tuple>
#include <filesystem>

std::tuple<std::vector<int>,std::vector<int>> extract_to_pair(std::string& reference, std::string& query, std::vector<int>& left, std::vector<int>& right)
{
    std::vector<int> left_ext, right_ext;
    int ref_cum=0, pair_now;
    bool search;
    if(ref_cum==right[0])
        search=true;
    else
        search=false;
    for(int i=0; i<reference.size(); ++i)
    {
        if(reference[i]!='-')
        {
            ++ref_cum;
            if(query[i]!='-')
            {
                pair_now=ref_cum;
                if(search)
                {
                    right_ext.push_back(pair_now-1);
                    search=false;
                }
            }
            if(ref_cum==left[left_ext.size()])
            {
                if(search || left[left_ext.size()]==right[left_ext.size()])
                {
                    right_ext.push_back((right[right_ext.size()]+left[left_ext.size()])/2);
                    search=false;
                    left_ext.push_back(right_ext.back());
                }
                else
                    left_ext.push_back(pair_now);
            }
            if(right_ext.size()<right.size() && ref_cum==right[right_ext.size()])
                search=true;
            if(left_ext.size()==left.size())
                break;
        }
    }
    return std::make_tuple(left_ext,right_ext);
}

std::tuple<std::vector<int>,std::vector<int>> normal_2_c(std::vector<int>& left, std::vector<int>& right, std::vector<int>& left_exp, std::vector<int>& right_exp, std::string& mode)
{
    std::vector<int> c_left, c_right;
    if(mode=="nonoverlapping")
    {
        for(int i=0; i<left.size(); ++i)
        {
            c_left.push_back(std::min(left[i],left_exp[i]));
            c_right.push_back(std::max(right[i],right_exp[i]));
        }
    }
    else if(mode=="overlapping")
    {
        c_right.push_back(right.front());
        for(int i=0; i<left.size()-1; ++i)
        {
            c_left.push_back(std::min(left[i],left_exp[i]+std::max(right[i+1]-right_exp[i+1],0)));
            c_right.push_back(std::max(right[i+1],right_exp[i+1]+std::min(left[i]-left_exp[i],0)));
        }
        c_left.push_back(left.back());
    }
    for(int i=0; i<left.size(); ++i)
    {
        if(c_right[i]>c_left[i])
        {
            if(c_right[i]==right[i])
                c_left[i]=c_right[i];
            else if(c_left[i]==left[i])
                c_right[i]=c_left[i];
        }
    }
    return std::make_tuple(c_left,c_right);
}

std::vector<std::string> get_mid(std::string& reference, std::string& query, std::vector<int>& left, std::vector<int>& right)
{
    std::vector<std::string> MID(1,std::string());
    int ref_cum=0;
    bool insert=true;
    for(int i=0; i<reference.size(); ++i)
    {
        if(reference[i]!='-')
        {
            ++ref_cum;
            if(MID.size()<=right.size() && ref_cum==right[MID.size()-1]+1)
                insert=false;
        }
        if(insert && query[i]!='-')
            MID.back().push_back(query[i]);
        if(MID.size()<=left.size() && ref_cum==left[MID.size()-1])
        {
            MID.push_back(std::string());
            insert=true;
        }
    }
    return MID;
}

std::tuple<std::vector<int>,std::vector<int>,std::vector<std::string>,std::vector<int>,std::vector<int>,std::vector<std::string>,std::vector<int>,std::vector<int>,std::vector<std::string>> search_left_right(std::string& reference, std::string& query, std::vector<int>& S, std::vector<int>& left_exp, std::vector<int>& right_exp, std::string& mode, std::string& x)
{
    std::vector<int> left, right, c_left, c_right, c_left_in, c_right_in;
    std::vector<int> right_S(S.begin(),S.end()-1);
    std::vector<int> left_S(S.begin()+1,S.end());
    std::tie(left,right)=extract_to_pair(reference, query, left_S, right_S);
    std::tie(c_left,c_right)=normal_2_c(left, right, left_exp, right_exp, mode);
    std::tie(c_left_in,c_right_in)=extract_to_pair(reference, query, c_left, c_right);
    std::vector<std::string> MID=get_mid(reference, query, left, right);
    std::vector<std::string> c_MID_in=get_mid(reference, query, c_left_in, c_right_in);
    std::vector<std::string> c_MID;
    for(int i=0; i<MID.size(); ++i)
    {
        if(i==0)
            c_MID.push_back(MID[i]+x.substr(right[i],c_right[i]-right[i]));
        else if(i==MID.size()-1)
            c_MID.push_back(x.substr(c_left[i-1],left[i-1]-c_left[i-1])+MID[i]);
        else
            c_MID.push_back(x.substr(c_left[i-1],left[i-1]-c_left[i-1])+MID[i]+x.substr(right[i],c_right[i]-right[i]));
    }
    return std::make_tuple(left, right, MID, c_left, c_right, c_MID, c_left_in, c_right_in, c_MID_in);
}


std::tuple<std::map<char,int>,TD_array<int>,std::vector<int>,std::vector<int>,std::vector<int>,std::vector<int>,TD_array<int>,TD_array<int>> initialization(std::vector<int>& S, std::string& alg_type, int u, int v, int max_len, int S0, int S1)
{
    TD_array<int> gamma(5,5,1,S0);
    gamma(1,1,0)=gamma(2,2,0)=gamma(3,3,0)=gamma(4,4,0)=S1;
    std::vector<int> ve(S.back()+1,v), ue(S.back()+1,u), tvf(S.size()-1,v), tuf(S.size()-1,u);
    for(int j=0; j<S.size(); ++j)
    {
        if(alg_type=="local" || alg_type=="contain")
        {
            ve[S[j]]=0;
            ue[S[j]]=0;
        }
        if(alg_type=="local_imbed")
        {
            if(j!=0 && j!=S.size()-1)
            {
                ve[S[j]]=0;
                ue[S[j]]=0;
            }
        }
        if(alg_type=="local" || alg_type=="local_imbed" || alg_type=="imbed")
        {
            if(j!=S.size()-1)
            {
                tvf[j]=0;
                tuf[j]=0;
            }
        }
    }
    return std::make_tuple(std::map<char,int>{{'N',0},{'A',1},{'C',2},{'G',3},{'T',4},{'n',0},{'a',1},{'c',2},{'g',3},{'t',4}},gamma,ve,ue,tvf,tuf,TD_array<int>(S.size()-1,max_len+1,1,v),TD_array<int>(S.size()-1,max_len+1,1,u));
}

void ori_cmp(Back& back, int new_val, Back* address)
{
    if(back.max_val<=new_val)
    {
        if(back.max_val<new_val)
        {
            back.max_val=new_val;
            back.biters.clear();
        }
        back.biters.push_back(address);
    }
}

std::tuple<std::vector<Align>,std::vector<std::string>,std::vector<std::string>> column_wise(std::vector<int>::iterator ve_b, std::vector<int>::iterator ue_b, int row, std::vector<int>::iterator vf_b, std::vector<int>::iterator uf_b, std::vector<int>::iterator tvf_b, std::vector<int>::iterator tuf_b, TD_array<int>& gamma, std::string::iterator x_b, std::string::iterator o_b, std::string::iterator o_e, std::vector<int>::iterator S_b, std::vector<int>::iterator S_e, bool head, bool tail, bool ceil, bool floor, int ALIGN_MAX, std::map<char,int>& nt2int, TD_array<Back>& EFG)
{
    EFG(0,0,2).max_val=0;
    for(int j=0; j<S_e-S_b; ++j)
        for(int s=*(S_b+j)+1; s<=*(S_b+j+1); ++s)
        {
            if(ceil)
                EFG(s,0,2).max_val=EFG(s-1,0,2).max_val+((s==*(S_b+j)+1)?(*(tvf_b+j)):(*(tuf_b+j)));
            else
                EFG(s,0,2).max_val=EFG(s-1,0,2).max_val+((s==*(S_b+j)+1)?(*(vf_b+j)):(*(uf_b+j)));
            EFG(s,0,0).max_val=std::numeric_limits<int>::min()/2;
        }
    for(int w=1; w<=o_e-o_b; ++w)
    {
        EFG(0,w,0).max_val=EFG(0,w,2).max_val=EFG(0,w-1,2).max_val+((w==1 && head)?(*ve_b):(*ue_b));
        for(int j=0; j<S_e-S_b; ++j)
        {
            if(floor)
                EFG(*(S_b+j+1),w,1).max_val=std::numeric_limits<int>::min()/2;
            for(int s=*(S_b+j)+1; s<=*(S_b+j+1); ++s)
            {
                if(floor)
                    ori_cmp(EFG(*(S_b+j+1),w,1), EFG(s-1,w,2).max_val+*(tvf_b+j)+(*(S_b+j+1)-s)*(*(tuf_b+j)), &EFG(s-1,w,2));
                if(s!=*(S_b+j+1) || !floor)
                {
                    EFG(s,w,1).max_val=std::numeric_limits<int>::min()/2;
                    ori_cmp(EFG(s,w,1), EFG(s-1,w,2).max_val+*(vf_b+j+w*row), &EFG(s-1,w,2));
                    if(s!=*(S_b+j)+1)
                        ori_cmp(EFG(s,w,1), EFG(s-1,w,1).max_val+*(uf_b+j+w*row), &EFG(s-1,w,1));
                    if(ceil)
                        ori_cmp(EFG(s,w,1), EFG(*(S_b+j),w,2).max_val+*(tvf_b+j)+(s-*(S_b+j)-1)*(*(tuf_b+j)), &EFG(*(S_b+j),w,2));
                }
                EFG(s,w,0).max_val=std::numeric_limits<int>::min()/2;
                if(!tail && s==*S_e)
                    ori_cmp(EFG(s,w,0), EFG(s,w-1,2).max_val+*(ue_b+s), &EFG(s,w-1,2));
                else
                {
                    ori_cmp(EFG(s,w,0), EFG(s,w-1,2).max_val+*(ve_b+s), &EFG(s,w-1,2));
                    if(EFG(s,w-1,2).max_val!=EFG(s,w-1,0).max_val || *(ve_b+s)!=*(ue_b+s))
                        ori_cmp(EFG(s,w,0), EFG(s,w-1,0).max_val+*(ue_b+s), &EFG(s,w-1,0));
                }
                EFG(s,w,2).max_val=std::numeric_limits<int>::min()/2;
                ori_cmp(EFG(s,w,2), EFG(s-1,w-1,2).max_val+gamma(nt2int[*(x_b+s-1)],nt2int[*(o_b+w-1)],0), &EFG(s-1,w-1,2));
                ori_cmp(EFG(s,w,2), EFG(s,w,1).max_val, &EFG(s,w,1));
                ori_cmp(EFG(s,w,2), EFG(s,w,0).max_val, &EFG(s,w,0));
            }
        }
    }
    std::vector<std::vector<std::pair<int,int>>> rivets;
    std::vector<Back*> biters;
    rivets.push_back(std::vector<std::pair<int,int>>());
    biters.push_back(&EFG(*S_e,o_e-o_b,2));
    int i=0;
    while(i<biters.size())
    {
        rivets[i].push_back(std::make_pair(EFG.grow(biters[i]),EFG.gcol(biters[i])));
        if(rivets[i].back().first==0 || rivets[i].back().second==0)
        {
            rivets[i].push_back(std::make_pair(0,0));
            ++i;
        }
        else
        {
            std::vector<Back*>::iterator iter=biters[i]->biters.begin();
            std::vector<Back*>::iterator end_iter=biters[i]->biters.end();
            biters[i]=*iter;
            ++iter;
            while(biters.size()<ALIGN_MAX && iter!=end_iter)
            {
                rivets.push_back(rivets[i]);
                biters.push_back(*iter);
                ++iter;
            }
        }
    }
    std::vector<Align> aligns;
    std::vector<std::string> references, queries;
    for(int i=0; i<rivets.size(); ++i)
    {
        aligns.emplace_back();
        references.emplace_back();
        queries.emplace_back();
        aligns.back().max_score=EFG(*S_e,o_e-o_b,2).max_val;
        for(std::vector<std::pair<int,int>>::reverse_iterator iter=rivets[i].rbegin(); iter!=rivets[i].rend()-1; ++iter)
        {
            if(iter->first==(iter+1)->first)
            {
                references.back().append(std::string((iter+1)->second-iter->second,'-'));
                queries.back().append(o_b+iter->second,o_b+(iter+1)->second);
            }
            else if(iter->second==(iter+1)->second)
            {
                references.back().append(x_b+iter->first,x_b+(iter+1)->first);
                queries.back().append(std::string((iter+1)->first-iter->first,'-'));
            }
            else
            {
                references.back().append(x_b+iter->first,x_b+(iter+1)->first);
                queries.back().append(o_b+iter->second,o_b+(iter+1)->second);
            }
        }
    }
    return std::make_tuple(aligns, references, queries);
}

std::vector<Align> wapper_column_wise(std::string x, std::vector<int> S, std::string alg_type, int u, int v, int ALIGN_MAX, std::vector<std::string> os, std::vector<int> index, std::vector<double> num, int max_len, int S0, int S1, std::string file, int blockindex, std::vector<int> left_exp, std::vector<int> right_exp, std::string mode)
{
    std::map<char,int> nt2int;
    TD_array<int> gamma, vf, uf;
    std::vector<int> ve, ue, tvf, tuf;
    std::tie(nt2int,gamma,ve, ue, tvf, tuf, vf, uf)=initialization(S, alg_type, u, v, max_len, S0, S1);
    TD_array<Back> EFG(S.back()+1,max_len+1,3);
    std::vector<Align> aligns;
    std::ofstream fout("tmp/" + std::filesystem::path(file).filename().string() + "." + std::to_string(blockindex) + ".alg");
    for(int i=0; i<os.size(); ++i)
    {
        std::vector<Align> tmp;
        std::vector<std::string> references, queries;
        std::tie(tmp, references, queries)=column_wise(ve.begin(), ue.begin(), S.size()-1, vf.p2i(&vf(0,0,0)), uf.p2i(&uf(0,0,0)), tvf.begin(), tuf.begin(), gamma, x.begin(), os[i].begin(), os[i].end(), S.begin(), S.end()-1, true, true, true, true, ALIGN_MAX, nt2int, EFG);
        for(int j=0; j<tmp.size(); ++j)
        {
            tmp[j].index=index[i];
            tmp[j].num=num[i];
            std::tie(tmp[j].left, tmp[j].right, tmp[j].MID, tmp[j].c_left, tmp[j].c_right, tmp[j].c_MID, tmp[j].c_left_in, tmp[j].c_right_in, tmp[j].c_MID_in)=search_left_right(references[j], queries[j], S, left_exp, right_exp, mode, x);

            fout << tmp[j].index << '\t' << tmp[j].num << '\t' << tmp[j].max_score << '\t';
            for(int k=0; k<tmp[j].left.size(); ++k)
                fout << tmp[j].MID[k] << '\t' << tmp[j].right[k] << '\t' << tmp[j].left[k] << '\t';
            fout << tmp[j].MID.back() << '\t';
            for(int k=0; k<tmp[j].c_left.size(); ++k)
                fout << tmp[j].c_MID[k] << '\t' << tmp[j].c_right[k] << '\t' << tmp[j].c_left[k] << '\t';
            fout << tmp[j].c_MID.back() << '\t';
            for(int k=0; k<tmp[j].c_left_in.size(); ++k)
                fout << tmp[j].c_MID_in[k] << '\t' << tmp[j].c_right_in[k] << '\t' << tmp[j].c_left_in[k] << '\t';
            fout << tmp[j].c_MID_in.back() << '\n' << references[j] << '\n' << queries[j] << '\n';
        }
        std::move(tmp.begin(), tmp.end(), std::back_inserter(aligns));
    }
    fout.close();
    

    return aligns;
}

#endif
