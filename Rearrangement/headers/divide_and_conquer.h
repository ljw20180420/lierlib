#ifndef DIVIDE_AND_CONQUER_H
#define DIVIDE_AND_CONQUER_H

#include "TD_array.h"
#include <map>
#include <string>
#include <filesystem>

template<typename vi, typename si>
void row_wise_linear(vi ve_b, vi ue_b, int row, vi vf_b, vi uf_b, vi tvf_b, vi tuf_b, TD_array<int>& gamma, si x_b, si o_b, si o_e, vi S_b, vi S_e, std::map<char,int>& nt2int, vi G_b)
{
    std::vector<int> tildeG(o_e-o_b+1);
    *G_b=0;
    for(int w=1; w<=o_e-o_b; ++w)
        tildeG[w]=*(G_b+w)=*(G_b+w-1)+((w==1)?(*ve_b):(*ue_b));
    for(int j=0; j<S_e-S_b; ++j)
    {
        std::vector<int> F(o_e-o_b+1);
        std::vector<int> tildeF(o_e-o_b+1,std::numeric_limits<int>::min()/2);
        for(int s=*(S_b+j)+1; s<=*(S_b+j+1); ++s)
        {
            int hatG=*G_b;
            *G_b=*G_b+((s==*(S_b+j)+1)?(*(tvf_b+j)):(*(tuf_b+j)));
            int E=std::numeric_limits<int>::min()/2;
            for(int w=1; w<=o_e-o_b; ++w)
            {
                tildeF[w]=std::max(tildeF[w],*(G_b+w)+*(tvf_b+j)+(*(S_b+j+1)-s)*(*(tuf_b+j)));
                E=std::max(*(G_b+w-1)+*(ve_b+s),E+*(ue_b+s));
                int Gp=*(G_b+w);
                *(G_b+w)=hatG+gamma(nt2int[*(x_b+s-1)],nt2int[*(o_b+w-1)],0);
                *(G_b+w)=std::max(*(G_b+w),E);
                if(s!=*(S_b+j+1))
                {
                    if(s!=*(S_b+j)+1)
                        F[w]=std::max(F[w]+*(uf_b+j+w*row),Gp+*(vf_b+j+w*row));
                    else
                        F[w]=Gp+*(vf_b+j+w*row);
                    F[w]=std::max(F[w],tildeG[w]+*(tvf_b+j)+(s-*(S_b+j)-1)*(*(tuf_b+j)));
                    *(G_b+w)=std::max(*(G_b+w),F[w]);
                }
                else
                {
                    *(G_b+w)=std::max(*(G_b+w),tildeF[w]);
                    tildeG[w]=*(G_b+w);
                }
                hatG=Gp;
            }
        }
    }
}

template<typename vi, typename si>
void column_wise_linear(vi ve_b, vi ue_b, int row, vi vf_b, vi uf_b, vi tvf_b, vi tuf_b, TD_array<int>& gamma, si x_b, si o_b, si o_e, vi S_b, vi S_e, bool head, bool ceil, bool floor, std::map<char,int>& nt2int, vi E_b, vi G_b)
{
    int F;
    *(G_b)=0;
    for(int j=0; j<S_e-S_b; ++j)
        for(int s=*(S_b+j)+1; s<=*(S_b+j+1); ++s)
        {
            if(ceil)
                *(G_b+s)=*(G_b+s-1)+((s==*(S_b+j)+1)?(*(tvf_b+j)):(*(tuf_b+j)));
            else
                *(G_b+s)=*(G_b+s-1)+((s==*(S_b+j)+1)?(*(vf_b+j)):(*(uf_b+j)));
            *(E_b+s)=std::numeric_limits<int>::min()/2;
        }
    for(int w=1; w<=o_e-o_b; ++w)
    {
        int hatG=*G_b;
        *E_b=*G_b=*G_b+((w==1 && head)?(*ve_b):(*ue_b));
        for(int j=0; j<S_e-S_b; ++j)
            for(int s=*(S_b+j)+1; s<=*(S_b+j+1); ++s)
            {
                if(floor && s==*(S_b+j+1))
                {
                    F=std::numeric_limits<int>::min()/2;
                    for(int s=*(S_b+j)+1; s<=*(S_b+j+1); ++s)
                        F=std::max(F,*(G_b+s-1)+*(tvf_b+j)+(*(S_b+j+1)-s)*(*(tuf_b+j)));
                }
                else
                {
                    if(s!=*(S_b+j)+1)
                        F=std::max(F+*(uf_b+j+w*row),*(G_b+s-1)+*(vf_b+j+w*row));
                    else
                        F=*(G_b+s-1)+*(vf_b+j+w*row);
                    if(ceil)
                        F=std::max(F,*(G_b+*(S_b+j))+*(tvf_b+j)+(s-*(S_b+j)-1)*(*(tuf_b+j)));
                }
                *(E_b+s)=std::max(*(G_b+s)+*(ve_b+s),*(E_b+s)+*(ue_b+s));
                int Gp=*(G_b+s);
                *(G_b+s)=hatG+gamma(nt2int[*(x_b+s-1)],nt2int[*(o_b+w-1)],0);
                *(G_b+s)=std::max(*(G_b+s),*(E_b+s));
                *(G_b+s)=std::max(*(G_b+s),F);
                hatG=Gp;
            }
    }
}


void divide_and_conquer(std::vector<int>::iterator ve_b, std::vector<int>::iterator ue_b, int row, std::vector<int>::iterator vf_b, std::vector<int>::iterator uf_b, std::vector<int>::iterator tvf_b, std::vector<int>::iterator tuf_b, TD_array<int>& gamma, std::string::iterator x_b, std::string::iterator o_b, std::string::iterator o_e, std::vector<int>::iterator S_b, std::vector<int>::iterator S_e, bool head, bool tail, bool ceil, bool floor, std::map<char,int>& nt2int, Align& align, std::string& reference, std::string& query, std::vector<int>::iterator GU_b, std::vector<int>::reverse_iterator GD_b, std::vector<int>::iterator GL_b, std::vector<int>::reverse_iterator GR_b, std::vector<int>::iterator EL_b, std::vector<int>::reverse_iterator ER_b)
{
    if(o_b==o_e)
    {
        reference.append(x_b,x_b+*S_e);
        query.append(std::string(*S_e,'-'));
        for(int i=0; i<S_e-S_b; ++i)
        {
            if(ceil || floor)
                align.max_score+=(*(S_b+i+1)-1)*(*(tuf_b+i))+*(tvf_b+i);
            else
                align.max_score+=(*(S_b+i+1)-1)*(*(uf_b+i))+*(vf_b+i);
        }
        return;
    }
    if(*S_e==0)
    {
        reference.append(std::string(o_e-o_b,'-'));
        query.append(o_b,o_e);
        align.max_score+=(o_e-o_b)*(*ue_b);
        if(head || tail)
            align.max_score+=*ve_b-*ue_b;
        return;
    }
    if(o_e-o_b==1)
    {
        TD_array<Back> EFG(*S_e+1,o_e-o_b+1,3);
        std::vector<Align> aligns;
        std::vector<std::string> references, queries;
        std::tie(aligns,references,queries)=column_wise(ve_b, ue_b, row, vf_b, uf_b, tvf_b, tuf_b, gamma, x_b, o_b, o_e, S_b, S_e, head, tail, ceil, floor, 1, nt2int, EFG);
        reference.append(references.front());
        query.append(queries.front());
        align.max_score+=aligns.front().max_score;
        return;
    }
    if(S_e-S_b>1)
    {
        int js=(S_e-S_b)/2;
        row_wise_linear(ve_b, ue_b, row, vf_b, uf_b, tvf_b, tuf_b, gamma, x_b, o_b, o_e, S_b, S_b+js, nt2int, GU_b);
        std::vector<int> S_tmp(S_b+js,S_e+1);
        for(std::vector<int>::iterator iter=S_tmp.begin(); iter!=S_tmp.end(); ++iter)
            *iter=*S_e-*iter;
        row_wise_linear(std::vector<int>::reverse_iterator(ve_b+*S_e+1), std::vector<int>::reverse_iterator(ue_b+*S_e+1), row, std::vector<int>::reverse_iterator(vf_b+(S_e-S_b)+(o_e-o_b)*row), std::vector<int>::reverse_iterator(uf_b+(S_e-S_b)+(o_e-o_b)*row), std::vector<int>::reverse_iterator(tvf_b+(S_e-S_b)), std::vector<int>::reverse_iterator(tuf_b+(S_e-S_b)), gamma, std::string::reverse_iterator(x_b+*S_e), std::string::reverse_iterator(o_e), std::string::reverse_iterator(o_b), S_tmp.rbegin(), S_tmp.rend()-1, nt2int, GD_b);
        int max_tmp=std::numeric_limits<int>::min();
        int ws;
        for(int w=0; w<=o_e-o_b; ++w)
        {
            int tmp=*(GU_b+w)+*(GD_b+(o_e-o_b)-w);
            if(tmp>max_tmp)
            {
                max_tmp=tmp;
                ws=w;
            }
        }
        divide_and_conquer(ve_b, ue_b, row, vf_b, uf_b, tvf_b, tuf_b, gamma, x_b, o_b, o_b+ws, S_b, S_b+js, true, true, true, true, nt2int, align, reference, query, GU_b, GD_b, GL_b, GR_b, EL_b, ER_b);
        for(int i=0; i<S_tmp.size(); ++i)
            S_tmp[i]=*(S_b+js+i)-*(S_b+js);
        divide_and_conquer(ve_b+*(S_b+js), ue_b+*(S_b+js), row, vf_b+js+ws*row, uf_b+js+ws*row, tvf_b+js, tuf_b+js, gamma, x_b+*(S_b+js), o_b+ws, o_e, S_tmp.begin(), S_tmp.end()-1, true, true, true, true, nt2int, align, reference, query, GU_b, GD_b, GL_b, GR_b, EL_b, ER_b);
    }
    else
    {
        int ws=(o_e-o_b)/2;
        column_wise_linear(ve_b, ue_b, row, vf_b, uf_b, tvf_b, tuf_b, gamma, x_b, o_b, o_b+ws, S_b, S_e, head, ceil, floor, nt2int, EL_b, GL_b);
        std::vector<int> S_tmp{*S_e,*S_b};
        column_wise_linear(std::vector<int>::reverse_iterator(ve_b+*S_e+1), std::vector<int>::reverse_iterator(ue_b+*S_e+1), row, std::vector<int>::reverse_iterator(vf_b+(S_e-S_b)+(o_e-o_b)*row), std::vector<int>::reverse_iterator(uf_b+(S_e-S_b)+(o_e-o_b)*row), std::vector<int>::reverse_iterator(tvf_b+(S_e-S_b)), std::vector<int>::reverse_iterator(tuf_b+(S_e-S_b)), gamma, std::string::reverse_iterator(x_b+*S_e), std::string::reverse_iterator(o_e), std::string::reverse_iterator(o_b+ws), S_tmp.rbegin(), S_tmp.rend()-1, tail, floor, ceil, nt2int, ER_b, GR_b);
        int max_tmp_e=std::numeric_limits<int>::min();
        int max_tmp_g=std::numeric_limits<int>::min();
        int se, sg;
        for(int s=0; s<=*S_e; ++s)
        {
            int tmp_e=*(EL_b+s)+*(ER_b+*S_e-s)-*(ve_b+s)+*(ue_b+s);
            int tmp_g=*(GL_b+s)+*(GR_b+*S_e-s);
            if(tmp_e>max_tmp_e)
            {
                max_tmp_e=tmp_e;
                se=s;
            }
            if(tmp_g>max_tmp_g)
            {
                max_tmp_g=tmp_g;
                sg=s;
            }
        }
        int ss, deltaws;
        bool bool_tmp;
        if(max_tmp_e>max_tmp_g)
        {
            ss=se;
            deltaws=1;
            bool_tmp=false;
        }
        else
        {
            ss=sg;
            deltaws=0;
            bool_tmp=true;
        }
        bool bool2_tmp;
        if(ss==0)
            bool2_tmp=ceil;
        else if(ss==*S_e)
            bool2_tmp=floor;
        else
            bool2_tmp=false;
        S_tmp[0]=0;
        S_tmp[1]=ss;
        divide_and_conquer(ve_b, ue_b, row, vf_b, uf_b, tvf_b, tuf_b, gamma, x_b, o_b, o_b+ws-deltaws, S_tmp.begin(), S_tmp.end()-1, head, bool_tmp, ceil, bool2_tmp, nt2int, align, reference, query, GU_b, GD_b, GL_b, GR_b, EL_b, ER_b);
        reference.append(std::string(2*deltaws,'-'));
        query.append(o_b+ws-deltaws,o_b+ws+deltaws);
        if(deltaws>0)
            align.max_score+=(2*deltaws-1)*(*(ue_b+ss))+*(ve_b+ss);
        S_tmp[0]=0;
        S_tmp[1]=*S_e-ss;
        divide_and_conquer(ve_b+ss, ue_b+ss, row, vf_b+ws+deltaws, uf_b+ws+deltaws, tvf_b, tuf_b, gamma, x_b+ss, o_b+ws+deltaws, o_e, S_tmp.begin(), S_tmp.end()-1, bool_tmp, tail, bool2_tmp, floor, nt2int, align, reference, query, GU_b, GD_b, GL_b, GR_b, EL_b, ER_b);
    }
}


std::vector<Align> wapper_divide_and_conquer(std::string x, std::vector<int> S, std::string alg_type, int u, int v, std::vector<std::string> os, std::vector<int> index, std::vector<double> num, int max_len, int S0, int S1, std::string file, int blockindex, std::vector<int> left_exp, std::vector<int> right_exp, std::string mode)
{
    std::map<char,int> nt2int;
    TD_array<int> gamma, vf, uf;
    std::vector<int> ve, ue, tvf, tuf;
    std::tie(nt2int,gamma,ve, ue, tvf, tuf, vf, uf)=initialization(S, alg_type, u, v, max_len, S0, S1);
    std::vector<int> GU(max_len+1), GD(max_len+1), GL(S.back()+1), GR(S.back()+1), EL(S.back()+1), ER(S.back()+1);
    std::vector<Align> aligns;
    std::ofstream fout("tmp/" + std::filesystem::path(file).filename().string() + "." + std::to_string(blockindex) + ".alg");
    for(int i=0; i<os.size(); ++i)
    {
        aligns.emplace_back();
        std::string reference, query;
        divide_and_conquer(ve.begin(), ue.begin(), S.size()-1, vf.p2i(&vf(0,0,0)), uf.p2i(&uf(0,0,0)), tvf.begin(), tuf.begin(), gamma, x.begin(), os[i].begin(), os[i].end(), S.begin(), S.end()-1, true, true, true, true, nt2int, aligns.back(), reference, query, GU.begin(), GD.rbegin(), GL.begin(), GR.rbegin(), EL.begin(), ER.rbegin());
        aligns.back().index=index[i];
        aligns.back().num=num[i];
        std::tie(aligns.back().left, aligns.back().right, aligns.back().MID, aligns.back().c_left, aligns.back().c_right, aligns.back().c_MID, aligns.back().c_left_in, aligns.back().c_right_in, aligns.back().c_MID_in)=search_left_right(reference, query, S, left_exp, right_exp, mode, x);
        fout << aligns.back().index << '\t' << aligns.back().num << '\t' << aligns.back().max_score << '\t';
        for(int j=0; j<aligns.back().left.size(); ++j)
            fout << aligns.back().MID[j] << '\t' << aligns.back().right[j] << '\t' << aligns.back().left[j] << '\t';
        fout << aligns.back().MID.back() << '\t';
        for(int j=0; j<aligns.back().c_left.size(); ++j)
            fout << aligns.back().c_MID[j] << '\t' << aligns.back().c_right[j] << '\t' << aligns.back().c_left[j] << '\t';
        fout << aligns.back().c_MID.back() << '\t';
        for(int j=0; j<aligns.back().c_left_in.size(); ++j)
            fout << aligns.back().c_MID_in[j] << '\t' << aligns.back().c_right_in[j] << '\t' << aligns.back().c_left_in[j] << '\t';
        fout << aligns.back().c_MID_in.back() << '\n' << reference << '\n' << query << '\n';
    }
    fout.close();
    
    return aligns;
}

#endif