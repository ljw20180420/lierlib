#ifndef EM_H
#define EM_H

#include "TD_array.h"
#include "tools.h"
#include <fstream>
#include <map>
#include <filesystem>

std::vector<std::vector<double>> initial_EM(std::string& ini, int dim, std::vector<TD_array<double>>& N_mats)
{
    std::vector<std::vector<double>> vecs;
    if(ini.empty())
        for(auto & N_mat : N_mats)
        {
            int tmp=N_mat.size(dim);
            vecs.push_back(std::vector<double>(tmp,1.0/tmp));
        }
    else
    {
        std::ifstream fin(ini);
        for(auto & N_mat : N_mats)
        {
            int tmp=N_mat.size(dim);
            vecs.push_back(std::vector<double>());
            std::string line;
            for(int rein=0; rein<tmp; ++rein)
            {
                getline(fin,line);
                vecs.back().push_back(std::stod(my_split(line," \t\n\v\f\r")[1]));
            }
        }
        fin.close();
    }
    return vecs;
}


void EM_predict(std::vector<Align>& aligns, std::string& x, std::vector<int>& left_exp, std::vector<int>& left_down, std::vector<int>& left_up, std::vector<int>& right_exp, std::vector<int>& right_down, std::vector<int>& right_up, int MID_MAX, std::string& ini_alpha, std::string& ini_beta, std::string& ini_pi, double thres, std::string& file)
{
    div_num_fun(aligns.begin(), aligns.end());
    
    std::vector<std::string> pos_2_MID;
        std::string ind_tmp("NACGT");
        for(int k=0;k<=MID_MAX;++k)
        {
            for(int w=0;w<std::pow(5,k);++w)
            {
                std::string MID;
                for(int i=0;i<k;i++)
                {
                    MID+=ind_tmp[w%5];
                    w/=5;
                }
                pos_2_MID.push_back(std::move(MID));
            }
        }
    
    std::vector<std::map<std::tuple<int,int,std::string>,std::pair<double,std::vector<double*>>,MY_LESS>> tuple_2_numpoi;
    std::vector<TD_array<double>> N_mats;
    for(int i=1; i<right_down.size(); ++i)
    {
        int n_row=left_up[i-1]-left_down[i-1]+1, n_col=right_up[i]-right_down[i]+1, n_slice=int((1-pow(5,MID_MAX+1))/(1-5)+0.5);
        N_mats.push_back(TD_array<double>(n_row, n_col, n_slice));
        tuple_2_numpoi.push_back(std::map<std::tuple<int,int,std::string>,std::pair<double,std::vector<double*>>,MY_LESS>(MY_LESS(x)));
        for(int r=0; r<n_row; ++r)
            for(int c=0; c<n_col; ++c)
                for(int s=0; s<n_slice; ++s)
                {
                    std::tuple<int,int,std::string> tmp=std::make_tuple(r+left_down[i-1],c+right_down[i],pos_2_MID[s]);
                    tuple_2_numpoi.back()[tmp].first=0;
                    tuple_2_numpoi.back()[tmp].second.push_back(&N_mats.back()(r,c,s));
                }
        
        for(auto& align : aligns)
        {
            std::tuple<int,int,std::string> tmp=std::make_tuple(align.left[i-1],align.right[i],align.MID[i]);
            if(tuple_2_numpoi.back().count(tmp)>0)
                tuple_2_numpoi.back()[tmp].first+=align.num;
        }
        
    }
    std::string pre_str;
    std::filesystem::path pathObj;
    pathObj=ini_alpha;
    pre_str.append(pathObj.filename().string());
    pathObj=ini_beta;
    if(!pre_str.empty() && !pathObj.filename().string().empty())
        pre_str.append("_");
    pre_str.append(pathObj.filename().string());
    pathObj=ini_pi;
    if(!pre_str.empty() && !pathObj.filename().string().empty())
        pre_str.append("_");
    pre_str.append(pathObj.filename().string());
    if(!pre_str.empty())
        pre_str.append("_");
    std::vector<std::vector<double>> alphas=initial_EM(ini_alpha, 0, N_mats);
    std::vector<std::vector<double>> betas=initial_EM(ini_beta, 1, N_mats);
    std::vector<std::vector<double>> pis=initial_EM(ini_pi, 2, N_mats);
    for(int i=0; i<left_down.size()-1; ++i)
    {
        double sum_tmp=0;
        for(auto & pair : tuple_2_numpoi[i])
            sum_tmp+=pair.second.first;
        if(sum_tmp==0)
        {
            std::fill(alphas[i].begin(),alphas[i].end(),0.0);
            std::fill(betas[i].begin(),betas[i].end(),0.0);
            std::fill(pis[i].begin(),pis[i].end(),0.0);
            N_mats[i].fill(0.0);
        }
        else
        {
            double tol=std::numeric_limits<double>::max();
            while(tol>thres)
            {
                // E-step
                for(int r=0;r<alphas[i].size();++r)
                    for(int c=0;c<betas[i].size();++c)
                        for(int s=0;s<pis[i].size();++s)
                            N_mats[i](r,c,s)=alphas[i][r]*betas[i][c]*pis[i][s];
                for(auto & pair : tuple_2_numpoi[i])
                {
                    double accum=0;
                    for(auto & add : pair.second.second)
                        accum+=*add;
                    if(accum>0)
                        for(auto & add : pair.second.second)
                            *add=pair.second.first*(*add)/accum;
                    else
                        for(auto & add : pair.second.second)
                            *add=pair.second.first/pair.second.second.size();
                }
                // M-step
                std::vector<double> alpha_old=alphas[i], beta_old=betas[i], pi_old=pis[i];
                std::tie(alphas[i],betas[i],pis[i])=N_mats[i].boundary_sum();
                for(auto & comp : alphas[i]) comp/=sum_tmp;
                for(auto & comp : betas[i]) comp/=sum_tmp;
                for(auto & comp : pis[i]) comp/=sum_tmp;
                
                tol=0;
                for(int r=0; r<alpha_old.size(); ++r)
                {
                    double tmp=std::abs(alpha_old[r]-alphas[i][r]);
                    if(tmp>tol) tol=tmp;
                }
                for(int c=0; c<beta_old.size(); ++c)
                {
                    double tmp=std::abs(beta_old[c]-betas[i][c]);
                    if(tmp>tol) tol=tmp;
                }
                for(int s; s<pi_old.size(); ++s)
                {
                    double tmp=std::abs(pi_old[s]-pis[i][s]);
                    if(tmp>tol) tol=tmp;
                }
            }
        }
        for(int r=0;r<alphas[i].size();++r)
            for(int c=0;c<betas[i].size();++c)
                for(int s=0;s<pis[i].size();++s)
                    N_mats[i](r,c,s)=alphas[i][r]*betas[i][c]*pis[i][s];
        std::ofstream fout(pre_str+file+std::to_string(i+1)+".box");
        for(auto & pair : tuple_2_numpoi[i])
        {
            double accum=0;
            for(auto & add : pair.second.second)
                accum+=*add;
            fout << pair.second.first << '\t' << accum*sum_tmp << '\t' << pair.second.second.size() << '\n';
        }
        fout.close();
    }
    std::ofstream fout;
    fout.open(pre_str+file+".alpha");
    for(int i=0; i<left_exp.size()-1; ++i)
        for(int j=left_down[i]; j<=left_up[i]; ++j)
            fout << j-left_exp[i] << '\t' << alphas[i][j-left_down[i]] << '\n';
    fout.close();
    fout.open(pre_str+file+".beta");
    for(int i=1; i<right_exp.size(); ++i)
        for(int j=right_down[i]; j<=right_up[i]; ++j)
            fout << j-right_exp[i] << '\t' << betas[i-1][j-right_down[i]] << '\n';
    fout.close();
    fout.open(pre_str+file+".pi");
    for(int i=0; i<left_exp.size()-1; ++i)
        for(int j=0; j<pos_2_MID.size(); ++j)
            fout << pos_2_MID[j] << '\t' << pis[i][j] << '\n';
    fout.close();
}

#endif