#include "headers/EM.h"
#include "headers/thread_pools.h"
#include "headers/output.h"
#include "headers/loader.h"
#include "headers/parser.h"
#include "headers/align.h"
#include "headers/divide_and_conquer.h"
#include <limits>

int main(int argc, char **argv)
{
    print_help(argc, argv);
    
    Command_content cc=command(argc, argv);
    
    std::string x;
    std::vector<int> S, left_exp, left_down, left_up, right_exp, right_down, right_up;
    load_ref(cc.ref_file, x, S, left_exp, left_down, left_up, right_exp, right_down, right_up);
    if(x.empty())
    {
        std::cerr << "empty reference loaded\n";
        return -1;
    }

    std::vector<std::pair<std::string,int>> sus;
    int total_read = 0, max_len = std::numeric_limits<int>::min();
    {
        std::ifstream fin(cc.file);
        std::string read_tmp;
        int num_tmp;
        while(fin >> read_tmp >> num_tmp)
        {
            sus.emplace_back(read_tmp, num_tmp);
            total_read += num_tmp;
            max_len=std::max(max_len,int(read_tmp.size()));
        }
        fin.close();
    }
    for (int i=0; i<sus.size(); ++i)
    {
        if (double(sus[i].second)/total_read < cc.per_thres)
        {
            sus.resize(i);
            break;
        }
    }

    std::filesystem::create_directory("tmp");
    thread_pool threads(cc.THR_NUM);
    std::vector<std::future<std::vector<Align>>> futures;
    int blockindex = 0;
    for(int so=0; so<sus.size(); so+=cc.SEQ_BATCH, ++blockindex)
    {
        std::vector<std::string> os;
        std::vector<double> num;
        std::vector<int> index;
        for (int i=so; i<so+cc.SEQ_BATCH && i<sus.size();)
        {
            os.push_back(sus[i].first);
            num.push_back(sus[i].second);
            index.push_back(++i);
        }
        if(cc.ALIGN_MAX>1 || !cc.DIVCON)
            futures.push_back(threads.submit(std::bind(wapper_column_wise, x, S, cc.alg_type, cc.u, cc.v, cc.ALIGN_MAX, std::move(os), std::move(index), std::move(num), max_len, cc.s0, cc.s1, cc.file, blockindex, left_exp, right_exp, cc.mode)));
        else
            futures.push_back(threads.submit(std::bind(wapper_divide_and_conquer, x, S, cc.alg_type, cc.u, cc.v, std::move(os), std::move(index), std::move(num), max_len, cc.s0, cc.s1, cc.file, blockindex, left_exp, right_exp, cc.mode)));
    }
    std::vector<Align> aligns;
    for(unsigned i=0; i<futures.size(); i++)
    {
        std::vector<Align> tmp=futures[i].get();
        std::move(tmp.begin(), tmp.end(), std::back_inserter(aligns));
    }
    if (cc.indel)
        indel_label_fun(aligns, left_exp, right_exp, cc.file);
    if (cc.EM)
        EM_predict(aligns, x, left_exp, left_down, left_up, right_exp, right_down, right_up, cc.MID_MAX, cc.ini_alpha, cc.ini_beta, cc.ini_pi, cc.thres, cc.file);

    return 0;
}

































