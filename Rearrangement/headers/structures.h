#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <string>
#include <vector>

struct Back
{
    int max_val;
    std::vector<Back*> biters;
};

struct Align
{
    int index, max_score;
    double num;
    std::vector<int> left, right, c_left, c_right, c_left_in, c_right_in;
    std::vector<std::string> MID, c_MID, c_MID_in;
};

struct JuncIndel
{
    double num;
    std::vector<std::string> ligation_indel;
    std::vector<std::vector<int>> insert; 
};

#endif