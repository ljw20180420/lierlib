#ifndef TD_ARRAY_H
#define TD_ARRAY_H

#include <functional>

template<typename T>
class TD_array
{
    std::vector<T> data;
    int row;
    int col;
    int rowcol;
    int slice;
public:
    TD_array(int row_=0, int col_=0, int slice_=0)
    {
        row=row_;
        col=col_;
        rowcol=row*col;
        slice=slice_;
        data.resize(row*col*slice);
    }
    TD_array(int row_, int col_, int slice_, int val_)
    {
        row=row_;
        col=col_;
        rowcol=row*col;
        slice=slice_;
        data.resize(row*col*slice, val_);
    }
    void resize(int row_, int col_, int slice_)
    {
        row=row_;
        col=col_;
        rowcol=row*col;
        slice=slice_;
        data.resize(row*col*slice);
    }
    void resize(int row_, int col_, int slice_, int val_)
    {
        row=row_;
        col=col_;
        rowcol=row*col;
        slice=slice_;
        data.resize(row*col*slice, val_);
    }
    T& operator()(int i, int j, int k)
    {
        return data[rowcol*k+row*j+i];
    }
    void fill(T val)
    {
        std::fill(data.begin(),data.end(),val);
    }
    double accumulate()
    {
        return accumulate(data.begin(), data.end(), 0.0);
    }
    int gsli(T* ptr)
    {
        return (ptr-data.data())/rowcol;
    }
    int gcol(T* ptr)
    {
        return ((ptr-data.data())%rowcol)/row;
    }
    int grow(T* ptr)
    {
        return (ptr-data.data())%row;
    }
    typename std::vector<T>::iterator p2i(T* ptr)
    {
        return data.begin()+(ptr-data.data());
    }
    int size(int dim)
    {
        if(dim==0)
            return row;
        else if(dim==1)
            return col;
        else if(dim==2)
            return slice;
        else
            return -1;
    }
    std::tuple<std::vector<T>,std::vector<T>,std::vector<T>> boundary_sum()
    {
        std::vector<T> alpha(row,0), beta(col,0), pi(slice,0);
        std::vector<T> tmp(col*slice,0);
        int index=0;
        int tind=0;
        for(int s=0;s<slice;++s)
            for(int c=0;c<col;++c)
            {
                for(int r=0;r<row;++r)
                {
                    alpha[r]+=data[index];
                    tmp[tind]+=data[index];
                    ++index;
                }
                beta[c]+=tmp[tind];
                pi[s]+=tmp[tind];
                ++tind;
            }
        return std::make_tuple(alpha, beta, pi);
    }
};

#endif