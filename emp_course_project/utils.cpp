// utils.cpp

#include "utils.h"

double scalar_product(
    const std::vector<double> &x, 
    const std::vector<double> &y
)
{
    double res = 0;
    
    for (int i = 0; i < x.size(); i++)
        res += x[i] * y[i];
    
    return res;
}

double vec_norm(
    const std::vector<double>& x
)
{
    return sqrt(scalar_product(x, x));
}

void vec_diff(
    std::vector<double>& res, 
    const std::vector<double>& x,
    const std::vector<double>& y
)
{
    for (int i = 0; i < res.size(); i++)
        res[i] = x[i] - y[i];
}

void vec_sum(
    std::vector<double>& res,
    const std::vector<double>& x,
    const std::vector<double>& y
)
{
    for (int i = 0; i < res.size(); i++)
        res[i] = x[i] + y[i];
}

void vec_mult_scalar(
    const double& a,
    std::vector<double>& x
)
{
    for (int i = 0; i < x.size(); i++)
        x[i] *= a;
}