// dense_matrix.cpp

#include "dense_matrix.h"

int SquareDenseMatrix::get_n()
{
    return n;
}

SquareDenseMatrix& SquareDenseMatrix::operator+=(
    const SquareDenseMatrix& op
)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            this->data[i][j] += op.data[i][j];
    
    return *this;
}

SquareDenseMatrix & SquareDenseMatrix::operator*=(
    const double& coef
)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            this->data[i][j] *= coef;
    
    return *this;
}

SquareDenseMatrix& SquareDenseMatrix::operator/=(
    const double& coef
)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            this->data[i][j] /= coef;
    
    return *this;
}

std::vector<double>& SquareDenseMatrix::operator[](
    int i
)
{
    return data[i];
}

void SquareDenseMatrix::dot_vector(
    const std::vector<double>& x,
    std::vector<double>& y
)
{
    for (int i = 0; i < n; i++)
    {
        y[i] = 0;
        
        for (int j = 0; j < n; j++)
            y[i] += data[i][j] * x[j];
    }
}