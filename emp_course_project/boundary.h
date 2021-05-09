// boundary.h

#pragma once

#include "node.h"
#include "dense_matrix.h"
#include "sparse_matrix.h"

// Учесть краевое условие 1-го рода
void first_kind_boundary_cond(
    SparseMatrix& A, 
    std::vector<double>& f_global, 
    const int& node1, 
    const double& f_val1,
    const int& node2, 
    const double& f_val2
);

// Учесть краевое условие 2-го рода
void second_kind_boundary_cond(
    std::vector<double>& b_global,
    const std::vector<int>& cur_cond,
    const int& cond_num,
    FiniteElement& el,
    const int& el_num,
    const std::vector<Node>& grid,
    const double& t,
    double(*bc2_func)(const double& x, const double& y, const double& t)
);

// Учесть краевое условие 3-го рода
void third_kind_boundary_cond(
    SparseMatrix& A,
    std::vector<double>& b_global,
    const std::vector<int>& cur_cond,
    const double& beta,
    const int& cond_num,
    FiniteElement& el,
    const int& el_num,
    const std::vector<Node>& grid,
    const double& t,
    double(*bc3_func)(const double& x, const double& y, const double& t),
    bool ignore_matrix = false
);
