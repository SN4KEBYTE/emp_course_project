// bi_cgm.h

#pragma once

#include <vector>

#include "utils.h"
#include "sparse_matrix.h"

// Максимальное число итераций по умолчанию
const int DEFAULT_MAX_ITER = 1e5;

// Относительная невязка по умолчанию
const double DEFAULT_EPS = 1e-13;

/*
Класс, описывающий решатель методом бисопряженных градиентов
*/

class BiCGM {
private:
    // Максимальное число итераций
    int max_iter = DEFAULT_MAX_ITER;

    // Величина относительной невязки
    double eps = DEFAULT_EPS;

    // Вспомогательные вектора
    std::vector<double> r, p, z, s, az, ats;

public:
    // Конструктор по умолчанию - использовать DEFAULT значения максимального числа итераций и невязки
    BiCGM() {}

    // Конструктор с параметрами - пользователь сам может задать максимальное число итераций и невязку
    BiCGM(
        const int& max_iter, 
        const double& eps
    ) : max_iter(max_iter), eps(eps) {}

    // Решить систему вида Ax = f, результат окажется в x
    void solve(
        SparseMatrix& A, 
        std::vector<double>& xk, 
        const std::vector<double>& f
    );
};