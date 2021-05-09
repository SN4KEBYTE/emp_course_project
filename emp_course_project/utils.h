// utils.h

#pragma once

#include <vector>
#include <fstream>

// Скалярное произведение векторов
double scalar_product(
    const std::vector<double>& x, 
    const std::vector<double>& y
);

// Евклидова норма вектора
double vec_norm(
    const std::vector<double>& x
);

// Разность векторов (res = x - y)
void vec_diff(
    std::vector<double>& res, 
    const std::vector<double>& x,
    const std::vector<double>& y
);

// Сумма векторов (res = x + y)
void vec_sum(
    std::vector<double>& res,
    const std::vector<double>& x,
    const std::vector<double>& y
);

// Умножить вектора на скаляр in-place
void vec_mult_scalar(
    const double& a,
    std::vector<double>& x
);

// Освобождение памяти, выделенной под вектор
template <typename T>
void vec_free(
    std::vector<T>& x
)
{
    std::vector<T>().swap(x);
}

// Вывод вектора
template <typename T>
void dump_vec(
    std::ostream& out, 
    const std::vector<T>& x, 
    std::string sep = " "
)
{
    for (const auto& el : x)
        out << el << sep;
}

// Освобождение памяти, выделенной под вектор векторов
template <typename T>
void vec_free(
    std::vector<std::vector<T>>& x
)
{
    for (auto& vec : x)
        vec_free(vec);
    
    std::vector<std::vector<T>>().swap(x);
}
