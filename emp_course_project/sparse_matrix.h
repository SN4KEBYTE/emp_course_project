// sparse_matrix.h

#pragma once

#include <vector>
#include "utils.h"
#include "dense_matrix.h"
#include "finite_element.h"

/*
Класс, описывающий разреженную матрицу, хранящуюся в строчно-столбцовом формате
*/
class SparseMatrix
{
private:
    // Размерность матрицы
    int n = 0;
    
    // ig - указатели на начала строк, jg - номера столбцов внедиагональных элементов
    std::vector<int> ig, jg;

public:
    // Диагональные элементы
    // ggu и ggl - внедиагональные элементы верхнего и нижнего треугольников соответственно
    std::vector<double> diag, ggu, ggl;
    
    // Конструкторы
    SparseMatrix(
        int n
    );

    SparseMatrix(
        int n,
        std::vector<int> ig,
        std::vector<int> jg,
        std::vector<double> ggu,
        std::vector<double> ggl,
        std::vector<double> diag
    );
    
    // Умножение матрицы на вектор (Ax = y)
    void dot_vector(
        const std::vector<double>& x, 
        std::vector<double>& y
    );

    SparseMatrix dot_scalar(
        const double& x
    );
    
    // Умножение транспонированной матрицы на вектор (A^(T)x = y)
    void T_dot_vector(
        const std::vector<double>& x, 
        std::vector<double>& y
    );
    
    // Генерация портрета конечноэлементной СЛАУ
    void build_portrait(
        const int& num_elems, 
        std::vector<FiniteElement>& elems
    );
    
    // Внесение локальной матрицы в глобальную
    void add_local(
        SquareDenseMatrix& mat, 
        const FiniteElement& el
    );
    
    // Обнулить строку
    void zero_out_row(
        const int& row_num
    );
    
    // In-place сложение матриц
    SparseMatrix& operator+=(
        const SparseMatrix& m
    );

    // Получить размерность матрицы
    int get_n();
};
