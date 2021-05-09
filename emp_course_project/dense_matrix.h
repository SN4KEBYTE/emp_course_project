// dense_matrix.h

#pragma once

#include <vector>

/*
Класс, описывающий плотную квадратную матрицу
*/
class SquareDenseMatrix
{
private:
    // Размерность
    int n = 0;
public:
    // Элементы
    std::vector<std::vector<double>> data;
    
    // Конструктор
    SquareDenseMatrix(
        int n
    ) : n(n)
    {
        data.resize(n);
        
        for (auto& row : data)
            row.resize(n);
    }
    
    // Getter для размерности
    int get_n();
    
    // Перегрузка in-place сложения
    SquareDenseMatrix& operator += (
        const SquareDenseMatrix& op
    );

    // Перегрзука in-place умножения на скаляр
    SquareDenseMatrix& operator *= (
        const double& coef
    );
   
    // Перегрзука in-place деления на скаляр
    SquareDenseMatrix& operator /= (
        const double& coef
    );
    
    // Перегрузка доступа по индексу
    std::vector<double>& operator[](
        int i
    );
    
    // Умножение матрицы на вектор (Ax = y)
    void dot_vector(
        const std::vector<double>& x,
        std::vector<double>& y
    );
};
