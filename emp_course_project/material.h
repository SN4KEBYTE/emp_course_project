#pragma once

#include <iostream>
#include <vector>

/*
Структура, хранящая значения коэффициентов на одном конечном элементе
*/
struct Material
{
    // Коэффициенты гамма в каждом узле конечного элемента
    std::vector<double> gamma = { 0, 0, 0 };
    
    // Коэффициент лямбда
    double lambda = 0;
    
    // Конструктор по умолчанию
    Material() {}
    
    // Конструктор с параметрами
    Material(
        double lambda, 
        const std::vector<double> gamma
    ) : lambda(lambda), gamma(gamma) {}
};

// Перегрузка оператора ввода (для удобства)
std::istream& operator >> (
    std::istream& in, 
    Material& m
);

// Перегрузка оператора вывода (для удобства)
std::ostream& operator << (
    std::ostream& out, 
    const Material& m
);