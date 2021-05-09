// node.h

#pragma once

#include <iostream>

/*
Структура, описывающая узел конечноэлементной сетки
*/
struct Node
{
    // x - x-координата узла, y - y-координата узла
    double x = 0, y = 0;
    
    // Конструктор по умолчанию
    Node() {}
    
    // Конструктор с параметрами
    Node(double x, double y) : x(x), y(y) {}
};

// Перегрузка оператора ввода (для удобства)
std::istream& operator >> (
    std::istream& in, 
    Node& n
);

// Перегрузка оператора вывода (для удобства)
std::ostream& operator << (
    std::ostream& out, 
    const Node& n
);

// Длина ребра
double mesG(
    const Node& n1, 
    const Node& n2
);