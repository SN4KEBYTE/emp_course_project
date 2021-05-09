// finite_element.h

#pragma once

#include <iostream>

/*
Структура, описывающая конечный элемент
*/
struct FiniteElement
{
    // ni - глобальный номер i-го узла элемента
    int n1 = 0, n2 = 0, n3 = 0;
    
    // Конструктор по умолчанию
    FiniteElement() {}
    
    // Конструктор с параметрами
    FiniteElement(
        int n1, 
        int n2, 
        int n3
    ) : n1(n1), n2(n2), n3(n3) {}
    
    // Доступ по индексу
    int& operator[](
        int i
    )
    {
        switch (i)
        {
        case 0:
            return n1;
            break;
        case 1:
            return n2;
            break;
        case 2:
            return n3;
            break;
        default:
            throw std::out_of_range("FiniteElement: index out of range");
            break;
        }
    }
};

// Перегрузка оператора ввода (для удобства)
std::istream& operator >> (
    std::istream& in, 
    FiniteElement& f_el
);

// Перегрузка оператора вывода (для удобства)
std::ostream & operator << (
    std::ostream & out, 
    const FiniteElement & f_el
);