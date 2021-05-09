#include "finite_element.h"

std::istream& operator >> (
    std::istream& in, 
    FiniteElement& f_el
)
{
    in >> f_el.n1;
    in >> f_el.n2;
    in >> f_el.n3;
    
    // Упорядочиваем номера узлов по возрастанию
    if (f_el.n1 > f_el.n2)
        std::swap(f_el.n1, f_el.n2);
    
    if (f_el.n2 > f_el.n3)
        std::swap(f_el.n2, f_el.n3);
    
    if (f_el.n1 > f_el.n2)
        std::swap(f_el.n1, f_el.n2);
    
    return in;
}

std::ostream& operator << (
    std::ostream& out, 
    const FiniteElement& f_el
)
{
    out << f_el.n1 << " " << f_el.n2 << " " << f_el.n3;
    
    return out;
}