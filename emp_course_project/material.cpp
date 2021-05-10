#include "material.h"

std::istream& operator >> (
    std::istream& in,
    Material& m
)
{
    in >> m.lambda;
    
    for (auto& g : m.gamma)
        in >> g;
    
    return in;
}

std::ostream& operator << (
    std::ostream& out, 
    const Material& m
)
{
    out << m.lambda << " ";
    
    for (const auto& g : m.gamma)
        out << g << " ";
    
    return out;
}
