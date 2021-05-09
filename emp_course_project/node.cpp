// node.cpp

#include "node.h"

std::istream& operator >> (
    std::istream& in, 
    Node& n
)
{
    in >> n.x >> n.y;
    
    return in;
}

std::ostream& operator << (
    std::ostream& out, 
    const Node& n
)
{
    out << n.x << " " << n.y;
    
    return out;
}

double mesG(
    const Node& n1, 
    const Node& n2
)
{
    double x_diff = n2.x - n1.x, y_diff = n2.y - n1.y;
    
    return sqrt(x_diff * x_diff + y_diff * y_diff);
}