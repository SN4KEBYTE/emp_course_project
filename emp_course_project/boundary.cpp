// boundary.cpp

#include "boundary.h"
#include "sparse_matrix.h"

void first_kind_boundary_cond(
    SparseMatrix& A,
    std::vector<double>& f_global,
    const int& node1,
    const double& f_val1,
    const int& node2,
    const double& f_val2
)
{
    // Ставим на соответствующие позиции в глобальном векторе правой части значения функции
    f_global[node1] = f_val1;
    f_global[node2] = f_val2;

    // На диагональ единицы
    A.diag[node1] = 1;
    A.diag[node2] = 1;

    // Обнуляем соответствующие строки
    A.zero_out_row(node1);
    A.zero_out_row(node2);
}

void second_kind_boundary_cond(
    std::vector<double>& b_global,
    const std::vector<int>& cur_cond,
    const int& cond_num,
    FiniteElement& el,
    const int& el_num,
    const std::vector<Node>& grid,
    const double& t,
    double(*bc2_func)(const double& x, const double& y, const double& t)
)
{
    int b1 = -1, b2 = -1, b3 = -1;

    // Последний элемент - номер функции
    for (int k = 0; k < cur_cond.size() - 1; k++)
    {
        if (cur_cond[k] == el.n1)
            b1 = cond_num;

        if (cur_cond[k] == el.n2)
            b2 = cond_num;

        if (cur_cond[k] == el.n3)
            b3 = cond_num;
    }

    // Если оказалось, что встретили 2 узла, принадлежащих одному ребру
    if (b1 != -1 && b2 == b1)
    {
        double coef = mesG(grid[el.n1], grid[el.n2]) / 6,
            teta1 = bc2_func(grid[el.n1].x, grid[el.n1].y, t),
            teta2 = bc2_func(grid[el.n2].x, grid[el.n2].y, t);
        b_global[el[0]] += coef * (2 * teta1 + teta2);
        b_global[el[1]] += coef * (teta1 + 2 * teta2);
    }

    if (b1 != -1 && b3 == b1)
    {
        double coef = mesG(grid[el.n1], grid[el.n3]) / 6,
            teta1 = bc2_func(grid[el.n1].x, grid[el.n1].y, t),
            teta2 = bc2_func(grid[el.n3].x, grid[el.n3].y, t);
        b_global[el[0]] += coef * (2 * teta1 + teta2);
        b_global[el[2]] += coef * (teta1 + 2 * teta2);
    }

    if (b3 != -1 && b2 == b3)
    {
        double coef = mesG(grid[el.n2], grid[el.n3]) / 6,
            teta1 = bc2_func(grid[el.n2].x, grid[el.n2].y, t),
            teta2 = bc2_func(grid[el.n3].x, grid[el.n3].y, t);
        b_global[el[1]] += coef * (2 * teta1 + teta2);
        b_global[el[2]] += coef * (teta1 + 2 * teta2);
    }
}

void third_kind_boundary_cond(
    SparseMatrix& A,
    std::vector<double>& b_global,
    const std::vector<int>& cur_cond,
    const double& beta,
    const int& cond_num,
    FiniteElement& el,
    const int& el_num,
    const std::vector<Node>& grid,
    const double& t,
    double(*bc3_func)(const double& x, const double& y, const double& t),
    bool ignore_matrix
)
{
    int b1 = -1, b2 = -1, b3 = -1;
    FiniteElement el_to_add;

    SquareDenseMatrix A_loc(2);
    A_loc.data = { {2, 1}, {1, 2} };

    // Последний элемент - номер функции
    for (int k = 0; k < cur_cond.size() - 1; k++)
    {
        // Ищем 2 номера узлов, которые принадлежат одному краевому условию
        if (cur_cond[k] == el.n1)
            b1 = cond_num;

        if (cur_cond[k] == el.n2)
            b2 = cond_num;

        if (cur_cond[k] == el.n3)
            b3 = cond_num;
    }

    if (b1 != -1 && b2 == b1)
    {
        double coef = beta * mesG(grid[el.n1], grid[el.n2]) / 6,
            u1 = bc3_func(grid[el.n1].x, grid[el.n1].y, t),
            u2 = bc3_func(grid[el.n2].x, grid[el.n2].y, t);

        A_loc *= coef;

        // Сразу добавляем нужные значения в глобальный вектор
        b_global[el[0]] += coef * (2 * u1 + u2);
        b_global[el[1]] += coef * (u1 + 2 * u2);

        el_to_add.n1 = el.n1;
        el_to_add.n2 = el.n2;

        // Добавляем локальную матрицу в глобальную
        if (!ignore_matrix)
            A.add_local(A_loc, el_to_add);
    }

    if (b1 != -1 && b3 == b1)
    {
        double coef = beta * mesG(grid[el.n1], grid[el.n3]) / 6,
            u1 = bc3_func(grid[el.n1].x, grid[el.n1].y, t),
            u2 = bc3_func(grid[el.n3].x, grid[el.n3].y, t);

        A_loc *= coef;

        // Сразу добавляем нужные значения в глобальный вектор
        b_global[el[0]] += coef * (2 * u1 + u2);
        b_global[el[2]] += coef * (u1 + 2 * u2);

        el_to_add.n1 = el.n1;
        el_to_add.n2 = el.n3;

        // Добавляем локальную матрицу в глобальную
        A.add_local(A_loc, el_to_add);
    }

    if (b3 != -1 && b2 == b3)
    {
        double coef = beta * mesG(grid[el.n2], grid[el.n3]) / 6,
            u1 = bc3_func(grid[el.n2].x, grid[el.n2].y, t),
            u2 = bc3_func(grid[el.n3].x, grid[el.n3].y, t);

        A_loc *= coef;

        // Сразу добавляем нужные значения в глобальный вектор
        b_global[el[1]] += coef * (2 * u1 + u2);
        b_global[el[2]] += coef * (u1 + 2 * u2);

        el_to_add.n1 = el.n2;
        el_to_add.n2 = el.n3;

        // Добавляем локальную матрицу в глобальную
        A.add_local(A_loc, el_to_add);
    }
}
