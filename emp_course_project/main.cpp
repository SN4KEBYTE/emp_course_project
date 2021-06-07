// main.cpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include "finite_element.h"
#include "node.h"
#include "material.h"
#include "bi_cgm.h"
#include "dense_matrix.h"
#include "boundary.h"

typedef std::vector<double(*)(const double& x, const double& y, const double& t)> f_funcs;
typedef std::vector<double(*)(const double& x, const double& y, const double& t)> bc1_funcs;
typedef std::vector<double(*)(const double& x, const double& y, const double& t)> bc2_funcs;
typedef std::vector<double(*)(const double& x, const double& y, const double& t)> bc3_funcs;
typedef std::vector<double(*)(const double& x, const double& y, const double& t)> es_funcs;
typedef std::vector<double(*)(const double& x, const double& y, const double& t)> u0_funcs;

const std::string BASE_PATH = "tests/";
const int TESTS_NUM = 11;

#pragma region Тест №0
// Правая часть
double test0_f(const double& x, const double& y, const double& t)
{
    return 0;
}

// Краевое условие 1-го рода
double test0_bc1_func(const double& x, const double& y, const double& t)
{
    return 1;
}

// Точное решение
double test0_exact_solution(const double& x, const double& y, const double& t)
{
    return 1;
}

// Начальное условие
double test0_u0(const double& x, const double& y, const double& t)
{
    return 1;
}
#pragma endregion

#pragma region Тест №1
// Правая часть
double test1_f(const double& x, const double& y, const double& t)
{
    return 6 * t * t + 4;
}

// Краевое условие 1-го рода
double test1_bc1_func(const double& x, const double& y, const double& t)
{
    return 1.5 * x + t * t * t + 2 * t;
}

// Краевое условие 2-го рода
double test1_bc2_func(const double& x, const double& y, const double& t)
{
    return 0;
}

// Краевое условие 3-го рода
double test1_bc3_func(const double& x, const double& y, const double& t)
{
    return 1.5 + 1.5 * x + t * t * t + 2 * t;
}

// Точное решение
double test1_exact_solution(const double& x, const double& y, const double& t)
{
    return 1.5 * x + t * t * t + 2 * t;
}

// Начальное условие
double test1_u0(const double& x, const double& y, const double& t)
{
    return 1.5 * x + t * t * t + 2 * t;
}
#pragma endregion

#pragma region Тест №2
// Правая часть
double test2_f(const double& x, const double& y, const double& t)
{
    return 0;
}

// Краевое условие 1-го рода
double test2_bc1_func(const double& x, const double& y, const double& t)
{
    return 1;
}

// Точное решение
double test2_exact_solution(const double& x, const double& y, const double& t)
{
    return 1;
}

// Начальное условие
double test2_u0(const double& x, const double& y, const double& t)
{
    return 1;
}
#pragma endregion

#pragma region Тест №3
// Правая часть
double test3_f(const double& x, const double& y, const double& t)
{
    return 2;
}

// Краевое условие 1-го рода
double test3_bc1_func(const double& x, const double& y, const double& t)
{
    return x + y + t;
}

// Точное решение
double test3_exact_solution(const double& x, const double& y, const double& t)
{
    return x + y + t;
}

// Начальное условие
double test3_u0(const double& x, const double& y, const double& t)
{
    return x + y + t;
}
#pragma endregion

#pragma region Тест №4
// Правая часть
double test4_f(const double& x, const double& y, const double& t)
{
    return 4 * t - 8;
}

// Краевое условие 1-го рода
double test4_bc1_func(const double& x, const double& y, const double& t)
{
    return x * x + y * y + t * t;
}

// Точное решение
double test4_exact_solution(const double& x, const double& y, const double& t)
{
    return x * x + y * y + t * t;
}

// Начальное условие
double test4_u0(const double& x, const double& y, const double& t)
{
    return x * x + y * y + t * t;
}
#pragma endregion

#pragma region Тест №5
// Правая часть
double test5_f(const double& x, const double& y, const double& t)
{
    return 6 * t * t - 12 * x - 12 * y;
}

// Краевое условие 1-го рода
double test5_bc1_func(const double& x, const double& y, const double& t)
{
    return x * x * x + y * y * y + t * t * t;
}

// Точное решение
double test5_exact_solution(const double& x, const double& y, const double& t)
{
    return x * x * x + y * y * y + t * t * t;
}

// Начальное условие
double test5_u0(const double& x, const double& y, const double& t)
{
    return x * x * x + y * y * y + t * t * t;
}
#pragma endregion

#pragma region Тест №6
// Правая часть
double test6_f(const double& x, const double& y, const double& t)
{
    return 8 * t * t * t - 24 * x * x - 24 * y * y;
}

// Краевое условие 1-го рода
double test6_bc1_func(const double& x, const double& y, const double& t)
{
    return x * x * x * x + y * y * y * y + t * t * t * t;
}

// Точное решение
double test6_exact_solution(const double& x, const double& y, const double& t)
{
    return x * x * x * x + y * y * y * y + t * t * t * t;
}

// Начальное условие
double test6_u0(const double& x, const double& y, const double& t)
{
    return x * x * x * x + y * y * y * y + t * t * t * t;
}
#pragma endregion

#pragma region Тест №7, 8, 9, 10
// Правая часть
double test7_f(const double& x, const double& y, const double& t)
{
    return 4 * t + 4;
}

// Краевое условие 1-го рода
double test7_bc1_func(const double& x, const double& y, const double& t)
{
    return t * t + 2 * t + 1;
}

// Точное решение
double test7_exact_solution(const double& x, const double& y, const double& t)
{
    return t * t + 2 * t + 1;
}

// Начальное условие
double test7_u0(const double& x, const double& y, const double& t)
{
    return t * t + 2 * t + 1;
}
#pragma endregion

void solve_parabolic_problem(
    const std::string& base_path,
    const int& num_test,
    double (*f)(const double& x, const double& y, const double& t),
    double (*u0)(const double& x, const double& y, const double& t),
    const bc1_funcs& bc1_f,
    const bc2_funcs& bc2_f,
    const bc3_funcs& bc3_f,
    double(*es_f)(const double& x, const double& y, const double& t)
);

double get_solution_in_point(
    const double& x,
    const double& y,
    const double& t,
    const std::vector<std::vector<double>>& q,
    const std::vector<Node>& grid,
    const std::vector<FiniteElement>& elems,
    const std::vector<double>& time_grid
);

int main()
{
#pragma region Функции правой части
    f_funcs all_f = { 
        test0_f,
        test1_f,
        test2_f, 
        test3_f,
        test4_f,
        test5_f,
        test6_f,
        test7_f,
        test7_f,
        test7_f,
        test7_f
    };
#pragma endregion

#pragma region Функции для краевых 1-го рода
    std::vector<bc1_funcs> all_bc1 = {
        {test0_bc1_func},
        {test1_bc1_func},
        {test2_bc1_func},
        {test3_bc1_func},
        {test4_bc1_func},
        {test5_bc1_func},
        {test6_bc1_func},
        {test7_bc1_func},
        {test7_bc1_func},
        {test7_bc1_func},
        {test7_bc1_func},
    };
#pragma endregion

#pragma region Функции для краевых 2-го рода
    std::vector<bc2_funcs> all_bc2 = {
        {},
        {test1_bc2_func},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
    };
#pragma endregion

#pragma region Функции для краевых 3-го рода
    std::vector<bc3_funcs> all_bc3 = {
        {},
        {test1_bc3_func},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
        {},
    };
#pragma endregion

#pragma region Точные решения
    es_funcs all_es = {
        test0_exact_solution,
        test1_exact_solution,
        test2_exact_solution,
        test3_exact_solution,
        test4_exact_solution,
        test5_exact_solution,
        test6_exact_solution,
        test7_exact_solution,
        test7_exact_solution,
        test7_exact_solution,
        test7_exact_solution,
    };
#pragma endregion

#pragma region Начальные условия
    u0_funcs all_u0 = {
        test0_u0,
        test1_u0,
        test2_u0,
        test3_u0,
        test4_u0,
        test5_u0,
        test6_u0,
        test7_u0,
        test7_u0,
        test7_u0,
        test7_u0,
    };
#pragma endregion

    for (int tn = 0; tn < TESTS_NUM; tn++)
    {
        std::cout << "RUNNING TEST CASE #" << tn << std::endl;

        solve_parabolic_problem(
            BASE_PATH,
            tn,
            all_f[tn],
            all_u0[tn],
            all_bc1[tn],
            all_bc2[tn],
            all_bc3[tn],
            all_es[tn]
        );
    }

    return 0;
}

void solve_parabolic_problem(
    const std::string& base_path,
    const int& num_test,
    double (*f)(const double& x, const double& y, const double& t),
    double (*u0)(const double& x, const double& y, const double& t),
    const bc1_funcs& bc1_f,
    const bc2_funcs& bc2_f,
    const bc3_funcs& bc3_f,
    double(*es_f)(const double& x, const double& y, const double& t)
)
{
    std::string path = base_path + std::to_string(num_test) + "/";

#pragma region Ввод пространственной сетки
    // Число узлов и число конечных элементов
    int num_nodes = 0, num_elems = 0;

    // Ввод глобальной информации - количество узлов, количество элементов
    std::ifstream in(path + "info.txt");

    if (!in.good())
    {
        std::cout << "error occurred while opening info.txt" << std::endl;

        return;
    }

    in >> num_nodes >> num_elems;
    in.close();

    // Ввод узлов конечноэлементной сетки
    in.open(path + "xy.txt");

    if (!in.good())
    {
        std::cout << "error occurred while opening xy.txt" << std::endl;
        return;
    }

    std::vector<Node> grid(num_nodes);

    for (auto& n : grid)
        in >> n;

    in.close();

    // Ввод конечных элементов
    in.open(path + "elem.txt");

    if (!in.good())
    {
        std::cout << "error occurred while opening elem.txt" << std::endl;
        return;
    }

    std::vector<FiniteElement> elems(num_elems);

    for (auto& el : elems)
        in >> el;

    in.close();

    // Ввод материалов
    in.open(path + "mat.txt");

    if (!in.good())
    {
        std::cout << "error occurred while opening mat.txt" << std::endl;
        return;
    }

    std::vector<Material> materials(num_elems);

    for (auto& m : materials)
        in >> m;

    in.close();
#pragma endregion

#pragma region Ввод временной сетки
    int num_layers = 0;

    in.open(path + "time.txt");

    if (!in.good())
    {
        std::cout << "error occurred while opening time.txt" << std::endl;

        return;
    }

    in >> num_layers;

    std::vector<double> time_grid(num_layers);

    for (auto& el : time_grid)
        in >> el;

    in.close();
#pragma endregion

#pragma region Формирование глобальных матриц масс, жесткости и вектора правой части
    // Глобальные матрицы масс и жесткости
    SparseMatrix M_global(num_nodes), G_global(num_nodes), M_S3(num_nodes);

    // Строим портреты
    M_global.build_portrait(num_elems, elems);
    G_global.build_portrait(num_elems, elems);
    M_S3.build_portrait(num_elems, elems);

    // Глобальный вектор правой части
    std::vector<double> b_global(num_nodes);

    SquareDenseMatrix D_inv(3), G_local(3), C(3), M_local(3);

    // В цикле по конечным элементам
    for (int i = 0; i < num_elems; i++)
    {
        Node n1 = grid[elems[i].n1]; // x1, y1
        Node n2 = grid[elems[i].n2]; // x2, y2
        Node n3 = grid[elems[i].n3]; // x3, y3

        double det_D = (n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) *
            (n2.y - n1.y);
        D_inv.data[0][0] = n2.x * n3.y - n3.x * n2.y; // a10
        D_inv.data[0][1] = n2.y - n3.y; // a11
        D_inv.data[0][2] = n3.x - n2.x; // a12
        D_inv.data[1][0] = n3.x * n1.y - n1.x * n3.y; // a20
        D_inv.data[1][1] = n3.y - n1.y; // a21
        D_inv.data[1][2] = n1.x - n3.x; // a22
        D_inv.data[2][0] = n1.x * n2.y - n2.x * n1.y; // a30
        D_inv.data[2][1] = n1.y - n2.y; // a31
        D_inv.data[2][2] = n2.x - n1.x; // a32
        D_inv /= det_D;

        // Локальная матрица жесткости
        double coef = materials[i].lambda * abs(det_D) / 2;

        for (int j = 0; j < G_local.data.size(); j++)
            for (int k = 0; k < G_local.data[j].size(); k++)
                G_local.data[j][k] = coef * (D_inv.data[j][1] * D_inv.data[k][1] + D_inv.data[j][2] * D_inv.data[k][2]);

        // Вспомогательная матрица C
        coef = abs(det_D) / 24;

        for (int j = 0; j < C.data.size(); j++)
            for (int k = 0; k < C.data[j].size(); k++)
                C.data[j][k] = coef * (j == k ? 2 : 1);

        double g1 = materials[i].gamma[0],
            g2 = materials[i].gamma[1],
            g3 = materials[i].gamma[2];

        // Локальная матрица масс
        M_local.data = {
            {6 * g1 + 2 * g2 + 2 * g3, 2 * g1 + 2 * g2 + g3, 2 * g1 + g2 + 2 * g3},
            {2 * g1 + 2 * g2 + g3, 2 * g1 + 6 * g2 + 2 * g3, g1 + 2 * g2 + 2 * g3},
            {2 * g1 + g2 + 2 * g3, g1 + 2 * g2 + 2 * g3, 2 * g1 + 2 * g2 + 6 * g3}
        };

        M_local *= abs(det_D) / 120;

        // Вносим локальные матрицы в глобальную
        M_global.add_local(M_local, elems[i]);
        G_global.add_local(G_local, elems[i]);

        // Локальный вектор правой части
        std::vector<double> f_loc(3), b_loc(3);
        f_loc[0] = f(n1.x, n1.y, time_grid[2]);
        f_loc[1] = f(n2.x, n2.y, time_grid[2]);
        f_loc[2] = f(n3.x, n3.y, time_grid[2]);

        C.dot_vector(f_loc, b_loc);

        // Вносим локальный вектор правой части в глобальный
        b_global[elems[i].n1] += b_loc[0];
        b_global[elems[i].n2] += b_loc[1];
        b_global[elems[i].n3] += b_loc[2];
    }
#pragma endregion

#pragma region Ввод краевых условий 1-го рода
    in.open(path + "bc1.txt");

    if (!in.good())
    {
        std::cout << "error occurred while opening bc1.txt" << std::endl;
        return;
    }

    int num_bc1 = 0;
    in >> num_bc1;

    std::vector<std::vector<int>> bc1(num_bc1, std::vector<int>(3));

    for (int i = 0; i < num_bc1; i++)
        // Краевое условие 1-го рода задается 2-мя глобальными номерами узлов и номером функции
        in >> bc1[i][0] >> bc1[i][1] >> bc1[i][2];

    in.close();
#pragma endregion

#pragma region Ввод краевых условий 2-го рода
    in.open(path + "bc2.txt");

    if (!in.good())
    {
        std::cout << "error occurred while opening bc2.txt" << std::endl;
        return;
    }

    int num_bc2 = 0;
    in >> num_bc2;

    std::vector<std::vector<int>> bc2(num_bc2);

    for (int i = 0; i < num_bc2; i++)
    {
        // Краевое условие задается ребром и номером функции
        bc2[i].resize(3);
        for (int j = 0; j < 3; j++)
            in >> bc2[i][j];
    }

    in.close();
#pragma endregion

#pragma region Ввод краевых условий 3-го рода
    in.open(path + "bc3.txt");

    if (!in.good())
    {
        std::cout << "error occurred while opening bc3.txt" << std::endl;
        return;
    }

    int num_bc3 = 0;
    in >> num_bc3;

    std::vector<std::vector<int>> bc3(num_bc3);
    std::vector<double> bc3_beta(num_bc3);

    for (int i = 0; i < num_bc3; i++)
    {
        // Краевое условие задается набором номеров узлов, номером функции и значением beta
        bc3[i].resize(3);
        for (int j = 0; j < 3; j++)
            in >> bc3[i][j];
        in >> bc3_beta[i];
    }

    in.close();
#pragma endregion

#pragma region Учет краевых условий 2-го и 3-го рода
    // В цикле по конечным элементам
    for (int i = 0; i < num_elems; i++)
    {
        // Обрабатываем каждое краевое условие 2-го рода
        for (int j = 0; j < num_bc2; j++)
            second_kind_boundary_cond(
                b_global,
                bc2[j],
                j,
                elems[i],
                i,
                grid,
                time_grid[2],
                bc2_f[bc2[j].back()]
            );

        // Обрабатываем каждое краевое условие 3-го рода
        for (int j = 0; j < num_bc3; j++)
            third_kind_boundary_cond(
                M_S3,
                b_global,
                bc3[j],
                bc3_beta[j],
                j,
                elems[i],
                i,
                grid,
                time_grid[2],
                bc3_f[bc3[j].back()]
            );
    }
#pragma endregion

#pragma region Решение СЛАУ
    std::vector<std::vector<double>> q(num_layers, std::vector<double>(num_nodes));

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < q[i].size(); j++)
            q[i][j] = u0(grid[j].x, grid[j].y, time_grid[i]);

    for (int i = 3; i < time_grid.size(); i++)
    {
        // Вычисляем коэффициенты
        double delta_t0 = time_grid[i] - time_grid[i - 1],
            delta_t1 = time_grid[i] - time_grid[i - 2],
            delta_t2 = time_grid[i] - time_grid[i - 3],
            delta_t3 = time_grid[i - 1] - time_grid[i - 2],
            delta_t4 = time_grid[i - 1] - time_grid[i - 3],
            delta_t5 = time_grid[i - 2] - time_grid[i - 3];

        double tetta0 = (delta_t1 * delta_t0 + delta_t0 * delta_t2 + delta_t1 * delta_t2) / delta_t1 / delta_t0 / delta_t2,
            tetta1 = delta_t0 * delta_t1 / delta_t2 / delta_t4 / delta_t5,
            tetta2 = delta_t0 * delta_t2 / delta_t1 / delta_t3 / delta_t5,
            tetta3 = delta_t1 * delta_t2 / delta_t0 / delta_t3 / delta_t4;

        std::vector<double> tmp(b_global.size()), right(b_global);

        SparseMatrix matM2 = M_global.dot_scalar(tetta1);
        matM2.dot_vector(q[i - 1], tmp);
        vec_sum(right, right, tmp);

        SparseMatrix matM3 = M_global.dot_scalar(tetta2);
        matM3.dot_vector(q[i - 2], tmp);
        vec_diff(right, right, tmp);

        matM2 = G_global.dot_scalar(-1.0);
        SparseMatrix matM4 = M_global.dot_scalar(tetta3);
        matM2 += matM4;
        matM2.dot_vector(q[i - 3], tmp);
        vec_sum(right, right, tmp);

        SparseMatrix A = M_global;
        A += M_S3;
        A = A.dot_scalar(tetta0);

        // Учет первых краевых
        for (int j = 0; j < num_bc1; j++)
        {
            int bc1_node1 = bc1[j][0], bc1_node2 = bc1[j][1], f_num = bc1[j][2];

            double val1 = bc1_f[f_num](grid[bc1_node1].x, grid[bc1_node1].y, time_grid[i]),
                val2 = bc1_f[f_num](grid[bc1_node2].x, grid[bc1_node2].y, time_grid[i]);

            first_kind_boundary_cond(
                A,
                right,
                bc1_node1,
                val1,
                bc1_node2,
                val2
            );
        }

        // Решаем СЛАУ
        BiCGM solver;
        solver.solve(
            A, 
            q[i], 
            right
        );

        // Обновляем глобальный вектор
        std::fill(
            b_global.begin(),
            b_global.end(), 
            0.
        );
        
        for (int j = 0; j < num_elems; j++)
        {
            Node n1 = grid[elems[j].n1]; // x1, y1
            Node n2 = grid[elems[j].n2]; // x2, y2
            Node n3 = grid[elems[j].n3]; // x3, y3

            double det_D = (n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y),
                coef = abs(det_D) / 24;

            for (int k = 0; k < C.data.size(); k++)
                for (int p = 0; p < C.data[k].size(); p++)
                    C.data[k][p] = coef * (k == p ? 2 : 1);

            // Локальный вектор правой части
            std::vector<double> f_loc(3), b_loc(3);
            f_loc[0] = f(n1.x, n1.y, time_grid[i]);
            f_loc[1] = f(n2.x, n2.y, time_grid[i]);
            f_loc[2] = f(n3.x, n3.y, time_grid[i]);

            C.dot_vector(f_loc, b_loc);

            // Вносим локальный вектор правой части в глобальный
            b_global[elems[j].n1] += b_loc[0];
            b_global[elems[j].n2] += b_loc[1];
            b_global[elems[j].n3] += b_loc[2];

            // Обрабатываем каждое краевое условие 2-го рода
            for (int k = 0; k < num_bc2; k++)
                second_kind_boundary_cond(
                    b_global,
                    bc2[k],
                    k,
                    elems[j],
                    j,
                    grid,
                    time_grid[i],
                    bc2_f[bc2[k].back()]
                );

            // Обрабатываем каждое краевое условие 3-го рода
            for (int k = 0; k < num_bc3; k++)
                third_kind_boundary_cond(
                    M_S3,
                    b_global,
                    bc3[k],
                    bc3_beta[k],
                    k,
                    elems[j],
                    j,
                    grid,
                    time_grid[i],
                    bc3_f[bc3[k].back()],
                    true
                );
        }
    }
#pragma endregion

#pragma region Вывод результата
    int num_test_points = 0;

    in.open(path + "test_points.txt");
    in >> num_test_points;

    std::ofstream out(path + "res.csv");
    out.precision(17);
    
    out << "x,y,t,u,u*,err\n";

    for (int i = 0; i < num_test_points; i++)
    {
        double test_x = 0, test_y = 0;
        in >> test_x >> test_y;

        for (int j = 0; j < time_grid.size(); j++)
        {
            double u = get_solution_in_point(test_x, test_y, time_grid[j], q, grid, elems, time_grid),
                u_es = es_f(test_x, test_y, time_grid[j]);

            if (j == 0)
                out << test_x << "," << test_y << ",";
            else
                out << ",,";

            out  << time_grid[j] << "," << u << "," << u_es << "," << abs(u - u_es) << "\n";
        }
    }
    
    in.close();
    out.close();
#pragma endregion
}

double get_solution_in_point(
    const double& x,
    const double& y,
    const double& t,
    const std::vector<std::vector<double>>& q,
    const std::vector<Node>& grid,
    const std::vector<FiniteElement>& elems,
    const std::vector<double>& time_grid
)
{
    double det_D = 0.;
    int i = 0;
    std::vector<double> s(3);
    FiniteElement cur_el;

    // Ищем конечный элемент, которому принадлежит точка
    for (i = 0; i < elems.size(); i++)
    {
        cur_el = elems[i];
        Node n1 = grid[elems[i].n1],
            n2 = grid[elems[i].n2],
            n3 = grid[elems[i].n3];
        det_D = abs((n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) *
            (n2.y - n1.y));

        s[0] = abs((n2.x - n1.x) * (y - n1.y) - (x - n1.x) * (n2.y - n1.y));
        s[1] = abs((n1.x - n3.x) * (y - n3.y) - (x - n3.x) * (n1.y - n3.y));
        s[2] = abs((n2.x - n3.x) * (y - n3.y) - (x - n3.x) * (n2.y - n3.y));

        if (abs(det_D - s[0] - s[1] - s[2]) < 1e-10)
            break;
    }

    if (i == elems.size())
        throw "Точка не принадлежит расчетной области";

    int j = 0;

    // Ищем временной отрезок, которому принадлежит точка
    for (j = 1; j < time_grid.size(); j++)
    {
        if (t >= time_grid[j - 1] && t <= time_grid[j])
            break;
    }

    if (j == time_grid.size())
        throw "Точка не принадлежит заданному временному отрезку";

    // Если не хватает слоев, аппроксимируем по времени линейным базисом
    if (j < 2)
    {
        double dt = time_grid[j] - time_grid[j - 1];
        auto q0 = q[j - 1],
            cur_q = q[j];

        double u0 = (q0[cur_el.n1] * s[2] + q0[cur_el.n2] * s[1] + q0[cur_el.n3] * s[0]) / det_D,
            u = (cur_q[cur_el.n1] * s[2] + cur_q[cur_el.n2] * s[1] + cur_q[cur_el.n3] * s[0]) / det_D;

        return u0 * (time_grid[j] - t) / dt + u * (t - time_grid[j - 1]) / dt;
    }
    // Если не хватает слоев, аппроксимируем по времени квадратичным базисом
    else if (j < 3)
    {
        double dt = time_grid[j] - time_grid[j - 2],
            dt1 = time_grid[j - 1] - time_grid[j - 2],
            dt0 = time_grid[j] - time_grid[j - 1];

        auto q0 = q[j - 2],
            q1 = q[j - 1],
            cur_q = q[j];

        double u1 = (q1[cur_el.n1] * s[2] + q1[cur_el.n2] * s[1] + q1[cur_el.n3] * s[0]) / det_D,
            u0 = (q0[cur_el.n1] * s[2] + q0[cur_el.n2] * s[1] + q0[cur_el.n3] * s[0]) / det_D,
            u = (cur_q[cur_el.n1] * s[2] + cur_q[cur_el.n2] * s[1] + cur_q[cur_el.n3] * s[0]) / det_D;

        return 1. / (dt1 * dt) * (t - time_grid[j - 1]) * (t - time_grid[j]) * u0 +
            (-1.) / (dt1 * dt0) * (t - time_grid[j - 2]) * (t - time_grid[j]) * u1 +
            1. / (dt * dt0) * (t - time_grid[j - 2]) * (t - time_grid[j - 1]) * u;
    }
    // Если хватает слоев, аппроксимируем по времени кубическим базисом
    else
    {
        double delta_t0 = time_grid[j] - time_grid[j - 1],
            delta_t1 = time_grid[j] - time_grid[j - 2],
            delta_t2 = time_grid[j] - time_grid[j - 3],
            delta_t3 = time_grid[j - 1] - time_grid[j - 2],
            delta_t4 = time_grid[j - 1] - time_grid[j - 3],
            delta_t5 = time_grid[j - 2] - time_grid[j - 3];

        // изначально
        double tetta0 = (t - time_grid[j - 3]) * (t - time_grid[j - 2]) * (t - time_grid[j - 1]) / (delta_t2 * delta_t1 * delta_t0),
            tetta1 = -(t - time_grid[j - 3]) * (t - time_grid[j - 2]) * (t - time_grid[j]) / (delta_t4 * delta_t3 * delta_t0),
            tetta2 = (t - time_grid[j - 3]) * (t - time_grid[j - 1]) * (t - time_grid[j]) / (delta_t5 * delta_t3 * delta_t1),
            tetta3 = -(t - time_grid[j - 2]) * (t - time_grid[j - 1]) * (t - time_grid[j]) / (delta_t5 * delta_t4 * delta_t2);

        /*double tetta0 = (t - time_grid[j - 3]) * (t - time_grid[j - 2]) * (t - time_grid[j - 1]) / (delta_t2 * delta_t1 * delta_t0),
            tetta1 = -(t - time_grid[j - 3]) * (t - time_grid[j - 2]) * (t - time_grid[j]) / (delta_t4 * delta_t3 * delta_t0),
            tetta2 = (t - time_grid[j - 3]) * (t - time_grid[j - 1]) * (t - time_grid[j]) / (delta_t5 * delta_t3 * delta_t1),
            tetta3 = -(t - time_grid[j - 2]) * (t - time_grid[j - 1]) * (t - time_grid[j]) / (delta_t5 * delta_t4 * delta_t2);*/

        auto q0 = q[j],
            q1 = q[j - 1],
            q2 = q[j - 2],
            q3 = q[j - 3];

        double u0 = (q0[cur_el.n1] * s[2] + q0[cur_el.n2] * s[1] + q0[cur_el.n3] * s[0]) / det_D,
            u1 = (q1[cur_el.n1] * s[2] + q1[cur_el.n2] * s[1] + q1[cur_el.n3] * s[0]) / det_D,
            u2 = (q2[cur_el.n1] * s[2] + q2[cur_el.n2] * s[1] + q2[cur_el.n3] * s[0]) / det_D,
            u3 = (q3[cur_el.n1] * s[2] + q3[cur_el.n2] * s[1] + q3[cur_el.n3] * s[0]) / det_D;

        return tetta3 * u3 + tetta2 * u2 + tetta1 * u1 + tetta0 * u0;
    }
}