#include "bi_cgm.h"

void BiCGM::solve(
    SparseMatrix& A,
    std::vector<double>& xk,
    const std::vector<double>& f
)
{
    // Выделяем память под вспомогательные вектора
    r.resize(A.get_n());
    p.resize(A.get_n());
    z.resize(A.get_n());
    s.resize(A.get_n());
    az.resize(A.get_n());
    ats.resize(A.get_n());

    // r = A * xk
    A.dot_vector(xk, r);

    // r = f - A * xk = f - r
    vec_diff(r, f, r);

    // Инициализируем вспомогательные вектора
    p = r;
    z = r;
    s = r;

    // Число итераций
    int iter = 0;

    // Вычисляем норму вектора правой части и начальное значение невязки
    double f_norm = vec_norm(f), r_norm = vec_norm(r), rr = r_norm / f_norm;

    // Пока не достигнуто заданное значение минимальной величины относительной невязки и не превышено максимальное число итераций
    while (rr >= eps && iter < max_iter)
    {
        // az = A * z
        A.dot_vector(z, az);

        // ats = A^(T) * s
        A.T_dot_vector(s, ats);

        // Скалярное произведение p и r с предыдщуей итерации необходимо сохранить,
        // оно понадобится для вычисления bk после изменения p и r
        double scal_pr = scalar_product(p, r), alpha_k = scal_pr / scalar_product(s, az);

        // Обновляем вектора x, r, p
        for (int i = 0; i < A.get_n(); i++)
        {
            xk[i] += alpha_k * z[i];
            r[i] -= alpha_k * az[i];
            p[i] -= alpha_k * ats[i];
        }

        // Вычисляем коэффициент Bk
        double beta_k = scalar_product(p, r) / scal_pr;

        // Обновляем вектора z, s
        for (int i = 0; i < A.get_n(); i++)
        {
            z[i] = r[i] + beta_k * z[i];
            s[i] = p[i] + beta_k * s[i];
        }

        // Вычисляем новое значение невязки
        r_norm = vec_norm(r);
        rr = r_norm / f_norm;

        // Увеличиваем число итераций
        iter++;
    }

    // После решения память можно освободить
    vec_free(r);
    vec_free(p);
    vec_free(z);
    vec_free(s);
    vec_free(az);
    vec_free(ats);
}
