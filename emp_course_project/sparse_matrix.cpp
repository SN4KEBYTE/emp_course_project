// sparse_matrix.cpp

#include "sparse_matrix.h"

SparseMatrix::SparseMatrix(
    int n
)
{
    this->n = n;
    diag.resize(n);
    ig.resize(n + 1);
}

SparseMatrix::SparseMatrix(
    int n,
    std::vector<int> ig,
    std::vector<int> jg,
    std::vector<double> ggu,
    std::vector<double> ggl,
    std::vector<double> diag
)
{
    this->n = n;
    this->ig = ig;
    this->jg = jg;
    this->ggu = ggu;
    this->ggl = ggl;
    this->diag = diag;
}

void SparseMatrix::dot_vector(
    const std::vector<double>& x, 
    std::vector<double>& y
)
{
    for (int i = 0; i < n; i++)
    {
        y[i] = diag[i] * x[i];
        
        for (int j = ig[i]; j < ig[i + 1]; j++)
        {
            y[i] += ggl[j] * x[jg[j]];
            y[jg[j]] += ggu[j] * x[i];
        }
    }
}

SparseMatrix SparseMatrix::dot_scalar(
    const double& x
)
{
    std::vector<double> tmp_ggu(ggu),
        tmp_ggl(ggl),
        tmp_diag(diag);

    for (int i = 0; i < tmp_ggu.size(); i++)
    {
        tmp_ggu[i] *= x;
        tmp_ggl[i] *= x;
    }

    for (int i = 0; i < diag.size(); i++)
        tmp_diag[i] *= x;

    return SparseMatrix(
        this->n,
        this->ig,
        this->jg,
        tmp_ggu,
        tmp_ggl,
        diag
    );
}

void SparseMatrix::T_dot_vector(
    const std::vector<double>& x, 
    std::vector<double>& y
)
{
    for (int i = 0; i < n; i++)
    {
        y[i] = diag[i] * x[i];
        
        for (int j = ig[i]; j < ig[i + 1]; j++)
        {
            y[i] += ggu[j] * x[jg[j]];
            y[jg[j]] += ggl[j] * x[i];
        }
    }
}

void SparseMatrix::build_portrait(
    const int& num_elems, 
    std::vector<FiniteElement>& elems
)
{
    // n - число глобальных базисных функций
    // num_elems - число конечных элементов
    // elems - конечные элементы
    // list - врменный массив для хранения списка
    /// для этого массива должно быть выделено вдвое больше памяти, чем для массива jg
    std::vector<std::vector<int>> list(2);
    
    for (auto& row : list)
        row.resize(n * (n - 1));
    
    // временный массив listbeg будет хранить номера элементов массива list,
    // с которых начинаются списки для каждой из базисных фукнций
    std::vector<int> listbeg(n, -1);
    
    // текущая длина массива list
    int list_size = -1;
    
    // цикл формирования списка одним проходом по элементам
    for (int i_elem = 0; i_elem < num_elems; i_elem++)
    {
        for (int i = 0; i < 3; i++)
        {
            int k = elems[i_elem][i];
            
            for (int j = i + 1; j < 3; j++)
            {
                int ind1 = k;
                int ind2 = elems[i_elem][j];
                
                if (ind2 < ind1)
                {
                    ind1 = ind2;
                    ind2 = k;
                }
                
                // заносить будем связь большего номера с меньшим, т.е ind2 с ind1
                int i_addr = listbeg[ind2];
                
                // списка еще не было, создаем его
                if (i_addr < 0)
                {
                    list_size++;
                    listbeg[ind2] = list_size;
                    list[0][list_size] = ind1;
                    list[1][list_size] = -1;
                }
                // список уже был
                else
                {
                    // ищем в списке ind1
                    while (list[0][i_addr] < ind1 && list[1][i_addr] > 0)
                        i_addr = list[1][i_addr];
                    
                    if (list[0][i_addr] > ind1)
                    {
                        // если найденный там элемент имеет больший номер,
                        // то нужно добавить перед ним, чтобы список был упорядоченным
                        list_size++;
                        list[0][list_size] = list[0][i_addr]; //перекладываем вперед найденный элемент
                        list[1][list_size] = list[1][i_addr];
                        list[0][i_addr] = ind1; // на его место помещаем новый
                        list[1][i_addr] = list_size;
                    }
                    else if (list[0][i_addr] < ind1)
                    {
                        //не нашли, а список закончился
                        //добавляем в конец списка
                        list_size++;
                        list[1][i_addr] = list_size;
                        list[0][list_size] = ind1;
                        list[1][list_size] = -1; //указываем, что это последний элемент списка
                    }
                }
            }
        }
    }
    
    // Создание портрета по списку
    jg.resize(list_size + 1);
    
    ig[0] = 0;
    
    for (int i = 0; i < n; i++)
    {
        // ig[i + 1] фактически указывает на ячейку массива jg, куда надо поместить следующий элемент
        ig[i + 1] = ig[i];
        
        int i_addr = listbeg[i];
        
        while (i_addr != -1)
        {
            jg[ig[i + 1]] = list[0][i_addr];
            ig[i + 1]++;
            i_addr = list[1][i_addr];
        }
    }
    
    // Теперь мы знаем, сколько элементов хранится в верхнем и нижнем треугольниках
    ggu.resize(ig.back());
    ggl.resize(ig.back());
    
    // Освободим память
    vec_free(listbeg);
    vec_free(list);
}

void SparseMatrix::add_local(
    SquareDenseMatrix& mat, 
    const FiniteElement& el
)
{
    std::vector<int> L = { el.n1, el.n2, el.n3 };
    
    for (int i = 0; i < mat.get_n(); i++)
        diag[L[i]] += mat.data[i][i];
    
    for (int i = 1; i < mat.get_n(); i++)
    {
        int m = ig[L[i]];

        for (int j = 0; j < i; j++)
            for (int s = m; s < ig[L[i] + 1]; s++)
                if (jg[s] == L[j])
                {
                    ggl[s] += mat.data[i][j];
                    ggu[s] += mat.data[j][i];
                    
                    m++;
                    
                    break;
                }
    }
}

void SparseMatrix::zero_out_row(
    const int& row_num
)
{
    // Обнуляем нужные элементы в нижнем треугольнике
    for (int i = ig[row_num]; i < ig[row_num + 1]; i++)
        ggl[i] = 0;
    
    // Обнуляем нужные элементы в верхнем треугольнике
    for (int i = row_num + 1; i < get_n(); i++)
        for (int j = ig[i]; j < ig[i + 1]; j++)
            if (jg[j] == row_num)
            {
                ggu[j] = 0;
                
                break;
            }
}

SparseMatrix& SparseMatrix::operator+=(const SparseMatrix& m)
{
    for (int i = 0; i < this->ggu.size(); i++)
    {
        this->ggu[i] += m.ggu[i];
        this->ggl[i] += m.ggl[i];
    }

    for (int i = 0; i < this->diag.size(); i++)
        this->diag[i] += m.diag[i];

    return *this;
}

int SparseMatrix::get_n()
{
    return n;
}