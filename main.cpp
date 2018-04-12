// Copyright rafirafi 2018
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <cstdio>
#include <ctime>
#include <cstring>

#include <algorithm>
#include <vector>

#include "simplex_solver_bounded.h"

#define NA -1

#ifndef D
#define D 3
#endif

#define N (D * D)
#define NN (D * D * D * D)

/**************************************************************************************************/

template <typename T>
bool is_zero(const T &value)
{
    return std::abs(value) <= T{1e-9};
}

template <typename T>
bool is_neg(const T &value)
{
    return value < -T{1e-9};
}

template <typename T>
bool is_pos(const T &value)
{
    return value > T{1e-9};
}


// http://www.math.ryerson.ca/~danziger/professor/MTH141/Handouts/Slides/gauss.pdf
template <typename T>
void gaussian_algo(std::vector<std::vector<T> > &mat)
{
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");
    assert(!mat.empty());
    int m = mat.size(), n = mat[0].size();

    // The Gaussian Algorithm

    // The  following  Algorithm  reduces  an m × n matrix to
    // REF  by  means  of  elementary  row  operations alone.

    //1.  For Each row i ( R i ) from 1 to m
    for (int i = 0; i < m; i++) {

        // (a)  If any row j below row i has non zero entries to the left of the
        //first non zero entry in row i exchange row i and j ( Ri ↔ Rj )
        //[Ensure We are working on the leftmost nonzero entry.]
        int first_not_zero_row = NA, first_not_zero_column = NA;
        T first_element{0};
        for (int j = i; j < n && first_not_zero_row == NA; j++) {
            for (int l = i; l < m; l++) {
                if (!is_zero(mat[l][j])) {
                    first_element = mat[l][j];
                    first_not_zero_row = l;
                    first_not_zero_column = j;
                    break;
                }
            }
        }
        if (first_not_zero_row == NA) { // done
            m = i;
            break;
        }

        if (first_not_zero_row != i) {
            for (int j = i; j < n; j++) {
                std::swap(mat[i][j], mat[first_not_zero_row][j]);
            }
        }

        // (b)  Preform Ri → 1/c*Ri where c =  the  first  non- zero entry
        // of row i .  [This ensures that row i starts with a one.]
        for (int j = first_not_zero_column; j < n; j++) {
            mat[i][j] /= first_element;
        }

        // (c)  For each row j ( Rj ) below row i (Each j > i )
        // Preform Rj → Rj − dRi where d =  the  entry  in  row j
        // which  is  directly below  the  pivot  in  row i .
        // [This  ensures that row j has a 0 below the pivot of row i .]
        for (int j = i + 1; j < m; j++) { // R
            if (is_zero(mat[j][first_not_zero_column])) {
                continue;
            }
            T d = mat[j][first_not_zero_column] / mat[i][first_not_zero_column];
            mat[j][first_not_zero_column] = 0;
            for (int l = first_not_zero_column + 1; l < n; l++) { // C
                mat[j][l] -= d * mat[i][l];
            }
        }

        // (d)  If any 0 rows have appeared exchange them  to the bottom  of the matrix.
        for (int j = i + 1; j < m; j++) { // R
            bool empty = true;
            for (int l = i + 1; l < n && empty; l++) {
                empty = is_zero(mat[j][l]);
            }
            if (empty) {
                m--;
                if (j < m) {
                    for (int l = i + 1; l < n; l++) {
                        std::swap(mat[j][l], mat[m][l]);
                    }
                    j--;
                }
            }
        }
    }

    mat.resize(m);

    fprintf(stderr, "%-30s m %d n %d\n", __func__, m , n);
}

template <typename T>
void gaussian_jordan_algo(std::vector<std::vector<T> > &mat)
{
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");
    assert(!mat.empty());
    int m = mat.size(), n = mat[0].size();

    // The  Gaussian-Jordan Algorithm
    // The  following  Algorithm  reduces  an n × m matrix to  RREF
    // by  means  of  elementary  row  operations alone.

    // 1.  Preform Gaussian elimination to get the matrixin REF

    // 2.  For each non zero row i ( Ri ) from n to 1 (bottom to top) (a)
    //  For each row j ( Rj ) above row i (Each j < i )
    // Preform Rj → Rj − bRi where b =  the  value  in  row j directly
    // above  the pivot  in  row i .   [This  ensures  that  row j
    // has a zero above the pivot in row i ]

    for (int i = m - 1; i >= 0; i--) {
        int pivot = NA;
        for (int k = 0; k < n; k++) {
            if (!is_zero(mat[i][k])) {
                pivot = k;
                break;
            }
        }
        assert(pivot != NA);
        for (int j = 0; j < i; j++) {
            T b = mat[j][pivot];
            for (int k = 0; k < n; k++) {
                mat[j][k] -= b * mat[i][k];
            }
        }
    }

    fprintf(stderr, "%-30s m %d n %d\n", __func__, m , n);
}

template <typename T>
void separate_unit_mat_from_free_vars(std::vector<std::vector<T> > &mat,
                                      std::vector<int> &col_order)
{
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");
    assert(!mat.empty());
    int m = mat.size(), n = mat[0].size();

    if (col_order.empty()) {
        col_order = std::vector<int>(n - 1); // from augmented mat
        std::iota(col_order.begin(), col_order.end(), 0);
    }
    assert((int)col_order.size() == n - 1);

    for (int i = 0; i < m; i++) {
        if (is_zero(mat[i][i])) {
            int pivot = NA;
            for (int j = i + 1; j < n - 1; j++) {
                if (!is_zero(mat[i][j])) {
                    pivot = j;
                    break;
                }
            }
            assert(pivot != NA);
            std::swap(col_order[i], col_order[pivot]);
            for (int j = 0; j < m; j++) {
                std::swap(mat[j][i], mat[j][pivot]);
            }
        }
    }
}

template <typename T>
std::vector<std::vector<T> > get_free_vars_augm_matrix(const std::vector<std::vector<T> > &mat)
{
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");
    assert(!mat.empty());
    int m = mat.size(), n = mat[0].size();

    std::vector<std::vector<T> > mat2;
    // remove the m*m unit matrix,
    // (n - 1) - m columns + 2 last columns : interval of possible value
    int nn = n - m + 1; // augment augmented matrix such that
    // mat[row][n- 2] : lower bound,  mat[row][n - 1] : upper bound
    // and lower bound >= 0
    for (int i = 0; i < m; i++) {
        mat2.push_back(std::vector<T>(nn));
        T mul = (!is_pos(mat[i][n - 1]) ? -1 : +1);
        for (int j = m, k = 0; j < n - 1; j++, k++) {
            mat2[i][k] = mul * mat[i][j];
        }
        if (is_neg(mul)) {
            mat2[i][nn - 2] = mul * mat[i][n - 1];
            mat2[i][nn - 1] = mul * mat[i][n - 1] + 1;
        } else {
            mat2[i][nn - 2] = mul * mat[i][n - 1] - 1;
            mat2[i][nn - 1] = mul * mat[i][n - 1];
        }
    }

    // remove duplicate rule
    std::sort(mat2.begin(), mat2.end());
    auto it2 = std::unique(mat2.begin(), mat2.end(), [](const std::vector<T> &va, const std::vector<T> &vb) {
        return std::equal(va.begin(), va.end(), vb.begin(), [](const T &a, const T &b){
            return is_zero(T{std::abs(a - b)});
        });
    } );
    mat2.erase(it2, mat2.end());

    fprintf(stderr, "%-30s m %d n %d\n", __func__, (int)mat2.size() , (int)mat2[0].size());

    // remove constraint without additional info as var is either 0 or 1
    auto it = std::remove_if(mat2.begin(), mat2.end(), [](const std::vector<T> &row) {
        if (!is_zero(row[row.size() - 2])) {
            return false;
        }
        int pos = NA;
        for (int i = 0, iend = row.size() - 2; i < iend; i++) {
            if (!is_zero(row[i])) {
                if (pos == NA) {
                    pos = i;
                } else {
                    pos = NA;
                    break;
                }
            }
        }
        return pos != NA && is_zero(row[row.size() - 1] - row[pos]);
    });
    mat2.erase(it, mat2.end());

    fprintf(stderr, "%-30s m %d n %d\n", __func__, (int)mat2.size() , (int)mat2[0].size());

    return mat2;
}

template <typename T>
void print_matrix(const std::vector<std::vector<T> > &mat)
{
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");
    assert(!mat.empty());
    int m = mat.size(), n = mat[0].size();

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(stderr, "%+.2f ", (float)mat[i][j]);
        } fprintf(stderr, "\n");
    }

    fprintf(stderr, "---------------------------------------\n");
}

template <typename T>
std::vector<T> gomory(const std::vector<std::vector<T> > &mat)
{
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");
    assert(!mat.empty());
    int m = mat.size(), n = mat[0].size();

    std::vector<std::vector<T> > A(2 * m, std::vector<T>(n - 2));
    std::vector<T> B(2 * m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n - 2; j++) {
            A[m + i][j] = mat[i][j];
            B[m + i] = mat[i][n - 1];

            A[i][j] = -mat[i][j];
            B[i] = -mat[i][n - 2];
        }
    }

    std::vector<T> C(n - 2, +1), X;
    simplex_solver_bounded<T> ssolver;

    T res = ssolver.optimize(A, B, C, X, false, 25, 1);
    if (std::isnan(res)) { // if integer search failed
        res = ssolver.optimize(A, B, C, X, true, 25, 1);
    }
    if (!std::isfinite(res)) {
        X.clear();
    }
#ifndef NDEBUG
    else {
        for (int i = 0; i < m; i++) {
            T val = 0.;
            for (int j = 0; j < n - 2; j++) {
                val += X[j] * mat[i][j];
            }
            if (!is_zero(val - mat[i][n - 2]) && !is_zero(val - mat[i][n - 1]))
            {
                fprintf(stderr, "val %f vs [%f, %f] | res %f\n", (float)val, (float)mat[i][n - 2], (float)mat[i][n - 1], (float)res);
            }
        }
    }
#endif
    return X;
}

/**************************************************************************************************/

template <typename T>
void lp_solve(std::vector<std::vector<T> > &mat, std::vector<int> &solution, std::vector<int> &col_order)
{
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");

    fprintf(stderr, "mat size %d sol size %d\n", (int)mat[0].size(), (int)solution.size());

    gaussian_algo(mat);
    gaussian_jordan_algo(mat);
    separate_unit_mat_from_free_vars(mat, col_order);

    auto mat2 = get_free_vars_augm_matrix(mat);
    auto mat2_res = gomory(mat2);
    if (mat2_res.empty()) {
        return;
    }

    int m = mat.size();
    int n = mat[0].size() - 1;
    // copy free vars values of mat
    for (int i = 0; i < n - m; i++) {
        assert(is_zero(mat2_res[i] - 0) || is_zero(mat2_res[i] - 1));
        int is_true = std::nearbyint(mat2_res[i]);
        if (is_true) {
            solution.push_back(col_order[m + i]);
        }
    }

    // retrieve other var values from mat
    for (int i = 0; i < m; i++) {
        auto val = mat[i][n];
        for (int j = m; j < n; j++) {
            val -= mat2_res[j - m] * mat[i][j];
        }
        // TODO => check dist from nearbyint
        assert(is_zero(val - 0) || is_zero(val - 1));
        int is_true = std::nearbyint(val);
        if (is_true) {
            solution.push_back(col_order[i]);
        }
    }

    fprintf(stderr, "mat size %d sol size %d\n", (int)mat[0].size(), (int)solution.size());
}

/**************************************************************************************************/

int grid_char_to_int(const char c)
{
    if (D == 3) {
        return c >= '1' && c <= '9' ? c - '1' : NA ;
    } else if (D == 4) {
        return (c >= '0' && c <= '9') ? c - '0'
                                      : ((c >= 'A' && c <= 'F') ? c - 'A' + 10
                                                                : ((c >= 'a' && c <= 'f')
                                                                   ? c - 'a' + 10 : NA));
    }
    return NA;
}

char int_to_grid_char(const int n)
{
    if (D == 3) {
        return '1' + n;
    } else if (D == 4) {
        return n < 10 ? '0' + n : 'A' + n - 10;
    }
    return NA;
}

void lp_solver_get_grid(const std::vector<int> &solution, char grid_str[])
{
    for (int i = 0; i < NN; i++) {
        grid_str[i] = '.';
    }
    for (auto v : solution) {
        grid_str[v / N] = int_to_grid_char(v % N);
    }
    grid_str[NN] = '\0';
}

template <typename T>
std::vector<std::vector<T> > empty_grid_to_mat()
{
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");

    // 4 type of constraints, NN * N vars + augmented
    std::vector<std::vector<T> > mat(4 * NN, std::vector<T>(NN * N + 1, 0));
    int rhs_col = NN * N;

    int cur_row = 0;
    // one var true per cell
    for (int u = 0; u < N * NN; u += N) {
        for (int i = 0; i < N; i++) {
            mat[cur_row][u + i] = 1;
        }
        mat[cur_row++][rhs_col] = 1;
    }
    // one var true per row, col, box
    for (int cand = 0; cand < N; cand++) {
        for (int col = 0; col < N; col++) {
            for (int row = 0; row < N; row++) {
                int u = (row * N + col) * N + cand;
                mat[cur_row][u] = 1;
            }
            mat[cur_row++][rhs_col] = 1;
        }
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                int u = (row * N + col) * N + cand;
                mat[cur_row][u] = 1;
            }
            mat[cur_row++][rhs_col] = 1;
        }
        for (int box = 0; box < N; box++) {
            int col_beg = (box % D) * D;
            int row_beg = (box / D) * D;
            for (int idx = 0; idx < N; idx++) {
                int col = col_beg + idx % D;
                int row = row_beg + idx / D;
                int u = (row * N + col) * N + cand;
                mat[cur_row][u] = 1;
            }
            mat[cur_row++][rhs_col] = 1;
        }
    }

    return mat;
}

std::vector<int> grid_to_given_lits(char grid_str[])
{
    assert(strlen(grid_str) == NN);

    std::vector<int> lits;
    for (int i = 0; i < NN; i++) {
        int n = grid_char_to_int(grid_str[i]);
        if (n != NA) {
            int u = i * N + n;
            lits.push_back(u);
        }
    }
    return lits;
}

template <typename T>
std::vector<int> mat_populate(const std::vector<int> &given_lits,
                         std::vector<std::vector<T> > &mat,
                         std::vector<int> &solution)
{
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");

    // sudoku rules =>
    // one and only one var is true by constraint
    // once one var is set to true we can
    // - remove every constraint where it appears
    // - set every vars in these constraints to false
    std::vector<bool> to_del_rows(mat.size(), false);
    std::vector<bool> to_del_cols(mat[0].size(), false);
    auto add_to_del = [&to_del_cols, &to_del_rows, &mat](int u) {
        for (int i = 0, iend = mat.size(); i < iend; i++) {
            if (!mat[i][u]) {
                continue;
            }
            to_del_rows[i] = true;
            for (int j = 0, jend = mat[0].size() - 1; j < jend; j++) {
                if (mat[i][j]) {
                    to_del_cols[j] = true;
                }
            }
        }
    };

    // add given vars to solution and remove them from the matrix
    for (int var_id : given_lits) {
        add_to_del(var_id);
        solution.push_back(var_id);
    }

    // here we remove deleted row and column from the matrix
    // and keep track of the name of the vars associated with a column index
    std::vector<int> col_order;
    int k = 0;
    for (int i = 0, iend = mat.size(), jend = mat[0].size(); i < iend; i++) {
        if (to_del_rows[i]) {
            continue;
        }
        int l = 0;
        for (int j = 0; j < jend; j++) {
            if (to_del_cols[j]) {
                continue;
            }
            if (k == 0 && j != jend - 1) { // don't keep augmented col
                col_order.push_back(j);
            }
            if (i != k || j != l) {
                mat[k][l] = mat[i][j];
            }
            l++;
        }
        mat[k].resize(l);
        k++;
    }
    mat.resize(k);

    return col_order;
}

bool solver_solve(char grid_str[])
{
    if (strlen(grid_str) != NN) {
        return false;
    }

    static auto empty_mat = empty_grid_to_mat<long double>();
    static std::vector<std::vector<long double> > mat;
    mat = empty_mat;

    std::vector<int> solution;
    auto col_order = mat_populate(grid_to_given_lits(grid_str), mat, solution);

    lp_solver_get_grid(solution, grid_str);
    fprintf(stderr, "%s\n", grid_str);

    if (solution.size() != NN) {
        lp_solve(mat, solution, col_order);
    }

    lp_solver_get_grid(solution, grid_str);
    fprintf(stderr, "%s\n", grid_str);

    return solution.size() == NN;
}

/**************************************************************************************************/

int main()
{
    static_assert(D == 3 || D == 4, "Invalid D value");

    char grid_str[NN * N + 1] = "";
    int grid_cnt = 0, solved_grid_cnt = 0;
    auto start = clock();

    while (scanf(" %s", grid_str) == 1)
    {
        grid_cnt++;
        solved_grid_cnt += solver_solve(grid_str);

        if (grid_cnt && (grid_cnt % 2000 == 0)) {
            auto end = clock();
            uint64_t us = ((end - start)/(double)CLOCKS_PER_SEC) * 1000000;
            fprintf(stderr, "solved %d / %d %3.3f%% time grid % 3.3f us time total %ld us\n",
                    solved_grid_cnt, grid_cnt, 100.f * solved_grid_cnt / (grid_cnt == 0 ? 1.f : (float)grid_cnt),
                    (float)us / (float)(grid_cnt == 0 ? 1 : grid_cnt), us);
            fflush(stderr);
        }
    }

    auto end = clock();
    uint64_t us = ((end - start)/(double)CLOCKS_PER_SEC) * 1000000;

    fprintf(stderr, "solved %d / %d %3.3f%% time grid % 3.3f us time total %ld us\n",
            solved_grid_cnt, grid_cnt, 100.f * solved_grid_cnt / (grid_cnt == 0 ? 1.f : (float)grid_cnt),
            (float)us / (float)(grid_cnt == 0 ? 1 : grid_cnt), us);
    fflush(stderr);

    return 0;
}
