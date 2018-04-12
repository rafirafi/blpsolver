#ifndef SIMPLEX_SOLVER_BOUNDED_H
#define SIMPLEX_SOLVER_BOUNDED_H

// Copyright rafirafi 2018, modified from :
//
// https://github.com/dieram3/competitive-programming-library/blob/master/include/cpl
//          Copyright Diego Ramirez 2015
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <cassert>
#include <cmath>

#ifndef NDEBUG
#include <cstdio>
#endif

#include <algorithm>
#include <limits>
#include <vector>

template <typename T>
class simplex_solver_bounded {
    static_assert(std::is_floating_point<T>::value, "Must be floating-point");

private:
    std::vector<std::vector<T> > tableau_objs_; // Part of the tableau with objectives functions
    std::vector<std::vector<T> > tableau_;  // Part of the tableau minus objectives functions, last col is constraint goal
    std::vector<int> basic_names_, not_basic_names_; // Basic and nonbasic variables

    std::vector<T> ub_;
    std::vector<bool> is_at_ub_;

    static const constexpr int kNA = -1;
    static const constexpr T epsilon = 1e-9; // use one epsilon only, +/- ok for primal/dual simplex and gomory

    bool is_neg(const T &val) const { return val < -epsilon; }
    bool is_pos(const T &val) const { return val >  epsilon; }
    bool is_zero(const T &val) const { return std::abs(val) <= epsilon; }
    bool approx(const T &v1, const T &v2) const { return std::abs(v1 - v2) <= epsilon; }

public:
    // input: Ax <= B , x >= 0
    // if 'minimize' false, maximize
    // if 'integer_search' true => search an integer solution, return nan if doesn't succeed
    //  -1 don't search integer solution, 0 no limit, 1-n max nb of gomory cuts added
    T optimize(const std::vector<std::vector<T> > &A,
               const std::vector<T> &B,
               const std::vector<T> &C,
               std::vector<T> &x,
               bool minimize = false,
               int integer_search_limit = -1,
               T global_upper_bound = std::numeric_limits<T>::infinity())
    {
        assert(A.size() == B.size());
        for (auto &row : A) { assert(row.size() == C.size()); }

        x.clear();

        T minimize_adjust = (minimize ? -1 : +1);

        int m = A.size(), n = A[0].size();
        tableau_.assign(m, std::vector<T>(n + 2, 0)); // + 2 : articial variable + Cj
        tableau_objs_.assign(2, std::vector<T>(n + 2, 0));

        // problem : phase 2 at 0
        for (int j = 0; j < n; j++) {
            if (!is_zero(C[j])) {
                tableau_objs_[0][j] = minimize_adjust * -C[j];
            }
        }

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (!is_zero(A[i][j])) {
                    tableau_[i][j] = A[i][j];
                }
            }
            if (!is_zero(B[i])) {
                tableau_[i][n + 1] = B[i];
            }
        }

        not_basic_names_.resize(n + 1);
        std::iota(not_basic_names_.begin(), not_basic_names_.end(), 0);
        basic_names_.resize(m);
        std::iota(basic_names_.begin(), basic_names_.end(), n + 1);

        is_at_ub_.assign(n + 1 + m, false);

        // NOTE: easy to transform to per var upper bound if needed
        ub_.assign(n + 1 + m, std::numeric_limits<T>::infinity());
        for (int j = 0; j < n; j++) {
            ub_[j] = global_upper_bound;
        }

        // artificial variable setup
        int min_b = kNA;
        for (int i = 0; i < m; i++) {
            if (is_neg(tableau_[i][n + 1])) {
                if (min_b == kNA || B[i] < B[min_b]) {
                    min_b = i;
                }
                tableau_[i][n] = -1;
            }
        }

        // phase I
        if (min_b != kNA) {
            tableau_objs_[1][n] = 1; // problem : phase 1 at 1

            pivot(min_b, n); // Pivot to make RHS positive.

            simplex(1); // Find optimal solution to phase I problem
            if (!is_zero(tableau_objs_[1][n + 1])) {
#ifndef NDEBUG
                fprintf(stderr, "%s Phase I failed val %f\n", __func__, (float)tableau_objs_[1][n + 1]);
#endif
                return std::numeric_limits<T>::quiet_NaN(); // Infeasible.
            }
            // Make 'n' a nonbasic variable if necessary
            const auto it = std::find(basic_names_.begin(), basic_names_.end(), n);
            if (it != basic_names_.end()) {
                int row = it - basic_names_.begin();
                assert(is_zero(tableau_[row][n + 1]));
                int col;
                for (col = 0; col <= n; ++col) {
                    if (!is_zero(tableau_[row][col])) {
                        break;
                    }
                }
                assert(col <= n);
                pivot(row, col);
            }
            // Nullify the artificial column.
            int art_pos = std::find(not_basic_names_.begin(), not_basic_names_.end(), n) - not_basic_names_.begin();

            assert(not_basic_names_.size() - art_pos > 0);
            for (int i = 0, iend = tableau_.size(); i < iend; i++) {
                tableau_[i][art_pos] = 0;
            }
            tableau_objs_[0][art_pos] = 0;
            // remove phase I problem
            tableau_objs_.resize(1);
        }

        // phase II
        if (!simplex(2)) {
#ifndef NDEBUG
            fprintf(stderr, "%s Phase II failed\n", __func__);
#endif
            return std::numeric_limits<T>::infinity(); // Unbounded.
        }

        if (integer_search_limit >= 0) {
            if (!gomory(integer_search_limit)) {
                x.clear();
                return std::numeric_limits<T>::quiet_NaN();
            }
            m = tableau_.size();
        }

        x.assign(n, 0);
        fprintf(stderr, "res %f\n", (float)(minimize_adjust * tableau_objs_[0][n + 1]));
        for (int i = 0; i < m; i++) {
            if (basic_names_[i] < n) {
                x[basic_names_[i]] = tableau_[i][n + 1];
                if (is_at_ub_[basic_names_[i]]) {
                    x[basic_names_[i]] = ub_[basic_names_[i]] - x[basic_names_[i]];
                }
                fprintf(stderr, "%f ", (float)x[basic_names_[i]]);
            }
        }
        for (int i = 0; i < n + 1; i++) {
            if (not_basic_names_[i] < n) {
                if (is_at_ub_[not_basic_names_[i]]) {
                    x[not_basic_names_[i]] = ub_[not_basic_names_[i]];
                    fprintf(stderr, "%f ", (float)x[not_basic_names_[i]]);
                }
            }
        }
        fprintf(stderr, "\n");

        // sanitize
        if (integer_search_limit >= 0) {
            bool ok = true;
            for (auto &val : x) {
                // TODO : use std::nearbyint instead
                T int_part;
                std::modf(val + 2 * epsilon, &int_part);
                if (approx(val, int_part)) {
                    val = int_part;
                } else {
                    ok = false;
                }
            }
            if (!ok) {
#ifndef NDEBUG
                fprintf(stderr, "%s failed, no integer result\n", __func__);
#endif
                return std::numeric_limits<T>::quiet_NaN();
            }
        }

        return minimize_adjust * tableau_objs_[0][n + 1];
    }

private:

    void multiply_row(int i, T mult)
    {
        for (int j = 0, jend = tableau_[0].size(); j < jend; j++) {
            tableau_[i][j] *= mult;
        }
    }

    void add_row(int row_from, T mult, int row_to)
    {
        assert(row_from != row_to);
        for (int j = 0, jend = tableau_[0].size(); j < jend; j++) {
            tableau_[row_to][j] += mult * tableau_[row_from][j];
        }
    }

    void add_row_to_objs(int row_from, T mult, int row_to)
    {
        for (int j = 0, jend = tableau_[0].size(); j < jend; j++) {
            tableau_objs_[row_to][j] += mult * tableau_[row_from][j];
        }
    }

    void pivot(int r, int c)
    {
        static std::vector<T> pivcol;
        int m = tableau_.size();
        pivcol.resize(m + tableau_objs_.size());

        for (int i = 0, iend = pivcol.size(); i < iend; i++) {
            if (i < m) {
                pivcol[i] = tableau_[i][c];
                tableau_[i][c] = 0;
            } else {
                pivcol[i] = tableau_objs_[i - m][c];
                tableau_objs_[i - m][c] = 0;
            }
        }
        tableau_[r][c] = 1;
        multiply_row(r, T{1} / pivcol[r]);
        for (int i = 0, iend = pivcol.size(); i < iend; ++i) {
            if (i == r) {
                continue;
            }
            if (i < m) {
                add_row(r, -pivcol[i], i);
            } else {
                add_row_to_objs(r, -pivcol[i], i - m);
            }
        }

        std::swap(basic_names_[r], not_basic_names_[c]);
    }

    /*******************************************************/

    // http://faculty.ndhu.edu.tw/~ywan/courses/network/notes/bounded_variable_new.pdf
    // http://web.mit.edu/15.053/www/AMP-Chapter-02.pdf
    void substitute_upper_bound(int r, int c)
    {
        assert(r != c && (r == kNA || c == kNA));
        int rhs_col = tableau_[0].size() - 1;
        if (c != kNA) {
            for (int i = 0, iend = tableau_.size(); i < iend; i++) {
                tableau_[i][rhs_col] -= tableau_[i][c] * ub_[not_basic_names_[c]];
                tableau_[i][c] = -tableau_[i][c];
            }
            for (int i = 0, iend = tableau_objs_.size(); i < iend; i++) {
                tableau_objs_[i][rhs_col] -= tableau_objs_[i][c] * ub_[not_basic_names_[c]];
                tableau_objs_[i][c] = -tableau_objs_[i][c];
            }
            is_at_ub_[not_basic_names_[c]] = !is_at_ub_[not_basic_names_[c]];
        } else {
            tableau_[r][rhs_col] -= ub_[basic_names_[r]];
            multiply_row(r, -1);
            is_at_ub_[basic_names_[r]] = !is_at_ub_[basic_names_[r]];
        }
    }

    /*******************************************************/

    enum search_result { must_pivot, must_flip_only, must_flip_and_pivot, optimized, unbounded };

    // Requires: A feasible solution must exist.
    bool simplex(int phase)
    {
        int r = 0, c = 0;
        search_result res;
        int max_loop_cnt = 10000, loop_cnt = 0; // arbitrary limit

        while (loop_cnt++ != max_loop_cnt && (res = simplex_find_pivot(r, c, phase)) != unbounded && res != optimized) {
            if (res == must_pivot) {
                pivot(r, c);
            } else if (res == must_flip_only) {
                substitute_upper_bound(kNA, c);
            } else if (res == must_flip_and_pivot) {
                substitute_upper_bound(r, kNA);
                pivot(r, c);
            }
        }
#ifndef NDEBUG
        fprintf(stderr, "primal pivot cnt %d\n", loop_cnt);
#endif
        assert(!(phase == 1 && res == unbounded)); // Phase 1 can't be unbounded.
        return res == optimized;
    }

    // Bland's rule pivot selecting.
    search_result simplex_find_pivot(int& r, int& c, int phase)
    {
        int m = tableau_.size();
        int n = tableau_[0].size();
        int rhs_col = n - 1;
        int objrow = (phase == 1 ? 1 : 0);

        // Find entering variable:
        c = kNA;
        for (int j = 0; j < rhs_col; j++) {
            if (!is_neg(tableau_objs_[objrow][j])) {
                continue;
            }
            if (c == kNA || not_basic_names_[j] < not_basic_names_[c]) {
                c = j;
            }
        }
        if (c == kNA) {
            return optimized;
        }

        // Find leaving variable:
        auto ratio_less = [&](const size_t r1, const T ratio1, const size_t r2, const T ratio2) {
            if (approx(ratio1, ratio2)) {
                return basic_names_[r1] < basic_names_[r2];
            }
            return ratio1 < ratio2;
        };

        r = kNA;
        T min_ratio{0}; // init to wathever
        search_result res = unbounded;
        for (int i = 0; i < m; i++) {
            if (is_pos(tableau_[i][c])) {
                T ratio = tableau_[i][rhs_col] / tableau_[i][c];
                if (r == kNA || ratio_less(i, ratio, r, min_ratio)) {
                    r = i;
                    min_ratio = ratio;
                    res = must_pivot;
                }
            }
            // upper bound
            else if (is_neg(tableau_[i][c])) {
                if (std::isfinite(ub_[basic_names_[i]])) {
                    T ratio = -(ub_[basic_names_[i]] - tableau_[i][rhs_col]) / tableau_[i][c];
                    if (r == kNA || ratio_less(i, ratio, r, min_ratio)) {
                        r = i;
                        min_ratio = ratio;
                        res = must_flip_and_pivot;
                    }
                }
            }
        }
        if (res != unbounded && ub_[not_basic_names_[c]] < min_ratio - epsilon) { // TODO => is_less(T a, T b)
            res = must_flip_only;
        }

        return res;
    }

    /*******************************************************/

    bool dual()
    {
        int r = 0, c = 0;
        search_result res = unbounded;
        int max_loop_cnt = 10000, loop_cnt = 0; // arbitrary limit
        while (loop_cnt++ != max_loop_cnt && (res = find_dual_pivot(r, c)) == must_pivot) {
            pivot(r, c);
        }
#ifndef NDEBUG
        fprintf(stderr, "dual pivot cnt %d\n", loop_cnt);
#endif
        return res == optimized;
    }

    // not really sure these pivot selection rules protect us from cycling
    search_result find_dual_pivot(int &r, int &c)
    {
        int m = tableau_.size();
        int n = tableau_[0].size();
        int objrow = 0;

        // Find leaving variable
        // also any oob variable is taken care of here
        r = kNA;
        for (int i = 0; i < m; i++) {
            // if basic var is > ub, flip it
            if (std::isfinite(ub_[basic_names_[i]]) && tableau_[i][n - 1] > ub_[basic_names_[i]] + epsilon) { // TODO => is_less(T a, T b)
                substitute_upper_bound(i, kNA);
            }

            if (!is_neg(tableau_[i][n - 1])) {
                continue;
            }
            if (r == kNA || basic_names_[i] < basic_names_[r]) {
                r = i;
            }
        }

        if (r == kNA) {
            return optimized;
        }

        // find entering variable:
        auto ratio_less = [&](const size_t c1, const size_t c2) {
            const auto ratio1 = tableau_objs_[objrow][c1] / tableau_[r][c1];
            const auto ratio2 = tableau_objs_[objrow][c2] / tableau_[r][c2];
            if (approx(ratio1, ratio2))
                return not_basic_names_[c1] < not_basic_names_[c2];
            return ratio1 > ratio2;
        };
        c = kNA;
        for (int j = 0; j < n - 1; j++) {
            if (!is_neg(tableau_[r][j])) {
                continue;
            }
            if (c == kNA || (ratio_less(j, c))) {
                c = j;
            }
        }

        return c == kNA ? unbounded : must_pivot;
    }

    /*******************************************************/

    T fractionnal_part(const T &value)
    {
        if (!is_zero(value)) {
            static T dummy;
            T res = std::modf(value, &dummy);
            if (is_neg(res)) {
                res = T{1} + res;
            }
            if (is_pos(res) && is_pos(T{1} - res)) {
                return res;
            }
        }
        return T{0};
    }

    int find_cut_row()
    {
        int m = tableau_.size();
        int n = tableau_[0].size();
        int objcol = n - 1;

        int r = kNA;
        T max = T{0};
        for (int i = 0; i < m; i++) {
            T f = fractionnal_part(tableau_[i][objcol]);
            if (!is_zero(f) && (r == kNA || (!approx(f, max) && f > max))) {
                r = i;
                max = f;
            }
        }
        return r;
    }

    bool add_frac_constraint()
    {
        int r = find_cut_row();
        if (r == kNA) {
            return false;
        }

        int m = tableau_.size();
        tableau_.resize(m + 1);
        int n = tableau_[0].size();
        tableau_[m].reserve(n);

        // raw gomory cuts : failure rate is not low.
        for (int j = 0; j < n; j++) {
            tableau_[m].push_back(-fractionnal_part(tableau_[r][j]));
        }

        // add the basic var associated to the added constraint
        basic_names_.push_back(basic_names_.size() + not_basic_names_.size());
        ub_.push_back(std::numeric_limits<T>::infinity());
        is_at_ub_.push_back(false);

        return true;
    }

    bool gomory(int limit)
    {
        int cnt = 0;
        while ((!limit || cnt++ < limit) && add_frac_constraint()) {
            if (!dual()) {
#ifndef NDEBUG
                fprintf(stderr, "%s dual failed at cnt %d\n", __func__, cnt);
#endif
                return false;
            }
        }
#ifndef NDEBUG
                fprintf(stderr, "%s cnt %d\n", __func__, cnt);
#endif
                return cnt <= limit;

    }

};

#endif // SIMPLEX_SOLVER_BOUNDED_H
