#pragma once

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace std;

class GramMethod
{
public:
    using Matrix = vector<vector<double>>;

private:
    static double dot(const vector<double>& a, const vector<double>& b)
    {
        double result(0);

        for (size_t i = 0; i < a.size(); ++i)
            result += a[i] * b[i];

        return result;
    }

    static double norm(const vector<double>& v) { return sqrt(dot(v, v)); }

public:
    static pair<Matrix, Matrix> decompose(const Matrix& A, const double epsilon)
    {
        const size_t m = A.size();
        const size_t n = A[0].size();

        Matrix Q(m, vector<double>(n, 0));
        Matrix R(n, vector<double>(n, 0));

        for (size_t k = 0; k < n; ++k)
        {
            vector<double> u(m);

            for (size_t i = 0; i < m; ++i)
                u[i] = A[i][k];

            for (size_t j = 0; j < k; ++j)
            {
                double r = 0;

                for (size_t i = 0; i < m; ++i)
                    r += Q[i][j] * A[i][k];

                R[j][k] = r;

                for (size_t i = 0; i < m; ++i)
                    u[i] -= r * Q[i][j];
            }

            R[k][k] = norm(u);

            if (R[k][k] <= epsilon)
                throw runtime_error("Грама-Шмидта прерван: столбцы линейно зависимы.");

            for (size_t i = 0; i < m; ++i)
                Q[i][k] = u[i] / R[k][k];
        }

        return {Q, R};
    }
};