#pragma once

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace std;

class GivensMethod
{
public:
    using Matrix = vector<vector<double>>;

private:
    static pair<double, double> givens(const double a, const double b)
    {
        if (b == 0) return {1, 0};
        if (a == 0) return {0, 1};

        const double r = hypot(a, b);

        return {a / r, b / r};
    }

    static void applyGivens(Matrix& R, Matrix& Qt,
                            const size_t i, const size_t j,
                            const double c, const double s)
    {
        for (size_t col = 0; col < R[0].size(); ++col)
        {
            const double a = R[i][col];
            const double b = R[j][col];

            R[i][col] = c * a + s * b;
            R[j][col] = -s * a + c * b;
        }

        for (size_t col = 0; col < Qt[0].size(); ++col)
        {
            const double a = Qt[i][col];
            const double b = Qt[j][col];

            Qt[i][col] = c * a + s * b;
            Qt[j][col] = -s * a + c * b;
        }
    }

public:
    static pair<Matrix, Matrix> decompose(const Matrix& A, const double epsilon)
    {
        const size_t m = A.size();
        const size_t n = A[0].size();

        Matrix R = A;
        Matrix Qt(m, vector<double>(m, 0));

        for (size_t i = 0; i < m; ++i)
            Qt[i][i] = 1;

        for (size_t col = 0; col < n; ++col)
        {
            for (size_t row = m; row-- > col + 1;)
            {
                const double b = R[row][col];

                if (abs(b) <= epsilon) continue;

                const auto [c, s] = givens(R[col][col], b);

                applyGivens(R, Qt, col, row, c, s);
            }
        }

        Matrix Qfull(m, vector<double>(m, 0));

        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < m; ++j)
                Qfull[i][j] = Qt[j][i];

        Matrix Q(m, vector<double>(n, 0));
        Matrix Rthin(n, vector<double>(n, 0));

        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                Q[i][j] = Qfull[i][j];

        for (size_t i = 0; i < n; ++i)
            for (size_t j = i; j < n; ++j)
                Rthin[i][j] = R[i][j];

        for (size_t i = 0; i < n; ++i)
        {
            if (abs(Rthin[i][i]) <= epsilon)
                throw runtime_error("Метод Гивенса выявил проблему ранга.");
        }

        return {Q, Rthin};
    }
};