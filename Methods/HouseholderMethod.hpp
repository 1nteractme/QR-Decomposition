#pragma once

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace std;

class HouseholderMethod
{
public:
    using Matrix = vector<vector<double>>;

private:
    static double norm(const vector<double>& v)
    {
        double result = 0;

        for (const double x : v)
            result += x * x;

        return sqrt(result);
    }

    static void applyHouseholder(Matrix& R, Matrix& Qt, const size_t k, const double epsilon)
    {
        const size_t m = R.size();
        const size_t n = R[0].size();

        vector<double> v(m - k);

        for (size_t i = k; i < m; ++i)
            v[i - k] = R[i][k];

        const double xNorm = norm(v);

        if (xNorm <= epsilon) return;

        const double alpha = v[0] >= 0 ? -xNorm : xNorm;

        v[0] -= alpha;

        const double vNorm = norm(v);

        if (vNorm <= epsilon) return;

        for (double& x : v)
            x /= vNorm;

        // R = H * R
        for (size_t j = k; j < n; ++j)
        {
            double projection = 0;

            for (size_t i = k; i < m; ++i)
                projection += v[i - k] * R[i][j];

            for (size_t i = k; i < m; ++i)
                R[i][j] -= 2 * v[i - k] * projection;
        }

        // Qt = H * Qt
        for (size_t j = 0; j < m; ++j)
        {
            double projection = 0;

            for (size_t i = k; i < m; ++i)
                projection += v[i - k] * Qt[i][j];

            for (size_t i = k; i < m; ++i)
                Qt[i][j] -= 2 * v[i - k] * projection;
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

        for (size_t k = 0; k < n; ++k)
            applyHouseholder(R, Qt, k, epsilon);

        // Q = Qt^T
        Matrix Q(m, vector<double>(n, 0));

        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                Q[i][j] = Qt[j][i];

        Matrix Rthin(n, vector<double>(n, 0));

        for (size_t i = 0; i < n; ++i)
            for (size_t j = i; j < n; ++j)
                Rthin[i][j] = R[i][j];

        for (size_t i = 0; i < n; ++i)
        {
            if (abs(Rthin[i][i]) <= epsilon)
                throw runtime_error("Метод Хаусхолдера выявил проблему ранга.");
        }

        return {Q, Rthin};
    }
};