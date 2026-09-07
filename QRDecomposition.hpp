#pragma once

#include <stdexcept>
#include <utility>
#include <vector>

#include "Methods/GivensMethod.hpp"
#include "Methods/GramMethod.hpp"
#include "Methods/HouseholderMethod.hpp"

using namespace std;

/// Класс QRDecomposition
/// @details QR-разложение превращает произвольные столбцы матрицы в ортонормированную систему
/// Q и верхнетреугольную систему коэффициентов R, позволяя заменить исходную матрицу A на произведение двух
/// более удобных матриц, сводя сложную задачу с A в более легкую с QR.
class QRDecomposition
{
public:
    using Matrix = vector<vector<double>>;
    using Vector = vector<double>;

    enum class Method
    {
        /// Метод Грама-Шмидта
        /// @brief Решает задачу QR-разложения методом Грама-Шмидта.
        GramSchmidt,

        /// Метод Хаусхолдера
        /// @brief Решает задачу QR-разложения методом Хаусхолдера.
        Householder,

        /// Метод Гивенса
        /// @brief Решает задачу QR-разложения методом Гивенса.
        Givens
    };

private:
    Matrix A_;
    size_t m_, n_;
    double epsilon_;

    /// Проверяет исходную матрицу на ошибки
    static void checkMatrix(const Matrix& A)
    {
        if (A.empty() || A[0].empty()) throw invalid_argument("Матрица не должна быть пустой.");

        const size_t n = A[0].size();

        for (const auto& row : A)
            if (row.size() != n)throw invalid_argument("Матрица должна быть прямоугольной.");
    }

public:
    /// Конструктор класса QRDecomposition
    explicit QRDecomposition(const Matrix& A, const double epsilon) : A_(A), m_(A.size()),
                                                                      n_(A.empty() ? 0 : A[0].size()), epsilon_(epsilon)
    {
        checkMatrix(A_);

        if (m_ < n_) throw invalid_argument("Ожидается m >= n.");

        if (epsilon_ <= 0) throw invalid_argument("Эпсилон должен быть положительным.");
    }

    /// Позволяет выбрать метод QR-разложения
    [[nodiscard]] pair<Matrix, Matrix> decompose(const Method method) const
    {
        switch (method)
        {
        case Method::GramSchmidt:
            return GramMethod::decompose(A_, epsilon_);

        case Method::Householder:
            return HouseholderMethod::decompose(A_, epsilon_);

        case Method::Givens:
            return GivensMethod::decompose(A_, epsilon_);
        }

        throw invalid_argument("Неизвестный QR-метод.");
    }

    [[nodiscard]] const Matrix& originalMatrix() const { return A_; }

    [[nodiscard]] size_t rows() const { return m_; }

    [[nodiscard]] size_t cols() const { return n_; }
};