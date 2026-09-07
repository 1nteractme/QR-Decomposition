#include <sstream>

#include "Data.hpp"
#include "QRDecomposition.hpp"
#include "Table.hpp"

using namespace std;

int main()
{
#pragma region ИНИЦИАЛИЗАЦИЯ МАТРИЦЫ

    constexpr double EPSILON(1e-15);
    const Data DATA(a);
    const auto& MATRIX = DATA.getMatrix();
    const QRDecomposition QR(MATRIX, EPSILON);

#pragma endregion

#pragma region ИНИЦИАЛИЗАЦИЯ МЕТОДОВ

    auto [Qgs, Rgs] = QR.decompose(QRDecomposition::Method::GramSchmidt); // ГРАММ-ШМИДТ
    auto [Qgv, Rgv] = QR.decompose(QRDecomposition::Method::Givens); // ГИВЕНС
    auto [Qhh, Rhh] = QR.decompose(QRDecomposition::Method::Householder); // ХАУСХОЛДЕР

#pragma endregion

#pragma region ВЫВОД ДАННЫХ

    Table table;

    table
        .addColumn("МАТРИЦА")
        .addColumn("ГРАММ-ШМИДТ")
        .addColumn("ГИВЕНС")
        .addColumn("ХАУСХОЛДЕР");

    table.section("Q")
         .text("Q")
         .matrix(Qgs)
         .matrix(Qgv)
         .matrix(Qhh);

    table.section("R")
         .text("R")
         .matrix(Rgs)
         .matrix(Rgv)
         .matrix(Rhh);

    table.print();

#pragma endregion

    return 0;
}