#pragma once

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

/// Класс для создания и вывода таблицы с текстовыми данными и матрицами
class Table
{
public:
    using Matrix = vector<vector<double>>;

private:
    struct Section
    {
        string label;
        vector<Matrix> matrices;
    };

    vector<string> columns;
    vector<Section> sections;

    /// Вычисляет длину строки в UTF-8 символах
    static size_t utf8Length(const string& str)
    {
        size_t count = 0;

        for (const unsigned char c : str)
            if ((c & 0xC0) != 0x80) ++count;

        return count;
    }

    /// Добавляет пробелы справа до указанной ширины
    static string padRight(const string& str, size_t width)
    {
        const size_t length = utf8Length(str);

        if (length >= width) return str;

        return str + string(width - length, ' ');
    }

    /// Форматирует матрицу для вывода в таблицу
    static vector<string> formatMatrix(const Matrix& matrix, const int precision = 6)
    {
        vector<string> result;

        ostringstream oss;
        oss << fixed << setprecision(precision);

        for (const auto& row : matrix)
        {
            oss.str("");
            oss.clear();

            oss << "[ ";

            for (const double x : row) oss << setw(12) << x << ' ';

            oss << "]";

            result.push_back(oss.str());
        }

        return result;
    }

    /// Вычисляет ширину столбца
    [[nodiscard]] size_t columnWidth(size_t index) const
    {
        size_t width = utf8Length(columns[index]);

        // Минимальная ширина для первой колонки
        if (index == 0)
            width = max(width, static_cast<size_t>(7));

        // Остальные колонки
        for (const auto& [label, matrices] : sections)
        {
            if (index == 0)
                width = max(width, utf8Length(label));
            else
            {
                if (const size_t matrixIndex = index - 1; matrixIndex < matrices.size())
                {
                    const auto lines =
                        formatMatrix(matrices[matrixIndex]);

                    for (const auto& line : lines)
                        width = max(width, line.length());
                }
            }
        }

        // Отступ внутри ячейки
        return width + 2;
    }

    /// Выводит горизонтальную линию таблицы
    void horizontalLine() const
    {
        cout << '+';

        for (size_t i = 0; i < columns.size(); ++i)
            cout << string(columnWidth(i), '-') << '+';

        cout << '\n';
    }

public:
    /// Добавляет колонку
    Table& addColumn(const string& name)
    {
        columns.push_back(name);
        return *this;
    }

#pragma region СОЗДАНИЕ СЕКЦИИ

    class SectionBuilder
    {
    private:
        Table& table;
        size_t sectionIndex;

    public:
        /// Конструктор построителя секции
        SectionBuilder(Table& table, const size_t sectionIndex) : table(table), sectionIndex(sectionIndex) { }

        /// Конструктор построителя секции
        SectionBuilder& text(const string& text)
        {
            table.sections[sectionIndex].label = text;
            return *this;
        }

        /// Добавляет матрицу в секцию
        SectionBuilder& matrix(const Matrix& matrix)
        {
            auto& matrices =
                table.sections[sectionIndex].matrices;

            if (matrices.size() + 1 >= table.columns.size())
                throw runtime_error("Недостаточно свободных колонок для матрицы.");

            matrices.push_back(matrix);

            return *this;
        }
    };

    /// Создаёт новую секцию
    SectionBuilder section(const string& name)
    {
        sections.push_back({name, {}});

        return SectionBuilder(*this, sections.size() - 1);
    }

#pragma endregion

    /// Выводит таблицу
    void print() const
    {
        if (columns.empty()) return;

        horizontalLine();

        cout << '|';

        for (size_t i = 0; i < columns.size(); ++i)
        {
            const size_t width = columnWidth(i);

            cout << ' ' << padRight(columns[i], width - 2) << ' ' << '|';
        }

        cout << '\n';

        horizontalLine();


        for (const auto& [label, matrices] : sections)
        {
            vector<vector<string>> matrixLines;

            size_t maxRows = 1;

            for (const auto& matrix : matrices)
            {
                auto lines = formatMatrix(matrix);

                maxRows = max(maxRows, lines.size());

                matrixLines.push_back(std::move(lines));
            }

            for (size_t row = 0; row < maxRows; ++row)
            {
                cout << '|';

                {
                    string value;

                    if (row == 0) value = label;

                    const size_t width = columnWidth(0);

                    cout << ' ' << padRight(value, width - 2) << ' ' << '|';
                }

                for (size_t column = 1; column < columns.size(); ++column)
                {
                    const size_t matrixIndex = column - 1;

                    string value;

                    if (matrixIndex < matrixLines.size() && row < matrixLines[matrixIndex].size())
                        value = matrixLines[matrixIndex][row];

                    const size_t width = columnWidth(column);

                    cout << ' ' << padRight(value, width - 2) << ' ' << '|';
                }

                cout << '\n';
            }

            horizontalLine();
        }
    }
};