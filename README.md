# QR-разложение

## Введение

Реализация QR-разложения вещественной матрицы на языке **C++** тремя классическими методами ортогонализации:

- **Грамма–Шмидта**
- **Отражениями Хаусхолдера**
- **Вращениями Гивенса**.

Проект демонстрирует, как одна и та же факторизация $A = QR$ может быть получена принципиально разными вычислительными
алгоритмами.

**Репозиторий**: [1nteractme/QR-Decomposition](https://github.com/1nteractme/QR-Decomposition)

---

## Содержание

- [1. Постановка задачи](#1-постановка-задачи)

- [2. Математическая основа](#2-математическая-основа)
    - [2.1. QR-разложение](#21-qr-разложение)
    - [2.2. Геометрический смысл](#22-геометрический-смысл)
    - [2.3. Почему R верхнетреугольная](#23-почему-r-верхнетреугольная)

- [3. Реализованные методы](#3-реализованные-методы)
    - [3.1. Грамм–Шмидт](#31-граммшмидт)
    - [3.2. Хаусхолдер](#32-хаусхолдер)
    - [3.3. Гивенс](#33-гивенс)

- [4. Архитектура проекта](#4-архитектура-проекта)

- [5. Описание исходного кода](#5-описание-исходного-кода)
    - [5.1. QRDecomposition.hpp](#51qrdecompositionhpp)
    - [5.2. GramMethod.hpp](#52-grammethodhpp)
    - [5.3. HouseholderMethod.hpp](#53-householdermethodhpp)
    - [5.4. GivensMethod.hpp](#54-givensmethodhpp)
    - [5.5. Data.hpp](#55-datahpp)
    - [5.6. Table.hpp](#56-tablehpp)
    - [5.7. main.c++](#57-mainc++)
    - [5.8. CMakeLists.txt](#58-cmakeliststxt)

- [6. API](#6-api)

- [7. Роль `epsilon`](#7-рольepsilon)

- [8. Сравнение методов](#8-сравнение-методов)

- [9. Тестовые данные](#9-тестовые-данные)

- [10. Сборка и запуск](#10-сборка-и-запуск)

- [11. Результат работы](#11-результат-работы)

- [12. Ограничения реализации](#12-ограничения-реализации)

- [Заключение](#заключение)

---

# 1. Постановка задачи

Пусть задана вещественная матрица $A\in Rm×n,\dots m>=n$.
Цель QR-разложения – представить её в виде произведения $A = QR$, где:

- $Q\in\mathbb{R}^{m\times n}$ – матрица с ортонормированными столбцами.
- $R\in\mathbb{R}^{n\times n}$ – верхнетреугольная матрица.

Для $Q$ выполняется $Q^TQ = In$. Следовательно, столбцы $Q$ удовлетворяют $qiTqj={1, i=j, 0, i≠j}$.

Проект реализует 3 способа получения такой факторизации:

| Метод       | Основная идея                                 |
|-------------|-----------------------------------------------|
| Грамм–Шмидт | Последовательное вычитание проекций           |
| Хаусхолдер  | Отражение вектора относительно гиперплоскости |
| Гивенс      | Поворот в двумерной координатной плоскости    |

Таким образом, математическая задача одна, но вычислительные подходы различаются.

---

# 2. Математическая основа

## 2.1. QR-разложение

Пусть $A=[a_1, a_2,\dots, a_n]$.

QR-разложение представляет каждый исходный столбец как линейную комбинацию ортонормированных
векторов: $a_j=r_1jq_1+r_2jq_2+\dots+r_jjq_j$.

Матрица $Q$ задаёт новую ортонормированную систему направлений, а $R$ хранит координаты исходных столбцов в этой
системе.

> Важно отличать **QR-разложение** от **QR-алгоритма**.
> QR-разложение – факторизация одной матрицы $A=QR$, тогда как QR-алгоритм является итерационным методом поиска
> собственных значений: $Ak=QkRk,\ldots Ak+1=RkQk$.

---

## 2.2. Геометрический смысл

Исходные столбцы $a_1,a_2,\dots ,a_n$ могут иметь произвольные длины и направления и не обязаны быть взаимно
перпендикулярными.

QR-разложение заменяет их системой $q_1, q_2, \dots, q_n$, которая:

1. Состоит из единичных векторов
2. Является взаимно ортогональной
3. Порождает то же линейное подпространство

Иными словами, происходит не изменение подпространства, а изменение системы координат внутри него.

Для двух векторов первый ортонормированный вектор имеет вид $q_1=a_1||a_1||$.

Для второго сначала удаляется компонента вдоль $q_1$: $u_2=a_2−(q_1Ta_2)q_1$, после чего выполняется
нормировка: $q_2=u_2||u_2||$.

В результате $q_1Tq_2=0,\qquad ||q1||=||q_2||=1$.

Обобщение этого процесса и лежит в основе метода Грамма–Шмидта.

---

## 2.3. Почему $R$ верхнетреугольная

При построении $q_k$ используются только первые $k$ исходных направлений.

Для 1 столбца: $a_1=r_{11}q_1$.
Для 2: $a_2=r_{12}q_1 + r_{22}q_2$.
Для 3: $a_3=r_{13}q_1 + r_{23}q_2 + r_{33}q_3$.
В общем случае: $a_i=\sum_{j=1}^{i}r_{ji}q_j$.

Следовательно, $r_{ij}=0,\qquad i>j$.

Именно поэтому $R$ является верхнетреугольной.

---

# 3. Реализованные методы

## 3.1. Грамм–Шмидт

Классический процесс Грамма–Шмидта последовательно ортонормирует столбцы матрицы.

Для $k$-го столбца:

$$
\boxed{
u_k =
a_k-\sum_{j=1}^{k-1}
(q_j^*a_k)q_j
}
$$

где $q_k=\frac{u_k}{\|u_k\|}$

Коэффициенты матрицы $R$: $r_{jk}=q_j^Ta_k,\qquad j<k $, а диагональный элемент: $ r_{kk}=||u_k||$.

---

### Реализация

Основная функция находится в [`Methods/GramMethod.hpp`](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GramMethod.hpp#L365-L415).
Она имеет интерфейс:

```c++
static pair<Matrix, Matrix> decompose(
    const Matrix& A,
    const double epsilon
);
```

#### Параметры

| Параметр  | Тип             | Назначение                                                        |
|-----------|-----------------|-------------------------------------------------------------------|
| `A`       | `const Matrix&` | Исходная матрица                                                  |
| `epsilon` | `double`        | Порог, ниже которого диагональный элемент $R$ считается нулевым |

Возвращаемое значение:

```c++
pair<Matrix, Matrix>
```

где:
- `first`– $Q$
- `second`– $R$

---

### Внутренние функции:

#### `dot`

[`GramMethod.hpp`, строки 347–357](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GramMethod.hpp#L347-L357)

```c++
static double dot(
    const vector<double>& a,
    const vector<double>& b
);
```

Параметры:
- `a`– первый вектор
- `b`– второй вектор

Вычисляет и возвращает скалярное произведение.

---

#### `norm`

[`GramMethod.hpp`, строки 361](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GramMethod.hpp#L361)

```c++
static double norm(const vector<double>& v);
```

Вычисляет евклидову норму.

---

#### `decompose`

[`GramMethod.hpp`, строки 365–415](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GramMethod.hpp#L365-L415)

Алгоритм последовательно:

1. Копирует $k$-й столбец $A$ в`u`
2. Вычитает из него проекции на уже найденные столбцы $Q$
3. Сохраняет соответствующие коэффициенты в $R$
4. Вычисляет $|u|$
5. Проверяет его относительно`epsilon`
6. Нормирует $u$ и записывает результат в $Q$

Ключевой участок:

```c++
R[j][k] = r;
u[i] -= r * Q[i][j];
```

соответствует математическому выражению $u_k=u_k−r_{jk}q_j$.

---

## 3.2. Хаусхолдер

Отражение Хаусхолдера задаётся матрицей $H=I-2uu^T$, где $\|u\|_2=1$.

Для такого преобразования: $H^TH=H^2=I$.

Следовательно, $H$ является ортогональной матрицей.

Главное свойство метода – возможность одним преобразованием занулить сразу все элементы ниже диагонального элемента
текущего столбца.

Для численной устойчивости знак выбирается так, чтобы избежать вычитания близких чисел.

### Реализация

#### Основная функция:

[`Methods/HouseholderMethod.hpp`, строки 515–562](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/HouseholderMethod.hpp#L515-L562)

```c++
static pair<Matrix, Matrix> decompose(
    const Matrix& A,
    const double epsilon
);
```

---

#### Внутренняя функция `norm`

[`HouseholderMethod.hpp`, строки 433–445](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/HouseholderMethod.hpp#L433-L445)

```c++
static double norm(const vector<double>& v);
```

---

#### `applyHouseholder`

[`HouseholderMethod.hpp`, строки 447–511](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/HouseholderMethod.hpp#L447-L511)

```c++
static void applyHouseholder(
    Matrix& R,
    Matrix& Qt,
    const size_t k,
    const double epsilon
);
```

Параметры:

| Параметр  | Тип       | Назначение                         |
|-----------|-----------|------------------------------------|
| `R`       | `Matrix&` | Текущая преобразуемая матрица      |
| `Qt`      | `Matrix&` | Накопленная матрица преобразований |
| `k`       | `size_t`  | Номер текущего шага                |
| `epsilon` | `double`  | Порог малой нормы                  |

На шаге $k$ из текущего столбца создаётся вектор и нормирует его. Само отражение применяется в компактной форме, что
позволяет не создавать полную матрицу $H$ явно.

### Формирование $Q$

После применения отражений вычисляется транспонированная накопленная матрица:

```c++
Q[i][j] = Qt[j][i];
```

Это реализовано в строках 535–542.

Затем из преобразованной матрицы извлекается верхняя треугольная часть $R$ – строки
544–550. ([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/HouseholderMethod.hpp "QR-Decomposition/Methods/HouseholderMethod.hpp at main · 1nteractme/QR-Decomposition · GitHub"))

---

## 3.3. Гивенс

Вращение Гивенса действует только на две координаты.

В выбранной двумерной плоскости оно имеет вид $G=(c_s - s_c)$, где $c_2 + s_2 = 1$.

Для пары $(a, b)$ необходимо получить $r_0$.

В отличие от Хаусхолдера, одно вращение Гивенса обычно уничтожает**один выбранный элемент**.

---

#### `givens`

[`Methods/GivensMethod.hpp`, строки 409–421](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GivensMethod.hpp#L409-L421)

```c++
static pair<double, double> givens(
    const double a,
    const double b
);
```

Параметры:

- `a`– диагональный элемент;
- `b`– элемент, который необходимо занулить.

Возвращает пару

```c++
{c, s}
```

В реализации используется`std::hypot(a, b)`, что предпочтительнее прямого вычисления $\sqrt{a^2+b^2}$ при работе с
числами существенно разных масштабов.

---

#### `applyGivens`

[`GivensMethod.hpp`, строки 423–458](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GivensMethod.hpp#L423-L458)

```c++
static void applyGivens(
    Matrix& R,
    Matrix& Qt,
    size_t i,
    size_t j,
    double c,
    double s
);
```

Параметры:

| Параметр | Назначение                               |
|----------|------------------------------------------|
| `R`      | Текущая матрица                          |
| `Qt`     | Накопленное ортогональное преобразование |
| `i`      | Первая строка вращения                   |
| `j`      | Вторая строка вращения                   |
| `c`      | Косинус вращения                         |
| `s`      | Синус вращения                           |

---

### `decompose`

[`GivensMethod.hpp`, строки 462–530](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GivensMethod.hpp#L462-L530)

```c++
static pair<Matrix, Matrix> decompose(
    const Matrix& A,
    const double epsilon
);
```

Алгоритм:

1. Создаётся копия $A$ в`R`
2. Создаётся единичная матрица`Qt`
3. Обходятся столбцы
4. Внутри каждого столбца обрабатываются элементы ниже диагонали
5. Для ненулевого элемента вычисляется вращение
6. Вращение применяется к $R$ и`Qt`
7. Из`Qt`формируется $Q$
8. Извлекается верхняя треугольная часть $R$
9. Проверяется диагональ $R$

Проверка выполняется:

```c++
if (abs(Rthin[i][i]) <= epsilon)
    throw runtime_error("Метод Гивенса выявил проблему ранга.");
```

Таким образом,`epsilon`одновременно выступает практическим критерием обнаружения вырожденности.

---

# 4. Архитектура проекта

Текущая структура репозитория:

```text
QR-Decomposition/
├── Methods/
│   ├── GivensMethod.hpp
│   ├── GramMethod.hpp
│   └── HouseholderMethod.hpp
├── .gitignore
├── CMakeLists.txt
├── Data.hpp
├── QRDecomposition.hpp
├── Table.hpp
└── main.c++
```

Структура непосредственно соответствует принципам ООП:

```text
                    QRDecomposition
                           |
          +----------------+----------------+
          |                |                |
          v                v                v
      GramMethod   HouseholderMethod   GivensMethod
          |                |                |
          +----------------+----------------+
                           |
                           v
                          Q,R
                           |
                           v
                         Table
```

Корневой файл`QRDecomposition.hpp`выступает фасадом над тремя
реализациями. ([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/QRDecomposition.hpp "QR-Decomposition/QRDecomposition.hpp at main · 1nteractme/QR-Decomposition · GitHub"))

---

# 5. Описание исходного кода

## 5.1. `QRDecomposition.hpp`

Файл:[`QRDecomposition.hpp`](https://github.com/1nteractme/QR-Decomposition/blob/main/QRDecomposition.hpp)

Класс:

```c++
class QRDecomposition
```

предоставляет единый интерфейс для всех трёх методов.

### `Matrix`

```c++
using Matrix = vector<vector<double>>;
```

Матрица представляется как динамический массив строк.

### `Vector`

```c++
using Vector = vector<double>;
```

### `Method`

[
`QRDecomposition.hpp`, строки 397–416](https://github.com/1nteractme/QR-Decomposition/blob/main/QRDecomposition.hpp#L397-L416)

```c++
enum class Method
{
    GramSchmidt,
    Householder,
    Givens
};
```

Позволяет выбрать алгоритм во время вызова`decompose`.

---

### `checkMatrix`

[
`QRDecomposition.hpp`, строки 430–441](https://github.com/1nteractme/QR-Decomposition/blob/main/QRDecomposition.hpp#L430-L441)

```c++
static void checkMatrix(const Matrix& A);
```

Проверяет:

1. Матрица не пустая;
2. Строки имеют одинаковую длину.

При нарушении условий выбрасывается`std::invalid_argument`.

---

### Конструктор

[`QRDecomposition.hpp`, строки 447–458](https://github.com/1nteractme/QR-Decomposition/blob/main/QRDecomposition.hpp#L447-L458)

```c++
explicit QRDecomposition(
    const Matrix& A,
    const double epsilon
);
```

Параметры:

| Параметр  | Назначение       |
|-----------|------------------|
| `A`       | исходная матрица |
| `epsilon` | численный порог  |

Проверяется условие $m >= n$. Если $m < n$, выбрасывается исключение:

```text
Ожидается m >= n.
```

Также проверяется $\epsilon > 0$.

---

### `decompose`

[`QRDecomposition.hpp`, строки 462–484](https://github.com/1nteractme/QR-Decomposition/blob/main/QRDecomposition.hpp#L462-L484)

```c++
[[nodiscard]]
pair<Matrix, Matrix> decompose(
    const Method method
) const;
```

Параметр:

- `method`– выбранный алгоритм.

Возвращает:

```text
{Q, R}
```

Внутри выполняется диспетчеризация:

```c++
case Method::GramSchmidt:
    return GramMethod::decompose(A_, epsilon_);

case Method::Householder:
    return HouseholderMethod::decompose(A_, epsilon_);

case Method::Givens:
    return GivensMethod::decompose(A_, epsilon_);
```

---

### `originalMatrix`

```c++
const Matrix& originalMatrix() const;
```

Возвращает исходную матрицу $A$.

### `rows`

```c++
size_t rows() const;
```

Возвращает количество строк $m$.

### `cols`

```c++
size_t cols() const;
```

Возвращает количество столбцов $n$.

---

## 5.2. `GramMethod.hpp`

Файл:[`Methods/GramMethod.hpp`](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GramMethod.hpp)

Реализует классический Грамм–Шмидт.

Основной вычислительный цикл находится в строках
376–413. ([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GramMethod.hpp "QR-Decomposition/Methods/GramMethod.hpp at main · 1nteractme/QR-Decomposition · GitHub"))

### Сложность

Для плотной матрицы $m\times n$: $O(m_n^2)$.
Память: $O(m^n+n^2)$.

Главная вычислительная нагрузка приходится на скалярные произведения и последовательное вычитание проекций.

---

## 5.3. `HouseholderMethod.hpp`

Файл: [`Methods/HouseholderMethod.hpp`](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/HouseholderMethod.hpp)

Алгоритм выполняет $n$ последовательных преобразований.
На каждом шаге $k$ рассматривается хвост текущего столбца: $x=Rk:m,k$.
В реализации полная матрица Хаусхолдера не формируется.

([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/HouseholderMethod.hpp "QR-Decomposition/Methods/HouseholderMethod.hpp at main · 1nteractme/QR-Decomposition · GitHub"))

### Сложность

Для плотной матрицы: $O(m_n^2)$.

При $m\ge n$ асимптотика совпадает с Граммом–Шмидтом, но константы и структура операций отличаются.

Для плотных матриц Хаусхолдер обычно предпочтительнее классического Грамма–Шмидта по численной устойчивости.

---

## 5.4. `GivensMethod.hpp`

Файл: [`Methods/GivensMethod.hpp`](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GivensMethod.hpp)

Внешний цикл:

```c++
for (size_t col = 0; col < n; ++col)
```

проходит по столбцам, а внутренний:

```c++
for (size_t row = m; row-- > col + 1;)
```

зануляет элементы под диагональю. ([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/Methods/GivensMethod.hpp "QR-Decomposition/Methods/GivensMethod.hpp at main · 1nteractme/QR-Decomposition · GitHub"))

Для плотной матрицы количество вращений имеет порядок $O(m_n)$.

Каждое вращение требует обработки соответствующей части строки/матрицы, поэтому общая вычислительная сложность для плотного случая также имеет порядок $O(m_n^2)$.

Однако Гивенс особенно интересен для разреженных матриц: можно выбирать конкретные элементы, которые требуется занулить, и не затрагивать структуру матрицы больше необходимого.

---

## 5.5. `Data.hpp`

Файл: [`Data.hpp`](https://github.com/1nteractme/QR-Decomposition/blob/main/Data.hpp)

Класс:

```c++
class Data
```

предназначен для хранения тестовых матриц.

Перечисление:

```c++
enum MatrixType
{
    a, b, c, d, e, f, g, h, i, j
};
```

выбирает одну из 10 заранее заданных матриц.

Конструктор:

```c++
explicit Data(const MatrixType type);
```

Сохраняет тип матрицы.

### `getMatrix`

[`Data.hpp`, строки 455–488](https://github.com/1nteractme/QR-Decomposition/blob/main/Data.hpp#L455-L488)

```c++
const vector<vector<double>>& getMatrix() const;
```

Возвращает ссылку на соответствующую тестовую матрицу.

Например,`MatrixA`:

```c++
{
    {1, 2},
    {3, 4}
};
```

Именно она используется в демонстрационном запуске
программы. ([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/Data.hpp "QR-Decomposition/Data.hpp at main · 1nteractme/QR-Decomposition · GitHub"))

---

## 5.6. `Table.hpp`

Файл: [`Table.hpp`](https://github.com/1nteractme/QR-Decomposition/blob/main/Table.hpp)

Класс`Table`отвечает не за математические вычисления, а за форматированный вывод результатов.

Основной класс:

```c++
class Table
```

использует:

```c++
using Matrix = vector<vector<double>>;
```

### `utf8Length`

[`Table.hpp`, строки 702–716](https://github.com/1nteractme/QR-Decomposition/blob/main/Table.hpp#L702-L716)

```c++
static size_t utf8Length(
    const string& str
);
```

Вычисляет количество UTF-8 символов. 
Это необходимо потому, что заголовки таблицы содержат кириллицу.

---

### `padRight`

[`Table.hpp`, строки 717–729](https://github.com/1nteractme/QR-Decomposition/blob/main/Table.hpp#L717-L729)

```c++
static string padRight(
    const string& str,
    size_t width
);
```

Добавляет пробелы справа до заданной ширины.

---

### `formatMatrix`

[`Table.hpp`, строки 731–760](https://github.com/1nteractme/QR-Decomposition/blob/main/Table.hpp#L731-L760)

```c++
static vector<string> formatMatrix(
    const Matrix& matrix,
    int precision = 6
);
```

Параметры:

| Параметр    | Назначение                               |
|-------------|------------------------------------------|
| `matrix`    | матрица для отображения                  |
| `precision` | количество знаков после десятичной точки |

По умолчанию: `precision=6`.

Для чисел используется формат:

```c++
fixed << setprecision(precision)
```

---

### `columnWidth`

[`Table.hpp`, строки 764–813](https://github.com/1nteractme/QR-Decomposition/blob/main/Table.hpp#L764-L813)

```c++
size_t columnWidth(size_t index) const;
```

Определяет необходимую ширину столбца таблицы.

---

### `SectionBuilder`

[`Table.hpp`, строки 845–905](https://github.com/1nteractme/QR-Decomposition/blob/main/Table.hpp#L845-L905)

Позволяет строить секции таблицы цепочкой вызовов.

Например:

```c++
table
    .section("Q")
    .text("Q")
    .matrix(Qgs)
    .matrix(Qgv)
    .matrix(Qhh);
```

---

### `section`

```c++
SectionBuilder section(
    const string& name
);
```

Создаёт новую секцию таблицы.

---

### `print`

[`Table.hpp`, строки 907–995](https://github.com/1nteractme/QR-Decomposition/blob/main/Table.hpp#L907-L995)

```c++
void print() const;
```

Выводит сформированную таблицу в консоль.

---

## 5.7. `main.c++`

Файл: [`main.c++`](https://github.com/1nteractme/QR-Decomposition/blob/main/main.c++)

`main.c++` объединяет все компоненты проекта.

## Настройка точности

```c++
constexpr double EPSILON(1e-15);
```

Используется $\epsilon=1e−15$.

Эта величина передаётся в `QRDecomposition`.
([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/main.c++ "QR-Decomposition/main.c++ at main · 1nteractme/QR-Decomposition · GitHub"))

---

## Выбор матрицы

```c++
const Data DATA(a);
const auto& MATRIX = DATA.getMatrix();
```

Выбирается`MatrixA`.

---

## Создание QR-разложения

```c++
const QRDecomposition QR(MATRIX, EPSILON);
```

Создаётся объект-фасад с исходной матрицей и численным порогом.

---

## Запуск Грамма–Шмидта

```c++
auto [Qgs, Rgs] =
    QR.decompose(QRDecomposition::Method::GramSchmidt);
```

Результаты:

- `Qgs`;
- `Rgs`.

---

## Запуск Гивенса

```c++
auto [Qgv, Rgv] =
    QR.decompose(QRDecomposition::Method::Givens);
```

---

## Запуск Хаусхолдера

```c++
auto [Qhh, Rhh] =
    QR.decompose(QRDecomposition::Method::Householder);
```

Все три вызова расположены в строках 325–329.
([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/main.c++ "QR-Decomposition/main.c++ at main · 1nteractme/QR-Decomposition · GitHub"))

---

## Формирование таблицы

В строках 334–364 создаются четыре столбца:

```text
МАТРИЦА
ГРАММ-ШМИДТ
ГИВЕНС
ХАУСХОЛДЕР
```

Затем отдельно выводятся секции $Q$ и $R$.
([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/main.c++ "QR-Decomposition/main.c++ at main · 1nteractme/QR-Decomposition · GitHub"))

---

## 5.8. `CMakeLists.txt`

Файл: [`CMakeLists.txt`](https://github.com/1nteractme/QR-Decomposition/blob/main/CMakeLists.txt)

Конфигурация содержит:

```cmake
cmake_minimum_required(VERSION 3.5)
project(QRDecomposition)
set(CMAKE_CXX_STANDARD 26)
add_executable(QRDecomposition main.c++)
```

Проект использует стандарт: C++26.

Исполняемый файл называется: `QRDecomposition`.

([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/CMakeLists.txt "QR-Decomposition/CMakeLists.txt at main · 1nteractme/QR-Decomposition · GitHub"))

---

# 6. API

## Создание объекта

```c++
QRDecomposition qr(A, epsilon);
```

где:

- `A`– матрица размера $m\times n$
- `epsilon`– положительный численный порог

Требования: $m>=n,\qquad \epsilon>0$

---

## Выполнение разложения

### Грамм–Шмидт

```c++
auto [Q, R] =
    qr.decompose(QRDecomposition::Method::GramSchmidt);
```

### Хаусхолдер

```c++
auto [Q, R] =
    qr.decompose(QRDecomposition::Method::Householder);
```

### Гивенс

```c++
auto [Q, R] =
    qr.decompose(QRDecomposition::Method::Givens);
```

Во всех случаях ожидается: $A\approx QR$.

---

# 7 Роль `epsilon`

Проект использует явный параметр $\epsilon$.

В демонстрационном запуске: `ε=1e−15`.

Если вычисленный диагональный элемент $R_{kk}$ удовлетворяет $|R_{kk}| <= \epsilon$, метод считает ситуацию проблемой ранга.

Это позволяет обнаруживать:

- Линейно зависимые столбцы
- Нулевые столбцы
- Практически зависимые направления

Однако`epsilon`не является универсальной оценкой обусловленности матрицы. Его значение зависит от:

- Масштаба входных данных
- Машинной точности
- Размера матрицы
- Структуры задачи

Поэтому фиксированное значение $10^{-15}$ следует рассматривать как параметр конкретного эксперимента, а не как
универсальную константу для любых матриц.

---

# 8. Сравнение методов

| Характеристика               | Грамм–Шмидт         | Хаусхолдер                    | Гивенс               |
|------------------------------|---------------------|-------------------------------|----------------------|
| Основная идея                | Проекции            | Отражения                     | Вращения             |
| Что строится непосредственно | $Q$                 | Преобразования                | Преобразования       |
| Что зануляется за один шаг   | Компоненты проекций | Целый хвост столбца           | Обычно один элемент  |
| Ортогональность              | Может деградировать | Высокая устойчивость          | Высокая устойчивость |
| Плотные матрицы              | Хорошо для обучения | Особенно эффективен           | Менее выгоден        |
| Разреженные матрицы          | Не всегда удобен    | Может разрушать разреженность | Особенно удобен      |
| Почти зависимые столбцы      | Проблематичен       | Надёжен                       | Надёжен              |
| Асимптотика                  | $O(mn^2)$           | $O(mn^2)$                     | $O(mn^2)$            |

---

# 9. Тестовые данные

Проект содержит десять матриц, специально подобранных для различных ситуаций.

Особенно важны:

### Линейно зависимые данные `MatrixC`
### Почти линейно зависимые данные `MatrixD`

Это хороший тест на чувствительность вычислений с плавающей точкой.

---

### Матрица с разными масштабами `MatrixG`

Cодержит одновременно очень большие и очень малые элементы. Такой тест позволяет увидеть, насколько алгоритм
чувствителен к масштабированию.

---

### Матрица Гильберта `MatrixI`

Матрицы Гильберта известны плохой обусловленностью, поэтому они являются полезным тестом для сравнения ортогонализации.

---

# 10. Сборка и запуск

## Требования

Необходимы:

- CMake версии не ниже 3.5
- Компилятор с поддержкой используемого стандарта C++
- Поддержка C++26 в соответствии с текущим`CMakeLists.txt`

---

## Клонирование

```bash
git clone https://github.com/1nteractme/QR-Decomposition.git
cd QR-Decomposition
```

---

## Сборка

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

После успешной сборки будет создан исполняемый файл:

```text
QRDecomposition
```

Запуск:

```bash
./QRDecomposition
```

---

# 11. Результат работы

```c++
const Data DATA(a);
```

[GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/Data.hpp "QR-Decomposition/Data.hpp at main · 1nteractme/QR-Decomposition · GitHub")

Программа выводит результаты трёх методов в единой таблице.

```text
+---------+-------------------------------+-------------------------------+-------------------------------+
| МАТРИЦА | ГРАММ-ШМИДТ                   | ГИВЕНС                        | ХАУСХОЛДЕР                    |
+---------+-------------------------------+-------------------------------+-------------------------------+
| Q       | [     0.316228     0.948683 ] | [     0.316228    -0.948683 ] | [    -0.316228     0.948683 ] |
|         | [     0.948683    -0.316228 ] | [     0.948683     0.316228 ] | [    -0.948683    -0.316228 ] |
+---------+-------------------------------+-------------------------------+-------------------------------+
| R       | [     3.162278     4.427189 ] | [     3.162278     4.427189 ] | [    -3.162278    -4.427189 ] |
|         | [     0.000000     0.632456 ] | [     0.000000    -0.632456 ] | [     0.000000     0.632456 ] |
+---------+-------------------------------+-------------------------------+-------------------------------+

Process finished with exit code 0
```

Формирование именно такой структуры вывода реализовано в`main.cpp`: создаются четыре колонки, а затем секции`Q`и`R`, в
которые последовательно передаются результаты Грамма–Шмидта, Гивенса и
Хаусхолдера. ([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/main.c++ "QR-Decomposition/main.c++ at main · 1nteractme/QR-Decomposition · GitHub"))

---

# 12. Ограничения реализации

Текущая версия проекта является учебной реализацией алгоритмов и имеет ряд осознанных ограничений.

## Только $m\ge n$

Класс проверяет:

```c++
if (m < n)
    throw invalid_argument("Ожидается m >= n.");
```

То есть поддерживаются прямоугольные матрицы с количеством строк не меньше количества
столбцов. 
([GitHub](https://github.com/1nteractme/QR-Decomposition/blob/main/QRDecomposition.hpp "QR-Decomposition/QRDecomposition.hpp at main · 1nteractme/QR-Decomposition · GitHub"))

---

## Использование`vector<vector<double>>`

Матрица представлена как:

```c++
vector<vector<double>>
```

`vector<vector<double>>`имеет:

- Дополнительную косвенность
- Отдельные аллокации строк
- Менее предсказуемую работу с кэшем
- Дополнительные накладные расходы

---

## Фиксированный `epsilon`

В демонстрационной программе используется:

```c++
constexpr double EPSILON(1e-15);
```

Универсального оптимального значения для всех матриц не существует.

Для практического численного пакета`epsilon`целесообразно выбирать с учётом:

- Нормы матрицы
- Размера матрицы
- Характера входных данных

---

# Заключение

Проект демонстрирует три фундаментальных способа построения QR-разложения: $A=QR$ с ортонормированной матрицей $Q$ и
верхнетреугольной матрицей $R$.

Для плотных матриц естественным выбором является Хаусхолдер, поскольку одно отражение уничтожает целый хвост столбца.

Гивенс особенно полезен, когда необходимо контролируемо занулять отдельные элементы, в том числе при работе с
разреженными структурами.

Грамм–Шмидт остаётся наиболее прозрачным с геометрической точки зрения и хорошо демонстрирует сам принцип построения
ортонормированного базиса.