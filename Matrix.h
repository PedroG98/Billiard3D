/* Matriz que soporta operaciones matriciales algebráicas
 * Implementa métodos relevantes como inversión, transposición, etc.
 * Construida como un vector de objetos Vector (Ver notas en Vector.h)
 * para facilitar el cálculo de operaciones elementales
 */
#ifndef MATRIX_H
#define MATRIX_H

#include <stdexcept>
#include <vector>
#include <iostream>
#include <cmath>

#include "Vector.h"



// Clase matriz -------------------------------------------------------
template <class T>
class Matrix
{
    private:
        // Array de N M-filas
        std::vector<Vector<T>> rows;
        unsigned int N;
        unsigned int M;

    public:
        // Constructores:
        // Vacío
        Matrix() {}
        // Con dimensiones y valor de relleno
        Matrix(const unsigned int N, const unsigned int M, T fill = T(0));
        // Por copia
        Matrix(const Matrix<T>& other);
        // Con vector de vectores
        Matrix(const std::vector<std::vector<T>>& entries);
        // Con vector y tamaño explícito
        Matrix(const unsigned int N, const unsigned int M, std::vector<T> entries);
        // Con lista y tamaño explícito
        Matrix(const unsigned int N, const unsigned int M, std::initializer_list<T> entries);

        // Asignación
        Matrix<T> operator=(const Matrix<T>& other);

        // Acceso a elemento (i,j) con operador ()
        // r
        T operator()(const unsigned int i, const unsigned int j) const;
        // rw
        T& operator()(const unsigned int i, const unsigned int j);
        // Acceso a fila i con operador []
        // r
        Vector<T> operator[](const unsigned int i) const;
        // rw
        Vector<T>& operator[](const unsigned int i);

        // Getters de dimensiones
        unsigned int n() const { return N; }
        unsigned int m() const { return M; }

        // Descomposiciones
        // Entrega par [Forma escalonada, Matriz de transformación]
        std::vector<Matrix<T>> rowEchelonForm() const;

        // Métodos
        // Determinante
        T det() const;
        // Inversa
        Matrix<T> inverse() const;
        // Transpuesta
        Matrix<T> transposed() const;
        // Aplanar
        Vector<T> flatten() const;
        // Aumentar tamaño
        void resize(unsigned int N, unsigned int M);
};

// Operaciones externas
// Extracción
template <class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& A);

// Aritmética de matrices
template <class T>
Matrix<T> operator+(const Matrix<T> A, const Matrix<T> B);
template <class T>
Matrix<T> operator-(const Matrix<T> A, const Matrix<T> B);
template <class T>
Matrix<T> operator*(const Matrix<T> A, const Matrix<T> B);

// Multiplicación con escalar/vector
template <class T>
Matrix<T> operator*(const T lambda, const Matrix<T> A);
template <class T>
Matrix<T> operator*(const Matrix<T> A, const T lambda);
template <class T>
Matrix<T> operator/(const Matrix<T> A, const T lambda);
template <class T>
Vector<T> operator*(const Matrix<T> A, const Vector<T> v);

// Comparación
template <class T>
bool operator==(const Matrix<T> A, const Matrix<T> B);



// Matrices especiales -------------------------------------------------
// Matriz identidad NxN
template <class T>
Matrix<T> Identity(unsigned int N)
{
    Matrix<T> id(N, N);
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int j = 0; j < N; j++)
            id(i, j) = (i == j) ? T(1) : T(0);
    return id;
}



// Definiciones de template --------------------------------------------
template <class T>
Matrix<T>::Matrix(const unsigned int N, const unsigned int M, T fill) : N(N), M(M)
{
    rows = std::vector<Vector<T>>(N);
    for (unsigned int i = 0; i < N; i++)
        rows[i] = Vector<T>(M);
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int j = 0; j < M; j++)
            (*this)(i, j) = fill;
}


template <class T>
Matrix<T>::Matrix(const Matrix<T>& other) : N(other.N), M(other.M)
{
    rows = std::vector<Vector<T>>(N);
    for (unsigned int i = 0; i < N; i++)
        rows[i] = other[i];

}


template <class T>
Matrix<T>::Matrix(const std::vector<std::vector<T>>& entries)
{
    // Si el input es vacío, error
    if (entries.empty())
        throw std::invalid_argument("Se intentó construir una matriz sin datos");

    // Definir tamaño
    N = entries.size();
    M = entries[0].size();
    rows = std::vector<Vector<T>>(N);
    for (unsigned int i = 0; i < N; i++)
        rows[i] = Vector<T>(M);

    // Revisar que dimensión es constante en todas las filas y escribir datos
    for (unsigned int i = 0; i < N; i++)
    {
        if (entries[i].size() != M)
            throw std::invalid_argument("Filas de datos no tienen igual dimension");

        for (unsigned int j = 0; j < M; j++)
            (*this)(i, j) = entries[i][j];
    }
}

template <class T>
Matrix<T>::Matrix(const unsigned int N, const unsigned int M, std::initializer_list<T> entries) : M(M), N(N)
{
    if (entries.size() == 0)
        throw std::invalid_argument("Se intentó construir una matriz sin datos");

    rows = std::vector<Vector<T>>(N);
    for (unsigned int i = 0; i < N; i++)
        rows[i] = Vector<T>(M);

    unsigned int i=0, j=0;
    for (const double entry : entries)
    {
        rows[j][i] = entry;
        i++;
        if (i >= M)
        {
            i = 0;
            j++;
        }
    }
}


template <class T>
Matrix<T>::Matrix(const unsigned int N, const unsigned int M, std::vector<T> entries) : N(N), M(M)
{
    // Si el input no es consistente con las dimensiones, error
    if (entries.size() != M*N)
        throw std::invalid_argument("El tamaño de los datos no corresponde a las dimensiones de la matriz");

    // Alocar memoria correspondientemente
    rows = std::vector<Vector<T>>(N);
    for (unsigned int i = 0; i < N; i++)
    {
        rows[i] = Vector<T>(M);
        for (unsigned int j = 0; j < M; j++)
            (*this)(i, j) = entries[i * M + j];
    }
}

template <class T>
Matrix<T> Matrix<T>::operator=(const Matrix<T>& other)
{
    if (this != &other)
    {
        rows = other.rows;
        N = other.N;
        M = other.M;
    }
    return *this;
}


template <class T>
T Matrix<T>::operator()(const unsigned int i, const unsigned int j) const
{
    return rows[i][j];
}

template <class T>
T& Matrix<T>::operator()(const unsigned int i, const unsigned int j)
{
    return rows[i][j];
}

template <class T>
Vector<T> Matrix<T>::operator[](const unsigned int i) const
{
    return rows[i];
}

template <class T>
Vector<T>& Matrix<T>::operator[](const unsigned int i)
{
    return rows[i];
}

// Extracción
template <class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& A)
{
    for (unsigned int i = 0; i < A.n(); i++)
        os << A[i] << std::endl;
    return os;
}

// Operadores aritméticos ------------------------------------------------
template <class T>
Matrix<T> operator+(const Matrix<T> A, const Matrix<T> B)
{
    Matrix<T> C(A.n(), A.m());
    for (unsigned int i = 0; i < A.n(); i++)
        C[i] = A[i] + B[i];
    return C;
}


template <class T>
Matrix<T> operator-(const Matrix<T> A, const Matrix<T> B)
{
    Matrix<T> C(A.n(), A.m());
    for (unsigned int i = 0; i < A.n(); i++)
        C[i] = A[i] - B[i];
    return C;
}


template <class T>
Matrix<T> operator*(const T lambda, const Matrix<T> A)
{
    Matrix<T> C(A.n(), A.m());
    for (unsigned int i = 0; i < A.n(); i++)
        C[i] = lambda * A[i];
    return C;
}


template <class T>
Matrix<T> operator*(const Matrix<T> A, const T lambda)
{
    Matrix<T> C(A.n(), A.m());
    for (unsigned int i = 0; i < A.n(); i++)
        C[i] = lambda * A[i];
    return C;
}


template <class T>
Matrix<T> operator/(const Matrix<T> A, const T lambda)
{
    Matrix<T> C(A.n(), A.m());
    for (unsigned int i = 0; i < A.n(); i++)
        C[i] = A[i] / lambda;
    return C;
}


template <class T>
Vector<T> operator*(const Matrix<T> A, const Vector<T> v) // Valores por copia pues copia es O(n2), lectura para multiplicar es O(n3)
{
    if (A.m() != v.length())
        throw std::domain_error("Matrices de tamaño invalido multiplicadas");

    Vector<T> C(A.n());
    T c;
    for (unsigned int i = 0; i < C.length(); i++)
    {
        c = T(0);
        for (unsigned int j = 0; j < v.length(); j++)
            c += A(i, j) * v[j];
        C[i] = c;
    }
    return C;
}



template <class T>
Matrix<T> operator*(const Matrix<T> A, const Matrix<T> B) // Valores por copia pues copia es O(n2), lectura para multiplicar es O(n3)
{
    if (A.m() != B.n())
        throw std::domain_error("Matrices de tamaño invalido multiplicadas");

    Matrix<T> C(A.n(), B.m());
    T c;
    for (unsigned int i = 0; i < C.n(); i++)
        for (unsigned int j = 0; j < C.m(); j++)
        {
            c = T(0);
            for (unsigned int k = 0; k < A.m(); k++)
                c += A(i, k) * B(k, j);
            C(i, j) = c;
        }
    return C;
}


// Comparación
template <class T>
bool operator==(const Matrix<T> A, const Matrix<T> B)
{
    if (A.n() != B.n() || A.m() != B.m())
        return false;

    for (unsigned int i = 0; i < A.n(); i++)
        if (A[i] != B[i])
            return false;
    return true;
}


template <class T>
std::vector<Matrix<T>> Matrix<T>::rowEchelonForm() const
{
    // La idea es aplicar operaciones inversas a I y M hasta que GM = U, G-1I = L, y M = LU
    auto L = Identity<T>(N);
    auto U = (*this);

    int pivot = -1;
    // Recorre columnas ordenando por pivotes
    for (unsigned int j = 0; j < M; j++)
    {
        int i;
        // Se recorre fila desde el pivote anterior+1 hasta encontrar elemento no nulo -> pivote
        for (i = pivot + 1; i < N; i++)
        {
            if (U(i,j) != 0)
                break;
        }
        // Si se llegó a N, la columna entera tras el pivote es 0 y se pasa a la siguiente
        if (i == N)
            continue;

        // Si hay un nuevo pivote, se intercambia la fila con la siguiente al pivote
        pivot++;
        Vector<T> tmp = U[i];
        U[i] = U[pivot];
        U[pivot] = tmp;
        tmp = L[pivot];
        L[i] = L[pivot];
        L[pivot] = tmp;
        // Se sustraen los elementos bajo el nuevo pivote
        for (unsigned int i = pivot + 1; i < N ; i++)
        {
            T factor = U(i, j) / U(pivot, j);
            U[i] -= U[pivot] * factor;
            L[i] -= L[pivot] * factor;
        }
    }
    std::vector<Matrix<T>> ret;
    ret.push_back(U);
    ret.push_back(L);
    return ret;
}


template <class T>
T Matrix<T>::det() const
{
    // Sólo en matrices cuadradas
    if (N != M)
        throw std::domain_error("Sólo matrices cuadradas tienen determinante");

    // Se escalona y se multiplican los elementos de la diagonal
    Matrix<T> U = rowEchelonForm()[0];
    T det = T(1);
    for (unsigned int i = 0; i < N; i++)
        det *= U(i, i);

    return det;
}


template <class T>
Matrix<T> Matrix<T>::inverse() const
{
    if (det() == 0)
        throw std::domain_error("La matriz no es invertible!");

    auto LU = rowEchelonForm();
    Matrix<T> U = LU[0];
    Matrix<T> L = LU[1];

    // Reducir forma escalonada de ambas hasta que U sea identidad -> L será inversa
    for (unsigned int i = 0; i < N; i++)
    {
        // Normalizar
        T factor = U(i, i);
        U[i] /= factor;
        L[i] /= factor;

        for (unsigned int k = 0; k < i; k++)
        {
            // Eliminar filas bajo diagonal
            factor = U(k, i);
            U[k] -= U[i] * factor;
            L[k] -= L[i] * factor;
        }
    }

    return L;
}


template <class T>
Matrix<T> Matrix<T>::transposed() const
{
    Matrix<T> ret(M, N);
    for (unsigned int i = 0; i < N; i++)
        for (unsigned int j = 0; j < M; j++)
            ret(j, i) = (*this)(i, j);
    return ret;
}


template <class T>
void Matrix<T>::resize(unsigned int N, unsigned int M)
{
    // Sólo si N y M son mayores a los parámetros actuales
    if (N < this->N || M < this->M)
        throw std::invalid_argument("No se puede reducir el tamaño de la matriz!");

    // Cambiar tamaño y alocar variables en nueva matriz
    Matrix<T> tmp = Matrix<T>(N, M);
    for (unsigned int i = 0; i < this->N; i++)
        for (unsigned int j = 0; j < this->M; j++)
            tmp(i, j) = (*this)(i, j);

    // Copiar datos a este objeto
    *this = tmp;
}


template <class T>
Vector<T> Matrix<T>::flatten() const
{
    Vector<double> ret(M*N);
    for (unsigned int i = 0; i < this->N; i++)
        for (unsigned int j = 0; j < this->M; j++)
            ret[i*N + j] = (*this)(i, j);
    return ret;
}


#endif // MATRIX_H
