/* Clase wrapper para std::vector
 * Incluye operaciones aritméticas vectoriales, extracción y otros métodos de utilidad numérica
 */
#ifndef VECTOR_H
#define VECTOR_H

#include <stdexcept>
#include <vector>
#include <iostream>
#include <initializer_list>
#include <cmath>


// Clase vector ------------------------------------
template <class T>
class Vector
{
    private:
        std::vector<T> data;
        unsigned int M;     // Dimensión

    public:
        // Constructores:
        // Vacío
        Vector() {}
        // Con tamaño y valor constante
        Vector(const unsigned int M, T val = T(0));
        // Por copia
        Vector(const Vector<T>& other);
        // Con un std::vector
        Vector(const std::vector<T>& vec);
        // Con initializer list (estilo C++11 con paréntesis de llave)
        Vector(std::initializer_list<T> list);
        // Concatenando vectores
        Vector(std::initializer_list<Vector<T>> vlist);

        // Asignacion
        Vector<T> operator=(const Vector<T>& other);
        Vector<T> operator+=(const Vector<T>& other);
        Vector<T> operator-=(const Vector<T>& other);
        Vector<T> operator*=(const T& other);
        Vector<T> operator/=(const T& other);

        // Acceso a componente i
        // r
        T operator[](const unsigned int i) const;
        // rw
        T& operator[](const unsigned int i);
        // Acceso a slice i-j
        Vector<double> operator()(const unsigned int i, const unsigned int j) const;

        // Comparación
        bool operator==(const Vector<T>& other);

        // Metodos relevantes
        unsigned int length() const { return M; }   // Getter dimensión
        T norm() const;     // Norma-2
};

// Operaciones aritméticas vectoriales (suma vectorial, multiplicación escalar)
template <class T>
Vector<T> operator+(const Vector<T>& vector1, const Vector<T>& vector2);
template <class T>
Vector<T> operator-(const Vector<T>& vector1, const Vector<T>& vector2);
template <class T>
Vector<T> operator*(const T lambda, const Vector<T>& vector);
template <class T>
Vector<T> operator*(const Vector<T>& vector, const T lambda);
template <class T>
T operator*(const Vector<T> v1, const Vector<T>& v2);
template <class T>
Vector<T> operator/(const Vector<T>& vector, const T lambda);
template <class T>
// Extracción
std::ostream& operator<<(std::ostream& os, const Vector<T>& vector);
// Producto punto
template <class T>
T dot(const Vector<T>& vector1, const Vector<T>& vector2);
// Producto cruz
template <class T>
Vector<T> cross(const Vector<T>& vector1, const Vector<T>& vector2);




// Definiciones -----------------------------------------------------------------
// Constructores
template <class T>
Vector<T>::Vector(const unsigned int M, T val) : M(M)
{
    data = std::vector<T>(M, val);
}


template <class T>
Vector<T>::Vector(const Vector<T>& other) : M(other.M)
{
    data = std::vector<T>(M);
    for (unsigned int i = 0; i < M; i++)
        data[i] = other[i];
}


template <class T>
Vector<T>::Vector(const std::vector<T>& vec)
{
    M = vec.size();
    data = vec;
}


template <class T>
Vector<T>::Vector(std::initializer_list<T> list) : M(list.size())
{
    // Iteración por rango estilo C++11
    for (T element : list)
        data.push_back(element);
}

template <class T>
Vector<T>::Vector(std::initializer_list<Vector<T>> vlist)
{
    unsigned int size = 0;
    for (auto &vec : vlist)
        size += vec.length();
    M = size;
    data = std::vector<T>(M);

    auto datait = data.begin();
    for (auto &vec : vlist)
        for (unsigned int i = 0; i < vec.length(); i++)
            *(datait++) = vec[i];
}


// Copia "profunda"
template <class T>
Vector<T> Vector<T>::operator=(const Vector<T>& other)
{
    if (this != &other)
    {
        M = other.M;
        data = other.data;
    }
    return (*this);
}


// Asignaciones llaman a operaciones aritméticas correspondientes
template <class T>
Vector<T> Vector<T>::operator+=(const Vector<T>& other)
{
    (*this) = (*this) + other;
    return *this;
}


template <class T>
Vector<T> Vector<T>::operator-=(const Vector<T>& other)
{
    (*this) = (*this) - other;
    return *this;
}

template <class T>
Vector<T> Vector<T>::operator*=(const T& other)
{
    (*this) = (*this) * other;
    return *this;
}

template <class T>
Vector<T> Vector<T>::operator/=(const T& other)
{
    (*this) = (*this) / other;
    return *this;
}


// Acceso
template <class T>
T Vector<T>::operator[](const unsigned int i) const
{
    return data[i];
}

template <class T>
T& Vector<T>::operator[](const unsigned int i)
{
    return data[i];
}

template <class T>
Vector<double> Vector<T>::operator()(const unsigned int i, const unsigned int j) const
{
    unsigned int N = j - i;
    Vector<T> ret(N);
    for (unsigned int k = 0; k < j - i; k++)
        ret[k] = (*this)[i + k];
    return ret;
}


// Comparación
template <class T>
bool Vector<T>::operator==(const Vector<T>& other)
{
    // Deben ser de igual largo
    if (M != other.M)
        return false;

    // Y tener todos sus elementos iguales
    for (unsigned int i = 0; i < M; i++)
        if ((*this)(i) != other(i))
            return false;
    return true;
}


// Operadores aritméticos
template <class T>
Vector<T> operator+(const Vector<T>& vector1, const Vector<T>& vector2)
{
    if (vector1.length() != vector2.length())
        throw std::domain_error("Se intentaron sumar vectores de distinto largo");

    // Por componentes
    Vector<T> ret(vector1.length());
    for (unsigned int i = 0; i < vector1.length(); i++)
        ret[i] = vector1[i] + vector2[i];

    return ret;
}

template <class T>
Vector<T> operator-(const Vector<T>& vector1, const Vector<T>& vector2)
{
    if (vector1.length() != vector2.length())
        throw std::domain_error("Se intentaron sumar vectores de distinto largo");

    // Por componentes
    Vector<T> ret(vector1.length());
    for (unsigned int i = 0; i < vector1.length(); i++)
        ret[i] = vector1[i] - vector2[i];

    return ret;
}


template <class T>
Vector<T> operator*(const T lambda, const Vector<T>& vector)
{
    // Por componentes
    Vector<T> ret(vector.length());
    for (unsigned int i = 0; i < vector.length(); i++)
        ret[i] = lambda * vector[i];

    return ret;
}

template <class T>
Vector<T> operator*(const Vector<T>& vector, const T lambda)
{
    // Simétrico
    return lambda * vector;
}

template <class T>
T operator*(const Vector<T> v1, const Vector<T>& v2)
{
    return dot(v1,v2);
}

template <class T>
Vector<T> operator/(const Vector<T>& vector, const T lambda)
{
    // Por componente
    Vector<T> ret(vector.length());
    for (unsigned int i = 0; i < vector.length(); i++)
        ret[i] = vector[i] / lambda;

    return ret;
}


// Extrae cada elemento
template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& vector)
{
    for (unsigned int i = 0; i < vector.length(); i++)
        os << vector[i] << '\t';
    return os;
}


// Norma-2
template <class T>
T Vector<T>::norm() const
{
    T ret = T(0);
    for (const T& element : data)
        ret += element * element;
    return std::sqrt(ret);
}


// Producto punto
template <class T>
T dot(const Vector<T>& vector1, const Vector<T>& vector2)
{
    if (vector1.length() != vector2.length())
        throw std::domain_error("El producto punto está definido para vectores de igual dimension");

    // Por componentes
    T ret = T(0);
    for (unsigned int i = 0; i < vector1.length(); i++)
        ret += vector1[i] * vector2[i];
    return ret;
}


// Cruz
template <class T>
Vector<T> cross(const Vector<T>& vector1, const Vector<T>& vector2)
{
    // Chequear tamaño?
    Vector<T> ret(3);
    ret[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
    ret[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
    ret[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];
    return ret;
}


#endif
