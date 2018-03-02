// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

// Reference: vtkTuple from VTK

#ifndef TUPLE_H_
#define TUPLE_H_

#include <cassert>
#include <iostream>
#include "thrust_common.h"

namespace edda {

// Compared to std::vector, Tuple is a fixed-length vector for space saving
template<typename T, int N>
class Tuple
{
protected:
    T vec[N];

public:
    enum { LENGTH = N };
    typedef int IsTuple;  // for compile-time generic type check

    __host__ __device__
    Tuple()
    {}

    // assign all values as scalar
    __host__ __device__
    explicit Tuple(const T& scalar)
    {
        for (int i = 0; i < N; ++i)
            this->vec[i] = scalar;
    }

    __host__ __device__
    explicit Tuple(const T* init)
    {
        for (int i = 0; i < N; ++i)
            this->vec[i] = init[i];
    }

    // assign Tuple
    template<typename U>
    __host__ __device__
    Tuple(const Tuple<U, N>& v)
    {
        for (int i = 0; i < N; ++i)
            this->vec[i] = v[i];
    }

    // Get the length of the tuple.
    __host__ __device__
    inline int length() const { return N; }

    // Get a pointer to the underlying data of the tuple.
    __host__ __device__
    inline T* getData() { return this->vec; }
    __host__ __device__
    inline const T* getData() const { return this->vec; }

    // Get a reference to the underlying data element of the tuple.
    // This works similarly to the way C++ STL containers work.  No
    // bounds checking is performed.
    __host__ __device__
    inline T& operator[](int i) {
        //assert("pre: index_in_bounds" && i >= 0 && i < N);
        return this->vec[i];
    }
    __host__ __device__
    inline const T& operator[](int i) const {
        //assert("pre: index_in_bounds" && i >= 0 && i < N);
        return this->vec[i];
    }

    // Cast the tuple to the specified type, returning the result.
    template<typename TR>
    __host__ __device__
    Tuple<TR, N> cast() const
    {
        Tuple<TR, N> result;
        for (int i = 0; i < N; ++i)
        {
            result[i] = static_cast<TR>(this->vec[i]);
        }
        return result;
    }

};

// Convenience use
template<typename T> struct Tuple1 : public Tuple<T, 1> {};
template<typename T> struct Tuple2 : public Tuple<T, 2> {
    Tuple2() {}
    explicit Tuple2(const T& x, const T& y)
        {this->vec[0] = x; this->vec[1] = y;}
};
template<typename T> struct Tuple3 : public Tuple<T, 3> {
    Tuple3() {}
    explicit Tuple3(const T& x, const T& y, const T& z)
        {this->vec[0] = x; this->vec[1] = y; this->vec[2] = z;}
};
template<typename T> struct Tuple4 : public Tuple<T, 4> {
    Tuple4() {}
    explicit Tuple4(const T& x, const T& y, const T& z, const T& w)
        {this->vec[0] = x; this->vec[1] = y; this->vec[2] = z; this->vec[3] = w;}
};
template<typename T> struct Tuple5 : public Tuple<T, 5> {};
template<typename T> struct Tuple6 : public Tuple<T, 6> {};
template<typename T> struct Tuple7 : public Tuple<T, 7> {};
template<typename T> struct Tuple8 : public Tuple<T, 8> {};
template<typename T> struct Tuple9 : public Tuple<T, 9> {};


// Output the contents of a tuple, mainly useful for debugging.
template<typename A, int N>
std::ostream& operator<<(std::ostream& out, const Tuple<A, N>& t)
{
    out << "(";
    bool first = true;
    for (int i = 0; i < N; ++i)
    {
        if (first)
        {
            first = false;
        }
        else
        {
            out << ", ";
        }
        out << t[i];
    }
    out << ")";
    return out;
}

// Equality operator performs an equality check on each component.
template<typename A, int N>
bool operator==(const Tuple<A, N>& t1, const Tuple<A, N>& t2)
{
    for (int i = 0; i < N; ++i)
    {
        if (t1[i] != t2[i])
        {
            return false;
        }
    }
    return true;
}


// Inequality for vector type.
template<typename A, int N>
inline bool operator!=(const Tuple<A, N>& t1, const Tuple<A, N>& t2)
{
    return !(t1 == t2);
}

}  // edda

#endif  // TUPLE_H_
