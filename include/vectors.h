#pragma once
#include<iostream>
#include<stdexcept>
#include<random>
#include <utility> //for swap
#include <complex>
#include <cmath>
using namespace std;

template<typename T>
class Vector {
	T* _data;
	size_t _size;
public:
	Vector(size_t size);
	Vector(size_t size, T value);
	Vector(size_t size, T lower_bound, T upper_bound);
	~Vector();
	Vector(const Vector<T>& other);
	void swap(Vector<T>& other) noexcept;
	Vector<T>& operator=(const Vector<T>& other);
	T& operator[](size_t index);

	Vector<T>& operator+(const Vector<T>& other);
	Vector<T>& operator-(const Vector<T>& other);

	T operator*(const Vector<T>& other);
	//complex<float> operator*(const Vector<complex<float>>& other);
	//complex<double> operator*(const Vector<complex<double>>& other);

	Vector<T> operator*(T scalar) const;
	Vector<T> operator/(T scalar) const;

	bool operator==(const Vector<T>& other) const;
	bool operator!=(const Vector<T>& other) const;

	size_t get_size();
	T* get_data();
};

template<typename T>
Vector<T> find_perpendiculator_unit_vector(const Vector<T>& vector);
