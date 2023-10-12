#include <vectors.h>
using namespace std;

template<typename T>
Vector<T>::Vector(size_t size) : _size(size) {
	_data = new T[size];
}

//конструктор с параметрами : длина вектора и значение для заполнения;
template<typename T>
Vector<T>::Vector(size_t size, T value) : _size(size) {
	if (size <= 0)
		throw invalid_argument("Vector size must be positive.");
	_data = new T[size];
	for (size_t i = 0; i < size; i++) {
		data[i] = value;
	}
}

//конструктор с параметрами, заполняющий вектор случайными значениями : 
	//длина вектора; нижняя граница, верхняя граница;
template<typename T>
Vector<T>::Vector(size_t size, T lower_bound, T upper_bound) : _size(size) {
	if (size <= 0)
		throw invalid_argument("Vector size must be positive.");
	if (lower_bound >= upper_bound)
		throw invalid_argument("Lower bound must be less than upper bound.");

	_data = new T[size];
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<T> distribution(lower_bound, upper_bound);
	for (size_t i = 0; i < size; i++) {
		_data[i] = distribution(gen);
	}
}

template<typename T>
Vector<T>::~Vector() {
	delete[] _data;
}

//copy constructor.
template<typename T>
Vector<T>:: Vector(const Vector& other) : _size(other._size) {
	delete[] _data;
	_data = new T[_size];
	for (size_t i = 0; i < _size; i++) {
		_data[i] = other._data[i];
	}
}

template<typename T>
void Vector<T>::swap(Vector<T>& other) noexcept {
	std::swap(_data, other._data);
	std::swap(_size, other._size);
}

template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& other) {
	//temporary
	Vector<T> tmp(other);
	swap(tmp);
	return *this;
}

//оператор[] для чтения / записи элемента вектора по указанному индексу;
template<typename T>
T& Vector<T>::operator[](size_t index) {
	if (index < 0 || index >= _size)
		throw invalid_argument("Index out of range.");
	return _data[index];
}

//операторы сложения и вычитания векторов;
template<typename T>
Vector<T>& Vector<T>::operator+(const Vector<T>& other) {
	if (_size != other._size)
		throw invalid_argument("Vectors must have the same size for addition.");
	Vector<T> result(_size);
	for (size_t i = 0; i < _size; i++) {
		result[i] = _data[i] + other._data[i];
	}
	return result;
}

template<typename T>
Vector<T>& Vector<T>::operator-(const Vector<T>& other) {
	if (_size != other._size)
		throw invalid_argument("Vectors must have the same size for suctraction.");
	Vector<T> result(_size);
	for (size_t i = 0; i < _size; i++) {
		result[i] = _data[i] - other._data[i];
	}
	return result;
}

//оператор умножения, выполняющий скалярное произведение векторов;
template<typename T>
T Vector<T>::operator*(const Vector<T>& other) {
	if (_size != other._size)
		throw invalid_argument("Vectors must have the same size for dot product.");
	T result = 0;
	for (size_t i = 0; i < _size; i++) {
		result += _data[i] * other._data[i];
	}
	return result;
}
// Specialization for std::complex<float>
template<>
complex<float> Vector<complex<float>>::operator*(const Vector<complex<float>>& other) {
	if (_size != other._size)
		throw invalid_argument("Vectors must have the same size for dot product.");
	complex<float> result = 0;
	for (size_t i = 0; i < _size; i++) {
		result += _data[i] * conj(other._data[i]);
	}
	return result;
}
// Specialization for std::complex<double> conjugate 
template<>
complex<double> Vector<complex<double>>::operator*(const Vector<complex<double>>& other) {
	if (_size != other._size)
		throw invalid_argument("Vectors must have the same size for dot product.");
	complex<double> result = 0;
	for (size_t i = 0; i < _size; i++) {
		result += _data[i] * std::conj(other._data[i]);
	}
	return result;
}

//оператор умножения вектора на скаляр(обеспечить коммутативность);
template<typename T>
Vector<T> Vector<T>::operator* (T scalar) const {
	Vector<T> result(_size);
	for (size_t i = 0; i < 0; i++) {
		result[i] = _data[i] * scalar;
	}
	return result;
}

//оператор деления вектора на скаляр.
template<typename T>
Vector<T> Vector<T>::operator/ (T scalar) const {
	Vector<T> result(_size);
	for (size_t i = 0; i < 0; i++) {
		result[i] = _data[i] / scalar;
	}
	return result;
}

// Comparison operators for vectors.
template<typename T>
bool Vector<T>::operator==(const Vector<T>& other) const {
	if (_size != other._size)
		return false;

	for (size_t i = 0; i < _size; i++) {
		if (_data[i] != other._data[i])
			return false;
	}

	return true;
}

template<typename T>
bool Vector<T>::operator!=(const Vector<T>& other) const {
	return !(*this == other);
}

template<typename T>
size_t Vector<T>::get_size() {
	return _size;
}

template<typename T>
T* Vector<T>::get_data() {
	return _data;
}

template<typename T>
Vector<T> find_perpendiculator_unit_vector(const Vector<T>& vector) {
	Vector<T> unitVector(vector.get_size(), 0);
	T length = 0;
	for (size_t i = 0; i < 2; i++) {
		length += vector[i] * vector[i];
	}

	length = sqrt(length);

	if (length == 0)
		throw invalid_argument("Cannot find perpendicular unit vector for a zero vector.");

	unitVector[0] = vector[1];
	unitVector[1] = -vector[0];
	unitVector = unitVector / length;
	return unitVector;
}
