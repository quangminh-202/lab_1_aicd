#include <random>
#include <complex>
#include <stdexcept>


namespace my_vector {
	template<class T>
	class vector {
		size_t _cnt_coords;
		T* _coords;
		static const inline double eps = 0.00000001;
	public:
		vector() : _cnt_coords(0), _coords(nullptr) {};

		vector(size_t cnt_coords, T* coords) : _cnt_coords(cnt_coords), _coords(new T[cnt_coords]) {
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] = coords[i];
			}
		}

		vector(const vector& copy) :_cnt_coords(copy._cnt_coords), _coords(new T[_cnt_coords]) {
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] = copy[i];
			}
		}

		vector(size_t cnt_coords, double min_value, double max_value) : _cnt_coords(cnt_coords), _coords(new T[_cnt_coords]) {
			std::mt19937 re;
			std::uniform_real_distribution<double> dist(min_value, max_value);
			for (size_t i = 0; i < cnt_coords; ++i) {
				_coords[i] = dist(re);
			}
		}
		//overloaded operators
		T& operator[] (size_t index)
		{
			if (index < 0 || index >= _cnt_coords) {
				throw std::out_of_range("invalid_index");
			}
			return _coords[index];
		}

		vector& operator=(const vector& second) {
			if (_cnt_coords != second._cnt_coords) {
				delete[] _coords;
				_coords = new T[second._cnt_coords];
			}
			_cnt_coords = second._cnt_coords;
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] = second[i];
			}
		}

		const T& operator[] (size_t index) const
		{
			if (index < 0 || index >= _cnt_coords) {
				throw std::out_of_range("invalid_index");
			}
			return _coords[index];
		}

		vector& operator+=(const vector& second) {
			if (second._cnt_coords != this->_cnt_coords) {
				throw std::invalid_argument("unequal_lengths");
			}
			for (size_t i = 0; i < _cnt_coords; ++i) {
				this->_coords[i] += second._coords[i];
			}
			return *this;
		}

		friend vector operator+(vector first, const vector& second) {
			return first += second;
		}

		vector& operator-=(const vector& second) {
			if (second._cnt_coords != this->_cnt_coords) {
				throw std::invalid_argument("unequal_lengths");
			}
			for (size_t i = 0; i < _cnt_coords; ++i) {
				this->_coords[i] -= second._coords[i];
			}
			return *this;
		}

		friend vector operator-(const vector& first, const vector& second) {
			auto copy(first);
			return copy -= second;
		}

		vector& operator*=(const double value) {
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] *= value;
			}
			return *this;
		}

		friend vector operator*(vector first, const double value) {
			return first *= value;
		}

		friend vector operator*(const double value, vector first) {
			return first *= value;
		}

		friend double operator*(const vector& first, const vector& second) {
			if (second._cnt_coords != first._cnt_coords) {
				throw std::invalid_argument("unequal_lengths");
			}
			double sum = 0.0;
			for (size_t i = 0; i < first._cnt_coords; ++i) {
				sum += first[i] * second[i];
			}
			return sum;
		}

		vector& operator/=(const double value) {
			if (value == 0) {
				throw std::invalid_argument("the_denominator_must_not_be_equal_to_0");
			}
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] /= value;
			}
			return *this;
		}

		friend vector operator/(vector first, const double value) {
			return first /= value;
		}

		friend bool operator==(const vector& first, const vector& second) {
			if (first._cnt_coords != second._cnt_coords) return false;
			for (size_t i = 0; i < first._cnt_coords; ++i) {
				if (fabs(first[i] - second[i]) > eps) return false;
			}
			return true;
		}

		friend bool operator!=(const vector& first, const vector& second) {
			return !(first == second);
		}

		friend std::ostream& operator<<(std::ostream& os, const vector& vect) {
			os << "{ ";
			for (size_t i = 0; i < vect._cnt_coords; ++i) {
				os << vect[i] << " ";
			}
			os << "}";
			return os;
		}
		//methods
		T get_coord(size_t i) const {
			if (i < 0 || i >= _cnt_coords) {
				throw std::out_of_range("invalid_index");
			}
			return _coords[i];
		}

		double get_length() const {
			double sum = 0;
			for (size_t i = 0; i < this->_cnt_coords; ++i) {
				sum += pow(this->_coords[i], 2);
			}
			return (double)(sqrtf(sum));
		}

		size_t get_cnt_coords() const {
			return _cnt_coords;
		}

		~vector() {
			delete[] _coords;
		}

	};
}


namespace my_vector {
	template<class T>
	class vector<std::complex<T>> {
		size_t _cnt_coords;
		std::complex<T>* _coords;
		static const inline double eps = 0.00000001;
	public:
		vector(size_t cnt_coords, const std::complex<T>* coords) : _cnt_coords(cnt_coords), _coords(new std::complex<T>[_cnt_coords]) {
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] = coords[i];
			}
		}

		vector(size_t cnt_coords, double min_value, double max_value) : _cnt_coords(cnt_coords), _coords(new std::complex<T>[_cnt_coords]) {
			std::mt19937 re;
			std::uniform_real_distribution<double> dist(min_value, max_value);
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] = std::complex<T>(dist(re), dist(re));
			}
		}

		vector(const vector<std::complex<T>>& copy) :_cnt_coords(copy._cnt_coords), _coords(new std::complex<T>[_cnt_coords]) {
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] = copy[i];
			}
		}

		//overloaded operators
		std::complex<T>& operator[] (size_t i)
		{
			if (i < 0 || i >= _cnt_coords) {
				throw std::out_of_range("invalid_index");
			}
			return _coords[i];
		}

		const std::complex <T>& operator[] (size_t i) const
		{
			if (i < 0 || i >= _cnt_coords) {
				throw std::out_of_range("invalid_index");
			}
			return _coords[i];
		}

		vector<std::complex<T>>& operator=(const vector<std::complex<T>>& second) {
			if (_cnt_coords != second._cnt_coords) {
				delete[] _coords;
				_coords = new std::complex<T>[second._cnt_coords];
			}
			_cnt_coords = second._cnt_coords;
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] = second[i];
			}
			return *this;
		}

		vector<std::complex<T>>& operator+=(const vector<std::complex<T>>& second) {
			if (second._cnt_coords != this->_cnt_coords) {
				throw std::invalid_argument("unequal_lengths");
			}
			for (size_t i = 0; i < _cnt_coords; ++i) {
				this->_coords[i] += second._coords[i];
			}
			return *this;
		}

		friend vector<std::complex<T>> operator+(vector<std::complex<T>> first, const vector<std::complex<T>>& second) {
			return first += second;
		}

		vector<std::complex<T>>& operator-=(const vector<std::complex<T>>& second) {
			if (second._cnt_coords != this->_cnt_coords) {
				throw std::invalid_argument("unequal_lengths");
			}
			for (size_t i = 0; i < _cnt_coords; ++i) {
				this->_coords[i] -= second._coords[i];
			}
			return *this;
		}

		friend vector<std::complex<T>> operator-(vector<std::complex<T>> first, const vector<std::complex<T>>& second) {
			return first -= second;
		}

		vector<std::complex<T>>& operator*=(const double value) {
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] *= value;
			}
			return *this;
		}

		friend vector<std::complex<T>> operator*(vector<std::complex<T>> first, const double value) {
			return first *= value;
		}

		friend vector<std::complex<T>> operator*(const double value, vector<std::complex<T>> first) {
			return first *= value;
		}

		friend std::complex<T> operator*(const vector<std::complex<T>>& first, const vector<std::complex<T>>& second) {
			if (second._cnt_coords != first._cnt_coords) {
				throw std::invalid_argument("unequal_lengths");
			}
			std::complex<T> sum = { 0.0,0.0 };
			for (size_t i = 0; i < first._cnt_coords; ++i) {
				sum += first[i] * (std::conj(second[i]));
			}
			return sum;
		}

		vector<std::complex<T>>& operator/=(const double value) {
			if (value == 0) {
				throw std::invalid_argument("the_denominator_must_not_be_equal_to_0");
			}
			for (size_t i = 0; i < _cnt_coords; ++i) {
				_coords[i] /= value;
			}
			return *this;
		}

		friend vector<std::complex<T>> operator/(vector first, const double value) {
			return first /= value;
		}

		friend bool operator==(const vector<std::complex<T>>& first, const vector<std::complex<T>>& second) {
			if (first._cnt_coords != second._cnt_coords) return false;
			for (size_t i = 0; i < first._cnt_coords; ++i) {
				if ((fabs(first[i].real() - second[i].real()) > eps) && (fabs(first[i].imag() - second[i].imag()) > eps)) return false;
			}
			return true;
		}

		friend bool operator!=(const vector<std::complex<T>>& first, const vector<std::complex<T>>& second) {
			return !(first == second);
		}

		friend std::ostream& operator<<(std::ostream& os, const vector<std::complex<T>>& vect) {
			os << "{ ";
			for (size_t i = 0; i < vect._cnt_coords; ++i) {
				os << vect[i] << " ";
			}
			os << "}";
			return os;
		}
		//methods
		std::complex<T> get_length() const {
			std::complex<T> sum = { 0.0, 0.0 };
			for (size_t i = 0; i < _cnt_coords; ++i) {
				sum += std::pow(_coords[i], 2);
			}
			return std::sqrt(sum);
		}

		std::complex<T> get_coord(int i) const {
			if (i < 0 || i >= _cnt_coords) {
				throw std::out_of_range("invalid_index");
			}
			return _coords[i];
		}

		size_t get_cnt_coords() const {
			return _cnt_coords;
		}
		~vector() {
			delete[] _coords;
		}
	};
}

template<typename T>
my_vector::vector<T> calculate_perpendicular_unit_vector(my_vector::vector<T> vector) {
	my_vector::vector<T> perpendicular_vector(vector.get_cnt_coords(), nullptr);
	if (std::is_same<T, std::complex<T>>::value) {
		perpendicular_vector[0] = -vector[1];
		perpendicular_vector[1] = vector[0];
	}
	else {
		perpendicular_vector[0] = -vector[1];
		perpendicular_vector[1] = vector[0];
	}
	

	T length = perpendicular_vector.get_length();
	for (size_t i = 0; i < perpendicular_vector.get_cnt_coords(); ++i) {
		perpendicular_vector[i] /= length;
	}

	return perpendicular_vector;
};
