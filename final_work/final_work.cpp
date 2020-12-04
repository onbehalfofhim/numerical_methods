#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;
class matrix {
public:
	matrix() = default;
	matrix(const matrix&) = default;
	~matrix() = default;
	matrix(const int& rows, const int& cols, const vector<double>& data_);

	double& at(const int& row_i, const int& col_i);
	const double& at(const int& row_i, const int& col_i) const;

	matrix& operator=(const matrix&);
	matrix& operator*=(const matrix& right);

	void print() const;

	double length();
	int cols() const { return n_cols; }
	int rows() const { return n_rows; }

	double left() const { return left_end; }
	double right() const  { return right_end; }
	void set_left_end(double a) { left_end = a; }
	void set_right_end(double b) { right_end = b; }

private:
	double left_end{ -10 };
	double right_end{ 10 };
	int n_rows{ 0 };
	int n_cols{ 0 };
	vector<double> data;
};

inline matrix operator*(const matrix& left, const matrix& right) {
	return matrix(left) *= right;
 }

inline matrix operator*(const matrix& left, double element) {
	vector<double> t;
	for (int i(0); i < left.rows(); ++i) {
		for (int j(0); j < left.cols(); ++j) {
			t.push_back(left.at(i + 1, j + 1)*element);
		}
	}

	return matrix(left.rows(), left.cols(), t);
}

inline matrix operator-(const matrix& left, const matrix& right) {
	vector<double> t;
	for (int i(0); i < left.rows(); ++i) {
		for (int j(0); j < left.cols(); ++j) {
			t.push_back(left.at(i + 1, j + 1)-right.at(i+1,j+1));
		}
	}

	return matrix(left.rows(), left.cols(), t);
}

inline matrix operator+(const matrix& left, const matrix& right) {
	vector<double> t;
	for (int i(0); i < left.rows(); ++i) {
		for (int j(0); j < left.cols(); ++j) {
			t.push_back(left.at(i + 1, j + 1) + right.at(i + 1, j + 1));
		}
	}

	return matrix(left.rows(), left.cols(), t);
}

matrix::matrix(const int& rows, const int& cols, const vector<double>& data_) : n_rows(rows), n_cols(cols), data(data_) {
	if (data.size() != n_rows * n_cols)
		throw invalid_argument("Change the input vector");
}

matrix& matrix::operator=(const matrix& rhs) {
	if (this == &rhs)
		return *this;
	n_cols = rhs.cols();
	n_rows = rhs.rows();
	data.clear();
	for (int i(0); i < n_rows; ++i) {
		for (int j(0); j < n_cols; ++j) {
			data.push_back(rhs.at(i + 1, j + 1));
		}
	}
	return *this;
}

double& matrix::at(const int& row_i, const int& col_i) {
	if ((row_i > n_rows) || (col_i > n_cols) || (row_i <= 0) || (col_i <= 0)) {
		throw out_of_range("Invalid index");
	}
	return data[(row_i - 1) * n_cols + (col_i - 1)];
}

const double& matrix::at(const int& row_i, const int& col_i) const {
	if ((row_i > n_rows) || (col_i > n_cols) || (row_i <= 0) || (col_i <= 0)) {
		throw out_of_range("Invalid index");
	}
	return data[(row_i - 1) * n_cols + (col_i - 1)];
}

void matrix::print() const {
	for (int i(1); i <= n_rows; ++i) {
		for (int j(1); j <= n_cols; ++j)
			cout << setw(9) << setprecision(5) << at(i, j);
		cout << "\n";
	}
	cout << "\n";
}

matrix& matrix::operator*=(const matrix& right) {
	vector<double> c;
	for (int i(0); i < n_rows; ++i) {
		for (int j(0); j < right.cols(); ++j) {
			c.push_back(0);
			for (int k(0); k < n_cols; ++k) {
				c.back() += at(i + 1, k + 1) * right.at(k + 1, j + 1);
			}
		}
	}
	matrix C{ n_rows, right.cols(), c };
	*this = C;
	return *this;
}

double matrix::length() {
	if (n_cols != 1) {
		throw invalid_argument("It is not a vector, but a matrix! Be careful!");
	}
	double len(0);
	for (int i(0); i < n_rows; ++i) {
		len += data[i] * data[i];
	}
	return sqrt(len);
}


matrix Froben_make(const matrix& Arg, vector<vector<double>>& F) {
	int k(2);
	matrix A(Arg);
	while (A.rows() - k >= 0) {
		vector<double> b;     //формирование матриц B1, B2, B3 и тд.
		for (int i(0); i < A.rows(); ++i) {
			for (int j(0); j < A.cols(); ++j) {
				if (i == A.rows() - k) {
					if (i != j)  b.push_back(-A.at(A.rows() - k + 2, j + 1) / A.at(A.rows() - k + 2, A.rows() - k + 1));
					else b.push_back(1 / A.at(A.rows() - k + 2, A.rows() - k + 1));
				}
				else {
					if (i == j) b.push_back(1);
					else b.push_back(0);
				}
			}
		}
		matrix B(A.rows(), A.cols(), b );
		F.push_back(b);

		vector<double> b_1;                       //формирование матриц B1^-1, B2^-1, B3^-1 и тд.
		for (int i(0); i < A.rows(); ++i) {
			for (int j(0); j < A.cols(); ++j) {
				if (i == A.rows() - k) {
					b_1.push_back(A.at(A.rows() - k + 2, j + 1));
				}
				else {
					if (i == j) b_1.push_back(1);
					else b_1.push_back(0);
				}
			}
		}
		matrix B_1(A.rows(), A.cols(), b_1);

		A = B_1*A*B;
		++k;
	}
	return A;
}

vector<double> find_roots(const matrix& A, vector<double> cof) {
	vector<double> x;
	double a(A.left());
	double b(a + 1);
	while (b < A.right()) {
		double f_a(0);
		double f_b(0);
		for (int i(0); i < cof.size(); ++i) {
			f_a += (-1)*cof[i] * pow(a, cof.size() - 1 - i);
			f_b += (-1)*cof[i] * pow(b, cof.size() - 1 - i);
		}
		if (f_a*f_b > 0) {
			++b;
			++a;
		}
		else {
			double left(a);
			double right(b);
			double c(0);
			while ((right - left) / 2 > 0.0001) {
				c = (left + right) / 2;
				double f_a(0);
				double f_c(0);
				for (int i(0); i < cof.size(); ++i) {
					f_a += (-1)*cof[i] * pow(right, cof.size() - 1 - i);
					f_c += (-1)*cof[i] * pow(c, cof.size() - 1 - i);
				}
				if (f_a*f_c > 0) right = c;
				else left = c;
			}
			x.push_back(c);
			++a;
			++b;
		}
	}

	return x;
}

void Danilevsky(const matrix &A) {
	vector<vector<double>> f;
	matrix P(Froben_make(A, f));

	vector<double> p; //коэффициенты характеристического полинома
	p.push_back(-1);
	for (int i(0); i < P.cols(); ++i) {
		p.push_back(P.at(1, i + 1));
	}

	vector<double> x = find_roots(A, p);		//матодом хорд ищем все корни характерисчтического полинома на заданном отрезке. Отрезок по-умолчанию [-10,10]
	cout << "Собственные значения матрицы A:" << '\n';
	for (int i(0); i < x.size()-1; ++i) {
		cout << "l_" << i + 1 << " = " << x[i] << ",   ";
	}
	cout << "l_" << x.size() << " = " << x[x.size() - 1] << "\n\n";
	
	//найдем матрицу подобия В = B_1*B_2*B_3...B_(n-1)
	matrix B(A.rows(), A.cols(), f[0]);
	for (int i(1); i < f.size(); ++i) {
		B *= matrix(A.rows(), A.cols(), f[i]);
	}
	
	//найдем собственный векторы матрицы P
	vector<vector<double>> y;
	for (int i(0); i < x.size(); ++i) {
		vector<double> t;
		for (int j(1); j <= A.rows(); ++j) {
			t.push_back(pow(x[i], A.rows() - j));
		}
		y.push_back(t);
	}

	//находим собственные векторы матрицы A: X_i = B*y_i
	vector<matrix> X;
	for (int i(0); i < y.size(); ++i) {
		X.push_back(matrix(B*matrix(y[i].size(), 1, y[i])));
	}

	cout << "Собственные векторы матрицы A:" << '\n';
	for (int j(0); j < X.size(); ++j) {
		cout << "x" << j + 1 << ":\n";
		X[j].print();
	}

	//нормируем полученные векторы
	vector<vector<double>> x_ort;
	for (int i(0); i < X.size(); ++i) {
		vector<double> t;
		for (int j(0); j < X[i].rows(); ++j) {
			t.push_back(X[i].at(j + 1, 1) / X[i].length());
		}
		x_ort.push_back(t);
	}

	vector<matrix> X_ort;		//ортонормированная система векторов матрица А
	for (int i(0); i < x_ort.size(); ++i) {
		X_ort.push_back(matrix(x_ort[i].size(), 1, x_ort[i]));
	}

	cout << "Ортонормированная система собственных векторов матрицы A:" << '\n';
	for (int j(0); j < X.size(); ++j) {
		cout << "x" << j + 1 << ":\n";
		X_ort[j].print();
	}

}

void Leverrie_Fadeev(const matrix& Ar) {
	vector<double> e;
	for (int i(0); i < Ar.rows(); ++i) {
		for (int j(0); j < Ar.cols(); ++j) {
			if (i != j) e.push_back(0);
			else e.push_back(1);
		}
	}

	vector<double> p;		//сюда собираем коэффициенты
	p.push_back(-1);
	vector<matrix> B;		//сюда собираем первые столбца матриц B

	int k(1);
	matrix A(Ar);
	while (k <= Ar.rows()) {
		double temp(0);
		matrix E(A.rows(), A.cols(), e);

		for (int i(0); i < A.rows(); ++i) {			//считаем сумму диагональных элементов матрицы A_i
			for (int j(0); j < A.cols(); ++j) {
				if (i == j)  temp += A.at(i + 1, j + 1);
			}
		}
		p.push_back(temp / k);			//записываем коэффициент
		matrix B_i(A - E * p.back());	//вычисляем B_i матрицу


		vector<double> b;			//сохраняем первые столбцы найденных матриц В, так как они понадобятся для вычисления собственных векторов
		for (int i(0); i < B_i.rows(); ++i)
			b.push_back(B_i.at(i + 1, 1));
		B.push_back(matrix(B_i.rows(), 1, b));

		A = Ar*B_i;		//вычисляем новую A_i матрицу
		++k;
	}

	B.pop_back();		//убираем столбец матрицы B_n,так как он не нужен в дальнейших вычислениях

	vector<double> x = find_roots(Ar, p);		//матодом хорд ищем все корни характерисчтического полинома на заданном отрезке. Отрезок по-умолчанию [-10,10]
	cout << "Собственные значения матрицы A:" << '\n';
	for (int i(0); i < x.size() - 1; ++i) {
		cout << "l_" << i + 1 << " = " << x[i] << ",   ";
	}
	cout << "l_" << x.size() << " = " << x[x.size() - 1] << "\n\n";

	vector<matrix> X;
	for (int i(0); i < x.size(); ++i) {
		vector<matrix> y;
		vector<double> t;		//формируем y_0
		t.push_back(1);
		for (int j(1); j < Ar.rows(); ++j)
			t.push_back(0);
		y.push_back(matrix(Ar.rows(), 1, t));

		for (int k(0); k < B.size(); ++k) {		//вычисляем y_n
			y.push_back(matrix(y.back()*x[i] + B[k]));
		}

		X.push_back(y.back());
	}

	cout << "Собственные векторы матрицы A:" << '\n';
	for (int j(0); j < X.size(); ++j) {
		cout << "x" << j + 1 << ":\n";
		X[j].print();
	}

	//нормируем полученные векторы
	vector<vector<double>> x_ort;
	for (int i(0); i < X.size(); ++i) {
		vector<double> t;
		for (int j(0); j < X[i].rows(); ++j) {
			t.push_back(X[i].at(j + 1, 1) / X[i].length());
		}
		x_ort.push_back(t);
	}

	vector<matrix> X_ort;		//ортонормированная система векторов матрица А
	for (int i(0); i < x_ort.size(); ++i) {
		X_ort.push_back(matrix(x_ort[i].size(), 1, x_ort[i]));
	}

	cout << "Ортонормированная система собственных векторов матрицы A:" << '\n';
	for (int j(0); j < X.size(); ++j) {
		cout << "x" << j + 1 << ":\n";
		X_ort[j].print();
	}
}

int main() {
	setlocale(LC_ALL, "Russian");

	matrix A( 4,4, vector<double>{
		2.2, 1, 0.5, 2,
		1, 1.3, 2, 1,
		0.5, 2, 0.5, 1.6,
		2, 1, 1.6, 2
	} );

	matrix A1(3, 3, vector<double>{
		5, 1, 2, 
		1, 4, 1,
		2, 1, 3
	});

	cout << "Матрица А:" << '\n';
	A.print();
	cout << "---------------------------------------------------------------------------------------------\n";

	cout << "Поиск собственных значений и собственных векторов методом А.М.Данилевского" << '\n';
	Danilevsky(A);
	cout << "---------------------------------------------------------------------------------------------\n";

	cout << "Поиск собственных значений и собственных векторов методом Леверрье-Фадеева" << '\n';
	Leverrie_Fadeev(A);

	cin.get();
	return 0;
}