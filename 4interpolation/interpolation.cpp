#include <iostream>
#include <vector>

using namespace std;

class point {
public:
	point() = default;
	point(const point&) = default;
	point& operator=(const point&) = default;
	~point() = default;

	point(double x_, double y_) { x = x_; y = y_; }

	double getX() { return x; }
	double getY() { return y; }

	std::ostream& write_to(std::ostream& ostrm) const {
		ostrm << '(' << x << ", " << y << ')';
		return ostrm;
	}

private:
	double x;
	double y;
};

inline std::ostream& operator<<(std::ostream& ostrm, const point& rhs) {
	return rhs.write_to(ostrm);
}

void function(vector<point>& data) {
	cout << "Сеточная функция:\n";
	for (size_t i(0); i < data.size(); ++i) {
		cout << data[i];
		(i != data.size()-1) ? (cout << ", ") : (cout << "\n\n\n");
	}
}

void Lagrange_output(vector<point>& data) {
	size_t n(data.size() - 1);

	vector<double> znam;

	cout << "Интерполяционный многочлен Лагранжа:\nL_" << n + 1 << "(x) = ";
	for (size_t j(0); j <= n; ++j) {
		znam.push_back(1);
		for (size_t i(0); i <= n; ++i) {
			if (j == i) continue;

			if(data[i].getX() > 0) cout << "(x - " << data[i].getX() << ") * ";
			else if (data[i].getX() == 0) cout << "x * ";
			else cout << "(x + " << abs(data[i].getX()) << ") * ";

			znam.back() *= data[j].getX() - data[i].getX();
		}
		(data[j].getY() > 0) ?
			(cout << data[j].getY()) :
			(cout << "(-" << abs(data[j].getY()) << ")");
		(znam.back() > 0) ?
			(cout << " / " << znam.back()) :
			(cout << " / (" << znam.back() << ")");
		if (j != n) cout << " + ";
	}
}

void Lagrange_count(vector<point>& data, double xp) {
	size_t n(data.size() - 1);

	double chesl;
	double znam;

	double result(0);

	for (size_t j(0); j <= n; ++j) {
		chesl = 1;
		znam = 1;
		for (size_t i(0); i <= n; ++i) {
			if (j == i) continue;

			chesl *= xp - data[i].getX();
			znam *= data[j].getX() - data[i].getX();
		}

		result += data[j].getY()*chesl / znam;
	}

	point r(xp, result);

	cout << "Значение интерполяционного многочлена Лагранжа в точке:";
	cout << r;
	
}

void Newton_output(vector<point>& data) {
	size_t n(data.size() - 1);

	vector<double> znam;
	cout << "Интерполяционный многочлен Ньютона:\nL_" << n + 1 << "(x) = ";
	cout << data[0].getY() << " + ";

	vector<double> delta;
	// Разделенные разности 0-го порядка совпадают со значениями функции в узлах
	for (auto& it : data)
		delta.push_back(it.getY());

	for (size_t m(1); m <= n; ++m) {
		// Разделенные разности k-го порядка определяются через разделенные разности порядка k-1
		vector<double> delta_new;
		for (size_t i(0); i <= n - m; ++i) {
			delta_new.push_back((delta[i + 1] - delta[i]) / (data[i + m].getX() - data[i].getX()));
		}
		delta = delta_new;

		(delta[0] > 0) ?
			(cout << delta[0] << " * ") :
			(cout << "(-" << abs(delta[0]) << ") * ");

		for (int i(0); i <= m - 1; ++i) {
			if (data[i].getX() > 0) cout << "(x - " << data[i].getX() << ")";
			else if (data[i].getX() == 0) cout << "x";
			else cout << "(x + " << abs(data[i].getX()) << ")";
			if (i != m - 1) cout << " * ";
		}
		if (m != n) cout << " + ";
	}
}

void Newton_count(vector<point>& data, double xp) {
	size_t n(data.size() - 1);

	double result(data[0].getY());
	double mnogitel;

	vector<double> delta;

	for (auto& it : data)
		delta.push_back(it.getY());

	for (size_t m(1); m <= n; ++m) {
		vector<double> delta_new;
		for (size_t i(0); i <= n - m; ++i) {
			delta_new.push_back((delta[i + 1] - delta[i]) / (data[i + m].getX() - data[i].getX()));
		}
		delta = delta_new;
		mnogitel = 1;
		for (int i(0); i <= m - 1; ++i) {
			mnogitel *= xp - data[i].getX();
		}
		result += delta[0] * mnogitel;
	}

	point r(xp, result);

	cout << "Значение интерполяционного многочлена Ньютона в точке:";
	cout << r;
}

int main() {
	setlocale(LC_ALL, "Russian");

	vector<point> data;
	data.push_back(point(-2, -5.68));
	data.push_back(point(-1, -5.89));
	data.push_back(point(0, 3.16));
	data.push_back(point(1, 4.11));
	data.push_back(point(2, -7.24));
	data.push_back(point(3, -5.14));
	data.push_back(point(4, 6.73));

	function(data);

	Lagrange_output(data);
	cout << "\n\n\n";
	Newton_output(data);
	cout << "\n\n\n";

	Lagrange_count(data, -1.3);
	cout << "\n";
	Lagrange_count(data, 1.75);
	cout << "\n\n\n";

	Newton_count(data, -1.3);
	cout << "\n";
	Newton_count(data, 1.75);

	cin.get();
	return 0;
}