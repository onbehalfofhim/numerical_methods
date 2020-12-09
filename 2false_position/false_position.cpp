#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <istream>


struct point {
	point() = default;
	point(const point&) = default;
	~point() = default;
	point& operator=(const point&) = default;

	double x{ 0 };
	double y{ 0 };

	double f() const {
		return x * x * x - 3 * x * y * y - x * x + y * y + x - 1;
	}

	double g() const {
		return 3 * x * x * y - y * y * y - 2 * x * y + y;
	}

	std::ostream& write_to(std::ostream& ostrm) const {
		std::string str = "(" + std::to_string(x) + ", " + std::to_string(y) + "),";
		ostrm << str;
		return ostrm;
	}
};

inline std::ostream& operator<<(std::ostream& ostrm, const point& rhs) {
	return rhs.write_to(ostrm);
}

int sign(const double& x) {
	if (x < 0) return -1;
	if (x == 0) return 0;
	if (x > 0) return 1;
}


double S(const point& P, const point& Q, const point& R) {
	double a1(Q.x - P.x);
	double a2(R.x - P.x);
	double b1(Q.y - P.y);
	double b2(R.y - P.y);
	return std::abs(0.5 * (a1 * b2 - a2 * b1));
}

bool same_z_signs(const point& L, const point& R) {
	if (sign(L.f()) != sign(R.f()))
		return false;
	if (sign(L.g()) != sign(R.g()))
		return false;
	return true;
}

bool same_quater_pos(const point& L, const point& R) {
	if (sign(L.x) != sign(R.x))
		return false;
	if (sign(L.y) != sign(R.y))
		return false;
	return true;
}

bool same_quadrant(const point& L, const point& R) {
	if (!same_z_signs(L, R))
		return false;
	if (!same_quater_pos(L, R))\
		return false;
	return true;
}

void delete_point(point& P, point& Q, point& R, point& M) {
	// Удалить точку, которая лежит в одном квадранте с M, если такая есть
	if (same_quadrant(P, M)) {
		P = Q;
		Q = R;
		R = M;
	}
	else if (same_quadrant(Q, M)) {
		Q = R;
		R = M;
	}
	else if (same_quadrant(R, M)) {
		R = M;
	}
	// Если все точки лежат в разных квадрантах, выбрать две точки, расстояние до которых от M минимально
	// То есть площадь треугольника, образованного двумя точками и M, должна быть минимальна
	else {
		double s0 = S(Q, R, M);
		double s1 = S(P, R, M);
		double s2 = S(P, Q, M);

		if ((s0 >= s1) && (s0 >= s2)) {
			P = M;
		}
		else if ((s1 >= s0) && (s0 >= s2)) {
			Q = M;
		}
		else if ((s2 >= s0) && (s2 >= s1)) {
			R = M;
		}
	}
}

int false_pos(point P, point Q, point R, const double& eps) {
	point M0;
	size_t epoch(0);
	bool flag(true);

	while (flag) {
		// Определить три точки, лежащие на поверхностях z=f(x,y) и z=g(x,y)
		double zfp = P.f();		double zfq = Q.f();		double zfr = R.f();
		double zgp = P.g();		double zgq = Q.g();		double zgr = R.g();

		// Коэффициенты уравнения плоскости πf
		double Af = (Q.y - P.y) * (zfr - zfp) - (R.y - P.y) * (zfq - zfp);
		double Bf = (R.x - P.x) * (zfq - zfp) - (Q.x - P.x) * (zfr - zfp);
		double Cf = (Q.x - P.x) * (R.y - P.y) - (R.x - P.x) * (Q.y - P.y);
		double Df = -P.x * Af - P.y * Bf - zfp * Cf;

		// Коэффициенты уравнения плоскости πg
		double Ag = (Q.y - P.y) * (zgr - zgp) - (R.y - P.y) * (zgq - zgp);
		double Bg = (R.x - P.x) * (zgq - zgp) - (Q.x - P.x) * (zgr - zgp);
		double Cg = (Q.x - P.x) * (R.y - P.y) - (R.x - P.x) * (Q.y - P.y);
		double Dg = -P.x * Ag - P.y * Bg - zgp * Cg;

		// Найти точку пересечения плоскостей πf, πg и плоскости z=0
		double dt = Af * Bg - Ag * Bf;
		point M{ (-Df * Bg + Dg * Bf) / dt, (-Af * Dg + Ag * Df) / dt };
		double zfm(M.f());
		double zgm(M.g());

		// Проверка выполнения критерия
		double dv = std::abs(zfm) + std::abs(zgm);
		M0 = M;

		flag = (std::abs(zfm) > eps || std::abs(zgm) > eps) && epoch < 3000;

		delete_point(P, Q, R, M);
		++epoch;
	}

	if ((std::abs(M0.x - 0) <= eps) && (std::abs(M0.y - 1) <= eps))
		return 1;
	else if ((std::abs(M0.x - 0) <= eps) && (std::abs(M0.y + 1) <= eps))
		return -1;
	else if ((std::abs(M0.x - 1) <= eps) && (std::abs(M0.y - 0) <= eps))
		return 0;
	return 2;
}

int main() {
	double eps = 0.01;
	point A{ -2,    2 };
	point B{ -0.4,  2.2 };
	point C{ 0.6,  1.4 };
	point D{ 1.6,  0.4 };
	point E{ 2.2, -0.6 };
	point F{ 0.4, -0.4 };
	point G{ 0.4, -1.4 };
	point H{ -0.4, -2 };
	point I{ -1,   -1.6 };
	point J{ -0.6,  0.4 };

	std::vector<std::pair<point, char>> dots;
	dots.push_back(std::make_pair(A, 'A'));
	dots.push_back(std::make_pair(B, 'B'));
	dots.push_back(std::make_pair(C, 'C'));
	dots.push_back(std::make_pair(D, 'D'));
	dots.push_back(std::make_pair(E, 'E'));
	dots.push_back(std::make_pair(F, 'F'));
	dots.push_back(std::make_pair(G, 'G'));
	dots.push_back(std::make_pair(H, 'H'));
	dots.push_back(std::make_pair(I, 'I'));
	dots.push_back(std::make_pair(J, 'J'));


	std::vector<char> a1, a2, a3;
	int i(0), j(0), k(0);
	for (i; i < 8; ++i)
		for (j = i + 1; j < 9; ++j)
			for (k = j + 1; k < 10; ++k) {
				int result = false_pos(dots.at(i).first, dots.at(j).first, dots.at(k).first, eps);
				if (result == -1) {
					a1.push_back(dots.at(i).second); a1.push_back(dots.at(j).second); a1.push_back(dots.at(k).second);
				}
				else if (result == 1) {
					a2.push_back(dots.at(i).second); a2.push_back(dots.at(j).second); a2.push_back(dots.at(k).second);
				}
				else if (result == 0) {
					a3.push_back(dots.at(i).second); a3.push_back(dots.at(j).second); a3.push_back(dots.at(k).second);
				}
			}



	std::cout << "(0, -1): ";
	for (int i(0); i < a1.size(); ++i) {
		std::cout << a1[i];
		if (i % 3 == 2)
			std::cout << ", ";
	}
	std::cout << "\n(0,  1): ";
	for (int i(0); i < a2.size(); ++i) {
		std::cout << a2[i];
		if (i % 3 == 2)
			std::cout << ", ";
	}
	std::cout << "\n(1,  0): ";
	for (int i(0); i < a3.size(); ++i) {
		std::cout << a3[i];
		if (i % 3 == 2)
			std::cout << ", ";
	}

	std::cin.get();
}