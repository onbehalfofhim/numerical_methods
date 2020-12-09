#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

int main() {
	setlocale(LC_ALL, "Russian");

	int a = 1;
	int b = 3;
	int n = 6;
	double h = double(b - a) / n;
	vector<double> x{};
	vector<double> y;
	cout << "h: " << h << "\n";
	x.push_back(a);
	for (int i = 1; i <= n; ++i) {
		x.push_back(h + x[i - 1]);
	}

	for (int i = 0; i < n+1; ++i) {
		y.push_back(1 / sqrt((1 + x[i] * x[i])*log(x[i] + sqrt(1 + x[i] * x[i]))));

		cout << "x[" << i << "] =" << setw(8) << x[i] << '\t' << setw(8) << "y[" << i << "] = " << y[i] << "\n";
	}
	cout << "\n";

	double Itrap = h * ((y[0] + y[6]) / 2 + y[1] + y[2] + y[3] + y[4] + y[5]);
	cout << "С помощью формулы трапеции: " << Itrap << "\n";
	cout << "\n";

	double Isimps = h / 3 * (y[0] + y[6] + 4 * (y[1] + y[3] + y[5]) + 2 * (y[2] + y[4]));
	cout << "С помощью формулы Симпсона: " << Isimps << "\n";
	cout << "\n";

	vector<double> t = { -0.906180, -0.538469, 0, 0.538469, 0.906180 };
	vector<double> A = { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 };
	for (int i = 0; i < 5; i++) {
		cout << "t[" << i << "]: " << setw(8) << t[i] << '\t' << "A[" << i << "]: " << setw(8) << A[i] << "\n";
	}
	cout << "\n";

	double Igauss = 0;
	vector<double> x1;
	for (int i = 0; i < 5; i++) {
		x1.push_back((a + b) / 2 + ((b - a)*t[i]) / 2);
		cout << "х1[" << i << "]= " << setw(8) << x1[i] << '\t' << setw(8) << "y1[" << i << "] = " << A[i] * 1 / sqrt((1 + x1[i] * x1[i])*log(x1[i] + sqrt(1 + x1[i] * x1[i]))) << "\n";
		Igauss = Igauss + A[i] * 1 / sqrt((1 + x1[i] * x1[i])*log(x1[i] + sqrt(1 + x1[i] * x1[i])));
	}
	cout << "\n";

	cout << "С помощью квадратурной формулы Гаусса: " << Igauss << "\n";

	cin.get();

	return 0;
}