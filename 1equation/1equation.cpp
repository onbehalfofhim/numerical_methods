#include <iostream>
#include <cmath>

using namespace std;
float func(float x) {
	return x * pow(x / (4 - x), 1 / 2) - 3 / tan(x);
}

void bisection(double a, double b) {
	int count(0);
	double e = 0.01;
	while (e > 0.00001) {
		cout << e << endl;
		while (b - a >= e) {
			if (func(a)*func((b + a) / 2) == 0)					//������ ������
				break;
			else if (func(a)*func((b + a) / 2) > 0) {			//��������� ���� ������� � �������� ���������
				a = (b + a) / 2;								//���� ������ ����, ������ ����� ������� ���������
				count++;
				cout << count << ": x(left) = " << a << "; x(right) = " << b << endl;
			}
			else {
				b = (b + a) / 2;								//���� ���� ������ 0, �� ������ ������ ������� ���������
				count++;
				cout << count << ": x(left) = " << a << "; x(right) = " << b << endl;
			}
		}
		//cout << "f(left) = " << func(a) << "; f(right) = " << func(b) << endl;
		e = e / 10;
	}
}

void secant(double a, double b) {
	int count(0);
	double e = 0.01;
	double res = 0;
	while (e > 0.00001) {
		cout << e << endl;
		while (fabs(a-b) >= e) {
			res = b * func(a) / (func(a) - func(b)) + a * func(b) / (func(b) - func(a)); //��������� ��������� �����������
			a = b;																//������������� �������� (x(i+1), x(i))
			b = res;
			count++;
			cout << count << ": x(left) = " << a << "; x(right) = " << b << endl;
		}
		//cout << "f(left) = " << func(a) << "; f(right) = " << func(b) << endl;
		e = e / 10;
	}
}

void false_position(double a, double b) {
	int count(0);
	double e = 0.01;
	double res;
	double delx = b - a;
	while (e > 0.00001) {
		cout << e << endl;
		while (delx >= e) {
			res = b * func(a) / (func(a) - func(b)) + a * func(b) / (func(b) - func(a));		//��������� ��������� �����������

			if (func(a)*func(res) > 0) {						//��������� ����� ����� � �������������� ������ ������� ���������
				delx = res - a;
				a = res;
			}else{
				delx = b - res;
				b = res;
			}
			count++;
			cout << count << ": x(left) = " << a << "; x(right) = " << b << ", delta = " << delx << endl;
		}
		//cout << "f(left) = " << func(a) << "; f(right) = " << func(b) << endl;
		e = e / 10;
	}
}

double df(double x) {				//���� �����������
	double h = 0.000001;
	return (func(x + h) - func(x)) / h;
}

void newton(double a) {
	int count(0);
	double b = a - func(a) / df(a);			//��������� ������ �������� ���������
	double e = 0.01;
	while (e > 0.00001) {
		cout << e << endl;
		while (fabs(b - a) > e) {  //��������� �������� ���������
			a = b;
			b = a - func(a) / df(a);		//��������� ��������� �����������
			count++;
			cout << count << ": x(left) = " << a << "; x(right) = " << b << endl;
		}
		//cout << "f(left) = " << func(a) << "; f(right) = " << func(b) << endl;
		e = e / 10;
	}
}

int main() {
	bisection(1, 2);
	cout << endl;
	secant(1, 2);
	cout << endl;
	false_position(1, 2);
	cout << endl;
	newton(1);
	cin.get();
	cin.get();
	return 0;
}