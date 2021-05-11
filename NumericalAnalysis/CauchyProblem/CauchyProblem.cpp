/*
* Semester 2. Laboratory work No.4
* The Cauchy problem for ordinary differential equations.
* Heun's method. Method Runge-Kutta (third order accuracy).
*/
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

const int NUMBER_VARIANT = 5;
const int QUANTITY_VALUE = 11;
const int QUANTITY_FUNCTION = 2;
const double ALPHA = 2.0 + 0.5 * NUMBER_VARIANT;
const double STEP = 0.1;

// First function u'_1(t)
double first_function_system(const double& value_t, const double& first_value, const double& second_value) {
	return sin(ALPHA * pow(first_value, 2)) + value_t + second_value;
}

// Second function u'_2(t)
double second_function_system(const double& value_t, const double& first_value, const double& second_value) {
	return value_t + first_value - ALPHA * pow(second_value, 2) + 1.0;
}

void print_solution(const vector<vector<double>>& vectors_solution) {
	cout << left << setw(8) << "t" << setw(10) << "u_1(t)" << setw(10) << "u_2(t)" << endl;
	for (int index_solution = 0; index_solution < QUANTITY_VALUE; ++index_solution) {
		cout << left << setw(8) << index_solution * STEP << setw(10) << vectors_solution[index_solution][0]
			<< setw(10) << vectors_solution[index_solution][1] << endl;
	}
	cout << endl;
}

vector<vector<double>> method_Heuns() {
	vector<vector<double>> solutions(QUANTITY_VALUE, vector<double>(QUANTITY_FUNCTION));
	// Initial conditions
	// u_1(0) = 1
	// u_2(0) = 0.5
	solutions[0][0] = 1.0;
	solutions[0][1] = 0.5;
	// Find the following solutions in step 0.1, 0.2, ... 1.0
	for (int index_solution = 1; index_solution < QUANTITY_VALUE; ++index_solution) {
		// Method Euler`s
		// y*_{i+1} = y_i + h * f(t_i, y_i)
		double value_t = (index_solution - 1.0) * STEP;
		solutions[index_solution][0] = solutions[index_solution - 1][0] +
			STEP * first_function_system(value_t, solutions[index_solution - 1][0], solutions[index_solution - 1][1]);
		solutions[index_solution][1] = solutions[index_solution - 1][1] +
			STEP * second_function_system(value_t, solutions[index_solution - 1][0], solutions[index_solution - 1][1]);
		// Method Heun`s based on Euler`s method
		// y_{i+1} = y_i + h / 2 * [f(t_i, y_i) + f(t_{i+1}, y*_{i+1})] 
		double value_t_next = index_solution * STEP;
		solutions[index_solution][0] = solutions[index_solution - 1][0] + STEP / 2.0 *
			(first_function_system(value_t, solutions[index_solution - 1][0], solutions[index_solution - 1][1]) +
				first_function_system(value_t_next, solutions[index_solution][0], solutions[index_solution][1]));
		solutions[index_solution][1] = solutions[index_solution - 1][1] + STEP / 2.0 *
			(second_function_system(value_t, solutions[index_solution - 1][0], solutions[index_solution - 1][1]) +
				second_function_system(value_t_next, solutions[index_solution][0], solutions[index_solution][1]));
	}
	return solutions;
}

vector<vector<double>> method_Runge_Kutta() {
	vector<vector<double>> solutions(QUANTITY_VALUE, vector<double>(QUANTITY_FUNCTION));
	// Initial conditions
	// u_1(0) = 1
	// u_2(0) = 0.5
	solutions[0][0] = 1.0;
	solutions[0][1] = 0.5;
	// Find the following solutions in step 0.1, 0.2, ... 1.0
	for (int index_solution = 1; index_solution < QUANTITY_VALUE; ++index_solution) {
		double value_t = (index_solution - 1.0) * STEP;
		// For every function u'_1(t) and u'_2(t)
		// k_1 = f(t_n, y_n)
		// k_2 = f(t_n + h / 2, y_n + h / 2 * k_1)
		// k_3 = f(t_n + h, y_n - h * k_1 + 2 * h * k_2)
		// y_{n+1} = y_n + 1 / 6 * (k_1 + 4 * k_2 + k_3)
		double k1_first = first_function_system(value_t, solutions[index_solution - 1][0], solutions[index_solution - 1][1]);
		double k1_second = second_function_system(value_t, solutions[index_solution - 1][0], solutions[index_solution - 1][1]);
		double k2_frist = first_function_system(value_t + STEP / 2.0, solutions[index_solution - 1][0] + STEP / 2.0 * k1_first,
			solutions[index_solution - 1][1] + STEP / 2.0 * k1_second);
		double k2_second = second_function_system(value_t + STEP / 2.0, solutions[index_solution - 1][0] + STEP / 2.0 * k1_first,
			solutions[index_solution - 1][1] + STEP / 2.0 * k1_second);
		double k3_frist = first_function_system(value_t + STEP, solutions[index_solution - 1][0] - STEP * k1_first + 2.0 * STEP * k2_frist,
			solutions[index_solution - 1][1] - STEP * k1_second + 2.0 * STEP * k2_second);
		double k3_second = second_function_system(value_t + STEP, solutions[index_solution - 1][0] - STEP * k1_first + 2.0 * STEP * k2_frist,
			solutions[index_solution - 1][1] - STEP * k1_second + 2.0 * STEP * k2_second);
		solutions[index_solution][0] = solutions[index_solution - 1][0] + STEP * (k1_first + 4.0 * k2_frist + k3_frist) / 6.0;
		solutions[index_solution][1] = solutions[index_solution - 1][1] + STEP * (k1_second + 4.0 * k2_second + k3_second) / 6.0;
	}
	return solutions;
}

int main() {
	cout << setprecision(4);
	cout << "Solution DE received by method Heun`s: " << endl;
	vector<vector<double>> solution_Heuns = method_Heuns();
	print_solution(solution_Heuns);
	cout << "Solution DE received by method Runge-Kutta: " << endl;
	vector<vector<double>> solution_Runge_Kutta = method_Runge_Kutta();
	print_solution(solution_Runge_Kutta);
	return 0;
}