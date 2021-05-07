/*
* Laboratory work No.3
* Solving systems of nonlinear equations.
* Simple-iteration method. Newton method.
*/
#include <iostream>
#include <vector>

using namespace std;

const double EPSILON_SIMPLE = 1e-3;
const double EPSILON_NEWTON = 1e-6;
const int QUANTITY_EQUATIONS = 3;

struct Solution {
	double value_x;
	double value_y;
	double value_z;
};

Solution operator-(const Solution& lhs_solution, const Solution& rhs_solution) {
	return { lhs_solution.value_x - rhs_solution.value_x,
			 lhs_solution.value_y - rhs_solution.value_y,
			 lhs_solution.value_z - rhs_solution.value_z };
}

void print_solution(const Solution& solution, const string& name_output) {
	char char_output = 'x';
	if (name_output == "Residual") {
		char_output = 'r';
	}
	cout << char_output << "_1: " << fixed << solution.value_x << endl;
	cout << char_output << "_2: " << fixed << solution.value_y << endl;
	cout << char_output << "_3: " << fixed << solution.value_z << endl;
	cout << endl;
}

double first_func(const Solution& solution) { return 0.1 - pow(solution.value_x, 2) + 2.0 * solution.value_y * solution.value_z; }

double second_func(const Solution& solution) { return -0.2 + pow(solution.value_y, 2) - 3.0 * solution.value_x * solution.value_z; }

double third_func(const Solution& solution) { return 0.3 - pow(solution.value_z, 2) - 2.0 * solution.value_x * solution.value_y; }

double second_vector_norm(const Solution& norm_solution) {
	return sqrt(pow(norm_solution.value_x, 2) + pow(norm_solution.value_y, 2) + pow(norm_solution.value_z, 2));
}

Solution find_residual_vector(Solution& solution) {
	return { -first_func(solution) + solution.value_x,
			 -second_func(solution) + solution.value_y,
			 -third_func(solution) + solution.value_z };
}

void simple_iteration_method(Solution& solution) {
	int count_iteration = 0;
	Solution new_solution = solution;
	do {
		solution = new_solution;
		new_solution = { first_func(solution), second_func(solution), third_func(solution) };
		++count_iteration;
	} while (second_vector_norm(new_solution - solution) > EPSILON_SIMPLE);
	cout << "Quantity iteration in simple-iteration method: " << count_iteration << endl;
	cout << endl << "Solution obtained by Simple-iteration method: " << endl;
	print_solution(new_solution, "Solution");
	cout << "Residual vector obtained by Simple-iteration method:  " << endl;
	print_solution(find_residual_vector(new_solution), "Residual");
}

int main() {
	Solution solution = { 0.0, -0.17, 0.24 };
	simple_iteration_method(solution);
	return 0;
}