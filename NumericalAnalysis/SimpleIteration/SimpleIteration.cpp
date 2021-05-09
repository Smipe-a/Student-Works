/*
* Laboratory work No.3
* Solving systems of nonlinear equations.
* Simple-iteration method. Newton method.
*/
#include <iostream>
#include <iomanip>
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

// Overloading the subtraction operator for vectors (3x1)
Solution operator-(const Solution& lhs_solution, const Solution& rhs_solution) {
	return { lhs_solution.value_x - rhs_solution.value_x,
			 lhs_solution.value_y - rhs_solution.value_y,
			 lhs_solution.value_z - rhs_solution.value_z };
}

// Overloading the multiplication operator for matrix (3x3) and vector (3x1)
Solution operator*(const vector<vector<double>>& lhs_matrix, const Solution& rhs_vector) {
	vector<double> multiply_matrix(QUANTITY_EQUATIONS);
	for (int index_row = 0; index_row < QUANTITY_EQUATIONS; ++index_row) {
		multiply_matrix[index_row] = lhs_matrix[index_row][0] * rhs_vector.value_x +
									 lhs_matrix[index_row][1] * rhs_vector.value_y +
									 lhs_matrix[index_row][2] * rhs_vector.value_z;
	}
	return { multiply_matrix[0], multiply_matrix[1], multiply_matrix[2] };
}

void print_solution(const Solution& solution, const string& name_output) {
	cout << setprecision(12);
	char char_output = 'x';
	if (name_output == "Residual") { char_output = 'r'; }
	cout << char_output << "_1: " << fixed << solution.value_x << endl;
	cout << char_output << "_2: " << fixed << solution.value_y << endl;
	cout << char_output << "_3: " << fixed << solution.value_z << endl;
	cout << endl;
}

void print_inverse_jacobian(const vector<vector<double>>& jacobian) {
	cout << setprecision(3);
	for (int index_row = 0; index_row < QUANTITY_EQUATIONS; ++index_row) {
		for (int index_col = 0; index_col < QUANTITY_EQUATIONS; ++index_col) {
			cout << setw(8) << left << fixed << jacobian[index_row][index_col];
		}
		cout << endl;
	}
	cout << endl;
}

vector<vector<double>> find_inverse_matrix(vector<vector<double>> lhs_matrix) {
	vector<vector<double>> unit_matrix = { {1.0, 0.0, 0.0},
										   {0.0, 1.0, 0.0},
										   {0.0, 0.0, 1.0} };
	for (int index_row = 0; index_row < QUANTITY_EQUATIONS; ++index_row) {
		double permitting_element = lhs_matrix[index_row][index_row];
		for (int index_col = 0; index_col < QUANTITY_EQUATIONS; ++index_col) {
			lhs_matrix[index_row][index_col] /= permitting_element;
			unit_matrix[index_row][index_col] /= permitting_element;
		}
		for (int index_subrow = index_row + 1; index_subrow < QUANTITY_EQUATIONS; ++index_subrow) {
			permitting_element = lhs_matrix[index_subrow][index_row];
			for (int index_col = 0; index_col < QUANTITY_EQUATIONS; ++index_col) {
				lhs_matrix[index_subrow][index_col] -= lhs_matrix[index_row][index_col] * permitting_element;
				unit_matrix[index_subrow][index_col] -= unit_matrix[index_row][index_col] * permitting_element;
			}
		}
	}
	for (int index_subrow = QUANTITY_EQUATIONS - 1; index_subrow > 0; --index_subrow) {
		for (int index_row = index_subrow - 1; index_row >= 0; --index_row) {
			double permitting_element = lhs_matrix[index_row][index_subrow];
			for (int index_col = 0; index_col < QUANTITY_EQUATIONS; ++index_col) {
				lhs_matrix[index_row][index_col] -= lhs_matrix[index_subrow][index_col] * permitting_element;
				unit_matrix[index_row][index_col] -= unit_matrix[index_subrow][index_col] * permitting_element;
			}
		}
	}
	return unit_matrix;
}

double first_func(const Solution& solution) { return 0.1 - pow(solution.value_x, 2) + 2.0 * solution.value_y * solution.value_z; }

double second_func(const Solution& solution) { return -0.2 + pow(solution.value_y, 2) - 3.0 * solution.value_x * solution.value_z; }

double third_func(const Solution& solution) { return 0.3 - pow(solution.value_z, 2) - 2.0 * solution.value_x * solution.value_y; }

double second_vector_norm(const Solution& norm_solution) {
	return sqrt(pow(norm_solution.value_x, 2) + pow(norm_solution.value_y, 2) + pow(norm_solution.value_z, 2));
}

Solution find_residual_vector(const Solution& solution) {
	return { solution.value_x - first_func(solution),
			 solution.value_y - second_func(solution),
			 solution.value_z - third_func(solution) };
}

vector<vector<double>> find_matrix_jacobi(const Solution& solution) {
	return { {1 + 2 * solution.value_x, -2 * solution.value_z, -2 * solution.value_y},
			 {3 * solution.value_z, 1 - 2 * solution.value_y, 3 * solution.value_x},
			 {2 * solution.value_y, 2 * solution.value_x, 1 + 2 * solution.value_z} };
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
	cout << endl << "Solution obtained by simple-iteration method: " << endl;
	print_solution(new_solution, "Solution");
	cout << "Residual vector obtained by simple-iteration method:  " << endl;
	print_solution(find_residual_vector(new_solution), "Residual");
}

void method_Newton(Solution& solution) {
	int count_iteration = 0;
	Solution new_solution = solution;
	do {
		solution = new_solution;
		vector<vector<double>> jacobian = find_matrix_jacobi(solution);
		vector<vector<double>> inverse_jacobian = find_inverse_matrix(jacobian);
		cout << "Matrix Jacobian has the form on " << count_iteration + 1 << " iteration:" << endl;
		print_inverse_jacobian(inverse_jacobian);
		new_solution = solution - inverse_jacobian * find_residual_vector(solution);
		++count_iteration;
	} while (second_vector_norm(new_solution - solution) > EPSILON_NEWTON);
	
	cout << "Quantity iteration in Newton method: " << count_iteration << endl;
	cout << endl << "Solution obtained by Newton method: " << endl;
	print_solution(new_solution, "Solution");
	cout << "Residual vector obtained by Newton method:  " << endl;
	print_solution(find_residual_vector(new_solution), "Residual");
}

int main() {
	Solution solution_for_simple = { 0.0, -0.17, 0.24 };
	Solution solution_for_newton = solution_for_simple;
	simple_iteration_method(solution_for_simple);
	method_Newton(solution_for_newton);
	return 0;
}