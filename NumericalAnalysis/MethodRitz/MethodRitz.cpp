/*
* Semester 2. Laboratory work No.1
* Solve Boundary value problems by Ritz Method.
* Ritz method. Simpson's rule. Method Gauss. Sylvester's criterion.
*/
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

// Approximation accuracy for finding integral
const double APPROXIMATION_ACCURACY = 1e-6;
// Left edge and right edge the segment [a; b]
const double LEFT_EDGE = 0.0;
const double RIGHT_EDGE = 1.3;
// Border conditions
const double MU_FIRST = 1.5;
const double MU_TWO = 0.0;
const int QUANTITY_FUNCTIONS = 3;

void print_matrix(const vector<vector<double>>& output_matrix, const vector<double>& output_vector) {
	for (int index_row = 0; index_row < QUANTITY_FUNCTIONS; ++index_row) {
		for (int index_col = 0; index_col < QUANTITY_FUNCTIONS; ++index_col) {
			if (index_col == QUANTITY_FUNCTIONS - 1) {
				cout << fixed << setw(12) << output_matrix[index_row][index_col] << " | " << output_vector[index_row];
			} else {
				cout << setw(12) << output_matrix[index_row][index_col] << " ";
			}
		}
		cout << endl;
	}
	cout << endl;
}

double function_k(const double& value_x) { return sqrt(value_x / 5) + 5.0 / 3; }

double function_q(const double& value_x) { return (value_x + 4) / (value_x * value_x + 5.0 / 3); }

double function_f(const double& value_x) { return value_x / 3; }

// Function phi_0
double function_phi_null(const double& value_x) { return 1.5 - 1.5 / 1.3 * value_x; }

// First derivative function phi_0
double function_derivative_phi_null() { return -1.5 / 1.3; }

// Function phi_k, where k = 1,3
double function_phi(const double& value_x, const int& exponent) { return (RIGHT_EDGE - value_x) * pow(value_x, exponent); }

// First derivative function phi_k, where k = 1,3
double function_derivative_phi(const double& value_x, const int& exponent) { return pow(value_x, exponent - 1) * (exponent * (RIGHT_EDGE - value_x) - value_x); }

// Method Gauss for solution SLEs
vector<double> method_Gauss(vector<vector<double>> left_matrix, vector<double> right_vector) {
	vector<double> help_vector(QUANTITY_FUNCTIONS);
	for (int index_row = 0; index_row < QUANTITY_FUNCTIONS; ++index_row) {
		// Matrix transformation. Triangularization
		for (int index_subrow = index_row + 1; index_subrow <= QUANTITY_FUNCTIONS - 1; ++index_subrow) {
			help_vector[index_subrow] = -(left_matrix[index_subrow][index_row] / left_matrix[index_row][index_row]);
			right_vector[index_subrow] = right_vector[index_subrow] + help_vector[index_subrow] * right_vector[index_row];
			for (int index_col = 0; index_col <= QUANTITY_FUNCTIONS - 1; ++index_col) {
				left_matrix[index_subrow][index_col] += help_vector[index_subrow] * left_matrix[index_row][index_col];
			}
		}
	}
	// Solution SLEs c_j, where j=1,3
	for (int index_row = QUANTITY_FUNCTIONS - 1; index_row >= 0; --index_row) {
		double sum_elements = 0.0;
		for (int index_col = index_row + 1; index_col < QUANTITY_FUNCTIONS; ++index_col)
			sum_elements += left_matrix[index_row][index_col] * right_vector[index_col];
		right_vector[index_row] = (right_vector[index_row] - sum_elements) / left_matrix[index_row][index_row];
	}
	return right_vector;
}

// Find integral by Simpson's rule
double method_Simpson(const string& check, const int& number_row, const int& number_col, const int& nodes) {
	// Segment splitting step
	double step = (RIGHT_EDGE - LEFT_EDGE) / (2.0 * nodes);
	double value_integral = 0;

	if (check == "FOR_VECTOR") {
		for (int index_point = 0; index_point <= (QUANTITY_FUNCTIONS - 1) * nodes; ++index_point) {
			double point = LEFT_EDGE + index_point * step;
			// b_k = integral {f * phi_k - k * phi'_0 * phi'_k - q * phi_0 * phi_k } 
			double value_function = -function_k(point) * function_derivative_phi_null() * function_derivative_phi(point, number_row) -
				function_q(point) * function_phi_null(point) * function_phi(point, number_row) +
				function_f(point) * function_phi(point, number_row);
			if ((index_point == 0) || (index_point == ((QUANTITY_FUNCTIONS - 1) * nodes))) {
				value_integral += value_function;
			} else if (index_point % 2 == 0) {
				value_integral += 2 * value_function;
			} else if (index_point % 2 != 0) {
				value_integral += 4 * value_function;
			}
		}
	} else if (check == "FOR_MATRIX") {
		for (int index_point = 0; index_point <= (QUANTITY_FUNCTIONS - 1) * nodes; ++index_point) {
			double point = LEFT_EDGE + index_point * step;
			// a_ik = integral {k * phi'_k * phi'_i + q * phi_k * phi_i}
			double value_function = function_k(point) * function_derivative_phi(point, number_row) * function_derivative_phi(point, number_col) +
				function_q(point) * function_phi(point, number_row) * function_phi(point, number_col);
			if ((index_point == 0) || (index_point == ((QUANTITY_FUNCTIONS - 1) * nodes))) {
				value_integral += value_function;
			} else if (index_point % 2 == 0) {
				value_integral += 2 * value_function;
			} else if (index_point % 2 != 0) {
				value_integral += 4 * value_function;
			}
		}
	}
	return step / 3 * value_integral;
}

// Find coefficients a_ij matrix A, where i,j = 1,3
vector<vector<double>> find_coefficient_matrix() {
	vector<vector<double>> matrix_a(QUANTITY_FUNCTIONS, vector<double>(QUANTITY_FUNCTIONS, 0.0));
	// Quantity nodes splitting segment
	int nodes = 3;
	// First integral with step h
	double integral_prev = 0;
	// Second integral with step h / 2
	double integral_next = 0;

	for (int index_row = 1; index_row <= QUANTITY_FUNCTIONS; ++index_row) {
		for (int index_col = 1; index_col <= QUANTITY_FUNCTIONS; ++index_col) {
			do {
				integral_prev = method_Simpson("FOR_MATRIX", index_row, index_col, nodes);
				integral_next = method_Simpson("FOR_MATRIX", index_row, index_col, nodes * 2);
				nodes *= 2;
			} while (abs(integral_prev - integral_next) > APPROXIMATION_ACCURACY);
			matrix_a[index_row - 1][index_col - 1] = integral_prev;
		}
	}
	return matrix_a;
}

// Find coefficients right side SLEs Ac=b
vector<double> find_coefficient_vector() {
	vector<double> vector_b(QUANTITY_FUNCTIONS, 0.0);
	// Quantity nodes splitting segment
	int nodes = 3;
	// First integral with step h
	double integral_prev = 0;
	// Second integral with step h / 2
	double integral_next = 0;

	for (int index_coefficient = 1; index_coefficient <= QUANTITY_FUNCTIONS; ++index_coefficient) {
		do {
			integral_prev = method_Simpson("FOR_VECTOR", index_coefficient, 0, nodes);
			integral_next = method_Simpson("FOR_VECTOR", index_coefficient, 0, nodes * 2);
			nodes *= 2;
		} while (abs(integral_prev - integral_next) > APPROXIMATION_ACCURACY);
		vector_b[index_coefficient - 1] = integral_prev;
	}
	return vector_b;
}

// First matrix norm defined as max_{j=1,n}sum_{i=1,m}|a_ij|
double first_matrix_norm(const vector<vector<double>>& left_matrix) {
	vector<double> maximum_sum(QUANTITY_FUNCTIONS);
	for (int index_row = 0; index_row < QUANTITY_FUNCTIONS; ++index_row) {
		double sum_elements = 0.0;
		for (int index_col = 0; index_col < QUANTITY_FUNCTIONS; ++index_col) {
			sum_elements += abs(left_matrix[index_row][index_col]);
		}
		maximum_sum[index_row] = sum_elements;
	}
	return *std::max_element(maximum_sum.begin(), maximum_sum.end());
}

// Find inverse matrix for condition number
vector<vector<double>> find_inverse_matrix(vector<vector<double>> left_matrix) {
	vector<vector<double>> unit_matrix = { {1.0, 0.0, 0.0},
										   {0.0, 1.0, 0.0},
										   {0.0, 0.0, 1.0} };
	for (int index_row = 0; index_row < QUANTITY_FUNCTIONS; ++index_row) {
		double permitting_element = left_matrix[index_row][index_row];
		for (int index_col = 0; index_col < QUANTITY_FUNCTIONS; ++index_col) {
			left_matrix[index_row][index_col] /= permitting_element;
			unit_matrix[index_row][index_col] /= permitting_element;
		}
		for (int index_subrow = index_row + 1; index_subrow < QUANTITY_FUNCTIONS; ++index_subrow) {
			permitting_element = left_matrix[index_subrow][index_row];
			for (int index_col = 0; index_col < QUANTITY_FUNCTIONS; ++index_col) {
				left_matrix[index_subrow][index_col] -= left_matrix[index_row][index_col] * permitting_element;
				unit_matrix[index_subrow][index_col] -= unit_matrix[index_row][index_col] * permitting_element;
			}
		}
	}
	for (int index_subrow = QUANTITY_FUNCTIONS - 1; index_subrow > 0; --index_subrow) {
		for (int index_row = index_subrow - 1; index_row >= 0; --index_row) {
			double permitting_element = left_matrix[index_row][index_subrow];
			for (int index_col = 0; index_col < QUANTITY_FUNCTIONS; ++index_col) {
				left_matrix[index_row][index_col] -= left_matrix[index_subrow][index_col] * permitting_element;
				unit_matrix[index_row][index_col] -= unit_matrix[index_subrow][index_col] * permitting_element;
			}
		}
	}
	return unit_matrix;
}

// Sylvester's criterion
void check_positive_definite(const vector<vector<double>>& left_matrix) {
	double corner_minor_frist = left_matrix[0][0];
	double corner_minor_two = left_matrix[0][0] * left_matrix[1][1] - left_matrix[0][1] * left_matrix[1][0];
	double corner_minor_third = left_matrix[0][0] * left_matrix[1][1] * left_matrix[2][2] +
		left_matrix[1][0] * left_matrix[2][1] * left_matrix[0][2] +
		left_matrix[0][1] * left_matrix[1][2] * left_matrix[2][0] -
		left_matrix[2][0] * left_matrix[1][1] * left_matrix[0][2] -
		left_matrix[2][1] * left_matrix[1][2] * left_matrix[0][0] -
		left_matrix[1][0] * left_matrix[0][1] * left_matrix[2][2];
	if ((corner_minor_frist > 0) && (corner_minor_two > 0) && (corner_minor_third > 0)) {
		cout << "Matrix A by Sylvester's criterion positively defined" << endl;
	} else if ((corner_minor_frist < 0) && (corner_minor_two > 0) && (corner_minor_third < 0)) {
		cout << "Matrix A by Sylvester's criterion negatively defined" << endl;
	} else {
		cout << "Matrix A by Sylvester's criterion alternating" << endl;
	}
	cout << "Angular minors equals: " << endl;
	cout << corner_minor_frist << "   " << corner_minor_two << "   " << corner_minor_third << endl;
	cout << endl;
}

int main() {
	cout << setprecision(7);

	// Task No.1
	vector<vector<double>> matrix = find_coefficient_matrix();
	vector<double> right_vector = find_coefficient_vector();
	cout << "System linear algebraic equations Ac=b has the form: " << endl;
	print_matrix(matrix, right_vector);

	// Task No.2
	vector<double> solution_linear_equations = method_Gauss(matrix, right_vector);
	cout << "Solution SLEs:" << endl;
	for (int index_element = 0; index_element < QUANTITY_FUNCTIONS; ++index_element) {
		cout << "c" << index_element + 1 << ": " << solution_linear_equations[index_element] << endl;
	}

	// Task No.3
	cout << endl << "The condition number of the matrix A with considering first matrix norm: " <<
		first_matrix_norm(matrix) * first_matrix_norm(find_inverse_matrix(matrix)) << endl << endl;
	check_positive_definite(matrix);

	// Task No.4
	cout << "Values a function in some segment points [a ,b]:" << endl;
	double step = (RIGHT_EDGE - LEFT_EDGE) / 5;
	for (int index_point = 0; index_point <= 5; ++index_point) {
		double value_x = index_point * step;
		double value_function = function_phi_null(value_x) + solution_linear_equations[0] * function_phi(value_x, 1) +
			solution_linear_equations[1] * function_phi(value_x, 2) + solution_linear_equations[2] * function_phi(value_x, 3);
		cout << "y(" << value_x << ") = " << value_function << endl;
	}
	return 0;
}