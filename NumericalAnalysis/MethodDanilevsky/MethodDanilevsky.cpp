/*
* Semester 1. Laboratory work No.4
* Approximation the function. Calculation of eigenvalues ​​and eigenvectors of matrices.
* Chebyshev polynomials. Danilevsky's method. Bisection method.
*/
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

// Approximate coordinates (x_i, y_i)
const vector<double> SAMPLE_X = { 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7 };
const vector<double> SAMPLE_Y = { 0.5913, 0.9241, 0.7162, 0.8731, 0.9574, 0.9015, 1.3561, 1.2738, 1.2724, 1.1672 };
// Step for sample X
const double STEP = 0.1;
// Quantity approximate values X and Y beginning with x_1 ... x_9. Denote 0.8 as x_0
const int QUANTITY_VALUE = 9;
const int POLYNOMIAL_DEGREE = 3;
const int SIZE_MATRIX = 4;

// Overloading the multiplication operator for matrices (4x4)
vector<vector<double>> operator*(const vector<vector<double>>& lhs_matrix, const vector<vector<double>>& rhs_matrix) {
	vector<vector<double>> multiply_matrix(SIZE_MATRIX, vector<double>(SIZE_MATRIX, 0.0));
	for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
		for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
			for (int index_subrow = 0; index_subrow < SIZE_MATRIX; ++index_subrow) {
				multiply_matrix[index_row][index_col] += lhs_matrix[index_row][index_subrow] * rhs_matrix[index_subrow][index_col];
			}
		}
	}
	return multiply_matrix;
}

// Overloading the multiplication operator for matrix (4x4) and vector (4x1)
vector<double> operator*(const vector<vector<double>>& lhs_matrix, const vector<double>& rhs_vector) {
	vector<double> multiply_matrix(SIZE_MATRIX, 0.0);
	for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
		double value_multiply = 0.0;
		for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
			value_multiply += lhs_matrix[index_row][index_col] * rhs_vector[index_col];
		}
		multiply_matrix[index_row] = value_multiply;
	}
	return multiply_matrix;
}

void print_matrix(const vector<vector<double>>& matrix) {
	for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
		for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
			cout << setw(10) << left << fixed << matrix[index_row][index_col];
		}
		cout << endl;
	}
	cout << endl;
}

void print_vector(const vector<double>& eigenvector) {
	for (const auto& value : eigenvector) {
		cout << setw(10) << fixed << left << value;
	}
	cout << endl;
}

void print_polynom(const vector<double>& coefficients) {
	cout << setprecision(4) << "Chebyshev polynomial has the form:" << endl;
	cout << "Q(t) = ";
	for (int index_polinom = 0; index_polinom <= POLYNOMIAL_DEGREE; ++index_polinom) {
		if (index_polinom == POLYNOMIAL_DEGREE) {
			if (coefficients[index_polinom] > 0) {
				cout << "+" << coefficients[index_polinom] << "*P_" << index_polinom << "(t)" << endl;
			} else {
				cout << coefficients[index_polinom] << "*P_" << index_polinom << "(t)" << endl;
			}
		} else {
			if (coefficients[index_polinom] < 0) {
				cout << coefficients[index_polinom] << "*P_" << index_polinom << "(t)";
			} else {
				if (index_polinom == 0) {
					cout << coefficients[index_polinom] << "*P_" << index_polinom << "(t)";
				} else {
					cout << "+" << coefficients[index_polinom] << "*P_" << index_polinom << "(t)";
				}
			}
		}
	}
}

// Finding values in Chebyshev polynomials
double polinom_chebyshev(const int& number_polinom, const double& value_x) {
	double value_t = (value_x - SAMPLE_X[0]) / STEP;
	switch (number_polinom) {
	case 0:
		return 1;
	case 1:
		return 1 - 2 * (value_t / QUANTITY_VALUE);
	case 2:
		return 1 - 6 * (value_t / QUANTITY_VALUE) + 6 * (value_t * (value_t - 1)) / (QUANTITY_VALUE * (QUANTITY_VALUE - 1));
	case 3:
		return 1 - 12 * (value_t / QUANTITY_VALUE) + 30 * (value_t * (value_t - 1)) / (QUANTITY_VALUE * (QUANTITY_VALUE - 1))
			- 20 * (value_t * (value_t - 1) * (value_t - 2)) / (QUANTITY_VALUE * (QUANTITY_VALUE - 1) * (QUANTITY_VALUE - 2));
	}
}

// Finding value of the characteristic equation of the matrix A
// (-1)^{k+1} * (Lambda^4 - p_1 * Lambda^3 - p_2 * Lamda^2 - p_3 * Lambda - p_4)
double function_polinom(const double& lambda, const vector<double>& coeff) {
	return coeff[0] * pow(lambda, 4) - coeff[1] * pow(lambda, 3) - coeff[2] * pow(lambda, 2) - coeff[3] * lambda - coeff[4];
}

// Power method
void power_method(const vector<vector<double>> matrix, vector<double> vector_eigenvalue) {
	// Initial approximation
	vector<double> vector_preview(SIZE_MATRIX, 1.0);
	// Accuracy finding spectral radius
	vector<double> vector_next = matrix * vector_preview;
	const double EPSILON = 1e-3;
	double lambda_first = vector_next[0] / vector_preview[0];
	double lambda_second = lambda_first;
	vector_preview = vector_next;
	do {
		lambda_first = lambda_second;
		vector_preview = vector_next;
		vector_next = matrix * vector_preview;
		lambda_second = vector_next[0] / vector_preview[0];
	} while (abs(lambda_second - lambda_first) > EPSILON);
	cout << "Seeking value of the spectral radius: " << fixed << abs(lambda_second);
}

// The bisection method as root-finding method
vector<double> bisection_method(const vector<double>& ceoff_polinom) {
	// Vector roots equation
	vector<double> roots;
	// Will look for values ​​on this interval D:[0, 10]
	double left_edge = 0.0;
	double right_edge = 10.0;
	// Accuracy finding roots
	const double ACCURACY_SOLUTION = 1e-3;
	vector<double> edges_segments(right_edge);
	for (int index_section = left_edge; index_section < right_edge; ++index_section) {
		edges_segments[index_section] = index_section;
	}
	for (int quantity_section = 0; quantity_section < edges_segments.size() - 1; ++quantity_section) {
		bool flag_root = true;
		left_edge = edges_segments[quantity_section];
		right_edge = edges_segments[quantity_section + 1];
		while ((right_edge - left_edge) > ACCURACY_SOLUTION) {
			double middle_edge = (left_edge + right_edge) / 2;
			if (function_polinom(left_edge, ceoff_polinom) * function_polinom(middle_edge, ceoff_polinom) <= 0) {
				right_edge = middle_edge;
			} else if (function_polinom(right_edge, ceoff_polinom) * function_polinom(middle_edge, ceoff_polinom) > 0) {
				flag_root = false;
				break;
			} else {
				left_edge = middle_edge;
			}
		}
		if (flag_root) {
			roots.push_back((left_edge + right_edge) / 2);
		}
	}
	return roots;
}

// Danilevsky's method P = B^-1 * A * B
vector<double> method_Danilevsky(const vector<vector<double>>& matrixA) {
	vector<vector<double>> matrixFrobenius = matrixA;
	vector<vector<double>> matrix_similarityB = { {1.0, 0.0, 0.0, 0.0},
												  {0.0, 1.0, 0.0, 0.0},
												  {0.0, 0.0, 1.0, 0.0},
												  {0.0, 0.0, 0.0, 1.0} };
	for (int index_row = 0; index_row < SIZE_MATRIX - 1; ++index_row) {
		vector<vector<double>> matrix_left_multiply(SIZE_MATRIX, vector<double>(SIZE_MATRIX, 0.0));
		vector<vector<double>> matrix_right_multiply(SIZE_MATRIX, vector<double>(SIZE_MATRIX, 0.0));
		for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
			matrix_left_multiply[index_col][index_col] = 1;
			matrix_right_multiply[index_col][index_col] = 1;
			if (index_col == SIZE_MATRIX - index_row - 2) {
				matrix_right_multiply[SIZE_MATRIX - 2 - index_row][index_col] = 1 / matrixFrobenius[SIZE_MATRIX - 1 - index_row][SIZE_MATRIX - 2 - index_row];
			} else {
				matrix_right_multiply[SIZE_MATRIX - 2 - index_row][index_col] = -matrixFrobenius[SIZE_MATRIX - 1 - index_row][index_col] / 
					matrixFrobenius[SIZE_MATRIX - 1 - index_row][SIZE_MATRIX - 2 - index_row];
			}
			matrix_left_multiply[SIZE_MATRIX - 2 - index_row][index_col] = matrixFrobenius[SIZE_MATRIX - 1 - index_row][index_col];
		}
		matrix_similarityB = matrix_similarityB * matrix_right_multiply;
		matrixFrobenius = (matrix_left_multiply * matrixFrobenius) * matrix_right_multiply;
	}
	cout << "Frobenius matrix has the form:" << endl;
	print_matrix(matrixFrobenius);
	vector<double> coeff_polynom = { 1, matrixFrobenius[0][0], matrixFrobenius[0][1], matrixFrobenius[0][2], matrixFrobenius[0][3] };
	vector<double> roots = bisection_method(coeff_polynom);
	cout << "Eigenvalues matrix A has the form:" << endl;
	for (int index_value_roots = 0; index_value_roots < roots.size(); ++index_value_roots) {
		cout << "Lambda" << index_value_roots + 1 << ": " << fixed << roots[index_value_roots] << endl;
	}
	cout << endl;
	for (int quantity_vectors = 0; quantity_vectors < roots.size(); ++quantity_vectors) {
		vector<double> y_vector;
		for (int index_element = 0; index_element < roots.size(); ++index_element) {
			y_vector.push_back(pow(roots[quantity_vectors], 3 - index_element));
		}
		vector<double> eigenvector = matrix_similarityB * y_vector;
		cout << quantity_vectors + 1 << " eigenvector matrix A has the form:" << endl;
		print_vector(eigenvector);
	}
	return roots;
}

int main() {
	// Task No.1
	// Vector coefficients Chebyshev polynomials
	vector<double> coeff(POLYNOMIAL_DEGREE + 1, 0.0);
	for (int index_polinom = 0; index_polinom <= POLYNOMIAL_DEGREE; ++index_polinom) {
		// Find numerator and denominator coefficients orthogonal Chebyshev polynomials
		double numerator = 0.0;
		double denominator = 0.0;
		for (int index_element = 0; index_element <= QUANTITY_VALUE; ++index_element) {
			numerator += SAMPLE_Y[index_element] * polinom_chebyshev(index_polinom, SAMPLE_X[index_element]);
			denominator += pow(polinom_chebyshev(index_polinom, SAMPLE_X[index_element]), 2);
		}
		coeff[index_polinom] = numerator / denominator;
	}
	print_polynom(coeff);
	for (int index_value = 0; index_value <= QUANTITY_VALUE; ++index_value) {
		double new_value_sample = SAMPLE_X[index_value] + STEP / 2;
		cout << "y(" << new_value_sample << ") = " << coeff[0] * polinom_chebyshev(0, new_value_sample) + coeff[1] * polinom_chebyshev(1, new_value_sample) +
			coeff[2] * polinom_chebyshev(2, new_value_sample) + coeff[3] * polinom_chebyshev(3, new_value_sample) << endl;
	}
	cout << endl;

	// Task No.2
	vector<vector<double>> matrixA = { {1.6, 1.6, 1.7, 1.8},
									   {1.6, 2.6, 1.3, 1.3},
									   {1.7, 1.3, 3.6, 1.4},
									   {1.8, 1.3, 1.4, 4.6} };
	vector<double> eigenvalues = method_Danilevsky(matrixA);
	cout << endl;

	// Task No.3
	power_method(matrixA, eigenvalues);
	return 0;
}