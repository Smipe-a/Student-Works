/*
* Semester 1. Laboratory work No.2
* Method LU-decomposition and simple-iteration method.
* LU-decomposition. Simple-iteration method.
* The solution is given for a 4x4 matrix.
*/
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>

using namespace std;

const int SIZE_MATRIX = 4;
const double APPROX_EPSILON = 1e-6;

// Overloading the multiplication operator for matrices (4x4)
vector<vector<double>> operator*(const vector<vector<double>>& lhs_matrix, const vector<vector<double>>& rhs_matrix) {
    vector<vector<double>> multiply_matrix(SIZE_MATRIX, vector<double>(SIZE_MATRIX));
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
    vector<double> multiply_matrix(SIZE_MATRIX);
    for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
        double value_multiply = 0.0;
        for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
            value_multiply += lhs_matrix[index_row][index_col] * rhs_vector[index_col];
        }
        multiply_matrix[index_row] = value_multiply;
    }
    return multiply_matrix;
}

// Overloading the multiplication operator for number and vector (4x1)
vector<double> operator* (const double& number, const vector<double>& rhs_vector) {
    vector<double> scale_vector(SIZE_MATRIX);
    for (int index_element = 0; index_element < SIZE_MATRIX; ++index_element) {
        scale_vector[index_element] = number * rhs_vector[index_element];
    }
    return scale_vector;
}

// Overloading the subtraction operator for vectors of dimension (4x1)
vector<double> operator-(const vector<double>& lhs_vector, const vector<double>& rhs_vector) {
    vector<double> subtraction_vector(SIZE_MATRIX);
    for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
        subtraction_vector[index_row] = lhs_vector[index_row] - rhs_vector[index_row];
    }
    return subtraction_vector;
}

void print_matrix(const vector<vector<double>>& output_matrix) {
    cout << setprecision(2);
    for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
        for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
            cout << setw(8) << left << fixed << output_matrix[index_row][index_col];
        }
        cout << endl;
    }
    cout << endl;
}

void print_vector(const vector<double>& output_vector, const string& name_vector) {
    for (int index_element = 0; index_element < SIZE_MATRIX; ++index_element) {
        if (name_vector == "Solution") {
            cout << setprecision(7);
            cout << "x_" << index_element + 1 << ": " << output_vector[index_element] << endl;
        }
        else if (name_vector == "Residual") {
            cout << setprecision(15);
            cout << "r_" << index_element + 1 << ": " << fixed << output_vector[index_element] << endl;
        }
    }
    cout << endl;
}

class Equation {
public:
    // Setters
    void setMatrix(const vector<vector<double>>& new_matrix) { matrix = new_matrix; }
    void setRHSVector(const vector<double>& new_rhs_vector) { rhs_vector = new_rhs_vector; }
    void setSolutionVector(const vector<double>& new_solution) { solution_vector = new_solution; }

    // Getters
    vector<vector<double>> getMatrix() const { return matrix; }
    vector<double> getRHSVector() const { return rhs_vector; }
    vector<double> getSolutionVector() const { return solution_vector; }

private:
    vector<vector<double>> matrix;
    vector<double> rhs_vector;
    vector<double> solution_vector;
};

// Task No.1
void lu_decomposition(Equation& sles, vector<vector<double>>& matrix_L, vector<vector<double>>& matrix_U) {
    vector<vector<double>> matrix_A = sles.getMatrix();
    vector<double> rhs_vector = sles.getRHSVector();
    matrix_U = matrix_A;
    for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
        for (int index_row = index_col; index_row < SIZE_MATRIX; ++index_row) {
            matrix_L[index_row][index_col] = matrix_U[index_row][index_col] / matrix_U[index_col][index_col];
        }
    }
    for (int index_subrow = 1; index_subrow < SIZE_MATRIX; ++index_subrow) {
        for (int index_col = index_subrow - 1; index_col < SIZE_MATRIX; ++index_col) {
            for (int index_row = index_col; index_row < SIZE_MATRIX; ++index_row) {
                matrix_L[index_row][index_col] = matrix_U[index_row][index_col] / matrix_U[index_col][index_col];
            }
        }
        for (int index_row = index_subrow; index_row < SIZE_MATRIX; ++index_row) {
            for (int index_col = index_subrow - 1; index_col < SIZE_MATRIX; ++index_col) {
                matrix_U[index_row][index_col] -= matrix_L[index_row][index_subrow - 1] * matrix_U[index_subrow - 1][index_col];
            }
        }
    }
    // Finding a solution 'x'
    // Top Down Ly=f
    for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
        double value_sum = 0.0;
        for (int index_col = 0; index_col < index_row; ++index_col)
            value_sum += matrix_L[index_row][index_col] * rhs_vector[index_col];
        rhs_vector[index_row] = (rhs_vector[index_row] - value_sum) / matrix_L[index_row][index_row];
    }
    // Bottom Up Ux=y
    for (int index_row = SIZE_MATRIX - 1; index_row >= 0; --index_row) {
        double sum = 0.0;
        for (int index_col = index_row + 1; index_col < SIZE_MATRIX; ++index_col)
            sum += matrix_U[index_row][index_col] * rhs_vector[index_col];
        rhs_vector[index_row] = (rhs_vector[index_row] - sum) / matrix_U[index_row][index_row];
    }
    sles.setSolutionVector(rhs_vector);
}

// Task No.2
vector<double> find_residual_vector(Equation& sles) {
    return sles.getMatrix() * sles.getSolutionVector() - sles.getRHSVector();
}

// Task No.3
double first_matrix_norm(const vector<vector<double>>& matrix) {
    vector<double> quantity_sum_coloumns(SIZE_MATRIX);
    for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
        double sum_elements = 0.0;
        for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
            sum_elements += abs(matrix[index_row][index_col]);
        }
        quantity_sum_coloumns[index_col] = sum_elements;
    }
    return *max_element(quantity_sum_coloumns.begin(), quantity_sum_coloumns.end());
}

double second_vector_norm(const vector<double>& vector_by_norm) {
    double sum_elements = 0.0;
    for (int index_element = 0; index_element < SIZE_MATRIX; ++index_element) {
        sum_elements += pow(vector_by_norm[index_element], 2);
    }
    return sqrt(sum_elements);
}

void simple_iteration_method(Equation& sles) {
    const double TAU = 2.0 / first_matrix_norm(sles.getMatrix());
    vector<double> prev_approximation(SIZE_MATRIX);
    vector<double> rhs_vector = sles.getRHSVector();
    int count_iteration = 0;
    vector<double> next_approximation(SIZE_MATRIX);
    for (int index_value = 0; index_value < SIZE_MATRIX; ++index_value) {
        next_approximation[index_value] = static_cast<int>(rhs_vector[index_value]);
    }
    sles.setSolutionVector(next_approximation);
    do {
        prev_approximation = next_approximation;
        next_approximation = prev_approximation - TAU * find_residual_vector(sles);
        sles.setSolutionVector(next_approximation);
        ++count_iteration;
    } while (second_vector_norm(next_approximation - prev_approximation) > APPROX_EPSILON);
    cout << "Quantity iteration: " << count_iteration << endl;
    cout << endl << "Solution obtained by Simple-iteration method: " << endl;
    print_vector(sles.getSolutionVector(), "Solution");
    cout << "Residual vector obtained by Simple-iteration method: " << endl;
    print_vector(find_residual_vector(sles), "Residual");
}

// Task No.4
vector<vector<double>> find_inverse_matrix(vector<vector<double>> lhs_matrix) {
    //   Find the inverse matrix
    //   [ 15.0, 1.0, -5.0,   3.0  |  1.0, 0.0, 0.0, 0.0]
    //   [ 1.0,  10.0, 2.0,  -4.0  |  0.0, 1.0, 0.0, 0.0]
    //   [-5.0,  2.0,  14.0, -6.0  |  0.0, 0.0, 1.0, 0.0]
    //   [ 3.0, -4.0, -6.0,   16.0 |  0.0, 0.0, 0.0, 1.0]
    vector<vector<double>> unit_matrix = { {1.0, 0.0, 0.0, 0.0},
                                           {0.0, 1.0, 0.0, 0.0},
                                           {0.0, 0.0, 1.0, 0.0},
                                           {0.0, 0.0, 0.0, 1.0} };
    for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
        double permitting_element = lhs_matrix[index_row][index_row];
        for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
            lhs_matrix[index_row][index_col] /= permitting_element;
            unit_matrix[index_row][index_col] /= permitting_element;
        }
        for (int index_subrow = index_row + 1; index_subrow < SIZE_MATRIX; ++index_subrow) {
            permitting_element = lhs_matrix[index_subrow][index_row];
            for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
                lhs_matrix[index_subrow][index_col] -= lhs_matrix[index_row][index_col] * permitting_element;
                unit_matrix[index_subrow][index_col] -= unit_matrix[index_row][index_col] * permitting_element;
            }
        }
    }
    for (int index_subrow = SIZE_MATRIX - 1; index_subrow > 0; --index_subrow) {
        for (int index_row = index_subrow - 1; index_row >= 0; --index_row) {
            double permitting_element = lhs_matrix[index_row][index_subrow];
            for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
                lhs_matrix[index_row][index_col] -= lhs_matrix[index_subrow][index_col] * permitting_element;
                unit_matrix[index_row][index_col] -= unit_matrix[index_subrow][index_col] * permitting_element;
            }
        }
    }
    return unit_matrix;
}

int main() {
    Equation sles;
    sles.setMatrix({ {15.0, 1.0, -5.0,   3.0},
                     {1.0,  10.0, 2.0,  -4.0},
                     {-5.0, 2.0,  14.0, -6.0},
                     {3.0, -4.0, -6.0,   16.0} });
    sles.setRHSVector({ -24.0,
                        -47.0,
                         28.0,
                        -50.0 });
    // Task No. 1
    vector<vector<double>> matrix_L(SIZE_MATRIX, { 0.0, 0.0, 0.0, 0.0 }), matrix_U(SIZE_MATRIX, { 0.0, 0.0, 0.0, 0.0 });
    lu_decomposition(sles, matrix_L, matrix_U);
    cout << "Matrix L has the form: " << endl;
    print_matrix(matrix_L);
    cout << "Matrix U has the form: " << endl;
    print_matrix(matrix_U);
    cout << "Solution obtained by LU-decomposition: " << endl;
    print_vector(sles.getSolutionVector(), "Solution");

    // Task No.2
    cout << "Residual vector:" << endl;
    print_vector(find_residual_vector(sles), "Residual");

    // Task No.3
    simple_iteration_method(sles);

    // Task No. 4
    cout << "Number condition matrix M_A: " << first_matrix_norm(sles.getMatrix()) * first_matrix_norm(find_inverse_matrix(sles.getMatrix()));
    return 0;
}