/*
* Laboratory work No.2
*
* LU-decomposition. Simple-iteration method.
* The solution is given for a 4x4 matrix.
*/
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

const int SIZE_MATRIX = 4;

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
            cout << setprecision(4);
            cout << "x_" << index_element + 1 << ": " << output_vector[index_element] << endl;
        }
        else if (name_vector == "Residual") {
            cout << setprecision(15);
            cout << "r_" << index_element + 1 << ": " << fixed << output_vector[index_element] << endl;
        }
    }
    cout << endl;
}

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

    // Task No. 4
    vector<vector<double>> inverse = find_inverse_matrix(sles.getMatrix());
    return 0;
}