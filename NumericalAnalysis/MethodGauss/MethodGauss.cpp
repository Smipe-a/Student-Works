/*
* Laboratory work No.1
* Gaussian method with main element selection by row.
* Method Gauss. Residual vector. Determinant a matrix. Inverse matrix.
* The solution is given for a 4x4 matrix.
*/
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>

using namespace std;

const int SIZE_MATRIX = 4;

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

vector<double> operator-(const vector<double>& lhs_vector, const vector<double>& rhs_vector) {
    vector<double> subtraction_vector(SIZE_MATRIX);
    for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
        subtraction_vector[index_row] = lhs_vector[index_row] - rhs_vector[index_row];
    }
    return subtraction_vector;
}

class Equation {
public:
    //
    void SetMatrix(const vector<vector<double>>& new_matrix) { matrix = new_matrix; }
    void SetRHSVector(const vector<double>& new_rhs_vector) { rhs_vector = new_rhs_vector; }
    void SetSolutionVector(const vector<double>& new_solution) { solution_vector = new_solution; }
    void SetCountPermutation() { ++permutation; }
    
    // Return left part equations Ax=f
    vector<vector<double>> GetMatrix() const { return matrix; }
    // Return right part equations Ax=f 
    vector<double> GetRHSVector() const { return rhs_vector; }
    vector<double> GetSolutionVector() const { return solution_vector; }
    int GetCountPermutation() const { return permutation; }
    
private:
    int permutation = 0;

    vector<vector<double>> matrix;
    vector<double> rhs_vector;
    vector<double> solution_vector;
};

void print_equations(Equation& sles) {
    vector<vector<double>> lhs_matrix = sles.GetMatrix();
    vector<double> rhs_vector = sles.GetRHSVector();

    cout << "System linear algebraic equations has the form:" << endl;
    for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
        for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
            cout << setw(7) << left << lhs_matrix[index_row][index_col];
        }
        cout << "| " << rhs_vector[index_row];
        cout << endl;
    }
    cout << endl;
}

void print_vector(const vector<double>& output_vector, const string& name_vector) {
    for (int index_element = 0; index_element < SIZE_MATRIX; ++index_element) {
        if (name_vector == "Solution") {
            cout << setprecision(4);
            cout << "x_" << index_element + 1 << ": " << output_vector[index_element] << endl;
        } else {
            cout << setprecision(15);
            cout << "r_" << index_element + 1 << ": " << fixed << output_vector[index_element] << endl;
        }
    }
    cout << endl;
}

// Task No.1
void method_Gauss(Equation& sles) {
    vector<vector<double>> lhs_matrix = sles.GetMatrix();
    vector<double> rhs_vector = sles.GetRHSVector();
    vector<int> index_solution = {0, 1, 2, 3};
    vector<double> solution_equation(SIZE_MATRIX, 0.0);

    for (int index_row = 0; index_row < SIZE_MATRIX; ++index_row) {
        // Initialization main element
        double main_element = lhs_matrix[index_row][index_row];
        int index_main_element = index_row;

        for (int index_col = index_row; index_col < SIZE_MATRIX; ++index_col) {
            if (abs(main_element) < abs(lhs_matrix[index_row][index_col])) {
                main_element = lhs_matrix[index_row][index_col];
                index_main_element = index_col;
            }
        }
        if (index_main_element != index_row) {
            sles.SetCountPermutation();
            // Swap columns matrix A
            for (int index_swap = 0; index_swap < SIZE_MATRIX; ++index_swap) {
                // first - Один из индексов главной диагонали
                // second - Индекс главного элемента
                swap(lhs_matrix[index_swap][index_row], lhs_matrix[index_swap][index_main_element]);
            }
            // Также меняем массив индексов решения СЛАУ, потому что столбцы матрицы(x_first <-> x_second) менялись
            swap(index_solution[index_row], index_solution[index_main_element]);
        }

        for (int index_subrow = index_row + 1; index_subrow < SIZE_MATRIX; ++index_subrow) {
            double coeff_multiply = -(lhs_matrix[index_subrow][index_row] / lhs_matrix[index_row][index_row]);
            rhs_vector[index_subrow] += coeff_multiply * rhs_vector[index_row];
            for (int index_col = 0; index_col < SIZE_MATRIX; ++index_col) {
                lhs_matrix[index_subrow][index_col] += coeff_multiply * lhs_matrix[index_row][index_col];
            }
        }
    }

    // Solution x linear algebraic equations Ax=B
    for (int index_row = SIZE_MATRIX - 1; index_row >= 0; --index_row) {
        double coeff_multiply = 0.0;
        for (int index_col = index_row; index_col < SIZE_MATRIX; ++index_col) {
            coeff_multiply += lhs_matrix[index_row][index_col] * solution_equation[index_col];
            solution_equation[index_row] = (rhs_vector[index_row] - coeff_multiply) / lhs_matrix[index_row][index_row];
        }
    }

    // Swap elements solution, 
    // Bubble sort - O(n^2)
    for (int index_left_element = 0; index_left_element < SIZE_MATRIX; ++index_left_element) {
        for (int index_right_element = 0; index_right_element < SIZE_MATRIX - 1; ++index_right_element) {
            if (index_solution[index_right_element] > index_solution[index_right_element + 1]) {
                swap(index_solution[index_right_element + 1], index_solution[index_right_element]);
                swap(solution_equation[index_right_element + 1], solution_equation[index_right_element]);
            }
        }
    }
    sles.SetSolutionVector(solution_equation);
}

// Task No.2
vector<double> find_residual_vector(Equation& sles) {
    return sles.GetMatrix() * sles.GetSolutionVector() - sles.GetRHSVector();
}

int main() {
    Equation sles;
    sles.SetMatrix({ {-2.0,  3.01,  0.12, -0.11},
                     {2.92, -0.17,  0.11,  0.22},
                     {0.66,  0.52,  3.17,  2.11},
                     {3.01,  0.42, -0.27,  0.15} });
    sles.SetRHSVector({ 4.13,
                        3.46,
                        2.79,
                        1.01 });
    // Print first Ax=B
    print_equations(sles);
    method_Gauss(sles);
    // Print solution x
    cout << "Solution system linear algebraic equations:" << endl;
    print_vector(sles.GetSolutionVector(), "Solution");
    cout << "Residual vector:" << endl;
    print_vector(find_residual_vector(sles), "Residual");
    return 0;
}