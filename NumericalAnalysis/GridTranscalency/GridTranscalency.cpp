/*
* Semester 2. Laboratory work No.3
* The grid method for solving the initial-boundary value problem for the heat conduction equation.
* Explicit and Implicit Finite Difference Schemes. Tridiagonal matrix algorithm.
*/
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

const int LIMIT_T = 6;
const int LIMIT_X = 11;
const double STEP_X = 0.1;
const double STEP_T = 0.01;
const double VALUE_S = STEP_T / pow(STEP_X, 2);
// Coefficients for TMA
const double VALUE_A = VALUE_S;
const double VALUE_B = VALUE_S;
const double VALUE_C = 1.0 + 2.0 * VALUE_S;

double edge_func(const double& value_t) { return 2.5 * value_t; }

double function_f(const double& value_x, const double& value_t) { return 2.5 * (pow(value_x, 2) - 2 * value_t); }

void nullify_grid(vector<vector<double>>& grid) {
    for (int index_row = 0; index_row < LIMIT_X; ++index_row) {
        for (int index_col = 0; index_col < LIMIT_T; ++index_col) {
            grid[index_row][index_col] = 0.0;
        }
    }
}

vector<vector<double>> transpose_matrix(const vector<vector<double>>& grid) {
    vector<vector<double>> transpose_function_grid(6, vector<double>(11));
    for (int index_row = 0; index_row < LIMIT_X; ++index_row) {
        for (int index_col = 0; index_col < LIMIT_T; ++index_col) {
            transpose_function_grid[index_col][index_row] = grid[index_row][index_col];
        }
    }
    return transpose_function_grid;
}

void print_solution(const vector<vector<double>>& out_solution, const vector<double>& out_nodes) {
    vector<vector<double>> transpose_function_grid(LIMIT_T, vector<double>(LIMIT_X));
    transpose_function_grid = transpose_matrix(out_solution);
    for (int index_t = 0; index_t < LIMIT_T; ++index_t) {
        cout << setprecision(2);
        cout << setw(9) << left << out_nodes[LIMIT_T - index_t - 1];
        cout << setprecision(4);
        for (int index_x = 0; index_x < LIMIT_X; ++index_x) {
            cout << fixed << left << setw(9) << transpose_function_grid[LIMIT_T - index_t - 1][index_x];
        }
        cout << endl;
    }
    cout << setprecision(2) << setw(9) << "t/x";
    for (int index_step = 0; index_step < LIMIT_X; ++index_step) {
        cout << setw(9) << left << STEP_X * index_step;
    }
    cout << endl;
}

int main() {
    vector<double> nodes_x(LIMIT_X);
    for (int index_x = 0; index_x < LIMIT_X; ++index_x) {
        nodes_x[index_x] = index_x * STEP_X;
    }
    vector<double> nodes_t(LIMIT_T);
    for (int index_t = 0; index_t < LIMIT_T; ++index_t) {
        nodes_t[index_t] = index_t * STEP_T;
    }
    vector<vector<double>> function_grid(LIMIT_X, vector<double>(LIMIT_T, 0.0));
    // Explicit Finite Difference Scheme
    for (int index_t = 0; index_t < LIMIT_T; ++index_t) {
        function_grid[LIMIT_X - 1][index_t] = edge_func(nodes_t[index_t]);
    }
    for (int index_col = 0; index_col < LIMIT_T - 1; ++index_col) {
        for (int index_row = 1; index_row < LIMIT_X - 1; ++index_row) {
            function_grid[index_row][index_col + 1] = VALUE_S * function_grid[index_row - 1][index_col] + (1 - 2 * VALUE_S) * function_grid[index_row][index_col] +
                VALUE_S * function_grid[index_row + 1][index_col] + STEP_T * function_f(nodes_x[index_row], nodes_t[index_col]);
        }
    }
    cout << "Solution obtained by explicit finite difference schemes:" << endl;
    print_solution(function_grid, nodes_t);

    // Implicit Finite Difference Scheme
    nullify_grid(function_grid);
    for (int index_t = 1; index_t < LIMIT_T; ++index_t) {
        // Find other TMA coefficients alpha and betta
        // Source alpha_0 and betta_0, which are boundary values, where alpha_0 = betta_0 = 0
        vector<double> alpha(LIMIT_X, 0.0);
        vector<double> betta(LIMIT_X, 0.0);
        betta[LIMIT_X - 1] = edge_func(nodes_t[index_t]);
        for (int index_x = 1; index_x < LIMIT_X; ++index_x) {
            // Find other alpha_n and betta_n, n > 1
            alpha[index_x] = VALUE_B / (VALUE_C - alpha[index_x - 1] * VALUE_B);
            betta[index_x] = (STEP_T * function_f(nodes_x[index_x], nodes_t[index_t]) + function_grid[index_x][index_t - 1] + VALUE_A * betta[index_x - 1]) /
                (VALUE_C - alpha[index_x - 1] * VALUE_A);
        }
        function_grid[LIMIT_X - 1][index_t] = edge_func(nodes_t[index_t]);
        for (int index_x = 9; index_x > 0; --index_x) {
            function_grid[index_x][index_t] = alpha[index_x] * function_grid[index_x + 1][index_t] + betta[index_x];
        }
    }
    cout << endl << "Solution obtained by implicit finite difference schemes:" << endl;
    print_solution(function_grid, nodes_t);
    return 0;
}