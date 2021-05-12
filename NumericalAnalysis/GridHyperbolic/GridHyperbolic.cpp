/*
* Semester 2. Laboratory work No.2
* The grid method for solving an initial-boundary value problem for an equation of hyperbolic type.
* 
*/
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

const double NUMBER_VARIANT = 5;
const double LENGTH_SEGMENT = 1.0;
const double ALPHA = 0.5 + 0.1 * NUMBER_VARIANT;
const int LIMIT_X = 10;
const int LIMIT_T = 10;
const double STEP_X = LENGTH_SEGMENT / static_cast<double>(LIMIT_X);
const double STEP_T = LENGTH_SEGMENT / static_cast<double>(LIMIT_T);
const double VALUE_S = pow(STEP_T, 2) / pow(STEP_X, 2);

double first_edge_func_t(const double& value_t) { return 1.0 / (ALPHA * value_t + 1.0); }

double second_edge_func_t(const double& value_t) { return 1.0 / (ALPHA * value_t + 2.0); }

double func_betta(const double& value_x) { return -ALPHA / pow((1.0 + value_x), 2); }

double func_f(const double& value_x, const double& value_t) { return 2.0 * (pow(ALPHA, 2) - 1.0) / pow((value_x + ALPHA * value_t + 1.0), 3); }

double first_edge_func_x(const double& value_x) { return 1.0 / (1.0 + value_x); }

double second_edge_func_x(const double& value_x) { return 2.0 / pow(1.0 + value_x, 3); }

int main() {
	// Building a grid from a collection of points {(x_i, t_j): x_i = i * h, t_j = j * tau}
	vector<double> value_x(LIMIT_X + 1);
	vector<double> value_t(LIMIT_T + 1);
	for (int index_point = 0; index_point <= LIMIT_X; ++index_point) {
		value_x[index_point] = index_point * STEP_X;
		value_t[index_point] = index_point * STEP_T;
	}
	// Our grid
	vector<vector<double>> grid_y(LIMIT_X + 1, vector<double>(LIMIT_T + 1, 0.0));
	// Calculate values on edges for function y
	for (int index_point = 0; index_point <= LIMIT_T; ++index_point) {
		grid_y[0][index_point] = first_edge_func_t(value_t[index_point]);
		grid_y[LIMIT_T][index_point] = second_edge_func_t(value_t[index_point]);
	}
	for (int index_point = 1; index_point < LIMIT_X; ++index_point) {
		grid_y[index_point][0] = first_edge_func_x(value_x[index_point]);
		grid_y[index_point][1] = first_edge_func_x(value_x[index_point]) + STEP_T * func_betta(value_x[index_point]) +
			(pow(STEP_T, 2) / 2.0) * (second_edge_func_x(value_x[index_point]) + func_f(value_x[index_point], 0));
	}
	// Calculate values grid function on other layers
	for (int index_grid_first = 1; index_grid_first < LIMIT_X; ++index_grid_first) {
		for (int index_grid_second = 1; index_grid_second < LIMIT_X; ++index_grid_second) {
			grid_y[index_grid_second][index_grid_first + 1] = VALUE_S * grid_y[index_grid_second + 1][index_grid_first] +
				2.0 * (1.0 - VALUE_S) * grid_y[index_grid_second][index_grid_first] +
				VALUE_S * grid_y[index_grid_second - 1][index_grid_first] - grid_y[index_grid_second][index_grid_first - 1] +
				pow(STEP_T, 2) * func_f(value_x[index_grid_second], value_t[index_grid_first]);
		}
	}

	cout << "Grid function values on 10 layers:" << endl;
	for (int index_point_x = 0; index_point_x <= LIMIT_X; ++index_point_x) {
		cout << setprecision(2);
		cout << setw(9) << left << value_x[LIMIT_X - index_point_x];
		cout << setprecision(4);
		for (int index_point_t = 0; index_point_t <= LIMIT_T; ++index_point_t) {
			cout << fixed << left << setw(9) << grid_y[index_point_t][LIMIT_X - index_point_x];
		}
		cout << endl;
	}
	cout << setprecision(2) << setw(9) << "t/x";
	double step = 0.0;
	for (int index_step = 0; index_step < 11; ++index_step) {
		cout << setw(9) << left << step;
		step += 0.1;
	}
	return 0;
}