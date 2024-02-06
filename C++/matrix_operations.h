#include "main.h"

vector<vector<double>> create_sigma_matrix(vector<double>& diagonal, int size);
vector<vector<double>> mul_matrix_diag(vector<vector<double>>& matrix, vector<vector<double>>& diagonal);
vector<vector<double>> calculate_rotation_matrix(vector<vector<double>>& matrix, vector<vector<double>>& rotation, int p, int q);
double vec_len(vector<double>& v);
vector<vector<double>> matrixMultiply(vector<vector<double>> matrix1, vector<vector<double>> matrix2);
double dot_product(vector<double>& v1, vector<double>& v2);
vector<vector<double>> dot_product(vector<vector<double>>& matrix1, vector<vector<double>> matrix2);
vector<double> dot_product(double scalar, vector<double>& v);
vector<vector<double>> transpose(vector<vector<double>>& matrix);
vector<double> get_matrix_column(vector<vector<double>>& matrix, int column);
tuple<double, int, int> find_matrix_max_value(vector<vector<double>>& matrix);
vector<vector<double>> create_identity_matrix(int n);
vector<double> get_matrix_diagonal(vector<vector<double>>& matrix);
void print_matrix(vector<vector<double>>& matrix);
void print_vec(vector<double>& v);