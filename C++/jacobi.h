#include "main.h"

tuple<vector<vector<double>>, double, double> jacobi_rotation_matrix_calculation(vector<vector<double>>& matrix, int p, int q);
void calculateNewDmatrix(vector<vector<double>>& mtx, double cos, double sin, int i, int j);
bool update_off_diagonal_values(vector<vector<double>>& matrix);
vector<tuple<double, vector<double>>>eigenvector_eigenvalues_tuples_vector(vector<vector<double>>& eigenvectors, vector<double>& lamdas);
tuple<vector<vector<double>>, vector<double>>Jacobi(vector<vector<double>> matrix, int max_iterations);