#include "main.h"

tuple<vector<vector<double>>, vector<double>, vector<vector<double>>> SVD(vector<vector<double>>& input_matrix, int max_iterations);
void decomposition_checker(vector<vector<double>>& input_matrix, vector<vector<double>> U, vector<double> Sigma, vector<vector<double>> V_T);