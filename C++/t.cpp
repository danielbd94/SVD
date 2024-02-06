#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <tuple>
#include <math.h>
#include <chrono>

using namespace std;

const double epsilon = 1.0e-9;
const int max_iterations = 100000;
#define PI						3.14159265358979323846

void print_matrix(vector<vector<double>> matrix) {
    
    for (auto row : matrix) {
        for (auto element : row) {
            cout << element << " ";
        }
        cout << endl;
    }
}

vector<vector<double>> mat_mult(vector<vector<double>> A, vector<vector<double>> B) {
    int n = A.size();
    int m = B.size();
    int l = B.size();
    vector<vector<double>> C(n, vector<double>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < l; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

vector<vector<double>> transpose(vector<vector<double>> A) {
    int n = A.size();
    int m = A.size();
    vector<vector<double>> AT(m, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            AT[j][i] = A[i][j];
        }
    }
    return AT;
}

// Function to find the maximum non-diagonal element of a matrix
tuple<double, int, int> find_max_num(vector<vector<double>>& matrix) {
    //get matrix sizes
    int row = matrix.size();
    int col = matrix.size();
    int p, q;
    double max_val;
    //find the max value and its index [i,j]
    for (int i = 0; i < row; i++)
        for (int j = i + 1; j < col; j++) {
            max_val = abs(matrix[i][j]);
            p = i;
            q = j;
        }

    return make_tuple(max_val, p, q);
}

// Jacobi's rotation algorithm
vector<vector<double>> Jacobi(vector<vector<double>> A, double epsilon = 1.0e-9, int max_iterations = 100000) {
    int n = A.size();
    double max_val;
    int p, q;
    vector<vector<double>> D = A;
    vector<vector<double>> S = vector<vector<double>>(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        S[i][i] = 1;
    }
    // Performing the Jacobi's rotations until D becomes diagonal
    for (int i = 0; i < max_iterations; i++) {
        auto tuple1 = find_max_num(D); //tuple1 = [max_val,i index,j index]
        max_val = get<0>(tuple1);
        p = get<1>(tuple1);
        q = get<2>(tuple1);
        if (max_val < epsilon) { // If the largest element is below epsilon - The matrix is already diagonal
            return S;
        }
        // Calculate the rotational angle - theta
        double theta;
        if (D[p][p] == D[q][q]) {
            theta = PI / 4;
        }
        else {
            double a = (2 * D[p][q]) / (D[q][q] - D[p][p]);
            theta = 0.5 * atan(a);
        }
        // Compute the matrix S1
        double sinus = sin(theta);
        double cosinus = cos(theta);
        vector<vector<double>> S1 = vector<vector<double>>(n, vector<double>(n));
        for (int i = 0; i < n; i++) {
            S1[i][i] = 1;
        }
        S1[p][p] = S1[q][q] = cosinus;
        S1[p][q] = -1 * sinus;
        S1[q][p] = sinus;
        // Update S and D
        S = mat_mult(S, S1);
        D = mat_mult(mat_mult(transpose(S1), D), S1);
    }
    return S;
}

vector<vector<double>> SVD(vector<vector<double>> A) {
    int r = A.size();
    int c = A.size();
    print_matrix(A);
    // Compute A.T * A
    vector<vector<double>> A_AT = mat_mult(A, transpose(A));
    print_matrix(A_AT);
    // Find the eigenvectors of A * A.T
    vector<vector<double>> eigenvectors = Jacobi(A_AT, epsilon, max_iterations);

    // Initialize S and V as zero matrices
    vector<vector<double>> S(r, vector<double>(c));
    vector<vector<double>> V(c, vector<double>(c));

    // Compute S and V
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            S[i][j] = sqrt(A_AT[i][i]);
            V[i][j] = eigenvectors[i][j] * S[i][j];
        }
    }
    return V;
}

vector<vector<double>> create_matrix_from_file(string file_name) {
    vector<vector<double>> matrix;

    ifstream file(file_name);
    string line;
    while (getline(file, line)) {
        vector<double> row;
        stringstream ss(line);
        double value;
        while (ss >> value) {
            row.push_back(value);
        }
        matrix.push_back(row);
    }

    return matrix;
}

//long Get_Time() {
//    using chrono::high_resolution_clock;
//    auto t = high_resolution_clock::now();
//    auto nanosec = t.time_since_epoch();
//    return nanosec.count() / 1000000000;
//}
//
//int main() {
//    // Get the file name from the user
//    // string file_name = "svd_matrix(1000x2000).txt";
//    string file_name = "example6.txt";
//    // Create the matrix from the input file
//    vector<vector<double>> matrix = create_matrix_from_file(file_name);
//
//    // Perform SVD on the matrix
//    vector<vector<double>> V = SVD(matrix);
//    cout << "The V matrix is :" << endl;
//    for (auto i : V) {
//        for (auto j : i)
//            cout << j << " ";
//        cout << endl;
//    }
//
//    return 0;
//}