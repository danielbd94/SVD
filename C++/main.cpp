#include "main.h"
#include "jacobi.h"
#include "svd.h"
#include "matrix_operations.h"

int main() {
	vector<vector<double>> matrix;
    string file_name = "svd_matrix(1000x2000).txt";
    //string file_name = "example8.txt";
	cout << "Creating the matrix from text file ";
	cout << file_name << "\n";
	matrix = read_file(file_name, ',');
	cout << "Calculating singular value decomposition\n";
    long s1 = Get_Time();
    auto tuple1 = SVD(matrix, 5500);
    cout << "Total duration of our's SVD calculation: " << Get_Time() - s1 << " seconds" << endl;
    vector<vector<double>> U = get<0>(tuple1);
    vector<double> Sigma = get<1>(tuple1);
    vector<vector<double>> V_T = get<2>(tuple1);
    decomposition_checker(matrix, U, Sigma, V_T);
    sort(Sigma.begin(), Sigma.end(), greater<double>());
    print_vec(Sigma);
	return 0;
}


long Get_Time() {
    using chrono::high_resolution_clock;
    auto t = high_resolution_clock::now();
    auto nanosec = t.time_since_epoch();
    return nanosec.count() / 1000000000;
}


vector<vector<double>> read_file(const string& file_name, const char& delimiter) {
    vector<vector<double>> matrix;
    ifstream file(file_name);
    if (!file) {
        throw runtime_error("Could not open file: " + file_name);
    }
    string line;
    while (getline(file, line)) {
        vector<double> row;
        size_t pos = 0;
        while ((pos = line.find(delimiter)) != string::npos) {
            string token = line.substr(0, pos);
            row.push_back(stod(token));
            line.erase(0, pos + 1);
        }
        row.push_back(stod(line));
        matrix.push_back(row);
    }
    return matrix;
}