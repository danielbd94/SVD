#include "main.h"
#include "matrix_operations.h"

/*used to create sigma matrix:
input: reference to vector, size
output: squared matrix with size: size_val with the valules of vector v on its diagonal (rest matrix values 0)
*/
vector<vector<double>> create_sigma_matrix(vector<double>& diagonal, int size_val) {
	int n = diagonal.size();

	if (n > size_val) {
		string str = "The desire matrix size is smaller than the amount of elements in diagonal vectors";
		throw invalid_argument(str);
	}

	vector<double> v(size_val, 0);
	vector<vector<double>> new_matrix(size_val, v);

	for (int i = 0; i < n; i++)
		new_matrix[i][i] = diagonal[i];

	return new_matrix;
}

/*used to calculate multiplication between matrix to another matrix diagonal:
input: reference to 2 matrices
output: matrix after multiplication by another matrix diagonal
*/
vector<vector<double>> mul_matrix_diag(vector<vector<double>>& matrix, vector<vector<double>>& diagonal) {
	//get sizes of matrices
	int row_matrix = matrix.size();
	int col_matrix = matrix[0].size();
	int row_diagonal = diagonal.size();
	int col_diagonal = diagonal[0].size();

	if (col_matrix != row_diagonal) { //check if sizes fit
		string str = "Can't multiply (" + to_string(row_matrix) + ", " + to_string(col_matrix) + ") by (" + to_string(row_diagonal) + ", " + to_string(col_diagonal) + ").";
		throw invalid_argument(str);
	}
	//multiplication of each value in matrix by the diagonal values
	for (int i = 0; i < row_matrix; i++)
		for (int j = 0; j < row_diagonal; j++)
			matrix[i][j] = matrix[i][j] * diagonal[j][j];
	return matrix;
}

/*used to calculate vector length
input: reference to vector
output: vector length
*/
double vec_len(vector<double>& v) {
	double sum = 0;
	for (double num : v)
		sum += pow(num, 2);
	return sqrt(sum);
}

/*used to calculate the rotation matrix for Jacobi
input: reference to 2 matrices and 2 int values for the positions
output: rotation matrix
*/
vector<vector<double>> calculate_rotation_matrix(vector<vector<double>>& matrix, vector<vector<double>>& rotation, int p, int q) {
	vector<vector<double>> matrix_copy = matrix; //copy of matrix
	//get sizes of matrices
	int row_matrix = matrix_copy.size();
	int col_matrix = matrix_copy[0].size();
	int row_rotation = rotation.size();
	int col_rotation = rotation[0].size();

	if (col_matrix != row_rotation) { //check if sizes fit
		string str = "Can't multiply (" + to_string(row_matrix) + ", " + to_string(col_matrix) + ") by (" + to_string(row_rotation) + ", " + to_string(col_rotation) + ").";
		throw invalid_argument(str);
	}

	rotation = transpose(rotation); //transpose matrix rotation

	for (int i = 0; i < row_matrix; i++) {
		double index_ip;
		double index_iq;
		index_ip = dot_product(matrix_copy[i], rotation[p]);
		index_iq = dot_product(matrix_copy[i], rotation[q]);
		matrix_copy[i][p] = index_ip;
		matrix_copy[i][q] = index_iq;
	}
	return matrix_copy;
}

/*used to calculate multification between two vectors
input: reference to 2 vectors v1,v2
output: vectors multiplication (scalar)
*/
double dot_product(vector<double>& v1, vector<double>& v2) {
	int n;

	if ((n = v1.size()) != v2.size()) {
		string str = "Vectors are not in the same size.";
		throw invalid_argument(str);
	}

	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += (v1[i] * v2[i]);

	return sum;
}

/*used to calculate multification between two matrices
input: reference to 2 matrices
output: matrix of multiplication between matrix1, matrix2
*/
vector<vector<double>> dot_product(vector<vector<double>>& matrix1, vector<vector<double>> matrix2) {
	//get matrices size
	int row_matrix1 = matrix1.size();
	int col_matrix1 = matrix1[0].size();
	int row_matrix2 = matrix2.size();
	int col_matrix2 = matrix2[0].size();
	//check if sizes fit
	if (col_matrix1 != row_matrix2) {
		string str = "Can't multiply (" + to_string(row_matrix1) + ", " + to_string(col_matrix1) + ") by (" + to_string(row_matrix2) + ", " + to_string(col_matrix2) + ").";
		cout << str << endl;
		throw invalid_argument(str);
	}

	matrix2 = transpose(matrix2); //transpose matrix 2
	int row_matrix2_T = matrix2.size(); //get the row size after transpose

	vector<double> v(row_matrix2_T); //vector constructor - initalize new vector to the size of row_matrix2_T
	vector<vector<double>> new_matrix(row_matrix1, v); //initalize new matrix with the needed sizes

	//multifaction of matrices in a loop mul row vector by col vector
	for (int i = 0; i < row_matrix1; i++)
		for (int j = 0; j < row_matrix2_T; j++) {
			double num = dot_product(matrix1[i], matrix2[j]);
			new_matrix[i][j] = num;
		}

	return new_matrix;
}

/*used to calculate multification of vector by scalar
input: scalar and reference to vector
output: vector after multification with the scalar
*/
vector<double> dot_product(double scalar, vector<double>& v) {
	for (int i = 0; i < v.size(); i++)
		v[i] = scalar * v[i];
	return v;
}

/*used to get a spacific column of a matix
input: reference to matrix and index of requested column
output: vector - the column requested
*/
vector<double> get_matrix_column(vector<vector<double>>& matrix, int column) {
	//check if the column exists
	if (matrix[0].size() < column) {
		string str = "The given matrix do not have column number " + to_string(column) + ".";
		throw invalid_argument(str);
	}

	int n = matrix.size(); //get the size of the matrix column
	vector<double> v(n); //intialize vector in the size of matrix column

	//intalize the vector values to the column values
	for (int i = 0; i < n; i++)
		v[i] = matrix[i][column];
	return v;
}

/*used to calculate the transpose of a matrix
input: reference to matrix
output: the matrix after transpose (changes the input matrix itself)
*/
vector<vector<double>> transpose(vector<vector<double>>& matrix) {
	//get matix sizes
	int row = matrix.size();
	int col = matrix[0].size();
	//initalize new matrix with the same dimensions
	vector<double> v(row);
	vector<vector<double>> new_matrix(col, v);

	//transpose matrix - each value[i,j] = value[j,i]
	for (int i = 0; i < col; i++)
		for (int j = 0; j < row; j++)
			new_matrix[i][j] = matrix[j][i];
	return new_matrix;
}

/*used to crate the identity matrix
input: integer n
output: squared identity matrix in size nxn 
*/
vector<vector<double>> create_identity_matrix(int n) {
	vector<double> v(n, 0);
	vector<vector<double>> I(n, v);
	//set the diagonal values to 1 and the rest stays 0
	for (int i = 0; i < n; i++)
		I[i][i] = 1;
	return I;
}

/*used to get the diagonal of a given squared matrix
input: reference to a matix
output: vector that represents the matrix's diagonal
*/
vector<double> get_matrix_diagonal(vector<vector<double>>& matrix) {
	//get matrix sizes
	int row = matrix.size();
	int col = matrix[0].size();
	//check if sizes fit
	if (row != col) {
		string str = "The given matrix is not square matrix.";
		throw invalid_argument(str);
	}

	vector<double> diagonal(row); //vector constructor - initalize vector in the size of row
	//intalize the vector to diagonal values
	for (int i = 0; i < row; i++)
		diagonal[i] = matrix[i][i];
	return diagonal;
}

/*used to find the max value in a given matrix and its index
input: reference to matrix
output: tuple with three values: max value, row, column
*/
tuple<double, int, int> find_matrix_max_value(vector<vector<double>>& matrix) {
	//get matrix sizes
	int row = matrix.size();
	int col = matrix[0].size();
	int p, q;
	double max_val = DBL_MIN;
	//find the max value and its index [i,j]
	for (int i = 0; i < row; i++)
		for (int j = i + 1; j < col; j++) {
			if (max_val < abs(matrix[i][j])) {
				max_val = abs(matrix[i][j]);
				p = i;
				q = j;
			}
		}
	return make_tuple(max_val, p, q);
}

/*
	print_matrix method print the given matrix to the console
*/
void print_matrix(vector<vector<double>>& matrix) {
	bool flag = false;
	for (int i = 0; i < matrix.size(); i++)	{
		if (matrix.size() == 1) {
			cout << "[\n ";
			print_vec(matrix[i]);
			cout << "]" << endl;
		}
		else if (i == 0) {
			cout << "[ ";
			print_vec(matrix[i]);
		}
		else if (i + 1 == matrix.size()) {
			cout << " ";
			print_vec(matrix[i]);
			cout << "] " << endl;
		}
		else {
			cout << " ";
			print_vec(matrix[i]);

			//The case the given matrix is bigger than 10 elements:
			if (i > 2 && matrix.size() > 10 && !flag) {
				cout << " ...	\n";
				i = matrix.size() - 4;
				flag = true;
			}
		}
	}
}

/*
	print_vec method print the given vector to the console.
*/
void print_vec(vector<double>& v) {
	bool flag = false;
	for (int i = 0; i < v.size(); i++) {
		if (v.size() == 1)
			cout << "[" << v[i] << "]" << endl;
		else if (i == 0)
			cout << "[" << v[i] << ",	";
		else if (i + 1 == v.size())
			cout << v[i] << "]" << endl;
		else {
			cout << v[i] << ",	";
			//The case the given vector is bigger than 10 elements:
			if (i > 2 && v.size() > 10 && !flag) {
				cout << "...	";
				i = v.size() - 4;
				flag = true;
			}
		}
	}
}