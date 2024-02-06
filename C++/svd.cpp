#include "main.h"
#include "jacobi.h"
#include "svd.h"
#include "matrix_operations.h"

/*used to calculate the singular value decompostion
input: input matrix from text file
output: tuple of three elements: matrix U, vector of eigenvalues, matrix V_T
*/
tuple<vector<vector<double>>, vector<double>, vector<vector<double>>> SVD(vector<vector<double>>& input_matrix, int max_iterations) {
	//vector<vector<double>> A_matrix = input_matrix; //matrix from text file
	vector<vector<double>> AT = transpose(input_matrix); //calculate A transpose
	vector<vector<double>> A_AT = dot_product(AT, input_matrix); //calculate A*A_T
	auto jacobi_result = Jacobi(A_AT, max_iterations); //calculate eigenvalues and eigenvectors using Jacobi itartions method
	vector<vector<double>> eigenvectors = get<0>(jacobi_result); //get the eigenvectors matrix from the result tuple
	vector<double> eigenvalues = get<1>(jacobi_result); //get the eigenvalues vector from the result tuple
	auto eigenvector_eigenvalue_tuples_vector = eigenvector_eigenvalues_tuples_vector(eigenvectors, eigenvalues); //create a vector of tuples [eigenvector, eigenvalue]
	vector<double> eigenvalues_vector(eigenvector_eigenvalue_tuples_vector.size()); //initalize vector for eigenvalues
	vector<vector<double>> U; //initalize matrix U
	vector<vector<double>> V_T; //initalize matrix V_T
	for (int i = 0; i < eigenvector_eigenvalue_tuples_vector.size(); i++) {
		auto& eigen_tuple = eigenvector_eigenvalue_tuples_vector[i]; //reference to eigen tuple in index i
		double eigenvalue = get<0>(eigen_tuple); //get the eigenvalue from the tuple
		vector<double> eigenvector = get<1>(eigen_tuple); //get the eigenvector from the tuple
		double singular_value = sqrt(eigenvalue); //calculate the singular value for the sigma matrix
		double eigenvector_len = vec_len(eigenvector); //get the length of the eigenvector
		vector<double> normalized_eigenvector = dot_product(1 / eigenvector_len, eigenvector); //calculate the normalized eigenvector for matrix V_T
		double scalar = 1 / singular_value;
		vector<vector<double>> temp_v_matrix = { normalized_eigenvector }; //add vector to matrix temp V matrix
		temp_v_matrix = transpose(temp_v_matrix); //transpose temp matrix 
		temp_v_matrix = dot_product(input_matrix, temp_v_matrix); //matrices multiplication
		temp_v_matrix = transpose(temp_v_matrix); //transpose after multiplication
		vector<double> u = dot_product(scalar, temp_v_matrix[0]); //calculate vector for matrix U
		eigenvalues_vector[i] = singular_value; //add the singular value to a result vector
		U.push_back(u); //add vector u to result matrix U
		V_T.push_back(normalized_eigenvector); //add vector v to result matrix V_T
	}
	U = transpose(U); //transpose matrix U
	return make_tuple(U, eigenvalues_vector, V_T); //result tuple
}

/*used to check if decomposition is correct
input: reference to input matrix, matrix U by value, vector of eigenvalues Sigma by value and matrix V_T
output: calculate and printing the singular value decompostion result
*/
void decomposition_checker(vector<vector<double>>& input_matrix, vector<vector<double>> U, vector<double> Sigma, vector<vector<double>> V_T) {
	cout << "Checking the decomposition: " << endl;
	cout << "Input matrix = " << endl;
	print_matrix(input_matrix);
	cout << endl;
	vector<vector<double>> Sigma_matrix = create_sigma_matrix(Sigma, U[0].size());
	vector<vector<double>> temp_mtx = mul_matrix_diag(U, Sigma_matrix);
	vector<vector<double>> res = dot_product(temp_mtx, V_T);
	cout << "Result of the 3 matrices multiplication is:\nU * S * V.T =" << endl;
	print_matrix(res);
}
