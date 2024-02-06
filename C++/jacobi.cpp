#include "jacobi.h"
#include "matrix_operations.h"

/*used to calculate the eigenvalues and eigenvectors (Jacobi iterations)
input: reference to matrix A*A_T
output: tuple with three values: matrix of eigenvector, vector of diagonal values (eigenvalues) and number of iterations done
*/
tuple<vector<vector<double>>, vector<double>>Jacobi(vector<vector<double>> D, int max_iterations) {
	if (D.size() != D[0].size()) {
		string str = "The given matrix is not square matrix.";
		throw invalid_argument(str);
	}
	int n = D.size(), p, q;
	double max_val;
	vector<vector<double>> S = create_identity_matrix(n);

	for (int i = 0; i < max_iterations; i++) {
		auto tuple1 = find_matrix_max_value(D); //tuple1 = [max_val,i index,j index]
		max_val = get<0>(tuple1);
		p = get<1>(tuple1);
		q = get<2>(tuple1);
		// cout << max_val << " " << p << " " << q << endl;
		if (max_val < JACOBI_ERROR)
			return make_tuple(S, get_matrix_diagonal(D));

		// Compute the matrix S1
		auto tuple2 = jacobi_rotation_matrix_calculation(D, p, q);
		vector<vector<double>> S1 = get<0>(tuple2);
		double cosinus = get<1>(tuple2);
		double sinus = get<2>(tuple2);

		//vector<vector<double>> S1_T = transpose(S1);
		//S = dot_product(S, S1);
		//S = dot_product(S1_T, S);
		S = calculate_rotation_matrix(S, S1, p, q);
		calculateNewDmatrix(D, cosinus, sinus, p, q);
		
		if (update_off_diagonal_values(D))
			break;
	}
	return make_tuple(S, get_matrix_diagonal(D));
}

/*used to calculate the Jacobi's rotation matrix
input: reference to matrix and two integers
output: tuple with three values: Jacobi rotation matrix, cos(theta) value, sin(theta) value
*/
tuple<vector<vector<double>>, double, double>jacobi_rotation_matrix_calculation(vector<vector<double>>& matrix, int p, int q) {
	//check if matrix is squered
	if (matrix.size() != matrix[0].size()) {
		string str = "The given matrix is not a square matrix.";
		throw invalid_argument(str);
	}
	double theta;
	int n = matrix.size(); //size of input matrix
	vector<vector<double>> S = create_identity_matrix(n); //initalize identity matrix in size nxn

	//choose the theta angle according to jacobi formulas
	if (matrix[q][q] == matrix[p][p])
		theta = PI / 4;
	else {
		double a = (2 * matrix[p][q]) / (matrix[q][q] - matrix[p][p]);
		theta = 0.5 * atan(a);
	}

	//set values of rotation to the Jacobi matrix
	double cosinus, sinus;
	S[p][p] = S[q][q] = cosinus = cos(theta);
	S[p][q] = -1 * sin(theta);
	S[q][p] = sinus = sin(theta);

	return make_tuple(S, cosinus, sinus);
}

/*used to calculate Jacobi's matrix according formulas
input: reference to matrix ,cos(theta) value, sin(theta) value and two integers
output: update the Jacobi's matrix (the method update the input matrix)
*/
void calculateNewDmatrix(vector<vector<double>>& mtx, double cos, double sin, int i, int j) {
	double a_ii = mtx[i][i];
	double a_ij = mtx[i][j];
	double a_jj = mtx[j][j];
	double a_ji = mtx[j][i];

	mtx[i][i] = cos * cos * a_ii - 2 * sin * cos * a_ij + sin * sin * a_jj;
	mtx[j][j] = sin * sin * a_ii + 2 * sin * cos * a_ij + cos * cos * a_jj;
	mtx[i][j] = mtx[j][i] = (cos * cos - sin * sin) * a_ij + sin * cos * (a_ii - a_jj);

	for (int k = 0; k < mtx.size(); k++)
		if (k != i && k != j) {
			double a_ik = mtx[i][k];
			double a_jk = mtx[j][k];
			mtx[i][k] = mtx[k][i] = cos * a_ik - sin * a_jk;
			mtx[j][k] = mtx[k][j] = sin * a_ik + cos * a_jk;
		}
}

/*used to check the matrix after Jacobi interations and update the off-diagonal values to 0
input: reference to matrix
output: true if update was a sucsses
*/
bool update_off_diagonal_values(vector<vector<double>>& matrix) {
	int rows_num = matrix.size();
	int cols_num = matrix[0].size();
	//check if matrix squerd
	if (rows_num != cols_num) {
		string str = "The given matrix is not square matrix.";
		throw invalid_argument(str);
	}
	//update off-diagonal values to 0
	for (int i = 0; i < rows_num; i++)
		for (int j = i + 1; j < cols_num; j++) {
			if (abs(matrix[i][j]) < JACOBI_ERROR)
				matrix[i][j] = matrix[j][i] = 0;
			else
				return false;
		}
	return true;
}

/*used to calculte a vector of tuples: [eigenvector,eigenvalue]
input: eigenvectors matrix and eigenvalues vector
output: vector with a tuple in each place
*/
vector<tuple<double, vector<double>>> eigenvector_eigenvalues_tuples_vector(vector<vector<double>>& eigenvectors, vector<double>& lamdas) {
	vector<tuple<double, vector<double>>> t_vecs;
	bool flag = false;
	for (int i = 0; i < lamdas.size(); i++)
		if (lamdas[i] > 0) {
			auto tuple = make_tuple(lamdas[i], get_matrix_column(eigenvectors, i));
			if (t_vecs.size() == 0)
				t_vecs.push_back(tuple);
			else {
				for (int j = 0; j < t_vecs.size(); j++) {
					auto tuple1 = t_vecs[j];
					if (get<0>(tuple1) <= lamdas[i]) {
						auto itPos = t_vecs.begin() + j;
						t_vecs.insert(itPos, tuple);
						flag = true;
						break;
					}
				}
				if (!flag)
					t_vecs.push_back(tuple);
				flag = false;
			}
		}
	return t_vecs;
}