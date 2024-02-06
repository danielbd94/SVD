"""
Created on Sat Dec 24 16:48:19 2022
@author: Ohad
@author: Daniel
"""

import numpy as np
from time import time

"""
Find the maximum absolute value in the matrix and its row and column index
Parameters: mat (np.ndarray): input matrix
Returns:    Tuple: (float, int, int) The maximum absolute value in the matrix, the row index and the column index of the maximum value
"""
def find_max_num(mat):
    input_matrix = mat.copy()
    input_matrix[np.eye(mat.shape[0], dtype=bool)] = 0 # Set diagonal element to 0
    row, col = np.unravel_index(np.argmax(np.abs(input_matrix)), input_matrix.shape)  # Get the indices of the maximum element
    return np.abs(input_matrix[row, col]), row, col

"""
Compute the new D matrix after a Jacobi rotation
Parameters: 
            matrix (np.ndarray): input matrix, 
            cos (float): cosine value, 
            sin (float): sine value, 
            p (int): row index, 
            q (int): column index
Returns:    np.ndarray: new D matrix
"""
def calculateNewDmatrix(matrix, cos, sin, p, q):
    D = matrix.copy()
    
    D_ii = matrix[p, p]
    D_ij = matrix[p, q]
    D_jj = matrix[q, q]
    
    D[p, p] = cos * cos * D_ii - 2 * sin * cos * D_ij + sin * sin * D_jj
    D[q, q] = sin * sin * D_ii + 2 * sin * cos * D_ij + cos * cos * D_jj
    D[p, q] = D[q, p] = (cos * cos - sin * sin) * D_ij + sin * cos * (D_ii - D_jj)
    
    for k in range(len(D)):
        if k != p and k != q:
            D[p, k] = D[k, p] = cos * matrix[p, k] - sin * matrix[q, k]
            D[q, k] = D[k, q] = sin * matrix[p, k] + cos * matrix[q, k]
            
    return D

"""
Compute the eigenvectors and eigenvalues of a matrix using the Jacobi algorithm
Parameters: 
            A (np.ndarray): input matrix
            epsilon (float): tolerance value, the algorithm stops when the largest off-diagonal element is below this value (default 1.0e-9)
            max_iterations (int): maximum number of iterations (default 2500)
Returns:    Tuple: (np.ndarray, np.ndarray) The eigenvectors matrix and the eigenvalues matrix
"""
def Jacobi(A, epsilon = 1.0e-9, max_iterations = 5500):
    n = len(A)
    D = A.copy()
    S = np.eye(n)
    # Performing the Jacobi's rotations until D becomes diagonal
    for i in range(max_iterations):
        max_val, p, q = find_max_num(D)
        if max_val < epsilon: # If the largest element is below epsilon - The matrix is already diagonal
            return S, np.diag(D)
        
        # Calculate the rotational angle - theta
        if D[p, p] == D[q, q]:
            theta = np.pi / 4
        else:
            a = (2 * D[p, q]) / (D[q, q] - D[p, p])
            theta = 0.5 * np.arctan(a)

        # Compute the matrix S1
        sin = np.sin(theta)
        cos = np.cos(theta)
        S1 = np.eye(n)
        S1[p,p] = S1[q,q] = cos
        S1[p,q] = -1 * sin
        S1[q,p] = sin
        
        S = np.dot(S, S1) # S = S * S1
        S = np.dot(S1.T, S)
        D = calculateNewDmatrix(D,cos,sin,p,q)
    
    return S, np.diag(D)

"""
Compute the singular value decomposition of a matrix using the Jacobi algorithm
Parameters: A (np.ndarray): input matrix       
Returns:    Tuple: (np.ndarray, np.ndarray, np.ndarray, np.ndarray) The left singular vectors matrix, the singular values matrix, the right singular vectors matrix, and the singular values vector
"""
def SVD(A):
    # Compute A.T * A
    A_AT = np.dot(A, A.T)
    
    # Find the eigenvectors and eigenvalues of A * A.T
    eigenvectors, eigenvalues = Jacobi(A_AT)

    # Get the number of rows and columns of A
    r, c = A.shape

    # Initialize S and V as zero matrices
    S = np.zeros((r, c))
    V = np.zeros((c, c))

    # Set the elements of S to the square roots of the eigenvalues
    singular_values = np.sqrt(eigenvalues)
    S = np.diag(singular_values)
    S = np.pad(S, ((0, 0), (0, c - r)), 'constant', constant_values=(0, 0))

    # Set the left singular vectors U to eigenvectors of A * A.T
    U = eigenvectors
    
    # Setting the V matrix
    for i in range(len(eigenvalues)):
        if i < r and i < c and S[i][i] != 0:
            temp = 1 / S[i][i]
            tempM = np.dot(temp,A.T)
            V[i] = tempM.dot(U[:,i])
    return U, S, V, singular_values
 
    
"""
Read matrix from a text file.
Parameters: MatrixFile (str): input file name
Returns:    np.ndarray: The matrix read from the file        
"""
def readMatrixFromTXT(MatrixFile):
    with open(MatrixFile, 'r') as f:
        matrix = [[np.float64(num) for num in line.split(',')] for line in f] # Loop lines in the file, then loop each number in the line
    matrix = np.asarray(matrix, dtype = np.float64) # Creating NumPy array from list of lists
    return matrix

if __name__ == '__main__':
    MatrixFileName = "svd_matrix(1000x2000).txt"
    print("Creating the matrix from text file {}".format(MatrixFileName))
    inputMatrix = readMatrixFromTXT(MatrixFileName)
    print("Calculating singular value decomposition")
    numPySVDstartTime = time()
    numpy_U, numpy_S, numpy_V_T = np.linalg.svd(inputMatrix)
    print("Total duration of NumPy's SVD calculation: {0:.15f} seconds".format(time() - numPySVDstartTime))
    ourSVDstartTime = time()
    u, s, vt, t = SVD(inputMatrix)
    print("Total duration of our's   SVD calculation: {0:.15f} seconds".format(time() - ourSVDstartTime))
    res = u @ s @ vt
    t = np.flip(np.sort(t))
    eigensCompare = np.column_stack((numpy_S, t, np.abs(numpy_S - t), [(abs(x - y) / max(x, y)) for x, y in zip(numpy_S, t)]))
    if np.allclose(inputMatrix, res):
        print("\nOur implementation worked as expected: ")
    else:
        print("Calculation error!")
        
    difference = np.linalg.norm(inputMatrix - res)

    print("Difference between A and reconstructed matrix:", difference)            
        
    # print("U:")
    # print(u)
    # print("\nS:")
    # print(s)
    # print("\nVT:")
    # print(vt)
    # print("\nA = U * S * V.T:")
    # print(res)       
