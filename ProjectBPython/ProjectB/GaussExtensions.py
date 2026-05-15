# -*- coding: utf-8 -*-
"""
@Project: LinalgDat2022
@File: GaussExtensions.py

@Description: Project B Gauss extensions

"""

import math

from sys import path
path.append('../Core')
from Vector import Vector
from Matrix import Matrix

def AugmentRight(A: Matrix, v: Vector) -> Matrix:
    """
    Create an augmented matrix from a matrix and a vector.

    This function creates a new matrix by concatenating matrix A
    and vector v. See page 12 in "Linear Algebra for Engineers and
    Scientists", K. Hardy.

    Parameters:
         A: M-by-N Matrix
         v: M-size Vector
    Returns:
        the M-by-(N+1) matrix (A|v)
    """
    M = A.M_Rows
    N = A.N_Cols
    if v.size() != M:
        raise ValueError("number of rows of A and length of v differ.")

    B = Matrix(M, N + 1)
    for i in range(M):
        for j in range(N):
            B[i, j] = A[i, j]
        B[i, N] = v[i]
    return B


def ElementaryRowReplacement(A: Matrix, i: int, m: float, j: int) -> Matrix:
    """
    Replace row i of A by row i of A + m times row j of A.

    Parameters:
        A: M-by-N Matrix
        i: int, index of the row to be replaced
        m: float, the multiple of row j to be added to row i
        j: int, index or replacing row.

    Returns:
        A modified in-place after row replacement.
    """
    #Total columns
    N = A.N_Cols

    #Loops through row we are replacing
    for k in range(N):
        #Replaces value
        A[i, k] = A[i, k] + m * A[j, k]

    return A


def ElementaryRowInterchange(A: Matrix, i: int, j : int) -> Matrix:
    """
    Interchange row i and row j of A.

    Parameters:
        A: M-by-N Matrix
        i: int, index of the first row to be interchanged
        j: int, index the second row to be interchanged.

    Returns:
        A modified in-place after row interchange
    """
    M = A.M_Rows
    N = A.N_Cols
    #loops through columns
    for col in range(N):
       A[i, col], A[j, col] = A[j, col], A[i, col]
    return A


def ElementaryRowScaling(A: Matrix, i: int, c: float) -> Matrix:
    """
    Replace row i of A c * row i of A.

    Parameters:
        A: M-by-N Matrix
        i: int, index of the row to be replaced
        c: float, the scaling factor

    Returns:
        A modified in-place after row scaling.
    """
    N = A.N_Cols 
    for j in range(N):
        A[i, j] = c * A[i, j]
    return A
# Vi har en matrix i som vi skal scalerer, det gør vi ved at scalerer matrixen med konstanten c
# og så returnere den modificerede matrix A.



def ForwardReduction(A: Matrix) -> Matrix:
    """
    Forward reduction of matrix A.

    This function performs the forward reduction of A provided in the
    assignment text to achieve row echelon form of a given (augmented)
    matrix.

    Parameters:
        A:  M-by-N augmented matrix
    returns
        M-by-N matrix which is the row-echelon form of A (performed in-place,
        i.e., A is modified directly).
    """
    M = A.M_Rows
    N = A.N_Cols
    i = 0
    j = 0

    # Finds pivot column and element
    for column in range(N):
        for row in range(M):
            if round(A[row, column], 10) != 0:  #Rounding to compensate for very small non-zero values which should be zero
                (i, j) = (row, column)
                break
        else:
            continue # If inner loop is not broken, continue
        break
        
    # Ensures pivot row is at top 
    A = ElementaryRowInterchange(A, i, 0)

    # Reduces all elements below our pivot element
    for row in range(1, M):
        frac = A[row, j] / A[0, j]
        A = ElementaryRowReplacement(A, row, -frac, 0)

    # Gets submatrix without first row
    B = Matrix(M-1, N)
    for row in range(M-1):
        for column in range(N):
            B[row, column] = A[row + 1, column]

    # Check if more rows left
    if B.Size <= 0:
        return A

    # Check if remaining submatrix are 0 rows
    for row in range(M-1):
        if round(sum(B.Row(row)), 10) != 0: #Rounding to compensate for very small non-zero values which should be zero
            firstRowArray = [A.asArray()[0]]
            return Matrix.fromArray(firstRowArray + ForwardReduction(B).asArray())
    return A


def BackwardReduction(A: Matrix) -> Matrix:
    """
    Backward reduction of matrix A.

    This function performs the forward reduction of A provided in the
    assignment text given a matrix A in row echelon form.

    Parameters:
        A:  M-by-N augmented matrix in row-echelon form
    returns
        M-by-N matrix which is the reduced form of A (performed in-place,
        i.e., A is modified directly).
    """
    M = A.M_Rows
    N = A.N_Cols

    for i in range(M - 1, -1, -1):
        # find pivot i rækken
        pivot = -1
        for j in range(N):
            if round(A[i, j], 10) != 0:  # Rounding to compensate for very small non-zero values which should be zero
                pivot = j
                break
        # spring nul rækkerne over
        if pivot == -1:
            continue
        # gør pivot = 1
        ElementaryRowScaling(A, i, 1 / A[i, pivot])
        # fjerner alt over pivot
        for k in range(i):
            ElementaryRowReplacement(A, k, -A[k, pivot], i)
    return A

        




def GaussElimination(A: Matrix, v: Vector) -> Vector:
    """
    Perform Gauss elimination to solve for Ax = v.

    This function performs Gauss elimination on a linear system given
    in matrix form by a coefficient matrix and a right-hand-side vector.
    It is assumed that the corresponding linear system is consistent and
    has exactly one solution.

    Hint: combine AugmentRight, ForwardReduction and BackwardReduction!

    Parameters:
         A: M-by_N coefficient matrix of the system
         v: N-size vector v, right-hand-side of the system.
    Return:
         M-size solution vector of the system.
    """
    
    M = A.M_Rows
    N = A.N_Cols

    # Laver udvidet matrix
    A = AugmentRight(A, v)

    # Lav Gauss elimination ved at gøre de tideligere trin
    A = ForwardReduction(A)
    A = BackwardReduction(A)

    # Få løsningen ved at udtrække den sidste kolonne
    x = Vector(N)
    for i in range(N):
        x[i] = A[i, N]
    return x