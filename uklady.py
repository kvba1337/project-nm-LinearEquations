import math
import time
import matplotlib.pyplot as plt
from copy import deepcopy

# zadanie A
def create_matrix_and_vector(N = 952, a1 = 7, a2 = -1, a3 = -1):
    A = [[0 for _ in range(N)] for _ in range(N)]

    for y in range(N):
        for x in range(N):
            if x == y:
                A[y][x] = a1
            elif x == y - 1 or x == y + 1:
                A[y][x] = a2
            elif x == y - 2 or x == y + 2:
                A[y][x] = a3

    b = [[math.sin(4*n)] for n in range(1, N + 1)]

    return A, b

def print_matrix_and_vector(A, b):
    print("Macierz A:")
    for idx, row in enumerate(A):
        if idx < 10:
            print(row)

    print("\nWektor b:")
    for idx, val in enumerate(b):
        if idx < 10:
            print(val)

# zadanie B
def norm(vector: list):
    return sum([(x[0])**2 for x in vector])**2

def mul(matrix1: list, matrix2: list):
    x1 = len(matrix1[0])
    x2 = len(matrix2[0])
    y1 = len(matrix1)

    matrix = [[0 for _ in range(x2)] for _ in range(y1)]

    for y in range(len(matrix)):
        for x in range(len(matrix[0])):
            matrix[y][x] = sum([matrix1[y][z] * matrix2[z][x] for z in range(x1)])

    return matrix

def sub(matrix1: list, matrix2: list) -> list:
    matrix = [[0 for _ in range(len(matrix1[0]))] for _ in range(len(matrix1))]

    for y in range(len(matrix1)):
        for x in range(len(matrix1[0])):
            matrix[y][x] = matrix1[y][x] - matrix2[y][x]
    return matrix

def residuum(A: list, b: list, x: list) -> list:
    return sub(mul(A, x), b)

def jacobi_method(A: list, b: list,  diverges = False, epsilon = 10**-9) -> (list, int):
    iterations = 0
    norms = []
    n = len(A)
    vector_x = [[0] for _ in range(n)]

    while norm(residuum(A, b, vector_x)) > epsilon:
        if diverges and iterations >= 50:
            break

        iterations += 1
        vector_x_2 = [[0] for _ in range(n)]
        norms.append(norm(residuum(A, b, vector_x)))
        
        for i in range(n):
            sig = sum([0 if i == j else vector_x[j][0]*A[i][j] for j in range(n)])
            vector_x_2[i][0] = (b[i][0] - sig)/A[i][i]
        
        for i in range(len(vector_x)):
            vector_x[i][0] = vector_x_2[i][0]

    norms.append(norm(residuum(A, b, vector_x)))
    return vector_x, iterations, norms

def gauss_seidel_method(A: list, b:list, diverges = False, epsilon = 10**-9) -> (list, int):
    iterations = 0
    norms = []
    n = len(A)
    vector_x = [[0] for _ in range(n)]

    while norm(residuum(A, b, vector_x)) > epsilon:
        if diverges and iterations >= 50:
            break
        
        iterations += 1
        vector_x_2 = [[0] for _ in range(n)]
        norms.append(norm(residuum(A, b, vector_x)))
        
        for i in range(n):
            sig1 = sum([A[i][j] * vector_x_2[j][0] for j in range(i)])
            sig2 = sum([A[i][j] * vector_x[j][0] for j in range(i+1, n)])
            vector_x_2[i][0] = (b[i][0] - sig1 - sig2)/A[i][i]
        
        for i in range(len(vector_x)):
            vector_x[i][0] = vector_x_2[i][0]

    norms.append(norm(residuum(A, b, vector_x)))
    return vector_x, iterations, norms

def solve_linear_system(N = 952, a1 = 7, a2 = -1, a3 = -1):
    A, b = create_matrix_and_vector(N, a1, a2, a3)

    # Metoda Jacobiego
    start_time = time.time()
    _, iterations_jacobi, residuals_jacobi = jacobi_method(A, b)
    end_time = time.time()
    time_jacobi = end_time - start_time

    # Metoda Gaussa-Seidla
    start_time = time.time()
    _, iterations_gauss_seidel, residuals_gauss_seidel = gauss_seidel_method(A, b)
    end_time = time.time()
    time_gauss_seidel = end_time - start_time

    print("Metoda Jacobiego:")
    print("Liczba iteracji:", iterations_jacobi)
    print(f"Czas wykonania: {time_jacobi} s")

    print("\nMetoda Gaussa-Seidla:")
    print("Liczba iteracji:", iterations_gauss_seidel)
    print(f"Czas wykonania: {time_gauss_seidel} s")

    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(residuals_jacobi) + 1), residuals_jacobi, label='Metoda Jacobiego')
    plt.plot(range(1, len(residuals_gauss_seidel) + 1), residuals_gauss_seidel, label='Metoda Gaussa-Seidla')
    plt.yscale('log')
    plt.xlabel('Iteracja')
    plt.ylabel('Norma residuum')
    plt.title('Zmiana normy residuum w kolejnych iteracjach')
    plt.legend()
    plt.grid(True)
    plt.show()

# zadanie C
def argmin_min(arr):
    m = arr[0]
    m_i = 0

    for i in range(len(arr)):
        if arr[i] < m:
            m = arr[i]
            m_i = i
    return m_i, m

def check_convergence():
    A, b = create_matrix_and_vector(952, 3, -1, -1)

    _, _, residuals_jacobi = jacobi_method(A, b, diverges = True)
    _, _, residuals_gauss_seidel = gauss_seidel_method(A, b, diverges = True)

    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(residuals_jacobi) + 1), residuals_jacobi, label='Metoda Jacobiego')
    plt.plot(range(1, len(residuals_gauss_seidel) + 1), residuals_gauss_seidel, label='Metoda Gaussa-Seidla')
    plt.yscale('log')
    plt.xlabel('Iteracja')
    plt.ylabel('Norma residuum')
    plt.title('Zmiana normy residuum w kolejnych iteracjach')
    plt.legend()
    plt.grid(True)
    plt.show()

# zadanie D
def identify_matrix(n):
    return [[0 if x != y else 1 for x in range(n)] for y in range(n)]

def lu_decomposition(matrix: list) -> (list, list):
    m = len(matrix)
    L = identify_matrix(m)
    U = deepcopy(matrix)

    for k in range(0, m-1):
        for j in range(k+1, m):
            L[j][k] = U[j][k]/U[k][k]
            for i in range(k, m):
                U[j][i] -= L[j][k] * U[k][i]
    return L, U

def solve_upper_triangle(matrix: list, b: list) -> list:
    n = len(matrix)
    vector_x = [[0] for _ in range(n)]
    vector_x[-1][0] = b[n-1][0]/matrix[n-1][n-1]
    
    for i in range(n-2, -1, -1):
        vector_x[i][0] = (b[i][0] - sum([matrix[i][j] * vector_x[j][0] for j in range(n-1, i, -1)]))/matrix[i][i]

    return vector_x

def solve_lower_triangle(matrix: list, b: list) -> list:
    n = len(matrix)
    vector_x = [[0] for _ in range(n)]
    vector_x[0][0] = b[0][0]/matrix[0][0]
    
    for i in range(1, n):
        vector_x[i][0] = (b[i][0] - sum([matrix[i][j] * vector_x[j][0] for j in range(0, i)]))/matrix[i][i]

    return vector_x

def lu(A: list, b: list) -> list:
    L, U = lu_decomposition(A)
    y = solve_lower_triangle(L, b)
    return solve_upper_triangle(U, y)

def solve_linear_lu():
    A, b = create_matrix_and_vector(952, 3, -1, -1)
    gauss_x = lu(A, b)
    print(f"LU residuum = {norm(residuum(A, b, gauss_x))}")

# zadanie E
def time_comparison():
    N_values = [100, 500, 1000, 1500, 2000]
    jacobi_times = []
    gauss_seidel_times = []
    lu_times = []

    for N in N_values:
        A, b = create_matrix_and_vector(N)

        # Metoda Jacobiego
        start_time = time.time()
        jacobi_method(A, b)
        end_time = time.time()
        jacobi_time = end_time - start_time
        jacobi_times.append(jacobi_time)

        # Metoda Gaussa-Seidla
        start_time = time.time()
        gauss_seidel_method(A, b)
        end_time = time.time()
        gauss_seidel_time = end_time - start_time
        gauss_seidel_times.append(gauss_seidel_time)

        # Metoda faktoryzacji LU
        start_time = time.time()
        lu(A, b)
        end_time = time.time()
        lu_time = end_time - start_time
        lu_times.append(lu_time)
        print(f"N={N}, LU = {lu_time} s, Jacobi = {jacobi_time} s, Gauss-Seidel = {gauss_seidel_time} s")

    plt.plot(N_values, jacobi_times, label='Metoda Jacobiego')
    plt.plot(N_values, gauss_seidel_times, label='Metoda Gaussa-Seidla')
    plt.plot(N_values, lu_times, label='Metoda faktoryzacji LU')
    plt.xlabel('Rozmiar macierzy (N)')
    plt.ylabel('Czas (s)')
    plt.title('Czas wyznaczenia rozwiązania w zależności od rozmiaru macierzy')
    plt.legend()
    plt.grid(True)
    plt.show()

    return jacobi_times, gauss_seidel_times, lu_times

#print_matrix_and_vector(*create_matrix_and_vector())
solve_linear_system()
check_convergence()
solve_linear_lu()
time_comparison()