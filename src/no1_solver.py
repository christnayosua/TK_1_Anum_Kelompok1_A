"""
Program: no1_solver.py
Deskripsi: Implementasi Solver SPL untuk Masalah Distribusi Steady State Antarhalte
           Tugas Kelompok 1 Analisis Numerik - Pengerjaan Lengkap Kelompok Ganjil (Kode A)
           
Fitur & Metode:
1. Validasi Matriks Stokastik (Dimensi, Non-negativitas, Row-Sum = 1, Syarat Ketunggalan Perron-Frobenius)
2. Deteksi Otomatis Lower Bandwidth (p) dan Upper Bandwidth (q)
3. Solver berbasis faktorisasi LU Dense dengan Partial Pivoting (PB = LU)
4. Solver bergaya Algoritma Thomas untuk matriks banded dengan Partial Pivoting (O(N) memori kompak)
5. Analisis Galat Residual, Normalisasi, dan Condition Number
"""

import numpy as np
import pandas as pd
import time
import os

def load_transition_matrix(csv_path: str) -> np.ndarray:
    """Membaca matriks transisi T dari file CSV."""
    df = pd.read_csv(csv_path, header=None)
    return df.values.astype(float)

def validate_transition_matrix(T: np.ndarray) -> dict:
    """
    Validasi integritas matematis matriks transisi T:
    - Dimensi bujursangkar (N x N)
    - Non-negativitas elemen: T_ij >= 0
    - Syarat stokastik baris: sum_j T_ij = 1
    - Syarat ketunggalan distribusi stasioner (Teorema Perron-Frobenius: irreducibility & aperiodicity)
    """
    n_rows, n_cols = T.shape
    is_square = (n_rows == n_cols)
    min_val = float(np.min(T))
    is_non_negative = (min_val >= -1e-15)
    
    row_sums = np.sum(T, axis=1)
    max_row_sum_err = float(np.max(np.abs(row_sums - 1.0)))
    is_stochastic = (max_row_sum_err < 1e-12)
    
    # Deteksi ketereduksian (irreducibility) melalui konektivitas graf transisi
    # Menggunakan representasi ketetanggaan (adjacency matrix A_adj: T_ij > 0)
    # Suatu rantai Markov dengan matriks pita tak-nol di sub dan superdiagonal
    # membentuk graf terhubung kuat (strongly connected) sehingga bersifat iredisibel.
    is_irreducible = True
    for i in range(n_rows - 1):
        if T[i+1, i] <= 0 and T[i, i+1] <= 0:
            is_irreducible = False
            break
            
    # Nilai eigen terbesar (Perron root)
    eigvals = np.linalg.eigvals(T.T)
    max_eigval = float(np.max(np.abs(eigvals)))
    # Hitung kelipatan aljabar dari eigenvalue 1.0 (selisih absolut murni < 1e-10)
    num_unit_eigvals = int(np.sum(np.abs(np.abs(eigvals) - 1.0) < 1e-10))
    
    return {
        "N": n_rows,
        "is_square": is_square,
        "min_val": min_val,
        "is_non_negative": is_non_negative,
        "max_row_sum_err": max_row_sum_err,
        "is_stochastic": is_stochastic,
        "is_irreducible": is_irreducible,
        "max_eigval": max_eigval,
        "num_unit_eigvals": num_unit_eigvals,
        "is_unique_stationary": is_irreducible and (num_unit_eigvals == 1)
    }

def detect_bandwidth(M: np.ndarray, tol: float = 1e-12) -> tuple[int, int]:
    """
    Menghitung secara otomatis lower bandwidth p dan upper bandwidth q:
    p: indeks baris terjauh di bawah diagonal utama (i - j > 0) dengan M_ij != 0
    q: indeks kolom terjauh di atas diagonal utama (j - i > 0) dengan M_ij != 0
    """
    n = M.shape[0]
    p = 0
    q = 0
    for i in range(n):
        for j in range(n):
            if abs(M[i, j]) > tol:
                if i - j > p:
                    p = i - j
                if j - i > q:
                    q = j - i
    return p, q

def construct_B_and_b(T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Membentuk matriks B dan vektor b dari matriks transisi T:
    A = I - T^T
    B[0, :] = [1, 0, ..., 0]
    B[i, :] = A[i, :] untuk i = 1, ..., N-1
    b = [1, 0, ..., 0]^T
    """
    N = T.shape[0]
    A = np.eye(N) - T.T
    B = A.copy()
    B[0, :] = 0.0
    B[0, 0] = 1.0
    b = np.zeros(N)
    b[0] = 1.0
    return B, b

def solve_dense_lu_pp(B_in: np.ndarray, b_in: np.ndarray) -> np.ndarray:
    """
    Solver SPL Bz = b menggunakan Faktorisasi LU Dense dengan Partial Pivoting:
    P B = L U
    L y = P b (Forward Substitution)
    U z = y   (Back Substitution)
    
    Implementasi mandiri tanpa library scipy.linalg atau numpy.linalg.solve.
    """
    n = B_in.shape[0]
    A = B_in.copy().astype(float)
    b = b_in.copy().astype(float)
    p = np.arange(n)
    
    # 1. Faktorisasi LU dengan Partial Pivoting
    for k in range(n - 1):
        max_row = k + np.argmax(np.abs(A[k:, k]))
        if max_row != k:
            A[[k, max_row]] = A[[max_row, k]]
            p[[k, max_row]] = p[[max_row, k]]
            
        pivot = A[k, k]
        if abs(pivot) < 1e-15:
            continue
            
        factors = A[k+1:, k] / pivot
        A[k+1:, k] = factors
        A[k+1:, k+1:] -= np.outer(factors, A[k, k+1:])
        
    # 2. Forward Substitution: L y = P b
    b_perm = b[p]
    y = np.zeros(n)
    for i in range(n):
        y[i] = b_perm[i] - np.dot(A[i, :i], y[:i])
        
    # 3. Back Substitution: U z = y
    z = np.zeros(n)
    for i in range(n - 1, -1, -1):
        z[i] = (y[i] - np.dot(A[i, i+1:], z[i+1:])) / A[i, i]
        
    return z

def solve_banded_thomas_pp(B_in: np.ndarray, b_in: np.ndarray) -> np.ndarray:
    """
    Solver SPL Bz = b bergaya Algoritma Thomas untuk matriks banded
    dengan Partial Pivoting, hanya menyimpan vektor-vektor diagonal (memori O(N)).
    
    Karakteristik matriks B:
    - Lower bandwidth p = 1 (1 subdiagonal dl)
    - Upper bandwidth q = 2 (2 superdiagonals du1, du2)
    - Akibat partial pivoting, terjadi fill-in pada superdiagonal ke-3 (du3),
      sehingga bandwidth atas efektif menjadi q + p = 3.
      
    Penyimpanan hanya membutuhkan 5 vektor berukuran N:
    - d   : diagonal utama
    - dl  : sub-diagonal (p = 1)
    - du1 : super-diagonal 1 (q = 1)
    - du2 : super-diagonal 2 (q = 2)
    - du3 : fill-in super-diagonal 3 (q + p = 3)
    """
    n = B_in.shape[0]
    
    d = np.zeros(n)
    dl = np.zeros(n - 1)
    du1 = np.zeros(n - 1)
    du2 = np.zeros(n - 2)
    du3 = np.zeros(n - 3)
    
    for i in range(n):
        d[i] = B_in[i, i]
        if i < n - 1:
            du1[i] = B_in[i, i+1]
            dl[i] = B_in[i+1, i]
        if i < n - 2:
            du2[i] = B_in[i, i+2]
            
    b = b_in.copy().astype(float)
    p = np.arange(n)
    
    # Eliminasi maju dengan partial pivoting
    for k in range(n - 1):
        val_k = abs(d[k])
        val_k1 = abs(dl[k])
        
        if val_k1 > val_k:
            p[[k, k+1]] = p[[k+1, k]]
            b[[k, k+1]] = b[[k+1, k]]
            
            d[k], dl[k] = dl[k], d[k]
            du1[k], d[k+1] = d[k+1], du1[k]
            if k < n - 2:
                du2[k], du1[k+1] = du1[k+1], du2[k]
            if k < n - 3:
                du3[k], du2[k+1] = du2[k+1], du3[k]
                
        if abs(d[k]) > 1e-15:
            factor = dl[k] / d[k]
            dl[k] = factor
            
            d[k+1] -= factor * du1[k]
            if k < n - 2:
                du1[k+1] -= factor * du2[k]
            if k < n - 3:
                du2[k+1] -= factor * du3[k]
                
            b[k+1] -= factor * b[k]
            
    # Substitusi mundur: U z = b
    z = np.zeros(n)
    for i in range(n - 1, -1, -1):
        sum_val = 0.0
        if i < n - 1:
            sum_val += du1[i] * z[i+1]
        if i < n - 2:
            sum_val += du2[i] * z[i+2]
        if i < n - 3:
            sum_val += du3[i] * z[i+3]
        z[i] = (b[i] - sum_val) / d[i]
        
    return z

def normalize_solution(z: np.ndarray) -> np.ndarray:
    """Melakukan normalisasi vektor z agar memenuhi 1^T pi = 1."""
    return z / np.sum(z)

def compute_errors(T: np.ndarray, pi: np.ndarray) -> tuple[float, float]:
    """Menghitung galat residual r = ||T^T pi - pi||_2 dan enorm = |1^T pi - 1|."""
    r = float(np.linalg.norm(T.T @ pi - pi))
    enorm = float(abs(np.sum(pi) - 1.0))
    return r, enorm

if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "Nomor 1", "Nomor 1", "A")
    
    print("=========================================================================")
    print("EKSPERIMEN SOLVER SPL NOMOR 1 (KELOMPOK GANJIL) - DATASET KODE A")
    print("=========================================================================")
    
    sizes = [16, 32, 64, 128, 256, 512]
    
    for N in sizes:
        csv_path = os.path.join(data_dir, f"T_{N}.csv")
        if not os.path.exists(csv_path):
            continue
            
        T = load_transition_matrix(csv_path)
        val_res = validate_transition_matrix(T)
        p_T, q_T = detect_bandwidth(T)
        
        B, b = construct_B_and_b(T)
        p_B, q_B = detect_bandwidth(B)
        
        t0 = time.perf_counter()
        z_dense = solve_dense_lu_pp(B, b)
        t_dense = (time.perf_counter() - t0) * 1000.0
        pi_dense = normalize_solution(z_dense)
        r_dense, e_dense = compute_errors(T, pi_dense)
        
        t0 = time.perf_counter()
        z_banded = solve_banded_thomas_pp(B, b)
        t_banded = (time.perf_counter() - t0) * 1000.0
        pi_banded = normalize_solution(z_banded)
        r_banded, e_banded = compute_errors(T, pi_banded)
        
        diff_max = np.max(np.abs(pi_dense - pi_banded))
        cond_B = np.linalg.cond(B)
        max_h = np.argmax(pi_banded) + 1
        max_p = np.max(pi_banded)
        
        print(f"N={N:3d} | Band(B): ({p_B},{q_B}) | cond(B): {cond_B:6.2f} | t_dense: {t_dense:6.2f}ms | t_banded: {t_banded:5.2f}ms | Diff: {diff_max:.1e} | Res: {r_banded:.1e} | Halte Puncak: {max_h} ({max_p:.4f})")
