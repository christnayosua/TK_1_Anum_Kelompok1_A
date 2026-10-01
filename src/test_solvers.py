"""
Program: test_solvers.py
Deskripsi: Unit Test & Validasi Presisi Tinggi untuk Solver Nomor 1 & Nomor 2
           Tugas Kelompok 1 Analisis Numerik - Fasilkom UI
"""

import numpy as np
import os
import sys

# Tambahkan path src
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from no1_solver import (
    load_transition_matrix,
    validate_transition_matrix,
    detect_bandwidth,
    construct_B_and_b,
    solve_dense_lu_pp,
    solve_banded_thomas_pp,
    normalize_solution,
    compute_errors
)
from no2_solver import (
    load_and_preprocess_stock,
    solve_normal_equations,
    demonstrate_gaussian_elimination_inconsistency,
    solve_iterative_refinement,
    householder_reflection_step1,
    solve_householder_qr,
    evaluate_out_of_sample
)

def test_nomor_1():
    print("==================================================")
    print("TEST SUITE: NOMOR 1 (DENSE LU VS BANDED THOMAS PP)")
    print("==================================================")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidates_no1 = [
        os.path.join(base_dir, "Nomor 1", "A"),
        os.path.join(base_dir, "Nomor 1", "Nomor 1", "A")
    ]
    data_dir = next((p for p in candidates_no1 if os.path.exists(os.path.join(p, "T_16.csv"))), candidates_no1[0])
    
    sizes = [16, 32, 64, 128, 256, 512]
    
    for N in sizes:
        csv_path = os.path.join(data_dir, f"T_{N}.csv")
        assert os.path.exists(csv_path), f"File {csv_path} not found!"
        
        T = load_transition_matrix(csv_path)
        
        # 1. Validasi Matriks Stokastik
        val = validate_transition_matrix(T)
        assert val['is_square'], f"T_{N} bukan bujursangkar"
        assert val['is_non_negative'], f"T_{N} memiliki elemen negatif"
        assert val['is_stochastic'], f"Jumlah baris T_{N} tidak sama dengan 1"
        assert val['is_irreducible'], f"Rantai Markov T_{N} tereduksi"
        assert val['is_unique_stationary'], f"Distribusi stasioner T_{N} tidak tunggal"
        
        # 2. Deteksi Bandwidth
        p_T, q_T = detect_bandwidth(T)
        assert p_T == 2 and q_T == 1, f"Bandwidth T_{N} tidak sesuai: p={p_T}, q={q_T}"
        
        B, b = construct_B_and_b(T)
        p_B, q_B = detect_bandwidth(B)
        assert p_B == 1 and q_B == 2, f"Bandwidth B_{N} tidak sesuai: p={p_B}, q={q_B}"
        
        # 3. Solver Eksekusi
        z_dense = solve_dense_lu_pp(B, b)
        pi_dense = normalize_solution(z_dense)
        
        z_banded = solve_banded_thomas_pp(B, b)
        pi_banded = normalize_solution(z_banded)
        
        # 4. Asersi Galat
        diff = np.max(np.abs(pi_dense - pi_banded))
        r_dense, e_dense = compute_errors(T, pi_dense)
        r_banded, e_banded = compute_errors(T, pi_banded)
        
        print(f"[TEST N={N:3d}] Diff: {diff:.2e} | Res: {r_banded:.2e} | Enorm: {e_banded:.2e} | Band(B): ({p_B},{q_B})")
        assert diff < 1e-12, f"Solusi berbeda pada N={N}: diff={diff}"
        assert r_banded < 1e-12, f"Residual tinggi pada N={N}: {r_banded}"
        assert e_banded < 1e-12, f"Enorm tinggi pada N={N}: {e_banded}"
        assert np.all(pi_banded >= -1e-15), f"Elemen negatif ditemukan pada N={N}"
        
    print(">>> SEMUA TEST NOMOR 1 PASSED! <<<\n")

def test_nomor_2():
    print("==================================================")
    print("TEST SUITE: NOMOR 2 (NORMAL EQ, QR & OUT-OF-SAMPLE)")
    print("==================================================")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidates_train = [
        os.path.join(base_dir, "Nomor 2", "stock_train.csv"),
        os.path.join(base_dir, "Nomor 2", "Nomor 2", "stock_train.csv")
    ]
    train_csv = next((p for p in candidates_train if os.path.exists(p)), candidates_train[0])
    
    candidates_test = [
        os.path.join(base_dir, "Nomor 2", "stock_test.csv"),
        os.path.join(base_dir, "Nomor 2", "Nomor 2", "stock_test.csv")
    ]
    test_csv = next((p for p in candidates_test if os.path.exists(p)), candidates_test[0])
    assert os.path.exists(train_csv), f"File {train_csv} not found!"
    assert os.path.exists(test_csv), f"File {test_csv} not found!"
    
    A, b, returns, dates = load_and_preprocess_stock(train_csv)
    assert A.shape == (300, 6), f"Bentuk A salah: {A.shape}"
    assert b.shape == (300,), f"Bentuk b salah: {b.shape}"
    
    # 0. Test Gaussian Elimination Inconsistency
    res_incons = demonstrate_gaussian_elimination_inconsistency(A, b)
    assert not res_incons['is_consistent'], "Sistem overdetermined seharusnya tidak konsisten!"
    print(f"[TEST N2 Inconsistency] Rank={res_incons['rank']}, Max Inconsistent Res={res_incons['max_inconsistent_residual']:.4f}")

    # 1. Test Normal Equations
    res_norm = solve_normal_equations(A, b)
    assert res_norm['cond_A'] > 100, "Condition number A tidak valid"
    assert abs(res_norm['cond_ATA'] - (res_norm['cond_A']**2)) / res_norm['cond_ATA'] < 1e-3, "cond(ATA) bukan kuadrat cond(A)"
    print(f"[TEST N2 Normal] cond(A)={res_norm['cond_A']:.2f}, cond(ATA)={res_norm['cond_ATA']:.2f}")
    
    # 1b. Test Iterative Error Refinement
    res_refine = solve_iterative_refinement(A, b, max_iter=2)
    delta_norm_1 = res_refine['history'][1]['correction_norm']
    print(f"[TEST N2 Iterative Refinement] Correction Norm Iter 1: {delta_norm_1:.2e}")
    assert delta_norm_1 < 1e-14, f"Koreksi refinement seharusnya mendekati presisi mesin: {delta_norm_1}"

    # 2. Test H1 Reflection
    res_h1 = householder_reflection_step1(A)
    print(f"[TEST N2 H1] Max subdiagonal error: {res_h1['max_subdiag_error']:.2e}")
    assert res_h1['max_subdiag_error'] < 1e-12, "Subdiagonal H1 A tidak tereliminasi menjadi nol!"
    
    # 3. Test QR vs Normal
    res_qr = solve_householder_qr(A, b)
    diff = np.linalg.norm(res_norm['x'] - res_qr['x_ls'])
    print(f"[TEST N2 QR vs Normal] Diff: {diff:.2e} | Res QR: {res_qr['residual_norm']:.6f}")
    assert diff < 1e-12, f"Solusi QR dan Normal berbeda: {diff}"
    assert abs(res_norm['residual_norm'] - res_qr['residual_norm']) < 1e-12, "Norm residual berbeda!"
    
    # 4. Test Out-of-Sample Evaluation
    eval_res = evaluate_out_of_sample(train_csv, test_csv, res_qr['x_ls'])
    print(f"[TEST N2 Out-of-Sample] RMSE Train: {eval_res['rmse_train']:.6f} | RMSE Test Ind: {eval_res['rmse_test_independent']:.6f}")
    assert eval_res['rmse_train'] < 0.015, "RMSE Train tidak wajar"
    assert eval_res['rmse_test_independent'] < 0.02, "RMSE Test tidak wajar"
    
    print(">>> SEMUA TEST NOMOR 2 PASSED! <<<\n")

if __name__ == "__main__":
    test_nomor_1()
    test_nomor_2()
    print("ALL UNIT TESTS PASSED SUCCESSFULLY!")
