"""
Program: generate_visualizations.py
Deskripsi: Skrip otomasi pembuatan visualisasi grafik untuk Technical Report Tugas Kelompok 1
           Analisis Numerik - Fasilkom UI
           
Grafik yang dihasilkan:
1. figures/fig1_distribusi_halte.png : Distribusi probabilitas stasioner pi_i terhadap halte
2. figures/fig2_performa_solver.png  : Perbandingan waktu komputasi Dense LU vs Banded Thomas (skala log-log)
3. figures/fig3_setar_train_overlay.png : Overlay deret return aktual vs estimasi model SETAR (Train)
4. figures/fig4_setar_continuous_overlay.png : Continuous time-series plot (Train + Test) dengan garis batas
"""

import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Konfigurasi gaya matplotlib akademik
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['grid.color'] = '#cccccc'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.alpha'] = 0.5

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fig_dir = os.path.join(base_dir, 'figures')
os.makedirs(fig_dir, exist_ok=True)

# Import solver mandiri
from no1_solver import (
    load_transition_matrix,
    construct_B_and_b,
    solve_dense_lu_pp,
    solve_banded_thomas_pp,
    normalize_solution
)
from no2_solver import (
    load_and_preprocess_stock,
    solve_householder_qr
)

# =========================================================================
# 1. GRAFIK DISTRIBUSI STASIONER HALTE (FIG 1)
# =========================================================================
def generate_fig1():
    print("Membuat Gambar 1: Distribusi Probabilitas Stasioner Halte...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=300)
    
    # Subplot A: Kasus N = 16 (Detail Halte 1 sampai 16)
    csv_16 = os.path.join(base_dir, 'Nomor 1', 'Nomor 1', 'A', 'T_16.csv')
    T16 = load_transition_matrix(csv_16)
    B16, b16 = construct_B_and_b(T16)
    z16 = solve_banded_thomas_pp(B16, b16)
    pi16 = normalize_solution(z16)
    
    halte_idx16 = np.arange(1, 17)
    max_idx16 = np.argmax(pi16)
    
    colors16 = ['#1f77b4' if i != max_idx16 else '#d62728' for i in range(16)]
    bars = ax1.bar(halte_idx16, pi16, color=colors16, alpha=0.85, edgecolor='black', linewidth=0.7)
    ax1.plot(halte_idx16, pi16, color='#2c3e50', linestyle='-', marker='o', markersize=4, alpha=0.7)
    
    ax1.annotate(f'Puncak: Halte {max_idx16+1}\n($\\pi = {pi16[max_idx16]:.4f}$)',
                 xy=(max_idx16+1, pi16[max_idx16]),
                 xytext=(max_idx16+2.5, pi16[max_idx16] + 0.003),
                 arrowprops=dict(facecolor='black', shrink=0.08, width=1, headwidth=6),
                 fontsize=9, fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.3", fc="#ffeaa7", ec="#b2bec3", lw=0.8))
    
    ax1.set_title('(a) Distribusi Probabilitas Stasioner Koridor ($N=16$ Halte)', fontsize=11, fontweight='bold', pad=10)
    ax1.set_xlabel('Nomor Halte ($i$)', fontsize=10)
    ax1.set_ylabel('Peluang Stasioner ($\\pi_i$)', fontsize=10)
    ax1.set_xticks(halte_idx16)
    ax1.grid(True)
    ax1.set_ylim(0, np.max(pi16) * 1.18)
    
    # Subplot B: Perbandingan Multi-Skala N in {16, 32, 64, 128, 256, 512} dengan Posisi Ternormalisasi
    palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    sizes = [16, 32, 64, 128, 256, 512]
    
    for idx, N in enumerate(sizes):
        csv_path = os.path.join(base_dir, 'Nomor 1', 'Nomor 1', 'A', f'T_{N}.csv')
        T = load_transition_matrix(csv_path)
        B, b = construct_B_and_b(T)
        z = solve_banded_thomas_pp(B, b)
        pi = normalize_solution(z)
        norm_pos = np.linspace(0, 1, N)
        # Tampilkan densitas proporsional (N * pi_i) agar kurva skala seimbang
        ax2.plot(norm_pos, N * pi, label=f'$N={N}$ (Halte Puncak: {np.argmax(pi)+1})',
                 color=palette[idx], linewidth=1.5 if N in [16, 512] else 1.0)
                 
    ax2.set_title('(b) Kerapatan Relatif Penumpang Sepanjang Koridor (Multi-Skala $N$)', fontsize=11, fontweight='bold', pad=10)
    ax2.set_xlabel('Posisi Koridor Ternormalisasi ($i/N$)', fontsize=10)
    ax2.set_ylabel('Kerapatan Relatif ($N \\times \\pi_i$)', fontsize=10)
    ax2.legend(fontsize=8, loc='upper right', framealpha=0.9)
    ax2.grid(True)
    
    plt.tight_layout()
    out_path = os.path.join(fig_dir, 'fig1_distribusi_halte.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Gambar 1 tersimpan di: {out_path}")

# =========================================================================
# 2. GRAFIK KINERJA KOMPUTASI DENSE VS BANDED (FIG 2)
# =========================================================================
def generate_fig2():
    print("Membuat Gambar 2: Perbandingan Kinerja Solver SPL...")
    sizes = [16, 32, 64, 128, 256, 512]
    t_dense_list = []
    t_banded_list = []
    
    for N in sizes:
        csv_path = os.path.join(base_dir, 'Nomor 1', 'Nomor 1', 'A', f'T_{N}.csv')
        T = load_transition_matrix(csv_path)
        B, b = construct_B_and_b(T)
        
        # Dense LU PP timing (multi-run average for precision)
        runs_dense = 5 if N <= 128 else 2
        t0 = time.perf_counter()
        for _ in range(runs_dense):
            _ = solve_dense_lu_pp(B, b)
        t_dense = (time.perf_counter() - t0) / runs_dense * 1000.0
        t_dense_list.append(t_dense)
        
        # Banded Thomas PP timing
        runs_banded = 20
        t0 = time.perf_counter()
        for _ in range(runs_banded):
            _ = solve_banded_thomas_pp(B, b)
        t_banded = (time.perf_counter() - t0) / runs_banded * 1000.0
        t_banded_list.append(t_banded)
        
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=300)
    
    # Subplot A: Skala Log-Log Waktu Eksekusi vs N
    ax1.plot(sizes, t_dense_list, 'o-', color='#d62728', linewidth=2, label='Faktorisasi Dense LU dengan PP ($O(N^3)$)')
    ax1.plot(sizes, t_banded_list, 's-', color='#2ca02c', linewidth=2, label='Banded Thomas-Style dengan PP ($O(N)$)')
    
    # Garis acuan teoretis
    c_dense = t_dense_list[0] / (sizes[0]**3)
    c_banded = t_banded_list[0] / sizes[0]
    ref_N = np.array(sizes)
    ax1.plot(ref_N, c_dense * (ref_N**3), '--', color='#e74c3c', alpha=0.5, label='Tren Teoretis $O(N^3)$')
    ax1.plot(ref_N, c_banded * ref_N, '--', color='#27ae60', alpha=0.5, label='Tren Teoretis $O(N)$')
    
    ax1.set_xscale('log', base=2)
    ax1.set_yscale('log')
    ax1.set_title('(a) Waktu Komputasi vs Ukuran Matriks $N$ (Skala Log-Log)', fontsize=11, fontweight='bold', pad=10)
    ax1.set_xlabel('Ukuran Matriks ($N$)', fontsize=10)
    ax1.set_ylabel('Waktu Eksekusi (milidetik)', fontsize=10)
    ax1.set_xticks(sizes)
    ax1.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax1.legend(fontsize=8.5, loc='upper left')
    ax1.grid(True, which="both", ls="--")
    
    # Subplot B: Rasio Percepatan (Speedup Factor: Dense / Banded)
    speedups = [d / b for d, b in zip(t_dense_list, t_banded_list)]
    bars = ax2.bar([str(s) for s in sizes], speedups, color='#3498db', edgecolor='#2980b9', width=0.55)
    for bar in bars:
        h = bar.get_height()
        ax2.annotate(f'{h:.1f}x',
                     xy=(bar.get_x() + bar.get_width() / 2, h),
                     xytext=(0, 4), textcoords="offset points",
                     ha='center', va='bottom', fontsize=9, fontweight='bold')
                     
    ax2.set_title('(b) Faktor Percepatan (*Speedup Factor*) Banded atas Dense', fontsize=11, fontweight='bold', pad=10)
    ax2.set_xlabel('Ukuran Matriks ($N$)', fontsize=10)
    ax2.set_ylabel('Rasio Percepatan ($t_{dense} / t_{banded}$)', fontsize=10)
    ax2.set_ylim(0, max(speedups) * 1.18)
    ax2.grid(True, axis='y')
    
    plt.tight_layout()
    out_path = os.path.join(fig_dir, 'fig2_performa_solver.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Gambar 2 tersimpan di: {out_path}")

# =========================================================================
# 3. GRAFIK OVERLAY RETURN SETAR DATA LATIH (FIG 3)
# =========================================================================
def generate_fig3():
    print("Membuat Gambar 3: Overlay Return Model SETAR (Data Latih)...")
    train_csv = os.path.join(base_dir, 'Nomor 2', 'Nomor 2', 'stock_train.csv')
    A, b, returns, dates = load_and_preprocess_stock(train_csv)
    
    res_qr = solve_householder_qr(A, b)
    x_ls = res_qr['x_ls']
    b_pred = A @ x_ls
    residuals = b - b_pred
    
    # Tanggal untuk baris 0 s.d. m-1 (t=3 s.d. 302)
    obs_dates = dates[2:]
    t_idx = np.arange(1, len(b) + 1)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7), dpi=300, sharex=True,
                                   gridspec_kw={'height_ratios': [2.5, 1]})
    
    # Panel Atas: Aktual vs Prediksi
    ax1.plot(t_idx, b * 100, color='#34495e', linewidth=1.2, alpha=0.8, label='Return Aktual ($R_t$)')
    ax1.plot(t_idx, b_pred * 100, color='#e74c3c', linewidth=1.4, linestyle='-', label='Estimasi SETAR 2-Rezim ($\\hat{R}_t$)')
    
    # Shading rezim: Bullish vs Bearish
    # Rezim ditentukan oleh R_{t-1} >= 0
    r_prev = returns[1:-1]  # R_{t-1} untuk t = 2 ... 301
    for i in range(len(b)):
        if r_prev[i] >= 0:
            ax1.axvspan(i + 0.5, i + 1.5, color='#2ecc71', alpha=0.08, lw=0)
        else:
            ax1.axvspan(i + 0.5, i + 1.5, color='#e67e22', alpha=0.08, lw=0)
            
    # Tampilkan box legenda khusus rezim
    ax1.plot([], [], color='#2ecc71', alpha=0.4, linewidth=8, label='Fase Bullish ($R_{t-1} \\geq 0$)')
    ax1.plot([], [], color='#e67e22', alpha=0.4, linewidth=8, label='Fase Bearish ($R_{t-1} < 0$)')
    
    rmse_tr = np.sqrt(np.mean(residuals**2))
    ax1.set_title(f'Prediksi Deret Waktu Model SETAR vs Return Aktual (Data Latih: 300 Observasi, RMSE = {rmse_tr*100:.3f}%)',
                  fontsize=11, fontweight='bold', pad=10)
    ax1.set_ylabel('Return Harian (%)', fontsize=10)
    ax1.legend(loc='upper right', fontsize=8.5, framealpha=0.9)
    ax1.grid(True)
    
    # Panel Bawah: Galat Residual e_t = R_t - R_hat_t
    ax2.plot(t_idx, residuals * 100, color='#7f8c8d', linewidth=0.8, label='Residual ($e_t = R_t - \\hat{R}_t$)')
    ax2.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax2.set_xlabel('Indeks Waktu Observasi ($t$)', fontsize=10)
    ax2.set_ylabel('Residual (%)', fontsize=10)
    ax2.set_ylim(-max(np.abs(residuals))*120, max(np.abs(residuals))*120)
    ax2.grid(True)
    ax2.legend(loc='upper right', fontsize=8)
    
    plt.tight_layout()
    out_path = os.path.join(fig_dir, 'fig3_setar_train_overlay.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Gambar 3 tersimpan di: {out_path}")

# =========================================================================
# 4. GRAFIK KONTINU GABUNGAN TRAIN + TEST (FIG 4)
# =========================================================================
def generate_fig4():
    print("Membuat Gambar 4: Grafik Kontinu Deret Waktu (Train vs Test Overlay)...")
    train_csv = os.path.join(base_dir, 'Nomor 2', 'Nomor 2', 'stock_train.csv')
    test_csv = os.path.join(base_dir, 'Nomor 2', 'Nomor 2', 'stock_test.csv')
    
    df_tr = pd.read_csv(train_csv)
    df_te = pd.read_csv(test_csv)
    
    # Solusi x_ls dari Train
    A_tr, b_tr, r_tr, d_tr = load_and_preprocess_stock(train_csv)
    res_qr = solve_householder_qr(A_tr, b_tr)
    x_ls = res_qr['x_ls']
    
    pred_tr = A_tr @ x_ls
    
    # Penggabungan kontinu deret harga
    p_all = np.concatenate([df_tr['Close'].values, df_te['Close'].values])
    r_all = (p_all[1:] - p_all[:-1]) / p_all[:-1]  # 405 return
    
    # Data test kontinu: dimulai dari indeks len(r_tr)
    start_test_idx = len(r_tr)
    b_te_cont = []
    A_te_cont = []
    
    for i in range(start_test_idx, len(r_all)):
        rt = r_all[i]
        rt1 = r_all[i-1]
        rt2 = r_all[i-2]
        b_te_cont.append(rt)
        if rt1 >= 0:
            A_te_cont.append([1.0, rt1, rt2, 0.0, 0.0, 0.0])
        else:
            A_te_cont.append([0.0, 0.0, 0.0, 1.0, rt1, rt2])
            
    A_te_cont = np.array(A_te_cont)
    b_te_cont = np.array(b_te_cont)
    pred_te_cont = A_te_cont @ x_ls
    
    # Gabungkan sumbu x dan kurva kontinu
    # Train: t = 1 ... 300
    # Test : t = 301 ... 403
    t_train = np.arange(1, len(b_tr) + 1)
    t_test = np.arange(len(b_tr) + 1, len(b_tr) + 1 + len(b_te_cont))
    
    fig, ax = plt.subplots(figsize=(15, 6), dpi=300)
    
    # Plot Aktual
    ax.plot(t_train, b_tr * 100, color='#2c3e50', linewidth=1.1, alpha=0.75, label='Return Aktual (Train)')
    ax.plot(t_test, b_te_cont * 100, color='#16a085', linewidth=1.1, alpha=0.85, label='Return Aktual (Test / Out-of-Sample)')
    
    # Plot Estimasi SETAR
    ax.plot(t_train, pred_tr * 100, color='#e74c3c', linewidth=1.3, label='Estimasi SETAR 2-Rezim (Train Fit)')
    ax.plot(t_test, pred_te_cont * 100, color='#8e44ad', linewidth=1.4, linestyle='-', label='Proyeksi SETAR 2-Rezim (Test Forecast)')
    
    # Garis pembatas Train - Test
    boundary_t = len(b_tr) + 0.5
    ax.axvline(x=boundary_t, color='#c0392b', linestyle='--', linewidth=2.0, zorder=5)
    
    # Shading latar belakang segmen
    ax.axvspan(1, boundary_t, color='#f5f6fa', alpha=0.6, label='Segmen Data Latih ($N=300$)')
    ax.axvspan(boundary_t, len(b_tr) + len(b_te_cont), color='#e8f8f5', alpha=0.6, label='Segmen Data Uji ($N=103$)')
    
    # Anotasi batas
    ax.text(boundary_t - 4, np.max(b_tr)*85, 'BATAS DATA LATIH\n(Hingga 30 Okt 2024)',
            ha='right', va='top', fontsize=9, fontweight='bold', color='#c0392b',
            bbox=dict(boxstyle="round,pad=0.3", fc="#fdfefe", ec="#c0392b", lw=1))
    ax.text(boundary_t + 4, np.max(b_tr)*85, 'SEGMEN DATA UJI\n(31 Okt 2024 - 10 Feb 2025)',
            ha='left', va='top', fontsize=9, fontweight='bold', color='#27ae60',
            bbox=dict(boxstyle="round,pad=0.3", fc="#fdfefe", ec="#27ae60", lw=1))
            
    rmse_train = np.sqrt(np.mean((b_tr - pred_tr)**2))
    rmse_test = np.sqrt(np.mean((b_te_cont - pred_te_cont)**2))
    
    ax.set_title(f'Evaluasi Deret Waktu Kontinu Model SETAR 2-Rezim: Train (RMSE = {rmse_train*100:.3f}%) vs Test (RMSE = {rmse_test*100:.3f}%)',
                 fontsize=11.5, fontweight='bold', pad=12)
    ax.set_xlabel('Kronologi Waktu Observasi Terpadu ($t$)', fontsize=10)
    ax.set_ylabel('Return Harian Saham (%)', fontsize=10)
    ax.legend(loc='lower left', fontsize=8.5, framealpha=0.92, ncol=3)
    ax.grid(True)
    
    plt.tight_layout()
    out_path = os.path.join(fig_dir, 'fig4_setar_continuous_overlay.png')
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Gambar 4 tersimpan di: {out_path}")

if __name__ == '__main__':
    print("Memulai pembuatan seluruh visualisasi akademik...")
    generate_fig1()
    generate_fig2()
    generate_fig3()
    generate_fig4()
    print("Semua visualisasi berhasil digenerasi 100%!")
