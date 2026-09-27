"""
Program: build_notebook.py
Deskripsi: Membuat berkas Jupyter Notebook (TK1_Bagian_Yosua_Solvers.ipynb) yang rapi,
           lengkap dengan markdown penjelasan, implementasi kode mandiri bagian Yosua,
           pengujian data riil, dan referensi bergaya APA 7th Edition untuk diunggah ke Google Colab.
"""

import json
import os

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
notebook_path = os.path.join(base_dir, 'TK1_Bagian_Yosua_Solvers.ipynb')

cells = []

def add_markdown(source):
    cells.append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [s + "\n" for s in source.split("\n")]
    })

def add_code(source):
    cleaned = source.replace(r'\"\"\"', '"""')
    cells.append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [s + "\n" for s in cleaned.split("\n")]
    })

# =========================================================================
# 1. JUDUL & INFORMASI UMUM
# =========================================================================
add_markdown(r"""# Tugas Kelompok 1 — Analisis Numerik (CSCM603117)
## Implementasi Solver SPL dan Dekomposisi QR Householder
**Kelompok:** Best Multi-Cohort Group (Kelompok Ganjil — Kode Data A)  
**Semester:** Gasal 2026/2027 — Fakultas Ilmu Komputer, Universitas Indonesia

---

### Cakupan Implementasi Komputasi:
1. **Nomor 1 (iv) — Implementasi Solver SPL Distribusi Stasioner**:
   - Solver 1: Faktorisasi LU Dense dengan *Partial Pivoting* ($PB = LU$) dengan kompleksitas $O(N^3)$ waktu dan $O(N^2)$ memori.
   - Solver 2: Solver Banded Thomas-Style dengan *Partial Pivoting* yang hanya menyimpan 5 vektor 1D ($O(N)$ waktu dan $O(N)$ memori).
2. **Nomor 2 (iii) — Penyelesaian Least Squares via Persamaan Normal**:
   - Pembentukan sistem normal $A^T A \vec{x} = A^T \vec{b}$.
   - Analisis pemburukan bilangan kondisi: pembuktian analitis dan empiris $\kappa_2(A^T A) = (\kappa_2(A))^2$.
3. **Nomor 2 (iv) — Pemilihan Strategi Dekomposisi QR & Penyelesaian LSP (Kelompok Ganjil: Householder)**:
   - Justifikasi keunggulan metode *Householder Reflections* atas *Givens Rotations* pada matriks *dense tall*.
   - Verifikasi analitis dan komputasi reflektor langkah awal $H_1$ untuk mengeliminasi seluruh sub-diagonal kolom pertama ($H_1 A$).
   - Penyelesaian kuadrat terkecil $A\vec{x} \approx \vec{b}$ via Householder QR dan substitusi mundur $R_1 \vec{x}_{LS} = \vec{c}_1$.
""")

# =========================================================================
# 2. SETUP LINGKUNGAN & AKUISISI DATA
# =========================================================================
add_markdown(r"""## 0. Setup Lingkungan dan Pemuatan Dataset
Notebook ini dirancang adaptif: dapat dijalankan di lingkungan **Google Colab** (dengan mengunggah folder atau clone repository) maupun pada lingkungan lokal. Seluruh algoritma diimplementasikan **from scratch** menggunakan struktur array dasar tanpa mengandalkan fungsi solver bawaan (`scipy.linalg` atau `numpy.linalg.solve`).
""")

add_code(r"""import os
import glob
import zipfile
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Konfigurasi presisi tampilan floating-point
np.set_printoptions(precision=8, suppress=True)

# 1. Otomatis ekstrak berkas zip jika diunggah ke Google Colab (/content/) atau direktori lokal
zip_files = glob.glob("/content/*.zip") + glob.glob("*.zip")
for zfile in zip_files:
    # Hindari mengekstrak zip yang bukan dataset
    if "Nomor" in zfile or "dataset" in zfile.lower():
        print(f"Mengekstrak berkas zip: {zfile}...")
        try:
            with zipfile.ZipFile(zfile, 'r') as zip_ref:
                dest = "/content" if os.path.exists("/content") else os.getcwd()
                zip_ref.extractall(dest)
        except Exception as e:
            print(f"Peringatan saat ekstrak {zfile}: {e}")

# 2. Deteksi direktori dataset adaptif (mendukung lokal, Colab, hierarki langsung atau bertingkat)
def get_data_paths():
    current_dir = os.getcwd()
    candidates_no1 = [
        os.path.join(current_dir, "Nomor 1", "A"),
        os.path.join(current_dir, "Nomor 1", "Nomor 1", "A"),
        "/content/Nomor 1/A",
        "/content/Nomor 1/Nomor 1/A",
        os.path.join(current_dir, "A")
    ]
    candidates_no2 = [
        os.path.join(current_dir, "Nomor 2"),
        os.path.join(current_dir, "Nomor 2", "Nomor 2"),
        "/content/Nomor 2",
        "/content/Nomor 2/Nomor 2",
        current_dir
    ]
    
    path_no1 = next((p for p in candidates_no1 if os.path.exists(p) and os.path.isdir(p) and os.path.exists(os.path.join(p, "T_16.csv"))), candidates_no1[0])
    path_no2 = next((p for p in candidates_no2 if os.path.exists(p) and os.path.isdir(p) and os.path.exists(os.path.join(p, "stock_train.csv"))), candidates_no2[0])
    return path_no1, path_no2

path_no1, path_no2 = get_data_paths()
print(f"Path Dataset Nomor 1: {path_no1} (Ditemukan: {os.path.exists(path_no1)})")
print(f"Path Dataset Nomor 2: {path_no2} (Ditemukan: {os.path.exists(path_no2)})")
""")

# =========================================================================
# 3. NOMOR 1: IMPLEMENTASI SOLVER SPL
# =========================================================================
add_markdown(r"""## 1. Nomor 1 (iv) — Implementasi Solver SPL Distribusi Stasioner

### Kerangka Matematis:
Diberikan matriks transisi stokastik $T \in \mathbb{R}^{N \times N}$ dengan struktur pita ($p_T=2, q_T=1$). Distribusi stasioner $\pi$ memenuhi:
$$T^T \pi = \pi \iff (I - T^T)\pi = \mathbf{0}$$
Karena matriks homogen $A = I - T^T$ singular ($\det(A) = 0$), baris pertama digantikan dengan penetapan non-homogen $B_{1, :} = [1, 0, \dots, 0]$ dan $b = [1, 0, \dots, 0]^T$ untuk menghasilkan sistem tak-singular:
$$Bz = b$$
Vektor solusi $z$ kemudian diproyeksikan ke ruang probabilitas terikat melalui normalisasi:
$$\pi = \frac{z}{\mathbf{1}^T z} = \frac{z}{\sum_{i=1}^N z_i}$$

Matriks $B$ memiliki *lower bandwidth* $p_B = 1$ dan *upper bandwidth* $q_B = 2$. Selama proses *partial pivoting*, pertukaran baris $k$ dan $k+1$ dapat memunculkan 1 super-diagonal tambahan (*fill-in*), sehingga $q_{eff} = q_B + p_B = 3$.
""")

add_code(r"""def load_transition_matrix(csv_path: str) -> np.ndarray:
    \"\"\"Membaca matriks transisi T dari file CSV.\"\"\"
    df = pd.read_csv(csv_path, header=None)
    return df.values.astype(float)

def construct_B_and_b(T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    \"\"\"
    Membentuk matriks B dan vektor b dari matriks transisi T:
    A = I - T^T
    B[0, :] = [1, 0, ..., 0]
    B[i, :] = A[i, :] untuk i = 1, ..., N-1
    b = [1, 0, ..., 0]^T
    \"\"\"
    N = T.shape[0]
    A = np.eye(N) - T.T
    B = A.copy()
    B[0, :] = 0.0
    B[0, 0] = 1.0
    b = np.zeros(N)
    b[0] = 1.0
    return B, b

def normalize_solution(z: np.ndarray) -> np.ndarray:
    \"\"\"Normalisasi vektor z agar memenuhi aksioma probabilitas total 1^T pi = 1.\"\"\"
    return z / np.sum(z)

def compute_errors(T: np.ndarray, pi: np.ndarray) -> tuple[float, float]:
    \"\"\"Menghitung galat residual stasioner ||T^T pi - pi||_2 dan galat normalisasi |1^T pi - 1|.\"\"\"
    r = float(np.linalg.norm(T.T @ pi - pi))
    enorm = float(abs(np.sum(pi) - 1.0))
    return r, enorm
""")

add_markdown(r"""### 1.1 Solver 1: Faktorisasi Dense LU dengan Partial Pivoting ($PB = LU$)
Memfaktorisasi matriks koefisien $B$ menjadi $PB = LU$ dengan pertukaran baris terarah untuk menjaga kestabilan numerik (menjamin $|L_{ij}| \le 1$).  
- **Kompleksitas FLOPs**: $\frac{2}{3}N^3 + 2N^2$  
- **Kebutuhan Memori**: $N \times N$ elemen float64 ($8N^2$ bytes)
""")

add_code(r"""def solve_dense_lu_pp(B_in: np.ndarray, b_in: np.ndarray) -> np.ndarray:
    \"\"\"
    Solver SPL Bz = b berbasis Faktorisasi LU Dense dengan Partial Pivoting.
    1. Faktorisasi: P B = L U
    2. Substitusi Maju: L y = P b
    3. Substitusi Mundur: U z = y
    \"\"\"
    n = B_in.shape[0]
    A = B_in.copy().astype(float)
    b = b_in.copy().astype(float)
    p = np.arange(n)
    
    # 1. Eliminasi Maju / Faktorisasi
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
        
    # 2. Substitusi Maju: L y = P b
    b_perm = b[p]
    y = np.zeros(n)
    for i in range(n):
        y[i] = b_perm[i] - np.dot(A[i, :i], y[:i])
        
    # 3. Substitusi Mundur: U z = y
    z = np.zeros(n)
    for i in range(n - 1, -1, -1):
        z[i] = (y[i] - np.dot(A[i, i+1:], z[i+1:])) / A[i, i]
        
    return z
""")

add_markdown(r"""### 1.2 Solver 2: Banded Thomas-Style Solver dengan Partial Pivoting
Solver khusus yang memanfaatkan struktur pita sempit ($p_B=1, q_B=2, q_{eff}=3$). Hanya menyimpan **5 vektor 1D berukuran $N$**:
1. $d \in \mathbb{R}^N$: Diagonal utama
2. $dl \in \mathbb{R}^{N-1}$: Sub-diagonal ($p=1$)
3. $du_1 \in \mathbb{R}^{N-1}$: Super-diagonal 1 ($q_1$)
4. $du_2 \in \mathbb{R}^{N-2}$: Super-diagonal 2 ($q_2$)
5. $du_3 \in \mathbb{R}^{N-3}$: Super-diagonal 3 (*fill-in* akibat pivoting)
- **Kompleksitas FLOPs**: $O(N \cdot p_B \cdot q_{eff}) = O(N)$ (linear)  
- **Kebutuhan Memori**: $5N$ elemen float64 ($40N$ bytes) alih-alih $8N^2$ bytes
""")

add_code(r"""def solve_banded_thomas_pp(B_in: np.ndarray, b_in: np.ndarray) -> np.ndarray:
    \"\"\"
    Solver SPL Bz = b bergaya Algoritma Thomas untuk matriks banded
    dengan Partial Pivoting (memori kompak O(N)).
    \"\"\"
    n = B_in.shape[0]
    
    # Ekstraksi komponen pita
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
    
    # 1. Eliminasi Maju dengan Partial Pivoting (antara baris k dan k+1)
    for k in range(n - 1):
        if abs(dl[k]) > abs(d[k]):
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
            
    # 2. Substitusi Mundur Terkompresi: U z = b
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
""")

add_markdown(r"""### 1.3 Eksperimen Komparasi Kinerja dan Hasil Uji pada Seluruh Ukuran ($N=16$ hingga $512$)
Kami mengeksekusi kedua solver pada dataset Kode A untuk ukuran $N \in \{16, 32, 64, 128, 256, 512\}$, mencatat waktu eksekusi, selisih solusi, galat residual, serta halte probabilitas puncak.
""")

add_code(r"""results_no1 = []
sizes = [16, 32, 64, 128, 256, 512]

print(f"{'N':>4} | {'t_dense':>10} | {'t_banded':>10} | {'Speedup':>8} | {'Diff Max':>10} | {'Residual':>10} | {'Halte Puncak':>12}")
print("-" * 78)

for N in sizes:
    csv_file = os.path.join(path_no1, f"T_{N}.csv")
    if not os.path.exists(csv_file):
        print(f"Berkas {csv_file} tidak ditemukan. Lewati N={N}.")
        continue
        
    T = load_transition_matrix(csv_file)
    B, b = construct_B_and_b(T)
    
    # Dense LU
    t0 = time.perf_counter()
    z_dense = solve_dense_lu_pp(B, b)
    t_dense = (time.perf_counter() - t0) * 1000.0
    pi_dense = normalize_solution(z_dense)
    r_dense, e_dense = compute_errors(T, pi_dense)
    
    # Banded Thomas
    t0 = time.perf_counter()
    z_banded = solve_banded_thomas_pp(B, b)
    t_banded = (time.perf_counter() - t0) * 1000.0
    pi_banded = normalize_solution(z_banded)
    r_banded, e_banded = compute_errors(T, pi_banded)
    
    diff_max = float(np.max(np.abs(pi_dense - pi_banded)))
    speedup = t_dense / t_banded if t_banded > 0 else 1.0
    max_h = int(np.argmax(pi_banded) + 1)
    max_val = float(np.max(pi_banded))
    
    results_no1.append({
        "N": N,
        "t_dense (ms)": t_dense,
        "t_banded (ms)": t_banded,
        "speedup": speedup,
        "diff_max": diff_max,
        "residual": r_banded,
        "max_halte": max_h,
        "max_pi": max_val
    })
    
    print(f"{N:>4} | {t_dense:>8.3f} ms | {t_banded:>8.3f} ms | {speedup:>7.1f}x | {diff_max:>10.2e} | {r_banded:>10.2e} | H{max_h} ({max_val:.4f})")
""")

add_markdown(r"""### 1.4 Detail Distribusi Stasioner Ukuran Kecil ($N=16$)
Berikut adalah vektor probabilitas stasioner $\pi$ untuk setiap halte pada koridor kecil $N=16$ (terlihat konsentrasi penumpang terbesar berada pada Halte 7 dengan probabilitas $\approx 0.0775$):
""")

add_code(r"""# Tampilkan nilai distribusi halte pada N=16
csv_16 = os.path.join(path_no1, "T_16.csv")
if os.path.exists(csv_16):
    T16 = load_transition_matrix(csv_16)
    B16, b16 = construct_B_and_b(T16)
    pi16 = normalize_solution(solve_banded_thomas_pp(B16, b16))
    
    df_pi16 = pd.DataFrame({
        "Nomor Halte": [f"Halte {i+1}" for i in range(16)],
        "Probabilitas Stasioner (pi_i)": pi16
    })
    print(df_pi16.to_string(index=False))
""")

add_markdown(r"""### Referensi Akademik (Nomor 1):
- Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.
- Heath, M. T. (2002). *Scientific Computing: An Introductory Survey* (2nd ed.). McGraw-Hill.
- Jiwanggi, M. A., Basaruddin, T., & Ibrohim, M. O. (2026). *Slide Perkuliahan Analisis Numerik: System of Linear Equations*. Fakultas Ilmu Komputer, Universitas Indonesia.
- Sauer, T. (2012). *Numerical Analysis* (2nd ed.). Pearson.
""")

# =========================================================================
# 4. NOMOR 2: MODEL SETAR & PERSAMAAN NORMAL VS QR HOUSEHOLDER
# =========================================================================
add_markdown(r"""---
## 2. Nomor 2 (iii, iv) — Model SETAR 2-Rezim, Persamaan Normal & QR Householder

### Kerangka Matematis:
Model SETAR 2-rezim dengan *lag order* $p=2$ dan threshold $c=0$:
$$R_t = \begin{cases} 
\alpha_1 + \phi_{1,1} R_{t-1} + \phi_{1,2} R_{t-2} + \varepsilon_t, & \text{jika } R_{t-1} \ge 0 \quad (\text{Rezim 1: Bullish}) \\
\alpha_2 + \phi_{2,1} R_{t-1} + \phi_{2,2} R_{t-2} + \varepsilon_t, & \text{jika } R_{t-1} < 0 \quad (\text{Rezim 2: Bearish})
\end{cases}$$
Dilinierkan ke sistem *overdetermined* $A\vec{x} \approx \vec{b}$ ($A \in \mathbb{R}^{300 \times 6}, \vec{x} \in \mathbb{R}^6, \vec{b} \in \mathbb{R}^{300}$).
""")

add_code(r"""def load_and_preprocess_stock(csv_path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    \"\"\"
    Membaca data harga saham dari CSV dan menyusun matriks overdetermined A dan vektor target b.
    \"\"\"
    df = pd.read_csv(csv_path)
    prices = df['Close'].values.astype(float)
    dates = df['Date'].values
    
    # 1. Deret return harian: R_t = (P_t - P_{t-1}) / P_{t-1}
    returns = (prices[1:] - prices[:-1]) / prices[:-1]
    
    # 2. Konstruksi matriks A (m x 6) dan b (m) dengan lag p = 2
    m = len(returns) - 2  # 300 baris
    A = np.zeros((m, 6))
    b = np.zeros(m)
    
    for i in range(m):
        t_idx = i + 2
        R_t = returns[t_idx]
        R_t_1 = returns[t_idx - 1]
        R_t_2 = returns[t_idx - 2]
        
        b[i] = R_t
        if R_t_1 >= 0.0:
            A[i, 0] = 1.0       # alpha_1
            A[i, 1] = R_t_1     # phi_1,1
            A[i, 2] = R_t_2     # phi_1,2
        else:
            A[i, 3] = 1.0       # alpha_2
            A[i, 4] = R_t_1     # phi_2,1
            A[i, 5] = R_t_2     # phi_2,2
            
    return A, b, returns, dates

train_csv = os.path.join(path_no2, "stock_train.csv")
A, b, returns, dates = load_and_preprocess_stock(train_csv)
print(f"Dimensi Matriks Desain A : {A.shape}")
print(f"Dimensi Vektor Target b  : {b.shape}")
""")

add_markdown(r"""### 2.1 Penyelesaian via Persamaan Normal (Butir iii)
Mencari solusi kuadrat terkecil melalui sistem normal:
$$A^T A \vec{x} = A^T \vec{b}$$
- Bilangan kondisi matriks $A$: $\kappa_2(A) = \frac{\sigma_{max}(A)}{\sigma_{min}(A)}$
- Bilangan kondisi matriks Gramian $A^T A$: $\kappa_2(A^T A) = \left(\frac{\sigma_{max}(A)}{\sigma_{min}(A)}\right)^2 = (\kappa_2(A))^2$
- **Bahaya Numerik**: Kehilangan presisi berbanding lurus dengan $\log_{10}(\kappa)$. Pada Persamaan Normal, digit presisi yang hilang adalah $\log_{10}(3.71 \times 10^4) \approx 4.57$ digit.
""")

add_code(r"""def solve_normal_equations(A: np.ndarray, b: np.ndarray) -> dict:
    \"\"\"
    Menyelesaikan Least Squares via Persamaan Normal mandiri: (A^T A) x = A^T b
    \"\"\"
    ATA = A.T @ A
    ATb = A.T @ b
    
    cond_A = float(np.linalg.cond(A))
    cond_ATA = float(np.linalg.cond(ATA))
    
    # Selesaikan SPL 6x6 menggunakan Faktorisasi LU PP mandiri
    x_normal = solve_dense_lu_pp(ATA, ATb)
    residual_norm = float(np.linalg.norm(A @ x_normal - b))
    
    return {
        "x": x_normal,
        "ATA": ATA,
        "ATb": ATb,
        "cond_A": cond_A,
        "cond_ATA": cond_ATA,
        "residual_norm": residual_norm
    }

res_normal = solve_normal_equations(A, b)
print("--- HASIL PERSAMAAN NORMAL ---")
print(f"Condition Number kappa_2(A)     : {res_normal['cond_A']:.4f}")
print(f"Condition Number kappa_2(A^T A) : {res_normal['cond_ATA']:.4f}")
print(f"Pembuktian Rasio kappa_2(A^T A) / kappa_2(A)^2 : {res_normal['cond_ATA'] / (res_normal['cond_A']**2):.6f} (Eksak 1.000)")
print(f"Norm Residual ||A x - b||_2     : {res_normal['residual_norm']:.8f}")
""")

add_markdown(r"""### 2.2 Verifikasi Langkah Eliminasi Awal Reflektor $H_1$ (Butir iv)
Untuk Kelompok Ganjil, ditugaskan metode **Householder Reflections**.  
Reflektor langkah awal $H_1 \in \mathbb{R}^{300 \times 300}$ dihitung dari kolom pertama $A$:
$$\vec{a}_1 = A_{:, 0}, \quad \alpha_1 = -\operatorname{sgn}(A_{0, 0}) \|\vec{a}_1\|_2, \quad \vec{u}_1 = \vec{a}_1 - \alpha_1 \mathbf{e}_1, \quad \vec{v}_1 = \frac{\vec{u}_1}{\|\vec{u}_1\|_2}$$
$$H_1 = I - 2 \vec{v}_1 \vec{v}_1^T$$
Diverifikasi bahwa perkalian $H_1 A$ mengeliminasi seluruh 299 elemen sub-diagonal pada kolom pertama menjadi nol.
""")

add_code(r"""def householder_reflection_step1(A: np.ndarray) -> dict:
    \"\"\"
    Menghitung vektor refleksi v1, matriks refleksi H1, dan verifikasi eliminasi H1 A.
    \"\"\"
    m = A.shape[0]
    a1 = A[:, 0].copy()
    norm_a1 = np.linalg.norm(a1)
    
    # Pemilihan tanda aman dari pembatalan katastropik
    sign_val = np.sign(a1[0]) if a1[0] != 0 else 1.0
    alpha = -sign_val * norm_a1
    
    v1 = a1.copy()
    v1[0] -= alpha
    v1 = v1 / np.linalg.norm(v1)
    
    H1 = np.eye(m) - 2.0 * np.outer(v1, v1)
    H1A = H1 @ A
    max_subdiag_error = float(np.max(np.abs(H1A[1:, 0])))
    
    return {
        "v1": v1,
        "alpha": alpha,
        "H1": H1,
        "H1A": H1A,
        "max_subdiag_error": max_subdiag_error
    }

res_h1 = householder_reflection_step1(A)
print("--- VERIFIKASI LANGKAH AWAL ELIMINASI HOUSEHOLDER (H1) ---")
print(f"Skalar alpha_1                            : {res_h1['alpha']:.8f}")
print(f"Elemen diagonal pertama (H1 A)[0, 0]      : {res_h1['H1A'][0, 0]:.8f}")
print(f"Maksimum absolut sub-diagonal (H1 A)[1:,0]: {res_h1['max_subdiag_error']:.2e}")
assert res_h1['max_subdiag_error'] < 1e-14, "Sub-diagonal tidak tereliminasi menjadi nol!"
print("Status Verifikasi: PASSED (Sub-diagonal tereliminasi hingga batas presisi mesin!)")
""")

add_markdown(r"""### 2.3 Penyelesaian Penuh LSP via Householder QR (Butir iv)
Transformasi ortogonal penuh diterapkan pada seluruh kolom $k = 1, \dots, n$:
$$Q^T A = R = \begin{bmatrix} R_1 \\ \mathbf{0} \end{bmatrix}, \quad Q^T \vec{b} = \begin{bmatrix} \vec{c}_1 \\ \vec{c}_2 \end{bmatrix}$$
Solusi kuadrat terkecil diperoleh via substitusi mundur sistem segitiga atas non-singular:
$$R_1 \vec{x}_{LS} = \vec{c}_1$$
""")

add_code(r"""def solve_householder_qr(A_in: np.ndarray, b_in: np.ndarray) -> dict:
    \"\"\"
    Penyelesaian Least Squares via Dekomposisi QR Householder (transformasi implisit).
    \"\"\"
    m, n = A_in.shape
    R = A_in.copy().astype(float)
    c = b_in.copy().astype(float)
    
    for k in range(n):
        x = R[k:, k]
        norm_x = np.linalg.norm(x)
        if norm_x < 1e-15:
            continue
            
        sign_val = np.sign(x[0]) if x[0] != 0 else 1.0
        alpha = -sign_val * norm_x
        
        u = x.copy()
        u[0] -= alpha
        v = u / np.linalg.norm(u)
        
        # Terapkan transformasi Householder tanpa membentuk matriks penuh
        v_dot_R = v @ R[k:, k:]
        R[k:, k:] -= 2.0 * np.outer(v, v_dot_R)
        
        v_dot_c = np.dot(v, c[k:])
        c[k:] -= 2.0 * v * v_dot_c
        
    # Substitusi Mundur Segitiga Atas R1 x = c1
    R1 = R[:n, :n]
    c1 = c[:n]
    x_ls = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x_ls[i] = (c1[i] - np.dot(R1[i, i+1:], x_ls[i+1:])) / R1[i, i]
        
    residual_norm = float(np.linalg.norm(A_in @ x_ls - b_in))
    
    return {
        "x_ls": x_ls,
        "residual_norm": residual_norm
    }

res_qr = solve_householder_qr(A, b)

# Komparasi Parameter
labels = ["alpha_1 (Drift Bull)", "phi_1,1 (Lag-1 Bull)", "phi_1,2 (Lag-2 Bull)",
          "alpha_2 (Drift Bear)", "phi_2,1 (Lag-1 Bear)", "phi_2,2 (Lag-2 Bear)"]

df_comp = pd.DataFrame({
    "Parameter": labels,
    "Persamaan Normal": res_normal['x'],
    "Householder QR": res_qr['x_ls'],
    "Selisih Absolut": np.abs(res_normal['x'] - res_qr['x_ls'])
})

print("--- PERBANDINGAN ESTIMASI PARAMETER ---")
print(df_comp.to_string(index=False))

diff_norm = np.linalg.norm(res_normal['x'] - res_qr['x_ls'])
print(f"\nSelisih Norma Antar-Metode ||x_normal - x_QR||_2: {diff_norm:.2e}")
""")

add_markdown(r"""### 2.4 Evaluasi Out-of-Sample pada Dataset Pengujian (`stock_test.csv`)
Parameter model terbaik $\vec{x}_{LS}$ diuji terhadap 103 observasi harga pengujian untuk menghitung *Root Mean Square Error* (RMSE):
$$\text{RMSE} = \sqrt{\frac{1}{N} \sum_{t=1}^N (R_t - \hat{R}_t)^2}$$
""")

add_code(r"""test_csv = os.path.join(path_no2, "stock_test.csv")
if os.path.exists(test_csv):
    df_test = pd.read_csv(test_csv)
    p_test = df_test['Close'].values.astype(float)
    r_test = (p_test[1:] - p_test[:-1]) / p_test[:-1]
    
    # 1. Evaluasi Data Latih (300 observasi)
    pred_train = A @ res_qr['x_ls']
    rmse_train = np.sqrt(np.mean((b - pred_train)**2))
    
    # 2. Evaluasi Data Uji Independen (100 observasi)
    m_test = len(r_test) - 2
    A_test = np.zeros((m_test, 6))
    b_test = np.zeros(m_test)
    for i in range(m_test):
        t_idx = i + 2
        rt = r_test[t_idx]
        rt1 = r_test[t_idx - 1]
        rt2 = r_test[t_idx - 2]
        b_test[i] = rt
        if rt1 >= 0:
            A_test[i, :3] = [1.0, rt1, rt2]
        else:
            A_test[i, 3:] = [1.0, rt1, rt2]
            
    pred_test = A_test @ res_qr['x_ls']
    rmse_test = np.sqrt(np.mean((b_test - pred_test)**2))
    
    print("--- EVALUASI AKURASI MODEL (RMSE) ---")
    print(f"RMSE Data Latih (Train, N=300) : {rmse_train:.6f} ({rmse_train*100:.3f}%)")
    print(f"RMSE Data Uji   (Test, N=100)  : {rmse_test:.6f} ({rmse_test*100:.3f}%)")
""")

add_markdown(r"""### 2.5 Visualisasi Deret Waktu Model SETAR (Prediksi vs. Aktual)
""")

add_code(r"""# Plot Prediksi vs Aktual Data Latih
plt.figure(figsize=(13, 5), dpi=150)
t_axis = np.arange(1, len(b) + 1)
plt.plot(t_axis, b * 100, label='Return Aktual (R_t)', color='#34495e', alpha=0.75, lw=1.2)
plt.plot(t_axis, pred_train * 100, label='Estimasi SETAR 2-Rezim', color='#e74c3c', lw=1.4)
plt.axhline(0, color='black', linestyle='--', lw=0.8, alpha=0.6)
plt.title(f'Prediksi Deret Waktu SETAR 2-Rezim vs Aktual (Data Latih, RMSE = {rmse_train*100:.3f}%)', fontsize=11, fontweight='bold')
plt.xlabel('Waktu Observasi (t)', fontsize=10)
plt.ylabel('Return Harian (%)', fontsize=10)
plt.legend(loc='upper right')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
""")

add_markdown(r"""### 2.6 Interpretasi Finansial Parameter Model
Persamaan empiris model SETAR 2-rezim yang diperoleh dari Householder QR:
$$\hat{R}_t = \begin{cases} 
+0.001513 + 0.172153 R_{t-1} + 0.166073 R_{t-2}, & \text{jika } R_{t-1} \ge 0 \quad (\text{Bullish}) \\
-0.001213 - 0.360377 R_{t-1} - 0.109808 R_{t-2}, & \text{jika } R_{t-1} < 0 \quad (\text{Bearish})
\end{cases}$$

1. **Rezim Bullish ($R_{t-1} \ge 0$)**:
   - Drift $\alpha_1 = +0.001513 > 0$: Pertumbuhan dasar $+0.15\%$ per hari.
   - Koefisien Lag $\phi_{1,1} = +0.1722 > 0$ dan $\phi_{1,2} = +0.1661 > 0$: Keduanya positif $\implies$ **Momentum Persistence** (kenaikan harga kemarin memicu tren kenaikan lanjutan).
2. **Rezim Bearish ($R_{t-1} < 0$)**:
   - Drift $\alpha_2 = -0.001213 < 0$: Tekanan depresiasi $-0.12\%$ per hari.
   - Koefisien Lag $\phi_{2,1} = -0.3604 < 0$ dan $\phi_{2,2} = -0.1098 < 0$: Keduanya negatif $\implies$ **Mean-Reversion / Technical Rebound** (penurunan tajam kemarin memicu aksi beli spekulatif yang menahan koreksi lanjutan).

---

### Referensi Akademik (Nomor 2):
- Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.
- Jiwanggi, M. A., Basaruddin, T., & Ibrohim, M. O. (2026). *Slide Perkuliahan Analisis Numerik: Least Square Problems*. Fakultas Ilmu Komputer, Universitas Indonesia.
- Tong, H. (1990). *Non-linear Time Series: A Dynamical System Approach*. Oxford University Press.
- Trefethen, L. N., & Bau, D. (1997). *Numerical Linear Algebra*. Society for Industrial and Applied Mathematics (SIAM).
""")

# =========================================================================
# 5. KOMPILASI STRUKTUR JSON NOTEBOOK
# =========================================================================
notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3"
        },
        "language_info": {
            "codemirror_mode": {
                "name": "ipython",
                "version": 3
            },
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "nbconvert_exporter": "python",
            "pygments_lexer": "ipython3",
            "version": "3.10.0"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 4
}

for fname in ['TK1_Bagian_Yosua_Solvers.ipynb', 'TK1_Anum_Kelompok1.ipynb']:
    out_path = os.path.join(base_dir, fname)
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, indent=2, ensure_ascii=False)
    print(f"Jupyter Notebook berhasil disimpan di: {out_path}")

print(f"Total cells: {len(cells)}")
