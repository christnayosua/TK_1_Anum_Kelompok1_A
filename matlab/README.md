# Panduan Eksekusi Implementasi MATLAB / GNU Octave

Folder ini berisi kode implementasi alternatif menggunakan bahasa **MATLAB / GNU Octave** untuk pengerjaan **Tugas Kelompok 1 (TK 1) Analisis Numerik Gasal 2026/2027**.

Seluruh algoritma utama diimplementasikan secara **mandiri (from scratch)** tanpa menggunakan operator black-box `\` (mldivide), fungsi `inv()`, `lu()`, atau `qr()`.

---

## 1. Struktur Berkas MATLAB

```
matlab/
├── no1_solver.m       # Solver Nomor 1: Dense LU PP vs Banded Thomas-Style PP
├── no2_solver.m       # Solver Nomor 2: Persamaan Normal, Householder QR, & Givens QR
├── run_all.m          # Master script untuk menjalankan seluruh eksperimen otomatis
└── README.md          # Petunjuk eksekusi ini
```

---

## 2. Cara Menjalankan Program

### Melalui Antarmuka MATLAB Desktop:
1. Buka aplikasi **MATLAB**.
2. Arahkan *Current Folder* ke direktori:
   `c:\Users\yosua\Documents\Archieve Kuliah\Archieve Semester 5\ALNUM\04_Tugas Kelompok (TK)\TK 1 (Rilis M03 - Deadline M06)\matlab`
3. Pada *Command Window*, ketikkan perintah:
   ```matlab
   run_all
   ```
   Atau jalankan masing-masing solver secara terpisah:
   ```matlab
   no1_solver  % Menjalankan validasi, solver banded vs dense, dan plot halte
   no2_solver  % Menjalankan estimasi SETAR via Normal Eq, QR Householder & Givens
   ```

### Melalui Terminal / GNU Octave:
Jika menggunakan GNU Octave melalui terminal shell:
```bash
cd "c:\Users\yosua\Documents\Archieve Kuliah\Archieve Semester 5\ALNUM\04_Tugas Kelompok (TK)\TK 1 (Rilis M03 - Deadline M06)\matlab"
octave run_all.m
```

---

## 3. Fitur Utama yang Diimplementasikan

1. **Nomor 1 (Analisis Perpindahan Penumpang Antarhalte)**:
   - Validasi matriks stokastik baris ($\sum_j T_{ij} = 1$), non-negativitas ($T_{ij} \ge 0$), dan ketunggalan distribusi stasioner via ketereduksian (*irreducibility*).
   - Deteksi otomatis *lower bandwidth* ($p$) dan *upper bandwidth* ($q$).
   - Konstruksi sistem tak-singular $B z = b$ dan normalisasi $\pi = z / (1^T z)$.
   - Faktorisasi Dense LU dengan *Partial Pivoting* ($P B = L U$).
   - Solver Banded Thomas-style dengan *Partial Pivoting* yang hanya memanfaatkan 5 vektor diagonal 1D (memori $O(N)$).
   - Evaluasi galat residual $\|T^T \pi - \pi\|_2$ dan galat normalisasi $|1^T \pi - 1|$ hingga batas presisi ganda ($O(10^{-16})$).

2. **Nomor 2 (Prediksi Return Saham Menggunakan SETAR 2-Rezim)**:
   - Perhitungan deret return harian $R_t = (P_t - P_{t-1})/P_{t-1}$.
   - Penyusunan matriks *overdetermined* $A \in \mathbb{R}^{300 \times 6}$ dan target $\vec{b} \in \mathbb{R}^{300}$.
   - Penyelesaian via Persamaan Normal ($A^T A \vec{x} = A^T \vec{b}$) dan pembuktian $\kappa_2(A^T A) = (\kappa_2(A))^2$.
   - Verifikasi refleksi Householder langkah awal ($H_1$ dan $H_1 A$) dengan eliminasi elemen sub-diagonal kolom pertama ke orde $10^{-15}$.
   - Formulasi rotasi Givens 2D ($G_1$ dan $G_1 A$) untuk mengeliminasi elemen sub-diagonal pertama.
   - Solver Householder QR mandiri dan substitusi mundur $R_1 \vec{x}_{LS} = \vec{c}_1$.
   - Solver Givens Rotations QR mandiri sebagai perbandingan komprehensif.
   - Evaluasi akurasi *out-of-sample* (RMSE) pada data uji `stock_test.csv`.
