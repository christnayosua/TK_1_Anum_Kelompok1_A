# Petunjuk Eksekusi Program (README)
## Tugas Kelompok 1 — Analisis Numerik Gasal 2026/2027

Dokumen ini menjelaskan struktur berkas kode dan panduan langkah demi langkah untuk menjalankan program solver numerik mandiri pada Tugas Kelompok 1 (TK 1), mencakup pengerjaan **Nomor 1 (iv)** dan **Nomor 2 (iii, iv)** untuk **Kelompok Ganjil** (Kode Data: A, Metode QR: Householder Reflections).

---

### 1. Struktur Direktori

```
04_Tugas Kelompok (TK)/TK 1 (Rilis M03 - Deadline M06)/
├── src/
│   ├── no1_solver.py          # Implementasi Dense LU PP & Banded Thomas-style PP
│   ├── no2_solver.py          # Implementasi Normal Equations & Householder QR + H1 Verifier
│   └── test_solvers.py        # Automated test suite (validasi akurasi numerik)
├── Nomor 1/
│   └── Nomor 1/
│       ├── A/                 # Dataset Kode A (T_16.csv, ..., T_512.csv)
│       └── B/                 # Dataset Kode B
├── Nomor 2/
│   └── Nomor 2/
│       ├── stock_train.csv    # Data deret harga latih (303 baris observasi)
│       └── stock_test.csv     # Data deret harga uji (103 baris observasi)
├── TK1_Bagian_Yosua_Technical_Report.md  # Naskah Technical Report Bagian Yosua
└── README.md                  # Berkas petunjuk ini
```

---

### 2. Prasyarat Lingkungan (Environment Requirements)

- **Python**: Versi 3.9 atau lebih baru.
- **Library yang Dibutuhkan**:
  - `numpy` (hanya digunakan untuk struktur data array dasar, manipulasi shape, dan perkalian matriks/vektor; **bukan** untuk fungsi solver/dekomposisi bawaan).
  - `pandas` (hanya digunakan untuk membaca file CSV ke dalam matriks numerik).

Instalasi dependensi dapat dilakukan dengan perintah:
```bash
pip install numpy pandas
```

---

### 3. Panduan Menjalankan Program

Semua perintah di bawah ini dijalankan dari direktori kerja utama tugas:
`c:\Users\yosua\Documents\Archieve Kuliah\Archieve Semester 5\ALNUM\04_Tugas Kelompok (TK)\TK 1 (Rilis M03 - Deadline M06)`

#### A. Menjalankan Solver Nomor 1 (Analisis Perpindahan Penumpang)
Program ini akan memuat seluruh matriks transisi $T$ pada folder `Nomor 1/Nomor 1/A` untuk ukuran $N \in \{16, 32, 64, 128, 256, 512\}$, membentuk matriks $B$ dan $b$, mengeksekusi kedua solver (Dense LU PP dan Banded Thomas-style PP), serta mencetak perbandingan waktu komputasi, galat residual, dan kondisi matriks.

```bash
python src/no1_solver.py
```

#### B. Menjalankan Solver Nomor 2 (Model SETAR 2-Rezim)
Program ini akan memuat data harga dari `Nomor 2/Nomor 2/stock_train.csv`, menghitung deret return harian $R_t$, menyusun matriks overdetermined $A \in \mathbb{R}^{300 \times 6}$ dan target $\vec{b} \in \mathbb{R}^{300}$, mengeksekusi penyelesaian Persamaan Normal, memverifikasi reflektor langkah awal $H_1$, dan menyelesaikan LSP melalui dekomposisi Householder QR.

```bash
python src/no2_solver.py
```

#### C. Menjalankan Seluruh Unit Test (Otomatis)
Untuk memverifikasi kebenaran matematis dan toleransi numerik secara komprehensif pada kedua nomor sekaligus:

```bash
python src/test_solvers.py
```

Hasil yang diharapkan adalah:
`ALL UNIT TESTS PASSED SUCCESSFULLY!`
dengan selisih perbedaan solusi dense vs banded bernilai $0.00$ dan residual di batas presisi mesin $O(10^{-16})$.

#### D. Menghasilkan Seluruh Grafik Visualisasi (Otomatis)
Untuk menghasilkan 4 grafik beresolusi tinggi (300 DPI) yang dicantumkan dalam laporan teknis:

```bash
python src/generate_visualizations.py
```

Seluruh gambar akan otomatis tersimpan dalam folder `figures/`:
1. `figures/fig1_distribusi_halte.png` : Kurva probabilitas stasioner $\pi_i$ dan kerapatan penumpang antarhalte.
2. `figures/fig2_performa_solver.png`  : Grafik komparasi runtime Dense LU vs Banded Thomas (skala log-log) serta faktor percepatan.
3. `figures/fig3_setar_train_overlay.png` : Time-series overlay return aktual vs estimasi SETAR 2-rezim pada data latih beserta deret residual.
4. `figures/fig4_setar_continuous_overlay.png` : Continuous time-series plot menggabungkan data latih dan uji dengan garis batas transisi.
