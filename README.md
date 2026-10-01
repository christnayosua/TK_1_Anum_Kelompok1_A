# Petunjuk Eksekusi Program (README)
## Tugas Kelompok 1 — Analisis Numerik Gasal 2026/2027
**Fakultas Ilmu Komputer, Universitas Indonesia**  
**Kelompok:** Kelompok Ganjil (Kode Data: A, Metode QR: Householder Reflections)  

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/christnayosua/TK_1_Anum_Kelompok1_A/blob/main/TK1_Anum_Kelompok1.ipynb)

---

### 1. Struktur Direktori dan Berkas

```
04_Tugas Kelompok (TK)/TK 1 (Rilis M03 - Deadline M06)/
├── TK1_Laporan_Lengkap_Technical_Report.md  # Naskah Technical Report Lengkap (Format Markdown)
├── TK1_Technical_Report.docx                # Naskah Laporan dalam format Microsoft Word / Google Docs
├── TK1_Anum_Kelompok1.ipynb                 # Jupyter Notebook lengkap dengan output komputasi & grafik (Root)
├── README.md                                # Berkas petunjuk ini
├── notebook/                                # Direktori Siap Pakai untuk Jupyter / Google Colab
│   ├── TK1_Anum_Kelompok1.ipynb             # Salinan notebook mandiri
│   ├── figures/                             # Visualisasi grafik resolusi tinggi (300 DPI)
│   ├── Nomor 1/                             # Dataset Matriks Transisi
│   └── Nomor 2/                             # Dataset Deret Waktu Harga Saham
├── figures/                                 # 4 Berkas Grafik Visualisasi Akademik Resolusi Tinggi (300 DPI)
│   ├── fig1_distribusi_halte.png            # Visualisasi sebaran probabilitas stasioner penumpang (Nomor 1)
│   ├── fig2_performa_solver.png             # Komparasi runtime log-log & rasio percepatan (Nomor 1)
│   ├── fig3_setar_train_overlay.png         # Overlay return aktual vs estimasi SETAR & residual (Nomor 2)
│   └── fig4_setar_continuous_overlay.png    # Continuous time-series plot Train + Test dengan pembatas (Nomor 2)
├── src/                                     # Modul Program Sumber Bahasa Python
│   ├── no1_solver.py                        # Solver SPL Nomor 1: Dense LU PP vs Banded Thomas-Style PP
│   ├── no2_solver.py                        # Solver LSP Nomor 2: Normal Eq, Householder QR, Givens QR, Iterative Refinement
│   ├── test_solvers.py                      # Automated Unit Test Suite (validasi presisi numerik)
│   ├── generate_visualizations.py           # Skrip otomasi pembuat seluruh grafik akademik
│   ├── build_notebook.py                    # Generator berkas Jupyter Notebook (.ipynb)
│   └── convert_to_docx.py                   # Generator berkas Microsoft Word (.docx)
├── matlab/                                  # Alternatif Implementasi Mandiri Bahasa MATLAB / GNU Octave
│   ├── no1_solver.m                         # Implementasi Solver Nomor 1 (Dense LU PP vs Banded Thomas PP)
│   ├── no2_solver.m                         # Implementasi Solver Nomor 2 (Normal Eq, Householder, Givens)
│   ├── run_all.m                            # Master script untuk menjalankan seluruh pengujian MATLAB
│   └── README.md                            # Panduan eksekusi program MATLAB / Octave
├── Nomor 1/                                 # Dataset Matriks Transisi
│   ├── A/                                   # Dataset Kode A (T_16.csv s.d. T_512.csv) [Kelompok Ganjil]
│   └── B/                                   # Dataset Kode B [Kelompok Genap]
└── Nomor 2/                                 # Dataset Deret Waktu Harga Saham
    ├── stock_train.csv                      # Data deret harga latih (303 observasi)
    └── stock_test.csv                       # Data deret harga uji (103 observasi)
```

---

### 2. Prasyarat Lingkungan (Environment Requirements)

- **Python**: Versi 3.9 atau lebih baru.
- **Library Python yang Digunakan**:
  - `numpy` (struktur data array dasar, manipulasi dimensi, dan perkalian matriks; **bukan** solver bawaan).
  - `pandas` (membaca berkas CSV).
  - `matplotlib` (menghasilkan visualisasi grafik 300 DPI).
  - `python-docx` (konversi laporan ke berkas docx).

Instalasi dependensi dapat dilakukan dengan perintah:
```bash
pip install numpy pandas matplotlib python-docx
```

---

### 3. Panduan Menjalankan Program Python

Semua perintah dijalankan dari direktori kerja utama tugas:
`c:\Users\yosua\Documents\Archieve Kuliah\Archieve Semester 5\ALNUM\04_Tugas Kelompok (TK)\TK 1 (Rilis M03 - Deadline M06)`

#### A. Menjalankan Solver Nomor 1 (Analisis Perpindahan Penumpang)
Mengeksekusi validasi stokastik matriks $T$, deteksi bandwidth $p$ dan $q$, konstruksi $B z = b$, serta membandingkan runtime dan galat antara solver *Dense LU PP* ($O(N^3)$) dan *Banded Thomas PP* ($O(N)$):
```bash
python src/no1_solver.py
```

#### B. Menjalankan Solver Nomor 2 (Model SETAR 2-Rezim)
Mengeksekusi pembentukan matriks *overdetermined* $A \in \mathbb{R}^{300 \times 6}$, penyelesaian Persamaan Normal ($\kappa_2(A^T A) = \kappa_2(A)^2$), verifikasi reflektor awal $H_1$, faktorisasi Householder QR, dan evaluasi *out-of-sample*:
```bash
python src/no2_solver.py
```

#### C. Menjalankan Unit Test Presisi Tinggi (Otomatis)
Untuk memverifikasi kebenaran matematis seluruh algoritma sekaligus:
```bash
python src/test_solvers.py
```
*Hasil yang diharapkan: `ALL UNIT TESTS PASSED SUCCESSFULLY!` dengan residual $O(10^{-16})$ dan selisih solusi $0.00$.*

#### D. Menghasilkan Ulang Seluruh Visualisasi Grafik
```bash
python src/generate_visualizations.py
```

#### E. Mengonversi Laporan Markdown ke Microsoft Word (.docx)
```bash
python src/convert_to_docx.py
```

#### F. Menjalankan Jupyter Notebook (Lokal atau Google Colab)
- **Eksekusi Lokal**:
  Buka terminal di root repositori dan jalankan:
  ```bash
  jupyter notebook TK1_Anum_Kelompok1.ipynb
  ```
  Atau masuk ke folder `notebook/`:
  ```bash
  cd notebook
  jupyter notebook TK1_Anum_Kelompok1.ipynb
  ```
- **Eksekusi di Google Colab**:
  1. Buka [Google Colab](https://colab.research.google.com/).
  2. Unggah berkas `TK1_Anum_Kelompok1.ipynb` dari folder `notebook/` (atau langsung dari root).
  3. Unggah folder dataset `Nomor 1/` dan `Nomor 2/` ke file explorer Colab (`/content/`).
  4. Jalankan seluruh *cells* secara berurutan (*Runtime -> Run all*). Seluruh visualisasi grafik dan tabel komparasi numerik akan langsung dirender secara interaktif.

---

### 4. Panduan Menjalankan Program MATLAB / GNU Octave

Seluruh algoritma mandiri juga disediakan dalam format **MATLAB / GNU Octave** di folder `matlab/`.
1. Buka aplikasi **MATLAB** atau jalankan terminal **GNU Octave**.
2. Arahkan *Current Directory* ke sub-folder `matlab/`:
   ```matlab
   cd 'matlab'
   ```
3. Jalankan master script:
   ```matlab
   run_all
   ```
   Atau jalankan masing-masing pengujian secara terpisah (`no1_solver` atau `no2_solver`).

---

### 5. Ringkasan Hasil Eksperimen Utama

1. **Nomor 1 (Solver SPL Distribusi Stasioner Halte):**
   - Struktur pita matriks koefisien: $p_B = 1$, $q_B = 2$.
   - Pada $N = 512$, solver *Banded Thomas PP* menghasilkan faktor percepatan (*speedup*) sebesar **$353.76\times$** dan memangkas konsumsi memori sebesar **$99.22\%$** dibandingkan *Dense LU PP*.
   - Norm residual stasioner konsisten berada di batas presisi mesin IEEE 754: $\|T^T\pi - \pi\|_2 \approx 1.33 \times 10^{-16}$.
   - Halte dengan kepadatan penumpang tertinggi pada koridor $N = 16$ adalah **Halte 7** ($\pi_7 = 7.75\%$).

2. **Nomor 2 (Estimasi Parameter SETAR 2-Rezim):**
   - Matriks desain berdimensi $300 \times 6$ dengan target 300 data return harian.
   - Bilangan kondisi: $\kappa_2(A) = 192.61 \implies \kappa_2(A^T A) = (\kappa_2(A))^2 = 37,096.84$ (terbukti eksak kuadrat).
   - Reflektor awal $H_1$ mengeliminasi elemen subdiagonal pertama dengan galat $1.44 \times 10^{-15}$.
   - Evaluasi akurasi peramalan: RMSE data latih = **$0.883\%$**, RMSE data uji independen = **$1.206\%$**.
   - Interpretasi finansial: Terdeteksi sifat *momentum persistence* pada fase Bullish ($\phi_{1,1} = +0.172, \phi_{1,2} = +0.166$) dan *technical rebound* yang kuat pada fase Bearish ($\phi_{2,1} = -0.360$).
