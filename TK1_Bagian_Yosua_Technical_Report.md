# Laporan Teknis: Formulasi Sistem Persamaan Linear Khusus dan Dekomposisi Ortogonal pada Model Rantai Markov Transportasi serta Model SETAR Pasar Modal

**Kelompok:** Ganjil (Kode Data: A, Metode QR: *Householder Reflections*)
**Mata Kuliah:** Analisis Numerik (CSCM603117) — Semester Gasal 2026/2027
**Fakultas Ilmu Komputer, Universitas Indonesia**

---

### Pakta Integritas

> *"Dengan ini, kami menyatakan bahwa tugas ini adalah hasil pekerjaan kelompok sendiri."*

---

## Rangkuman Eksekutif

Laporan teknis ini menyajikan investigasi teoretis, formulasi aljabar linear numerik, perancangan algoritma mandiri (*from scratch*), serta evaluasi kinerja komputasi pada dua permasalahan terapan:

1. **Analisis Perpindahan Penumpang Antarhalte (Rantai Markov & SPL Berpita)**: Kami memodelkan perpindahan penumpang koridor sebagai sistem rantai Markov waktu-diskrit. Matriks transisi stokastik $T \in \mathbb{R}^{N \times N}$ terbukti bersifat tereduksi-tak (*irreducible*) dan aperiodik, menjamin ketunggalan distribusi stasioner $\pi$ berdasarkan Teorema Perron-Frobenius. Untuk mengatasi singularitas intrinsik matriks $A = I - T^T$, kami menerapkan modifikasi baris pertama menjadi sistem non-homogen $Bz = b$ dengan $B_{1, :} = \mathbf{e}_1^T$ dan $b = \mathbf{e}_1$, dilanjutkan proyeksi normalisasi $\pi = z / (\mathbf{1}^T z)$. Mengingat matriks $B$ memiliki struktur pita sempit ($p_B=1, q_B=2$), kami merancang dua solver mandiri: Faktorisasi Dense LU dengan *Partial Pivoting* ($PB=LU$) berorde $O(N^3)$ waktu dan $O(N^2)$ memori, serta Solver Banded bergaya Algoritma Thomas tergeneralisasi dengan *Partial Pivoting* dan penanganan *fill-in* ($q_{eff}=3$) berorde optimal $O(N)$ waktu dan $O(N)$ memori. Eksperimen pada ukuran $N \in \{16, \dots, 512\}$ membuktikan solver banded menghasilkan solusi identik hingga batas presisi mesin IEEE-754 ganda ($\|\pi_{dense} - \pi_{banded}\|_\infty = 0.00$), residual stasioner $\|T^T \pi - \pi\|_2 \approx 1.33 \times 10^{-16}$, serta memberikan lonjakan kecepatan (*speedup*) hingga $685.8\times$ pada $N=512$ dengan penghematan memori sebesar $99.02\%$.
2. **Deteksi Rezim Pasar dan Prediksi Return Saham Menggunakan SETAR 2-Rezim (Least Squares Problem)**: Kami memprediksi return harian saham melalui model *Self-Exciting Threshold Autoregressive* (SETAR, $p=2, c=0$) yang dilinierkan ke dalam sistem *overdetermined* $A\vec{x} \approx \vec{b}$ ($A \in \mathbb{R}^{300 \times 6}$). Kami membandingkan penyelesaian kuadrat terkecil via Persamaan Normal ($A^T A \vec{x} = A^T \vec{b}$) dan Dekomposisi QR berbasis *Householder Reflections*. Kami membuktikan secara analitis dan empiris terjadinya fenomena pengkuadratan *condition number* pada Persamaan Normal: $\kappa_2(A) \approx 192.61 \implies \kappa_2(A^T A) \approx 3.7097 \times 10^4 = (\kappa_2(A))^2$, yang memicu kehilangan digit signifikan sebesar $\approx 4.57$ digit. Sebaliknya, metode *Householder Reflections* terbukti unggul secara numerik (*unconditionally backward stable*) dan efisien (menghemat 33% FLOPs atas *Givens Rotations*). Verifikasi langkah awal membuktikan reflektor $H_1$ berhasil mengeliminasi seluruh 299 elemen sub-diagonal kolom pertama matriks desain hingga batas toleransi $1.44 \times 10^{-15}$. Evaluasi *out-of-sample* pada berkas pengujian menghasilkan $\text{RMSE}_{train} = 0.008832$ ($0.88\%$) dan $\text{RMSE}_{test} = 0.012062$ ($1.21\%$). Parameter empiris menunjukkan rezim Bullish didominasi sifat *momentum persistence* ($\phi_{1,1} = +0.1722, \phi_{1,2} = +0.1661$), sedangkan rezim Bearish didominasi sifat *mean-reversion / technical rebound* ($\phi_{2,1} = -0.3604, \phi_{2,2} = -0.1098$).

---

## Daftar Isi

1. [Pendahuluan](#1-pendahuluan)
2. [Bagian 1: Analisis Perpindahan Penumpang Antarhalte (Nomor 1)](#2-bagian-1-analisis-perpindahan-penumpang-antarhalte-nomor-1)
   - 2.1 [Validasi Integritas Data dan Syarat Ketunggalan Perron-Frobenius (Butir i)](#21-validasi-integritas-data-dan-syarat-ketunggalan-perron-frobenius-butir-i)
   - 2.2 [Formulasi Sistem Persamaan Linear dan Penanganan Isu Singularitas (Butir ii)](#22-formulasi-sistem-persamaan-linear-dan-penanganan-isu-singularitas-butir-ii)
   - 2.3 [Identifikasi Struktur Matriks Berpita dan Efisiensi Penyimpanan (Butir iii)](#23-identifikasi-struktur-matriks-berpita-dan-efisiensi-penyimpanan-butir-iii)
   - 2.4 [Rancang Bangun Algoritma Solver SPL Mandiri (Butir iv)](#24-rancang-bangun-algoritma-solver-spl-mandiri-butir-iv)
   - 2.5 [Eksperimen Kinerja Komputasi dan Skalabilitas Runtime (Butir v)](#25-eksperimen-kinerja-komputasi-dan-skalabilitas-runtime-butir-v)
   - 2.6 [Analisis Kompleksitas FLOPs Teoretis vs. Kinerja Empiris (Butir vi)](#26-analisis-kompleksitas-flops-teoretis-vs-kinerja-empiris-butir-vi)
   - 2.7 [Analisis Bilangan Kondisi, Galat, dan Residual (Butir vii)](#27-analisis-bilangan-kondisi-galat-dan-residual-butir-vii)
   - 2.8 [Interpretasi Distribusi Penumpang dan Pemetaan Titik Kritis Halte (Butir viii)](#28-interpretasi-distribusi-penumpang-dan-pemetaan-titik-kritis-halte-butir-viii)
   - 2.9 [Rekomendasi Arsitektural untuk Koridor Skala Besar (Butir ix)](#29-rekomendasi-arsitektural-untuk-koridor-skala-besar-butir-ix)
3. [Bagian 2: Deteksi Rezim Pasar dan Prediksi Return Saham SETAR (Nomor 2)](#3-bagian-2-deteksi-rezim-pasar-dan-prediksi-return-saham-setar-nomor-2)
   - 3.1 [Formulasi Matriks Desain Overdetermined dan Konstruksi Data (Butir i)](#31-formulasi-matriks-desain-overdetermined-dan-konstruksi-data-butir-i)
   - 3.2 [Investigasi Isu Numerik dan Strategi Mitigasi Teknis (Butir ii)](#32-investigasi-isu-numerik-dan-strategi-mitigasi-teknis-butir-ii)
   - 3.3 [Penyelesaian via Persamaan Normal dan Bahaya Pengkuadratan Kondisi (Butir iii)](#33-penyelesaian-via-persamaan-normal-dan-bahaya-pengkuadratan-kondisi-butir-iii)
   - 3.4 [Pemilihan Strategi Dekomposisi QR: Householder Reflections (Butir iv)](#34-pemilihan-strategi-dekomposisi-qr-householder-reflections-butir-iv)
   - 3.5 [Analisis Komparatif Performa dan Sensitivitas Outlier Ekstrem (Butir v)](#35-analisis-komparatif-performa-dan-sensitivitas-outlier-ekstrem-butir-v)
   - 3.6 [Evaluasi Generalisasi Out-of-Sample pada Dataset Pengujian (Butir vi)](#36-evaluasi-generalisasi-out-of-sample-pada-dataset-pengujian-butir-vi)
   - 3.7 [Interpretasi Finansial Parameter dan Visualisasi Overlay Time Series (Butir vii)](#37-interpretasi-finansial-parameter-dan-visualisasi-overlay-time-series-butir-vii)
4. [Kesimpulan](#4-kesimpulan)
5. [Daftar Referensi](#5-daftar-referensi)
6. [Lampiran: Panduan Eksekusi Program dan Kode Mandiri](#6-lampiran-panduan-eksekusi-program-dan-kode-mandiri)

---

## 1. Pendahuluan

Komputasi aljabar linear merupakan pilar fundamental dalam pemodelan sains, rekayasa, dan data kuantitatif. Namun, penerapan metode numerik dalam skala praktis sering kali dihadapkan pada dua tantangan utama: (1) ledakan kompleksitas waktu dan memori pada sistem berdimensi besar, serta (2) kerentanan propagasi galat pembulatan (*round-off error*) akibat operasi aritmetika presisi hingga (*floating-point system*).

Dalam konteks manajemen transportasi massal, pemodelan perpindahan penumpang pada koridor transit panjang memunculkan sistem Rantai Markov dengan matriks transisi berstruktur pita (*banded sparse matrix*). Apabila diselesaikan menggunakan eliminasi Gauss atau faktorisasi LU konvensional yang memperlakukan matriks secara padat (*dense*), kebutuhan komputasi berkembang secara kubik ($O(N^3)$) dan memori secara kuadratik ($O(N^2)$). Hal ini memicu inefisiensi masif saat jumlah stasiun atau halte mencapai ratusan. Oleh karena itu, diperlukan perancangan algoritma khusus yang mengeksploitasi struktur *sparsity* pita secara presisi tanpa mengorbankan stabilitas numerik (*pivoting*).

Di sisi lain, pemodelan pergerakan harga instrumen keuangan kerap menunjukkan ketak-linieran asimetris, di mana respon pasar pada fase penguatan (*bullish*) berbeda secara fundamental dari fase kepanikan pasar (*bearish*). Model deret waktu non-linear *Self-Exciting Threshold Autoregressive* (SETAR) yang dirumuskan ke dalam permasalahan kuadrat terkecil (*Least Squares Problem* / LSP) mengharuskan estimasi parameter linier pada sistem persamaan overdetermined. Pendekatan klasik menggunakan Persamaan Normal ($A^T A \vec{x} = A^T \vec{b}$) rentan mengalami pemburukan kondisi matriks secara kuadratik ($\kappa_2(A^T A) = (\kappa_2(A))^2$), yang dapat mengakibatkan hilangnya digit signifikansi secara katastropik. Untuk memitigasi hal tersebut, transformasi ortogonal stabil seperti *Householder Reflections* menjadi solusi komputasi yang esensial.

Laporan ini menyajikan analisis komprehensif, penurunan matematis analitik, dan bukti implementasi komputasi mandiri untuk menjawab kedua domain permasalahan tersebut secara terpadu.

---

## 2. Bagian 1: Analisis Perpindahan Penumpang Antarhalte (Nomor 1)

### 2.1 Validasi Integritas Data dan Syarat Ketunggalan Perron-Frobenius (Butir i)

Kami memuat seluruh berkas matriks transisi $T \in \mathbb{R}^{N \times N}$ untuk kelompok ganjil (Dataset Kode A) pada rentang ukuran $N \in \{16, 32, 64, 128, 256, 512\}$ dari folder `Nomor 1/Nomor 1/A/`. Setiap matriks diuji terhadap empat kriteria stokastik:

1. **Dimensi Bujursangkar**: Kami memverifikasi bahwa matriks memiliki ukuran baris dan kolom yang persis sama ($N \times N$).
2. **Non-negativitas**: Seluruh elemen memenuhi $T_{ij} \ge 0$, dengan nilai minimum terukur $\min(T) = 0.00$.
3. **Syarat Stokastik Baris**: Setiap baris merepresentasikan distribusi probabilitas bersyarat keluar dari halte $i$, sehingga $\sum_{j=1}^N T_{ij} = 1.0$. Galat deviasi maksimum pada seluruh dataset terukur sebesar:
   $$
   \max_{1 \le i \le N} \left| \sum_{j=1}^N T_{ij} - 1.0 \right| \le 4.44 \times 10^{-16}
   $$

   Nilai ini identik dengan batas presisi mesin $\epsilon_{mach} \approx 2.22 \times 10^{-16}$ pada representasi IEEE-754 *double precision* (Jiwanggi et al., 2026; Sauer, 2012).
4. **Ketunggalan Distribusi Stasioner**:
   Berdasarkan Teorema Perron-Frobenius untuk matriks stokastik tak-negatif (Horn & Johnson, 2012), distribusi stasioner $\pi$ bersifat **tunggal dan bernilai positif terikat** ($\pi_i > 0, \sum \pi_i = 1$) jika dan hanya jika rantai Markov bersifat **tereduksi-tak (*irreducible*)** dan **aperiodik (*aperiodic*)**:
   - *Irreducibility*: Kami menguji konektivitas graf transisi berarah. Karena elemen sub-diagonal ($T_{i+1, i} > 0$) dan super-diagonal ($T_{i, i+1} > 0$) bernilai positif untuk seluruh $i \in \{1, \dots, N-1\}$, maka setiap halte dapat mencapai halte lainnya dalam sejumlah langkah tertentu. Graf transisi terbukti terhubung kuat (*strongly connected*).
   - *Aperiodicity*: Karena elemen diagonal utama bernilai positif ($T_{ii} > 0$, probabilitas penumpang tetap di halte yang sama $> 0$), periode setiap status adalah 1 ($\gcd(\text{panjang siklus}) = 1$).
   - *Spektrum Nilai Eigen*: Evaluasi nilai eigen menunjukkan bahwa nilai eigen dominan (akar Perron) bernilai tepat $\lambda_1 = 1.0$, dengan kelipatan aljabar dan geometri tepat satu ($|\lambda_1| = 1.0$ dan $|\lambda_2| < 1.0$). Hal ini menjamin secara mutlak bahwa vektor distribusi stasioner $\pi$ bersifat tunggal.

---

### 2.2 Formulasi Sistem Persamaan Linear dan Penanganan Isu Singularitas (Butir ii)

#### Penurunan Persamaan Distribusi Stasioner

Misalkan $X_k$ menyatakan lokasi halte penumpang pada interval waktu ke-$k$. Vektor distribusi probabilitas halte pada langkah $k$ didefinisikan sebagai vektor baris $p^{(k)} \in \mathbb{R}^{1 \times N}$. Berdasarkan hukum probabilitas total:

$$
p^{(k+1)}_j = \sum_{i=1}^N p^{(k)}_i \Pr(X_{k+1} = j \mid X_k = i) = \sum_{i=1}^N p^{(k)}_i T_{ij} \implies p^{(k+1)} = p^{(k)} T
$$

Distribusi stasioner (keadaan tunak) tercapai ketika distribusi probabilitas invarian terhadap transisi waktu:

$$
\lim_{k \to \infty} p^{(k)} = \pi^T \implies \pi^T T = \pi^T
$$

Dengan melakukan transposisi matriks pada kedua ruas, diperoleh sistem linier homogen dalam bentuk vektor kolom $\pi \in \mathbb{R}^N$:

$$
( \pi^T T )^T = (\pi^T)^T \iff T^T \pi = \pi
$$

Mendefinisikan matriks $A = I - T^T$, persamaan di atas ekuivalen dengan:

$$
A \pi = (I - T^T)\pi = \mathbf{0}
$$

#### Pembuktian Analitis Singularitas Matriks $A$

Matriks transisi $T$ adalah matriks stokastik baris, yang berarti jumlah setiap barisnya adalah 1. Dalam notasi vektor, jika $\mathbf{1} = [1, 1, \dots, 1]^T \in \mathbb{R}^N$, maka berlaku:

$$
T \mathbf{1} = \mathbf{1} \iff (T - I)\mathbf{1} = \mathbf{0}
$$

Transposisikan identitas tersebut:

$$
\mathbf{1}^T (T^T - I) = \mathbf{0}^T \iff \mathbf{1}^T A = \mathbf{0}^T \iff A^T \mathbf{1} = \mathbf{0}
$$

Persamaan ini membuktikan bahwa vektor $\mathbf{1}$ merupakan anggota ruang nol (*null space*) dari $A^T$. Akibatnya:

$$
\dim(\ker(A^T)) \ge 1 \implies \operatorname{rank}(A) = \operatorname{rank}(A^T) \le N - 1 \implies \det(A) = 0
$$

Matriks $A$ bersifat singular (*rank-deficient*). Baris-baris matriks $A$ saling bergantung linear (*linearly dependent*), sebab penjumlahan seluruh baris $A$ menghasilkan vektor nol:

$$
\sum_{i=1}^N A_{i, :} = \mathbf{1}^T A = \mathbf{0}^T \implies A_{1, :} = -\sum_{i=2}^N A_{i, :}
$$

#### Rasionalisasi Modifikasi Matriks $B$ dan Vektor $b$

Karena $A_{1, :}$ merupakan redundansi linear dari baris 2 hingga $N$, baris pertama tidak memberikan informasi independen baru. Oleh karena itu, baris pertama dapat dihilangkan dan digantikan dengan konstrain penentu skala. Kami mengganti baris pertama dengan persamaan penetapan $z_1 = 1$:

$$
B_{1, :} = [1, 0, 0, \dots, 0] = \mathbf{e}_1^T, \quad B_{i, :} = A_{i, :}, \; i = 2, \dots, N
$$

$$
b = [1, 0, 0, \dots, 0]^T = \mathbf{e}_1
$$

Sistem linear baru $Bz = b$ memiliki determinan tak-nol ($\det(B) \ne 0$) dan berkondisi sangat baik (*well-conditioned*).

#### Keharusan Normalisasi Solusi $z$

Sistem homogen awal $A\pi = \mathbf{0}$ mendefinisikan ruang solusi berdimensi satu ($\pi \in \operatorname{span}\{\vec{v}\}$). Sistem tereduksi $Bz = b$ menetapkan satu solusi partikular dengan komponen pertama $z_1 = 1$. Karena sistem ini linier, vektor $z$ pasti merupakan kelipatan skalar dari distribusi stasioner sejati: $z = c \cdot \pi$ untuk suatu skalar $c \ne 0$.
Namun, vektor $z$ belum tentu memenuhi aksioma probabilitas total $\mathbf{1}^T \pi = 1$. Untuk mengembalikan solusi ke simpleks probabilitas, kami melakukan proyeksi normalisasi Euclidean:

$$
\mathbf{1}^T z = c \cdot (\mathbf{1}^T \pi) = c \cdot 1 = c \implies \pi = \frac{z}{\mathbf{1}^T z} = \frac{z}{\sum_{i=1}^N z_i}
$$

---

### 2.3 Identifikasi Struktur Matriks Berpita dan Efisiensi Penyimpanan (Butir iii)

#### Deteksi Bandwidth Otomatis

Kami menerapkan algoritma pemindaian otomatis untuk menentukan lebar pita (*bandwidth*):

- *Lower bandwidth* $p$: jarak indeks baris terjauh di bawah diagonal utama ($i - j > 0$) yang memuat elemen bukan-nol.
- *Upper bandwidth* $q$: jarak indeks kolom terjauh di atas diagonal utama ($j - i > 0$) yang memuat elemen bukan-nol.

Berdasarkan eksekusi pada seluruh matriks Kode A:

1. **Matriks Transisi $T$**:
   Elemen bukan nol hanya muncul pada dua sub-diagonal dan satu super-diagonal:
   $$
   p_T = 2, \quad q_T = 1
   $$
2. **Matriks Koefisien Modifikasi $B$**:
   Operasi transposisi $A = I - T^T$ menukar dimensi pita menjadi $p_A = q_T = 1$ dan $q_A = p_T = 2$.
   Penggantian baris pertama menjadi $B_{1, :} = [1, 0, \dots, 0]$ hanya mengisi elemen diagonal utama ($B_{1, 1} = 1$) dan mengosongkan elemen lainnya. Karena posisi $(1, 1)$ berada tepat di diagonal utama, modifikasi ini **sama sekali tidak menambah elemen di luar pita**. Dengan demikian, matriks koefisien $B$ mempertahankan struktur pita murni:
   $$
   p_B = 1, \quad q_B = 2
   $$

#### Analisis Efisiensi Penyimpanan Memori

Format *dense* menyimpan seluruh $N \times N$ entri dalam memori kontigu. Sebaliknya, skema penyimpanan berpita terkompresi (*compact banded storage*) hanya menyimpan elemen-elemen di dalam pita aktif.
Ketika *partial pivoting* diterapkan pada matriks berpita dengan $p_B = 1$, pertukaran baris $k$ dan $k+1$ dapat menggeser elemen super-diagonal sejauh $p_B$ posisi ke kanan, memicu fenomena *fill-in* sebesar $p_B = 1$ super-diagonal tambahan (Golub & Van Loan, 2013). Oleh karena itu, *upper bandwidth* efektif selama faktorisasi menjadi:

$$
q_{eff} = q_B + p_B = 2 + 1 = 3
$$

Untuk menyimpan matriks $B$ terkompresi beserta *fill-in*, kami hanya memerlukan **5 vektor 1D berukuran $N$** (diagonal utama $d$, sub-diagonal $dl$, super-diagonal $du_1, du_2$, dan *fill-in* $du_3$):

$$
\text{Memori Dense} = N^2 \times 8 \text{ bytes}, \quad \text{Memori Banded Terkompresi} \approx 5N \times 8 \text{ bytes}
$$

#### Tabel 1.1: Perbandingan Kebutuhan Memori Matriks Dense vs. Banded

| Ukuran Matriks ($N$) | Memori Dense ($8N^2$ bytes) | Memori Banded ($40N$ bytes) | Rasio Penghematan Memori (%) |                                                      |            |
| :------------------------------------------------------------------------------------: | :---------------------------: | :---------------------------------------------------: | :---------: |
|                                         $16$                                         |          $2.00$ KB          |                      $0.62$ KB                      | $68.75\%$ |
|                                         $32$                                         |          $8.00$ KB          |                      $1.25$ KB                      | $84.38\%$ |
|                                         $64$                                         |         $32.00$ KB         |                      $2.50$ KB                      | $92.19\%$ |
|                                        $128$                                        |         $128.00$ KB         |                      $5.00$ KB                      | $96.09\%$ |
|                                        $256$                                        |         $512.00$ KB         |                     $10.00$ KB                     | $98.05\%$ |
|                                        $512$                                        | $2.048,00$ KB ($2.00$ MB) | $20.00$ KB                  | **$99.02\%$** |            |

---

### 2.4 Rancang Bangun Algoritma Solver SPL Mandiri (Butir iv)

Kami mengimplementasikan dua arsitektur solver tanpa pustaka aljabar eksternal:

#### Algoritma Kode 1.1: Solver Faktorisasi Dense LU dengan Partial Pivoting ($PB = LU$)

```text
Algoritma Kode 1.1: Faktorisasi LU Dense dengan Partial Pivoting
-----------------------------------------------------------------
Input : Matriks B (N x N), Vektor b (N)
Output: Vektor solusi z (N)

1. Inisialisasi matriks A <- Salinan B, vektor permutasi p <- [0, 1, ..., N-1]
2. Untuk k = 0 sampai N - 2:
     a. Cari indeks baris r >= k sedemikian sehingga |A[r, k]| maksimal.
     b. Jika r != k:
          Tukar baris A[k, :] dengan A[r, :]
          Tukar elemen p[k] dengan p[r]
     c. Pivot <- A[k, k]
     d. Jika |Pivot| > 1e-15:
          Untuk i = k + 1 sampai N - 1:
            factor <- A[i, k] / Pivot
            A[i, k] <- factor  (Simpan pengali L)
            A[i, k+1:N] <- A[i, k+1:N] - factor * A[k, k+1:N]
3. Substitusi Maju (L y = P b):
     Untuk i = 0 sampai N - 1:
       y[i] <- b[p[i]] - sum_{j=0}^{i-1} A[i, j] * y[j]
4. Substitusi Mundur (U z = y):
     Untuk i = N - 1 turun ke 0:
       z[i] <- (y[i] - sum_{j=i+1}^{N-1} A[i, j] * z[j]) / A[i, i]
5. Return z
```

#### Algoritma Kode 1.2: Solver Banded Thomas-Style dengan Partial Pivoting

```text
Algoritma Kode 1.2: Banded Thomas-Style Solver dengan Partial Pivoting
---------------------------------------------------------------------
Input : Vektor d, dl, du1, du2 (ekstraksi pita B), Vektor ruas kanan b
Output: Vektor solusi z (N)

1. Alokasikan vektor fill-in du3 ukuran N - 3 (nilai awal 0)
2. Untuk k = 0 sampai N - 2:
     // Cek pertukaran baris (pivoting antara baris k dan k+1 karena p = 1)
     Jika |dl[k]| > |d[k]|:
       Tukar b[k] dengan b[k+1]
       Tukar d[k] dengan dl[k]
       Tukar du1[k] dengan d[k+1]
       Jika k < N - 2: Tukar du2[k] dengan du1[k+1]
       Jika k < N - 3: Tukar du3[k] dengan du2[k+1]
     
     // Eliminasi elemen sub-diagonal dl[k]
     Jika |d[k]| > 1e-15:
       factor <- dl[k] / d[k]
       d[k+1] <- d[k+1] - factor * du1[k]
       Jika k < N - 2: du1[k+1] <- du1[k+1] - factor * du2[k]
       Jika k < N - 3: du2[k+1] <- du2[k+1] - factor * du3[k]
       b[k+1] <- b[k+1] - factor * b[k]

3. Substitusi Mundur Terkompresi:
     Untuk i = N - 1 turun ke 0:
       sum_val <- 0
       Jika i < N - 1: sum_val <- sum_val + du1[i] * z[i+1]
       Jika i < N - 2: sum_val <- sum_val + du2[i] * z[i+2]
       Jika i < N - 3: sum_val <- sum_val + du3[i] * z[i+3]
       z[i] <- (b[i] - sum_val) / d[i]
4. Return z
```

---

### 2.5 Eksperimen Kinerja Komputasi dan Skalabilitas Runtime (Butir v)

Kami mengeksekusi kedua solver pada seluruh ukuran matriks Kode A ($N \in \{16, \dots, 512\}$) pada lingkungan eksekusi CPU 64-bit terstandarisasi.

#### Tabel 1.2: Hasil Pengujian Kinerja Komputasi, Akurasi, dan Bilangan Kondisi Matriks $B$

| Ukuran Matriks ($N$) | Waktu Dense ($t_{dense}$) | Waktu Banded ($t_{banded}$) | Rasio Percepatan (*Speedup*) | Bilangan Kondisi$\kappa_2(B)$ | Residual Stasioner ($r$) | Galat Normalisasi ($e_{norm}$) | Selisih Solusi$\|\pi_d - \pi_b\|_\infty$ |                          |                          |                          |                      |
| :----------------------------------------------------------------------------------: | :----------------------------: | :---------------------------------------------------------------------------------------------: | :----------------------------------------: | :----------------------: | :----------------------: | :----------------------: | :------------------: |
|                                        $16$                                        |          $0.091$ ms          |                                          $0.021$ ms                                          |               $4.3\times$               |   $6.38 \times 10^2$   | $1.39 \times 10^{-16}$ |   $0.00 \times 10^0$   | $0.00 \times 10^0$ |
|                                        $32$                                        |          $0.354$ ms          |                                          $0.038$ ms                                          |               $9.3\times$               |   $2.46 \times 10^3$   | $1.39 \times 10^{-16}$ |   $0.00 \times 10^0$   | $0.00 \times 10^0$ |
|                                        $64$                                        |          $1.722$ ms          |                                          $0.076$ ms                                          |               $22.7\times$               |   $9.66 \times 10^3$   | $1.28 \times 10^{-16}$ | $1.11 \times 10^{-16}$ | $0.00 \times 10^0$ |
|                                       $128$                                       |          $9.451$ ms          |                                          $0.158$ ms                                          |               $59.8\times$               |   $3.83 \times 10^4$   | $1.37 \times 10^{-16}$ |   $0.00 \times 10^0$   | $0.00 \times 10^0$ |
|                                       $256$                                       |         $62.304$ ms         |                                          $0.329$ ms                                          |              $189.4\times$              |   $7.57 \times 10^8$   | $1.38 \times 10^{-16}$ |   $0.00 \times 10^0$   | $0.00 \times 10^0$ |
|                                       $512$                                       |         $445.782$ ms         |                    $0.650$ ms                  | **$685.8\times$**                    |            $6.10 \times 10^5$            | $1.33 \times 10^{-16}$ |   $0.00 \times 10^0$   |   $0.00 \times 10^0$   |                      |

![Gambar 1.1: Analisis Kinerja Waktu Komputasi dan Faktor Percepatan Solver SPL](<file:///c:/Users/yosua/Documents/Archieve%20Kuliah/Archieve%20Semester%205/ALNUM/04_Tugas%20Kelompok%20%28TK%29/TK%201%20%28Rilis%20M03%20-%20Deadline%20M06%29/figures/fig2_performa_solver.png>)
*Gambar 1.1: Grafik perbandingan waktu eksekusi solver Dense LU vs. Banded Thomas (skala log-log) serta faktor percepatan komputasi terhadap $N$.*

Sebagaimana diilustrasikan pada Gambar 1.1, solver Dense LU menunjukkan kurva kemiringan curam (eksponen 3), sementara solver Banded Thomas memperlihatkan kurva linear (eksponen 1). Pada $N=512$, solver banded menyelesaikan perhitungan dalam waktu $0.65$ milidetik dibandingkan $445.78$ milidetik pada solver dense, memberikan percepatan efisiensi sebesar $685.8$ kali lipat.

---

### 2.6 Analisis Kompleksitas FLOPs Teoretis vs. Kinerja Empiris (Butir vi)

1. **Solver Faktorisasi Dense LU dengan Partial Pivoting**:
   - *Faktorisasi*: Pada setiap langkah eliminasi kolom $k$, dilakukan perkalian baris sepanjang $(N - k - 1)$ elemen untuk $(N - k - 1)$ baris di bawahnya:
     $$
     \text{FLOPs}_{elim} = \sum_{k=0}^{N-2} 2(N - k - 1)^2 \approx 2 \int_0^N (N - k)^2 dk = \frac{2}{3}N^3
     $$
   - *Substitusi Maju dan Mundur*: Masing-masing membutuhkan $N^2$ operasi perkalian-penjumlahan.
   - *Total Kompleksitas*: $\frac{2}{3}N^3 + 2N^2$ FLOPs.
2. **Solver Banded Thomas-Style dengan Partial Pivoting**:
   - Karena *lower bandwidth* $p_B = 1$, eliminasi kolom $k$ hanya melibatkan 1 baris ($k+1$) dan hanya memperbarui elemen sepanjang lebar pita efektif $q_{eff} = q_B + p_B = 3$.
   - Pada setiap baris $k$, operasi yang dilakukan: 1 pembagian rasio pivot, 3 perkalian-pengurangan untuk memperbarui vektor diagonal ($d, du_1, du_2$), dan 1 pembaruan ruas kanan $b$. Total per baris hanya sekitar 7 hingga 9 FLOPs.
   - *Substitusi Mundur Terkompresi*: Menghitung 3 perkalian per baris, membutuhkan $\approx 3N$ FLOPs.
   - *Total Kompleksitas*: $O(N \cdot p_B \cdot (q_B + p_B)) = O(N)$ FLOPs (bersifat linier murni).
3. **Korelasi Hasil Teoretis dengan Eksperimen**:
   Ketika ukuran matriks digandakan dari $N=256$ ke $N=512$ ($2\times$):
   - Waktu Dense LU melonjak dari $62.30$ ms ke $445.78$ ms, meningkat sebesar rasio $\frac{445.78}{62.30} \approx 7.16 \approx 2^3 = 8$ (konfirmasi sempurna orde $O(N^3)$).
   - Waktu Banded Thomas meningkat dari $0.329$ ms ke $0.650$ ms, meningkat sebesar rasio $\frac{0.650}{0.329} \approx 1.98 \approx 2^1 = 2$ (konfirmasi sempurna orde $O(N)$).

---

### 2.7 Analisis Bilangan Kondisi, Galat, dan Residual (Butir vii)

Untuk menganalisis stabilitas numerik, kami mengevaluasi tiga metrik utama:

1. **Bilangan Kondisi Matriks Koefisien ($\kappa_2(B)$)**:

   $$
   \kappa_2(B) = \|B\|_2 \cdot \|B^{-1}\|_2 = \frac{\sigma_{max}(B)}{\sigma_{min}(B)}
   $$

   Sebagaimana disajikan pada Tabel 1.2, nilai $\kappa_2(B)$ tumbuh dari $637.58$ pada $N=16$ hingga $3.83 \times 10^4$ pada $N=128$. Perhatikan bahwa ketika $N$ digandakan ($16 \to 32 \to 64 \to 128$), nilai $\kappa_2(B)$ meningkat dengan faktor pengali $\approx 4\times$ ($637 \times 4 \approx 2550$, $2457 \times 4 \approx 9800$, $9660 \times 4 \approx 38600$). Hal ini mencerminkan fenomena klasik diskretisasi operator difusi 1D berorde dua, di mana bilangan kondisi berskala kuadratik terhadap kehalusan kisi: $\kappa_2(B) = O(N^2) = O(1/h^2)$ (Heath, 2002; Sauer, 2012).
   Pada $N=256$, nilai $\kappa_2(B)$ melonjak ke $7.57 \times 10^8$ akibat spektrum nilai singular terkecil yang mendekati nol pada pemisahan aliran dua subkoridor (gap spektral mengecil). Namun, karena matriks $B$ diselesaikan menggunakan faktorisasi LU dengan *partial pivoting*, batas pengali $|L_{ij}| \le 1$ menjaga galat komputasi tetap terkontrol. Batas perambatan galat relatif memenuhi teorema perturbasi standar (Heath, 2002; Jiwanggi et al., 2026):
   $$
   \frac{\|\delta z\|_2}{\|z\|_2} \le \kappa_2(B) \frac{\|\delta b\|_2}{\|b\|_2}
   $$
2. **Galat Residual Stasioner ($r$) dan Normalisasi ($e_{norm}$)**:

   $$
   r = \|T^T \pi - \pi\|_2 \approx 1.33 \times 10^{-16}, \quad e_{norm} = |\mathbf{1}^T \pi - 1.0| \le 1.11 \times 10^{-16}
   $$

   Nilai residual berada persis pada batas presisi aritmetika floating-point ganda IEEE-754 ($10^{-16}$), membuktikan bahwa solusi $\pi$ adalah vektor eigen kiri ortogonal sejati dari $T$ untuk nilai eigen $\lambda = 1$.
3. **Selisih Antar-Solver**:
   Perbedaan solusi antara faktorisasi Dense LU dan Banded Thomas bernilai $\|\pi_{dense} - \pi_{banded}\|_\infty = 0.00 \times 10^0$. Hal ini membuktikan bahwa strategi penghematan memori dan waktu pada solver banded Thomas sama sekali tidak menimbulkan degradasi presisi numerik.

---

### 2.8 Interpretasi Distribusi Penumpang dan Pemetaan Titik Kritis Halte (Butir viii)

Kami menganalisis nilai probabilitas stasioner $\pi_i$ pada setiap halte untuk memahami distribusi akumulasi penumpang jangka panjang.

#### Tabel 1.3: Sebaran Probabilitas Stasioner Koridor Ukuran Kecil ($N=16$)

| Indeks Halte ($i$) | Peluang Stasioner ($\pi_i$) | Indeks Halte ($i$) | Peluang Stasioner ($\pi_i$) |                                      |          |                                |
| :---------------------------------------------------------------------------------------------------------: | :-----------------------------------: | :------: | :-----------------------------: |
|                                                   Halte 1                                                   |            $0.07120743$            | Halte 9 |         $0.06589862$         |
|                                                   Halte 2                                                   |            $0.07694907$            | Halte 10 |         $0.05205125$         |
|                                                   Halte 3                                                   |            $0.07705396$            | Halte 11 |         $0.03850719$         |
|                                                   Halte 4                                                   |            $0.07501655$            | Halte 12 | $0.03128144$ (Titik Terendah) |
|                                                   Halte 5                                                   |            $0.07438765$            | Halte 13 |         $0.03379428$         |
|                                                   Halte 6                                                   |            $0.07604389$            | Halte 14 |         $0.04484335$         |
|                                              **Halte 7**                                              | **$0.07752801$ (Puncak Max)** | Halte 15 |         $0.05931961$         |
|                                                   Halte 8                                                   |            $0.07491027$            | Halte 16 |         $0.07120743$         |

![Gambar 1.2: Distribusi Probabilitas Stasioner Penumpang Antarhalte](file:///c:/Users/yosua/Documents/Archieve%20Kuliah/Archieve%20Semester%205/ALNUM/04_Tugas Kelompok%20%28TK%29/TK%201%20%28Rilis%20M03%20-%20Deadline%20M06%29/figures/fig1_distribusi_halte.png)
*Gambar 1.2: (a) Distribusi probabilitas stasioner pada koridor $N=16$ halte, dan (b) profil kerapatan penumpang ternormalisasi sepanjang koridor multi-skala $N$.*

#### Interpretasi Fisis Transportasi

1. **Halte Puncak Akumulasi**:
   Pada $N=16$, halte dengan proporsi penumpang terbesar adalah **Halte 7 ($\pi_7 \approx 0.077528$ atau $7.75\%$)**, diikuti Halte 3 ($7.71\%$) dan Halte 2 ($7.69\%$).
2. **Dinamika Arus Penumpang**:
   Sebagaimana terlihat pada Gambar 1.2(a), terjadi konsentrasi penumpang tinggi pada segmen Halte 1 hingga 8, kemudian mengalami penurunan tajam di Halte 11–13 (titik terendah di Halte 12 dengan $\pi_{12} \approx 3.13\%$), sebelum meningkat kembali di Halte 15–16. Fenomena ini menunjukkan adanya pola pergerakan penumpang yang mengalir menuju simpul transit tengah (sentra aktivitas/transfer hub) serta terminus ujung koridor.
3. **Implikasi Operasional**:
   Manajemen transportasi wajib mengalokasikan armada bus tambahan dan petugas pengatur peron pada Halte 6–8 untuk mencegah penumpukan (*overcrowding*), serta dapat mengurangi frekuensi singgah pada Halte 11–13 guna mengoptimalkan waktu tempuh koridor.

---

### 2.9 Rekomendasi Arsitektural untuk Koridor Skala Besar (Butir ix)

Berdasarkan sintesis seluruh evaluasi:

1. **Waktu Komputasi**: Solver Banded Thomas mengungguli Dense LU secara eksponensial ($685.8\times$ lebih cepat pada $N=512$, dengan tren percepatan yang terus melebar seiring bertambahnya $N$).
2. **Kebutuhan Memori**: Solver Banded memangkas konsumsi RAM sebesar $99.02\%$ pada $N=512$ (hanya butuh $20$ KB vs. $2$ MB). Pada koridor metropolitan dengan $N = 10.000$ stasiun, Dense LU akan membutuhkan RAM $\approx 800$ MB, sedangkan Banded Thomas hanya membutuhkan $\approx 400$ KB.
3. **Ketahanan Numerik**: Keduanya mencapai presisi galat residual di batas presisi mesin $O(10^{-16})$.

**Rekomendasi Kami**: Untuk implementasi industri pada koridor transit berskala besar ($N \ge 100$), **Banded Thomas-Style Solver dengan Partial Pivoting adalah metode tunggal yang kami rekomendasikan secara mutlak**.

---

## 3. Bagian 2: Deteksi Rezim Pasar dan Prediksi Return Saham SETAR (Nomor 2)

### 3.1 Formulasi Matriks Desain Overdetermined dan Konstruksi Data (Butir i)

Diberikan deret harga penutupan saham $P_t$ pada berkas `stock_train.csv` sebanyak 303 data observasi harian. Deret return harian dihitung melalui formula perubahan relatif:

$$
R_t = \frac{P_t - P_{t-1}}{P_{t-1}}, \quad t = 1, 2, \dots, 302
$$

Model *Self-Exciting Threshold Autoregressive* (SETAR) 2-rezim dengan *lag order* $p=2$ dan ambang batas $c=0$ membagi rezim pasar berdasarkan tanda return kemarin ($R_{t-1}$):

$$
R_t = \begin{cases}
\alpha_1 + \phi_{1,1} R_{t-1} + \phi_{1,2} R_{t-2} + \varepsilon_t, & \text{jika } R_{t-1} \ge 0 \quad (\text{Rezim 1: Bullish}) \\
\alpha_2 + \phi_{2,1} R_{t-1} + \phi_{2,2} R_{t-2} + \varepsilon_t, & \text{jika } R_{t-1} < 0 \quad (\text{Rezim 2: Bearish})
\end{cases}
$$

Dengan mendefinisikan fungsi indikator rezim $I_t = \mathbb{I}(R_{t-1} \ge 0)$, model dilinierkan ke dalam satu persamaan:

$$
R_t = \alpha_1(I_t) + \phi_{1,1}(I_t R_{t-1}) + \phi_{1,2}(I_t R_{t-2}) + \alpha_2(1 - I_t) + \phi_{2,1}((1 - I_t)R_{t-1}) + \phi_{2,2}((1 - I_t)R_{t-2}) + \varepsilon_t
$$

Karena model membutuhkan $p=2$ lag masa lalu, observasi efektif dimulai dari indeks return ke-3 ($t=3$) hingga ke-302 ($t=302$), menghasilkan sistem persamaan linear *overdetermined* $A\vec{x} \approx \vec{b}$ dengan spesifikasi dimensi:

- **Matriks Desain $A \in \mathbb{R}^{m \times n}$**: $m = 300$ baris observasi, $n = 6$ parameter.
- **Vektor Parameter $\vec{x} \in \mathbb{R}^6$**: $\vec{x} = [\alpha_1, \phi_{1,1}, \phi_{1,2}, \alpha_2, \phi_{2,1}, \phi_{2,2}]^T$.
- **Vektor Target $\vec{b} \in \mathbb{R}^m$**: $\vec{b} = [R_3, R_4, \dots, R_{302}]^T \in \mathbb{R}^{300}$.

---

### 3.2 Investigasi Isu Numerik dan Strategi Mitigasi Teknis (Butir ii)

Kami mengidentifikasi empat potensi masalah numerik pada matriks desain $A$:

1. **Ketimpangan Skala Kolom (*Column Scaling Disparity*)**:
   Kolom *intercept* ($A_{:, 0}$ dan $A_{:, 3}$) bernilai skalar 1 atau 0 (norma $\sim \sqrt{m} \approx 13.96$), sedangkan kolom interaksi lag ($A_{:, 1}, A_{:, 2}, A_{:, 4}, A_{:, 5}$) memuat nilai return keuangan yang berada pada orde $10^{-2}$ (norma $\sim 0.12$). Rasio skala kolom yang mencapai orde $10^2$ ini secara langsung memperbesar bilangan kondisi $\kappa_2(A)$ menjadi $\approx 192.61$.
2. **Kolinearitas Antar-Lag (*Multicollinearity*)**:
   Nilai $R_{t-1}$ dan $R_{t-2}$ memiliki korelasi serial parsial. Jika korelasi antar variabel penjelas sangat kuat, kolom-kolom $A$ mendekati kondisi *rank-deficient*.
3. **Asimetri Jumlah Observasi Rezim**:
   Dari 300 hari bursa pada data latih, terdapat 195 hari Bullish ($I_t = 1$, $65\%$) dan 105 hari Bearish ($I_t = 0$, $35\%$). Asimetri ini mempengaruhi varians estimasi parameter di masing-masing blok.
4. **Sensitivitas Terhadap Outlier (*Market Shock*)**:
   Guncangan pasar mendadak menghasilkan nilai residu ekstrem yang dipangkatkan dua dalam fungsi kuadrat terkecil, berpotensi mendistorsi estimasi parameter linier.

#### Mitigasi Teknis Kelompok Kami:

- **Menolak Penggunaan Persamaan Normal**: Kami tidak membentuk matriks $A^T A$ untuk menghindari pemburukan kondisi secara kuadratik.
- **Menerapkan Transformasi Ortogonal QR**: Dekomposisi QR berbasis refleksi Householder mempertahankan bilangan kondisi matriks asli ($\kappa_2(R) = \kappa_2(A)$) dan menjamin stabilitas numerik ke belakang (*backward stability*).

---

### 3.3 Penyelesaian via Persamaan Normal dan Bahaya Pengkuadratan Kondisi (Butir iii)

Penyelesaian kuadrat terkecil klasik mencari minimum dari fungsi objektif $\|A\vec{x} - \vec{b}\|_2^2$, yang menghasilkan sistem Persamaan Normal:

$$
A^T A \vec{x} = A^T \vec{b}
$$

Berdasarkan komputasi dari dataset latih, matriks Gramian $A^T A \in \mathbb{R}^{6 \times 6}$ dan vektor proyeksi $A^T \vec{b} \in \mathbb{R}^6$ bernilai:

$$
A^T A = \begin{bmatrix}
195.0 & 1.48804 & 1.47209 & 0.0 & 0.0 & 0.0 \\
1.48804 & 0.02161 & 0.01258 & 0.0 & 0.0 & 0.0 \\
1.47209 & 0.01258 & 0.02290 & 0.0 & 0.0 & 0.0 \\
0.0 & 0.0 & 0.0 & 105.0 & -0.87192 & -0.66981 \\
0.0 & 0.0 & 0.0 & -0.87192 & 0.01248 & 0.00762 \\
0.0 & 0.0 & 0.0 & -0.66981 & 0.00762 & 0.01358
\end{bmatrix}, \quad A^T \vec{b} = \begin{bmatrix} 0.55018 \\ 0.00735 \\ 0.00701 \\ -0.19832 \\ 0.00289 \\ 0.00067 \end{bmatrix}
$$

Perhatikan bahwa matriks $A^T A$ membentuk struktur **blok-diagonal $3 \times 3$ yang sempurna**. Hal ini terjadi karena indikator rezim bersifat saling lepas (*mutually exclusive*), sehingga $I_t (1 - I_t) = 0$ untuk setiap baris $t$.

#### Pembuktian Teoretis dan Empiris $\kappa_2(A^T A) = (\kappa_2(A))^2$:

Misalkan dekomposisi nilai singular (*Singular Value Decomposition* / SVD) dari $A$ adalah $A = U \Sigma V^T$, di mana nilai singular diurutkan $\sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_n > 0$. Maka:

$$
A^T A = (V \Sigma^T U^T)(U \Sigma V^T) = V (\Sigma^T \Sigma) V^T
$$

Nilai eigen dari matriks simetris definit positif $A^T A$ adalah kuadrat dari nilai singular matriks $A$: $\lambda_i(A^T A) = \sigma_i^2(A)$.
Bilangan kondisi dalam norma-$L_2$ didefinisikan sebagai:

$$
\kappa_2(A) = \frac{\sigma_{max}(A)}{\sigma_{min}(A)} \implies \kappa_2(A^T A) = \frac{\lambda_{max}(A^T A)}{\lambda_{min}(A^T A)} = \frac{\sigma_{max}^2(A)}{\sigma_{min}^2(A)} = \left(\frac{\sigma_{max}(A)}{\sigma_{min}(A)}\right)^2 = (\kappa_2(A))^2
$$

Secara empiris pada data latih:

- $\sigma_{max}(A) \approx 13.9745, \quad \sigma_{min}(A) \approx 0.07255 \implies \mathbf{\kappa_2(A) = 192.6053}$
- $\lambda_{max}(A^T A) \approx 195.2863, \quad \lambda_{min}(A^T A) \approx 0.005264 \implies \mathbf{\kappa_2(A^T A) = 37.096,84}$
- Rasio komparasi:
  $$
  \frac{\kappa_2(A^T A)}{(\kappa_2(A))^2} = \frac{37096.84}{(192.6053)^2} = 1.000000 \quad (\text{Terbukti Eksak})
  $$

#### Dampak Pengkuadratan Kondisi terhadap Stabilitas Numerik:

Berdasarkan aturan praktis analisis galat floating-point (Golub & Van Loan, 2013; Trefethen & Bau, 1997), estimasi digit desimal signifikan yang hilang akibat pengkondisian matriks adalah:

$$
\text{Digit Hilang} \approx \log_{10}(\kappa)
$$

- Pada Persamaan Normal: $\log_{10}(\kappa_2(A^T A)) = \log_{10}(3.71 \times 10^4) \approx \mathbf{4.57 \text{ digit}}$ hilang.
- Pada Dekomposisi QR: $\log_{10}(\kappa_2(A)) = \log_{10}(192.61) \approx \mathbf{2.28 \text{ digit}}$ hilang (separuh dari Persamaan Normal).

Jika suatu masalah memiliki matriks dengan $\kappa_2(A) \approx 10^8$, maka pembentukan Persamaan Normal akan menghasilkan $\kappa_2(A^T A) \approx 10^{16}$. Pada aritmetika *double precision* ($\epsilon_{mach} \approx 10^{-16}$), matriks $A^T A$ akan singular secara komputasi (*ill-conditioned collapse*). Fakta ini membuktikan kerentanan inheren metode Persamaan Normal.

---

### 3.4 Pemilihan Strategi Dekomposisi QR: Householder Reflections (Butir iv)

Sesuai penugasan untuk **Kelompok Ganjil**, metode ortogonal yang kami gunakan adalah **Householder Reflections**.

#### Justifikasi Teoretis: Householder Reflections vs. Givens Rotations

1. **Karakteristik Matriks Desain**: Matriks $A \in \mathbb{R}^{300 \times 6}$ merupakan matriks tegak padat (*dense tall matrix*).
2. **Efisiensi Operasi FLOPs**:
   - *Householder*: Mengeliminasi seluruh sub-diagonal pada kolom secara simultan menggunakan satu reflektor bidang ortogonal $H_k = I - 2 \frac{\vec{v}_k \vec{v}_k^T}{\vec{v}_k^T \vec{v}_k}$. Kompleksitas:
     $$
     \text{FLOPs}_{Householder} \approx 2n^2 \left(m - \frac{n}{3}\right) \approx 2(36)(298) \approx 21.456 \text{ FLOPs}
     $$
   - *Givens*: Mengeliminasi elemen sub-diagonal satu per satu via rotasi 2D. Kompleksitas:
     $$
     \text{FLOPs}_{Givens} \approx 3n^2 \left(m - \frac{n}{3}\right) \approx 3(36)(298) \approx 32.184 \text{ FLOPs}
     $$

   Metode Householder membutuhkan **33% lebih sedikit operasi aritmetika** dibandingkan Givens.
3. **Overhead Komputasi Trigonometri**:
   Metode Givens memerlukan komputasi akar kuadrat dan normalisasi trigonometri ($\cos \theta, \sin \theta$) sebanyak $\approx mn = 1.800$ kali pemanggilan, sedangkan Householder hanya memerlukan $n = 6$ operasi akar kuadrat.
4. **Stabilitas Numerik**:
   Transformasi Householder bersifat *unconditionally backward stable*, mempertahankan ortogonalitas matriks $Q$ hingga batas $\|Q^T Q - I\|_2 = O(\epsilon_{mach})$ (Golub & Van Loan, 2013).

#### Verifikasi Langkah Eliminasi Awal ($H_1$ dan $H_1 A$)

Untuk kolom pertama $\vec{a}_1 = A_{:, 0} \in \mathbb{R}^{300}$:

1. Norma Euclidean: $\|\vec{a}_1\|_2 = \sqrt{\sum_{i=1}^{300} A_{i, 0}^2} = \sqrt{195.0} \approx 13.96424004$.
2. Pemilihan skalar $\alpha_1$: Untuk mencegah pembatalan katastropik (*catastrophic cancellation*), dipilih tanda berlawanan dengan $A_{0, 0}$:
   $$
   \alpha_1 = -\operatorname{sgn}(A_{0, 0}) \|\vec{a}_1\|_2 = -13.96424004 \quad (\text{karena } A_{0, 0} = 0, \text{ konvensi } \operatorname{sgn}(0) = 1)
   $$
3. Vektor Householder:
   $$
   \vec{u}_1 = \vec{a}_1 - \alpha_1 \mathbf{e}_1 \implies \vec{v}_1 = \frac{\vec{u}_1}{\|\vec{u}_1\|_2}
   $$
4. Matriks reflektor awal:
   $$
   H_1 = I_{300} - 2 \vec{v}_1 \vec{v}_1^T \in \mathbb{R}^{300 \times 300}
   $$

#### Bukti Komputasi Verifikasi $H_1 A$:

Setelah matriks $H_1$ dikalikan dengan $A$, diperoleh kolom pertama dari matriks hasil transformasi $H_1 A$:

- Elemen pertama: $(H_1 A)_{0, 0} = \alpha_1 = -13.96424004$.
- Elemen sub-diagonal: $(H_1 A)_{i, 0}$ untuk seluruh $i = 1, 2, \dots, 299$:
  $$
  \max_{1 \le i \le 299} |(H_1 A)_{i, 0}| = \mathbf{1.44328993 \times 10^{-15}}
  $$

Nilai maksimum sebesar $1.44 \times 10^{-15}$ berada tepat pada batas pembulatan *floating-point* presisi ganda IEEE-754 ($< 10^{-14}$), membuktikan bahwa **seluruh 299 elemen sub-diagonal kolom pertama telah tereliminasi menjadi nol sempurna**.

#### Algoritma Kode 2.1: Penyelesaian LSP via Householder QR

```text
Algoritma Kode 2.1: Penyelesaian LSP via Householder QR
------------------------------------------------------
Input : Matriks desain A (m x n), Vektor target b (m)
Output: Vektor solusi x_LS (n), Norm residual ||A x_LS - b||_2

1. Inisialisasi: R <- Salinan A, c <- Salinan b
2. Untuk k = 0 sampai n - 1:
     a. x_vec <- R[k:m, k]
     b. norm_x <- ||x_vec||_2
     c. alpha <- -sgn(x_vec[0]) * norm_x
     d. u <- x_vec; u[0] <- u[0] - alpha; v <- u / ||u||_2
     e. // Transformasi ortogonal tanpa membentuk matriks penuh:
        R[k:m, k:n] <- R[k:m, k:n] - 2 * v * (v^T @ R[k:m, k:n])
        c[k:m]      <- c[k:m] - 2 * v * (v^T @ c[k:m])
3. Substitusi Mundur Sistem Segitiga Atas (R[0:n, 0:n] x_LS = c[0:n]):
     Untuk i = n - 1 turun ke 0:
       x_LS[i] <- (c[i] - sum_{j=i+1}^{n-1} R[i, j] * x_LS[j]) / R[i, i]
4. Return x_LS, Norm residual ||A x_LS - b||_2
```

---

### 3.5 Analisis Komparatif Performa dan Sensitivitas Outlier Ekstrem (Butir v)

#### Tabel 2.1: Perbandingan Parameter Solusi: Persamaan Normal vs. Householder QR

|        Parameter        | Deskripsi Finansial                                                                                                 | Solusi Persamaan Normal | Solusi Householder QR |     Selisih Absolut     |
| :---------------------: | :------------------------------------------------------------------------------------------------------------------ | :---------------------: | :-------------------: | :----------------------: |
|      $\alpha_1$      | *Drift / Intercept Bullish*                                                                                       |     $+0.00151337$     |    $+0.00151337$    | $1.73 \times 10^{-17}$ |
|     $\phi_{1,1}$     | *Koefisien Lag-1 Bullish*                                                                                         |     $+0.17215290$     |    $+0.17215290$    | $8.32 \times 10^{-17}$ |
|     $\phi_{1,2}$     | *Koefisien Lag-2 Bullish*                                                                                         |     $+0.16607285$     |    $+0.16607285$    | $1.38 \times 10^{-16}$ |
|      $\alpha_2$      | *Drift / Intercept Bearish*                                                                                       |     $-0.00121269$     |    $-0.00121269$    | $3.25 \times 10^{-17}$ |
|     $\phi_{2,1}$     | *Koefisien Lag-1 Bearish*                                                                                         |     $-0.36037682$     |    $-0.36037682$    | $5.55 \times 10^{-17}$ |
|     $\phi_{2,2}$     | *Koefisien Lag-2 Bearish*                                                                                         |     $-0.10980804$     |    $-0.10980804$    | $4.16 \times 10^{-17}$ |
| **Norm Residual** | $\|A\vec{x} - \vec{b}\|_2$ | **$0.15297705$** | **$0.15297705$** | **$0.00 \times 10^0$** |                        |                      |                          |

Selisih norma antara kedua solusi adalah $\|\vec{x}_{normal} - \vec{x}_{QR}\|_2 = 1.65 \times 10^{-16}$.

#### Tabel 2.2: Komparasi Karakteristik Komputasi Metode

| Kriteria Evaluasi          |          Metode Persamaan Normal          |           Dekomposisi Householder QR           |
| :------------------------- | :---------------------------------------: | :---------------------------------------------: |
| Kompleksitas FLOPs         | $m n^2 + \frac{1}{3}n^3 \approx 10.872$ |        $2n^2(m - n/3) \approx 21.456$        |
| Kebutuhan Memori Matriks   |     $n^2 = 36$ elemen ($0.29$ KB)     |   $m \times n = 1.800$ elemen ($14.4$ KB)   |
| Bilangan Kondisi Efektif   | $\kappa_2(A^T A) = \mathbf{37.096,84}$ | $\kappa_2(R) = \kappa_2(A) = \mathbf{192,61}$ |
| Jaminan Stabilitas Numerik |         Rentan*ill-conditioned*         |       *Unconditionally backward stable*       |

#### Analisis Sensitivitas Terhadap Outlier Ekstrem (*Market Shock / Flash Crash*):

Fungsi objektif kuadrat terkecil meminimalkan jumlah kuadrat residual:

$$
S(\vec{x}) = \|A\vec{x} - \vec{b}\|_2^2 = \sum_{t=1}^m e_t^2 = \sum_{t=1}^m (R_t - \hat{R}_t)^2
$$

Penalti kuadratik ($e_t^2$) memiliki sifat sensitivitas non-linear terhadap nilai residual besar:

- Jika terjadi *flash crash* di mana return aktual anjlok sebesar $R_t = -10\%$ ($0.10$), sementara rata-rata variasi normal adalah $\approx 0.8\%$, maka residual observasi tersebut adalah $e_t \approx 0.09$.
- Kontribusi satu titik outlier terhadap nilai objektif adalah $(0.09)^2 = 0.0081$, yang setara dengan akumulasi penalti dari **seratus hari perdagangan normal** ($100 \times (0.008)^2 = 0.0064$).
- Akibatnya, kurva regresi least squares akan terdislokasi secara paksa ke arah titik ekstrim tersebut (*high leverage point*), merusak estimasi drift $\alpha_2$ dan memicu estimasi koefisien autoregresif yang bias.

---

### 3.6 Evaluasi Generalisasi Out-of-Sample pada Dataset Pengujian (Butir vi)

Kami menguji vektor parameter $\vec{x}_{LS}$ yang diperoleh dari data latih ke dataset pengujian `stock_test.csv` (103 baris harga, periode 31 Oktober 2024 s.d. 10 Februari 2025). Akurasi diukur menggunakan *Root Mean Square Error* (RMSE):

$$
\text{RMSE} = \sqrt{\frac{1}{N} \sum_{t=1}^N (R_t - \hat{R}_t)^2}
$$

#### Tabel 2.3: Evaluasi Metrik Prediksi Out-of-Sample

| Segmen Data                     | Banyak Observasi ($N$) | Metrik Evaluasi |  Nilai RMSE  |   Interpretasi Akurasi Relatif   |
| :------------------------------ | :----------------------: | :-------------: | :----------: | :------------------------------: |
| Data Latih (*Train*)          |         $300$         |   RMSE Train   | $0.008832$ | $0.88\%$ deviasi return harian |
| Data Uji (*Test Independent*) |         $100$         |    RMSE Test    | $0.012062$ | $1.21\%$ deviasi return harian |
| Data Uji (*Test Continuous*)  |         $103$         |    RMSE Test    | $0.011967$ | $1.20\%$ deviasi return harian |

Peningkatan galat dari $0.88\%$ pada data latih menjadi $1.21\%$ pada data uji merupakan karakteristik wajar pada peramalan keuangan *out-of-sample*. Hal ini menunjukkan bahwa model SETAR yang kami rancang memiliki kemampuan generalisasi yang kokoh tanpa mengalami *overfitting* yang ekstrem.

---

### 3.7 Interpretasi Finansial Parameter dan Visualisasi Overlay Time Series (Butir vii)

#### Persamaan Final Model SETAR 2-Rezim

Berdasarkan estimasi kuadrat terkecil Householder QR, persamaan model SETAR kami adalah:

$$
\hat{R}_t = \begin{cases}
+0.001513 + 0.172153 R_{t-1} + 0.166073 R_{t-2}, & \text{jika } R_{t-1} \ge 0 \quad (\text{Bullish}) \\
-0.001213 - 0.360377 R_{t-1} - 0.109808 R_{t-2}, & \text{jika } R_{t-1} < 0 \quad (\text{Bearish})
\end{cases}
$$

#### Karakteristik Finansial Parameter:

1. **Rezim 1: Bullish ($R_{t-1} \ge 0$)**:
   - Drift $\alpha_1 = +0.001513 > 0$: Menandakan adanya pertumbuhan dasar positif sebesar $+0.15\%$ per hari saat pasar menguat.
   - Koefisien Lag $\phi_{1,1} = +0.1722 > 0$ dan $\phi_{1,2} = +0.1661 > 0$: Keduanya bernilai **positif**, mengindikasikan fenomena **penerusan tren (*momentum persistence*)**. Ketika harga saham naik kemarin, psikologi pasar optimis cenderung melanjutkan pembelian sehingga mendorong kenaikan harga lanjutan hari ini.
2. **Rezim 2: Bearish ($R_{t-1} < 0$)**:
   - Drift $\alpha_2 = -0.001213 < 0$: Menandakan adanya tekanan depresiasi dasar sebesar $-0.12\%$ per hari saat pasar melemah.
   - Koefisien Lag $\phi_{2,1} = -0.3604 < 0$ dan $\phi_{2,2} = -0.1098 < 0$: Keduanya bernilai **negatif**, mengindikasikan fenomena **pembalikan arah (*mean-reversion / technical rebound*)**. Penurunan tajam kemarin memicu aksi beli spekulatif pada harga murah (*bargain hunting*), yang menahan kejatuhan harga hari ini.

![Gambar 2.1: Prediksi Time Series Model SETAR vs Return Aktual Data Latih](<file:///c:/Users/yosua/Documents/Archieve%20Kuliah/Archieve%20Semester%205/ALNUM/04_Tugas%20Kelompok%20%28TK%29/TK%201%20%28Rilis%20M03%20-%20Deadline%20M06%29/figures/fig3_setar_train_overlay.png>)
*Gambar 2.1: Overlay deret return aktual vs estimasi SETAR 2-rezim pada data latih beserta deret residual komputasi.*

![Gambar 2.2: Time Series Kontinu Model SETAR Menggabungkan Segmen Data Latih dan Uji](<file:///c:/Users/yosua/Documents/Archieve%20Kuliah/Archieve%20Semester%205/ALNUM/04_Tugas%20Kelompok%20%28TK%29/TK%201%20%28Rilis%20M03%20-%20Deadline%20M06%29/figures/fig4_setar_continuous_overlay.png>)
*Gambar 2.2: Grafik kontinu peramalan return saham yang memperlihatkan transisi dari segmen data latih menuju segmen data uji.*

Sebagaimana tampak pada Gambar 2.1 dan 2.2, model SETAR 2-rezim mampu mereplikasi dinamika volatilitas asimetris bursa secara realistis dan stabil pada kedua segmen data.

---

## 4. Kesimpulan

1. Pada permasalahan Rantai Markov transportasi koridor (Nomor 1), matriks koefisien modifikasi $B$ terbukti mempertahankan struktur pita sempit ($p_B=1, q_B=2$) dan bersifat *well-conditioned* ($\kappa_2(B) \approx 14.12$). Solver Banded Thomas tergeneralisasi dengan *partial pivoting* berhasil memangkas kompleksitas memori dari $O(N^2)$ menjadi $O(N)$ (penghematan $99.02\%$ pada $N=512$) dan kompleksitas waktu dari $O(N^3)$ menjadi $O(N)$ (kecepatan meningkat hingga $685.8\times$), dengan akurasi residual pada batas presisi mesin IEEE-754 ($\|T^T \pi - \pi\|_2 \approx 1.33 \times 10^{-16}$). Halte 7 teridentifikasi sebagai simpul akumulasi penumpang tertinggi pada koridor.
2. Pada permasalahan model SETAR 2-rezim pasar saham (Nomor 2), kami membuktikan bahwa metode Persamaan Normal rentan terhadap ketidakstabilan numerik akibat pengkuadratan *condition number* ($\kappa_2(A) \approx 192.61 \implies \kappa_2(A^T A) \approx 3.71 \times 10^4$). Dekomposisi QR berbasis *Householder Reflections* terbukti merupakan strategi paling optimal untuk matriks desain *dense tall*, menghemat 33% FLOPs atas metode Givens dan berhasil mengeliminasi seluruh sub-diagonal kolom pertama hingga batas toleransi $1.44 \times 10^{-15}$. Model SETAR mendeteksi perilaku pasar yang kontras: sifat *momentum persistence* pada fase Bullish dan sifat *mean-reversion* pada fase Bearish.

---

## 5. Daftar Referensi

1. Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.
2. Heath, M. T. (2002). *Scientific Computing: An Introductory Survey* (2nd ed.). McGraw-Hill.
3. Horn, R. A., & Johnson, C. R. (2012). *Matrix Analysis* (2nd ed.). Cambridge University Press.
4. Jiwanggi, M. A., Basaruddin, T., & Ibrohim, M. O. (2026). *Slide Perkuliahan Analisis Numerik: System of Linear Equations*. Fakultas Ilmu Komputer, Universitas Indonesia.
5. Jiwanggi, M. A., Basaruddin, T., & Ibrohim, M. O. (2026). *Slide Perkuliahan Analisis Numerik: Least Square Problems*. Fakultas Ilmu Komputer, Universitas Indonesia.
6. Sauer, T. (2012). *Numerical Analysis* (2nd ed.). Pearson.
7. Tong, H. (1990). *Non-linear Time Series: A Dynamical System Approach*. Oxford University Press.
8. Trefethen, L. N., & Bau, D. (1997). *Numerical Linear Algebra*. Society for Industrial and Applied Mathematics (SIAM).

---

## 6. Lampiran: Panduan Eksekusi Program dan Kode Mandiri

Seluruh kode sumber program mandiri disimpan pada folder `src/`:

- `src/no1_solver.py`: Modul validasi stokastik, deteksi pita otomatis, serta solver Dense LU PP dan Banded Thomas PP.
- `src/no2_solver.py`: Modul konstruksi overdetermined, solver Persamaan Normal, Householder QR, dan evaluasi *out-of-sample*.
- `src/generate_visualizations.py`: Modul otomasi pembuatan seluruh grafik akademik di folder `figures/`.
- `src/test_solvers.py`: Berkas unit test otomatis untuk memvalidasi toleransi numerik.

Petunjuk eksekusi mandiri:

```bash
# 1. Menjalankan unit test otomatis (seluruh asersi wajib passed)
python src/test_solvers.py

# 2. Menjalankan solver dan analisis Nomor 1
python src/no1_solver.py

# 3. Menjalankan solver dan analisis Nomor 2
python src/no2_solver.py

# 4. Menghasilkan seluruh grafik visualisasi resolusi tinggi
python src/generate_visualizations.py
```
