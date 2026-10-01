% =========================================================================
% Program: no1_solver.m
% Deskripsi: Implementasi Solver SPL untuk Masalah Distribusi Steady State Antarhalte
%            Tugas Kelompok 1 Analisis Numerik - Fasilkom UI
%            Metode: Dense LU PP vs Banded Thomas-Style PP
% =========================================================================

function no1_solver()
    fprintf('=========================================================================\n');
    fprintf('  NOMOR 1: ANALISIS PERPINDAHAN PENUMPANG ANTARHALTE (SOLVER SPL)\n');
    fprintf('=========================================================================\n\n');

    % 1. Deteksi direktori dataset
    base_dir = fileparts(pwd);
    candidates = {
        fullfile(pwd, 'Nomor 1', 'A'), ...
        fullfile(pwd, 'Nomor 1', 'Nomor 1', 'A'), ...
        fullfile(base_dir, 'Nomor 1', 'A'), ...
        fullfile(base_dir, 'Nomor 1', 'Nomor 1', 'A'), ...
        fullfile(pwd, 'A')
    };
    
    data_dir = '';
    for k = 1:length(candidates)
        if exist(fullfile(candidates{k}, 'T_16.csv'), 'file')
            data_dir = candidates{k};
            break;
        end
    end
    
    if isempty(data_dir)
        error('Direktori dataset Nomor 1 (Folder A) tidak ditemukan. Pastikan file T_16.csv tersedia.');
    end
    fprintf('Dataset ditemukan di: %s\n\n', data_dir);

    sizes = [16, 32, 64, 128, 256, 512];
    results = struct('N', {}, 'cond_B', {}, 't_dense', {}, 't_banded', {}, ...
                     'res_dense', {}, 'res_banded', {}, 'enorm', {}, 'max_halte', {}, 'max_prob', {});

    for idx = 1:length(sizes)
        N = sizes(idx);
        csv_file = fullfile(data_dir, sprintf('T_%d.csv', N));
        T = readmatrix(csv_file);

        % i) Validasi Data & Sifat Stokastik
        val = validate_matrix(T);
        [p_T, q_T] = detect_bandwidth(T);

        % ii) Formulasi Matriks B dan vektor b
        [B, b] = construct_B_and_b(T);
        [p_B, q_B] = detect_bandwidth(B);

        % Hitung Condition Number kappa_2(B)
        cond_B = cond(B);

        % iv-a) Solver Dense LU dengan Partial Pivoting
        tic;
        z_dense = solve_dense_lu_pp(B, b);
        t_dense = toc * 1000; % ms
        pi_dense = normalize_solution(z_dense);

        % iv-b) Solver Banded Thomas-Style dengan Partial Pivoting
        tic;
        z_banded = solve_banded_thomas_pp(B, b);
        t_banded = toc * 1000; % ms
        pi_banded = normalize_solution(z_banded);

        % vii) Analisis Galat Residual dan Normalisasi
        [r_dense, e_dense] = compute_errors(T, pi_dense);
        [r_banded, e_banded] = compute_errors(T, pi_banded);
        diff_sol = max(abs(pi_dense - pi_banded));

        [max_prob, max_halte] = max(pi_banded);

        results(idx).N = N;
        results(idx).cond_B = cond_B;
        results(idx).t_dense = t_dense;
        results(idx).t_banded = t_banded;
        results(idx).res_dense = r_dense;
        results(idx).res_banded = r_banded;
        results(idx).enorm = e_banded;
        results(idx).max_halte = max_halte;
        results(idx).max_prob = max_prob;

        fprintf('-------------------------------------------------------------------------\n');
        fprintf('Ukuran Halte N = %d | Bandwidth T: (%d,%d) -> B: (%d,%d)\n', N, p_T, q_T, p_B, q_B);
        fprintf('Validasi: Stokastik=%d | Non-negatif=%d | Iredisibel=%d\n', ...
                val.is_stochastic, val.is_non_negative, val.is_irreducible);
        fprintf('Condition Number kappa(B)  : %12.2f\n', cond_B);
        fprintf('Runtime Dense LU PP        : %10.4f ms\n', t_dense);
        fprintf('Runtime Banded Thomas PP   : %10.4f ms (Percepatan: %.1fx)\n', t_banded, t_dense / t_banded);
        fprintf('Residual ||T^T*pi - pi||_2 : %12.4e\n', r_banded);
        fprintf('Galat Normalisasi |1^T*pi-1|: %12.4e\n', e_banded);
        fprintf('Perbedaan Solusi Dense-Band: %12.4e\n', diff_sol);
        fprintf('Halte dengan Proporsi Maks : Halte %d (pi = %.6f)\n', max_halte, max_prob);
    end

    fprintf('=========================================================================\n\n');

    % Plotting Distribusi Probabilitas Stasioner N = 16
    figure('Name', 'Distribusi Probabilitas Stasioner Halte', 'Position', [100, 100, 900, 400]);
    csv_16 = fullfile(data_dir, 'T_16.csv');
    T16 = readmatrix(csv_16);
    [B16, b16] = construct_B_and_b(T16);
    z16 = solve_banded_thomas_pp(B16, b16);
    pi16 = normalize_solution(z16);

    subplot(1, 2, 1);
    bar(1:16, pi16, 'FaceColor', [0.2, 0.5, 0.8], 'EdgeColor', 'k');
    hold on;
    [val_max, idx_max] = max(pi16);
    bar(idx_max, val_max, 'FaceColor', [0.85, 0.2, 0.2], 'EdgeColor', 'k');
    plot(1:16, pi16, '-k', 'LineWidth', 1.2);
    title('(a) Distribusi Probabilitas Stasioner (N=16 Halte)');
    xlabel('Nomor Halte (i)');
    ylabel('Peluang Stasioner (\pi_i)');
    grid on;

    subplot(1, 2, 2);
    plot([results.N], [results.t_dense], 'r-o', 'LineWidth', 1.8, 'DisplayName', 'Dense LU PP O(N^3)');
    hold on;
    plot([results.N], [results.t_banded], 'g-s', 'LineWidth', 1.8, 'DisplayName', 'Banded Thomas PP O(N)');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    title('(b) Runtime Komputasi (Skala Log-Log)');
    xlabel('Ukuran Matriks (N)');
    ylabel('Waktu Eksekusi (ms)');
    legend('Location', 'northwest');
    grid on;
end

% =========================================================================
% FUNGSI-FUNGSI SUB-ALGORITMA MANDIRI
% =========================================================================

function val = validate_matrix(T)
    [m, n] = size(T);
    val.is_square = (m == n);
    val.is_non_negative = all(T(:) >= -1e-15);
    row_sums = sum(T, 2);
    val.is_stochastic = max(abs(row_sums - 1.0)) < 1e-12;
    
    % Ketereduksian melalui konektivitas pita
    val.is_irreducible = true;
    for i = 1:(n - 1)
        if T(i+1, i) <= 0 && T(i, i+1) <= 0
            val.is_irreducible = false;
            break;
        end
    end
end

function [p, q] = detect_bandwidth(M, tol)
    if nargin < 2, tol = 1e-12; end
    n = size(M, 1);
    p = 0; q = 0;
    for i = 1:n
        for j = 1:n
            if abs(M(i, j)) > tol
                if (i - j) > p, p = i - j; end
                if (j - i) > q, q = j - i; end
            end
        end
    end
end

function [B, b] = construct_B_and_b(T)
    n = size(T, 1);
    A = eye(n) - T';
    B = A;
    B(1, :) = 0;
    B(1, 1) = 1.0;
    b = zeros(n, 1);
    b(1) = 1.0;
end

function z = solve_dense_lu_pp(B_in, b_in)
    % Faktorisasi Dense LU dengan Partial Pivoting (P B = L U)
    n = size(B_in, 1);
    A = B_in;
    b = b_in;
    p = 1:n;

    % 1. Faktorisasi / Eliminasi Maju
    for k = 1:(n - 1)
        [~, max_idx] = max(abs(A(k:n, k)));
        pivot_row = k + max_idx - 1;
        if pivot_row ~= k
            A([k, pivot_row], :) = A([pivot_row, k], :);
            p([k, pivot_row]) = p([pivot_row, k]);
        end

        pivot = A(k, k);
        if abs(pivot) < 1e-15, continue; end

        for i = (k + 1):n
            factor = A(i, k) / pivot;
            A(i, k) = factor;
            A(i, (k+1):n) = A(i, (k+1):n) - factor * A(k, (k+1):n);
        end
    end

    % 2. Forward Substitution: L y = P b
    b_perm = b(p);
    y = zeros(n, 1);
    for i = 1:n
        sum_val = 0;
        for j = 1:(i - 1)
            sum_val = sum_val + A(i, j) * y(j);
        end
        y(i) = b_perm(i) - sum_val;
    end

    % 3. Back Substitution: U z = y
    z = zeros(n, 1);
    for i = n:-1:1
        sum_val = 0;
        for j = (i + 1):n
            sum_val = sum_val + A(i, j) * z(j);
        end
        z(i) = (y(i) - sum_val) / A(i, i);
    end
end

function z = solve_banded_thomas_pp(B_in, b_in)
    % Solver Banded Thomas-Style dengan Partial Pivoting
    % Struktur B: Lower bandwidth p=1, Upper bandwidth q=2, fill-in q_eff=3
    % Menggunakan 5 vektor 1D (memori O(N))
    n = size(B_in, 1);
    d = diag(B_in);
    dl = zeros(n-1, 1);
    du1 = zeros(n-1, 1);
    du2 = zeros(n-2, 1);
    du3 = zeros(n-3, 1);

    for i = 1:(n - 1)
        dl(i) = B_in(i+1, i);
        du1(i) = B_in(i, i+1);
    end
    for i = 1:(n - 2)
        du2(i) = B_in(i, i+2);
    end

    b = b_in;

    % Eliminasi Maju dengan Partial Pivoting antar baris k dan k+1
    for k = 1:(n - 1)
        if abs(dl(k)) > abs(d(k))
            % Tukar baris b
            tmp_b = b(k); b(k) = b(k+1); b(k+1) = tmp_b;
            % Tukar elemen matriks
            tmp = d(k); d(k) = dl(k); dl(k) = tmp;
            tmp = du1(k); du1(k) = d(k+1); d(k+1) = tmp;
            if k <= n - 2
                tmp = du2(k); du2(k) = du1(k+1); du1(k+1) = tmp;
            end
            if k <= n - 3
                tmp = du3(k); du3(k) = du2(k+1); du2(k+1) = tmp;
            end
        end

        if abs(d(k)) > 1e-15
            factor = dl(k) / d(k);
            dl(k) = factor;

            d(k+1) = d(k+1) - factor * du1(k);
            if k <= n - 2
                du1(k+1) = du1(k+1) - factor * du2(k);
            end
            if k <= n - 3
                du2(k+1) = du2(k+1) - factor * du3(k);
            end

            b(k+1) = b(k+1) - factor * b(k);
        end
    end

    % Substitusi Mundur: U z = b
    z = zeros(n, 1);
    for i = n:-1:1
        sum_val = 0;
        if i <= n - 1, sum_val = sum_val + du1(i) * z(i+1); end
        if i <= n - 2, sum_val = sum_val + du2(i) * z(i+2); end
        if i <= n - 3, sum_val = sum_val + du3(i) * z(i+3); end
        z(i) = (b(i) - sum_val) / d(i);
    end
end

function pi = normalize_solution(z)
    pi = z / sum(z);
end

function [r, enorm] = compute_errors(T, pi)
    r = norm(T' * pi - pi, 2);
    enorm = abs(sum(pi) - 1.0);
end
