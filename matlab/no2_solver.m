% =========================================================================
% Program: no2_solver.m
% Deskripsi: Implementasi Solver Least Squares Problem (LSP) untuk Model SETAR 2-Rezim
%            Tugas Kelompok 1 Analisis Numerik - Fasilkom UI
%            Metode: Persamaan Normal, Householder Reflections, dan Givens Rotations
% =========================================================================

function no2_solver()
    fprintf('=========================================================================\n');
    fprintf('  NOMOR 2: DETEKSI REZIM PASAR & PREDIKSI RETURN SAHAM (SETAR)\n');
    fprintf('=========================================================================\n\n');

    % 1. Deteksi file dataset
    base_dir = fileparts(pwd);
    candidates = {
        pwd, ...
        fullfile(pwd, 'Nomor 2'), ...
        fullfile(pwd, 'Nomor 2', 'Nomor 2'), ...
        fullfile(base_dir, 'Nomor 2'), ...
        fullfile(base_dir, 'Nomor 2', 'Nomor 2')
    };

    train_file = '';
    test_file = '';
    for k = 1:length(candidates)
        f_tr = fullfile(candidates{k}, 'stock_train.csv');
        f_te = fullfile(candidates{k}, 'stock_test.csv');
        if exist(f_tr, 'file') && exist(f_te, 'file')
            train_file = f_tr;
            test_file = f_te;
            break;
        end
    end

    if isempty(train_file)
        error('Berkas dataset Nomor 2 (stock_train.csv dan stock_test.csv) tidak ditemukan.');
    end
    fprintf('Dataset Latih: %s\n', train_file);
    fprintf('Dataset Uji  : %s\n\n', test_file);

    % i) Formulasi Matriks Desain Overdetermined A dan Vektor Target b
    [A, b, returns_tr, dates_tr] = load_and_preprocess_stock(train_file);
    [m, n] = size(A);
    fprintf('Dimensi Matriks Desain A : %d x %d\n', m, n);
    fprintf('Dimensi Vektor Solusi x  : %d x 1\n', n);
    fprintf('Dimensi Vektor Target b  : %d x 1\n\n', m);

    % iii) Penyelesaian via Persamaan Normal
    ATA = A' * A;
    ATb = A' * b;
    cond_A = cond(A);
    cond_ATA = cond(ATA);
    x_normal = custom_lu_solve(ATA, ATb);
    res_normal = norm(A * x_normal - b, 2);

    fprintf('--- PENYELESAIAN VIA PERSAMAAN NORMAL ---\n');
    fprintf('Condition Number kappa_2(A)   : %12.6f\n', cond_A);
    fprintf('Condition Number kappa_2(A^T A): %12.6f\n', cond_ATA);
    fprintf('Kuadrat kappa_2(A)^2          : %12.6f\n', cond_A^2);
    fprintf('Norm Residual ||Ax - b||_2    : %12.6f\n\n', res_normal);

    % iv-a) Verifikasi Langkah Awal Householder Reflection (Kelompok Ganjil)
    [H1, v1, alpha_h1, H1A, max_subdiag_err_h1] = verify_householder_step1(A);
    fprintf('--- VERIFIKASI REFLEKSI HOUSEHOLDER LANGKAH 1 (H1) ---\n');
    fprintf('Skalar alpha H1               : %12.6f\n', alpha_h1);
    fprintf('Maks Galat Subdiagonal H1*A   : %12.4e (Presisi Mesin)\n\n', max_subdiag_err_h1);

    % iv-b) Formulasi Rotasi Givens Pertama (Kelompok Genap)
    [G1, c_g1, s_g1, G1A, err_subdiag_g1] = verify_givens_step1(A);
    fprintf('--- FORMULASI ROTASI GIVENS LANGKAH 1 (G1) ---\n');
    fprintf('Parameter Rotasi: c = %.8f, s = %.8f\n', c_g1, s_g1);
    fprintf('Galat Elemen Subdiagonal (2,1): %12.4e (Presisi Mesin)\n\n', err_subdiag_g1);

    % iv-c) Penyelesaian LSP via Faktorisasi QR Householder Mandiri
    [x_qr, R_qr, res_qr] = solve_householder_qr(A, b);
    fprintf('--- PENYELESAIAN VIA FAKTORISASI QR HOUSEHOLDER ---\n');
    fprintf('Norm Residual ||A*x_LS - b||_2: %12.6f\n', res_qr);
    fprintf('Selisih Solusi QR vs Normal   : %12.4e\n\n', norm(x_qr - x_normal, 2));

    % iv-d) Penyelesaian LSP via Faktorisasi QR Givens Rotations Mandiri
    [x_givens, res_givens] = solve_givens_qr(A, b);
    fprintf('--- PENYELESAIAN VIA FAKTORISASI QR GIVENS ROTATIONS ---\n');
    fprintf('Norm Residual ||A*x_LS - b||_2: %12.6f\n', res_givens);
    fprintf('Selisih Solusi Givens vs QR   : %12.4e\n\n', norm(x_givens - x_qr, 2));

    % Cetak Koefisien Parameter Model SETAR
    param_names = {
        'alpha_1 (Drift Bullish)', ...
        'phi_1,1 (Lag-1 Bullish)', ...
        'phi_1,2 (Lag-2 Bullish)', ...
        'alpha_2 (Drift Bearish)', ...
        'phi_2,1 (Lag-1 Bearish)', ...
        'phi_2,2 (Lag-2 Bearish)'
    };
    fprintf('=========================================================================\n');
    fprintf('                  ESTIMASI PARAMETER MODEL SETAR 2-REZIM\n');
    fprintf('=========================================================================\n');
    fprintf('%-25s | %-15s | %-15s\n', 'Parameter', 'Estimasi QR', 'Estimasi Normal');
    fprintf('-------------------------------------------------------------------------\n');
    for i = 1:6
        fprintf('%-25s | %15.8f | %15.8f\n', param_names{i}, x_qr(i), x_normal(i));
    end
    fprintf('=========================================================================\n\n');

    % vi) Evaluasi Out-of-Sample pada stock_test.csv
    eval_res = evaluate_out_of_sample(train_file, test_file, x_qr);
    fprintf('--- EVALUASI AKURASI OUT-OF-SAMPLE (RMSE) ---\n');
    fprintf('RMSE Data Latih (Train, 300 obs)          : %10.6f (%.3f%%)\n', ...
            eval_res.rmse_train, eval_res.rmse_train * 100);
    fprintf('RMSE Data Uji (Test Independen, 100 obs)  : %10.6f (%.3f%%)\n', ...
            eval_res.rmse_test_ind, eval_res.rmse_test_ind * 100);
    fprintf('RMSE Data Uji (Test Berkesinambungan, 103): %10.6f (%.3f%%)\n\n', ...
            eval_res.rmse_test_cont, eval_res.rmse_test_cont * 100);

    % vii) Plotting Visualisasi Overlay Time Series
    figure('Name', 'Model SETAR 2-Rezim - Prediksi vs Aktual', 'Position', [100, 100, 1000, 500]);
    b_pred = A * x_qr;
    
    subplot(2, 1, 1);
    plot(1:m, b, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 1.0, 'DisplayName', 'Return Aktual R_t');
    hold on;
    plot(1:m, b_pred, 'Color', [0.9, 0.3, 0.1], 'LineWidth', 1.2, 'DisplayName', 'Estimasi SETAR');
    title('Overlay Return Aktual vs Estimasi SETAR 2-Rezim (Data Latih)');
    xlabel('Indeks Observasi (t)');
    ylabel('Return Harian');
    legend('Location', 'northeast');
    grid on;

    subplot(2, 1, 2);
    residuals = b - b_pred;
    stem(1:m, residuals, 'Marker', 'none', 'Color', [0.4, 0.4, 0.4]);
    yline(0, 'r--', 'LineWidth', 1.0);
    title('Deret Residual \epsilon_t = R_t - \hat{R}_t');
    xlabel('Indeks Observasi (t)');
    ylabel('Residual');
    grid on;
end

% =========================================================================
% FUNGSI-FUNGSI SUB-ALGORITMA MANDIRI NOMOR 2
% =========================================================================

function [A, b, returns, dates] = load_and_preprocess_stock(csv_file)
    T = readtable(csv_file);
    prices = T.Close;
    dates = T.Date;
    
    % Return harian: R_t = (P_t - P_{t-1}) / P_{t-1}
    returns = (prices(2:end) - prices(1:end-1)) ./ prices(1:end-1);
    
    % Konstruksi matriks overdetermined (lag p = 2)
    m = length(returns) - 2;
    A = zeros(m, 6);
    b = zeros(m, 1);
    
    for i = 1:m
        t_idx = i + 2;
        rt = returns(t_idx);
        rt1 = returns(t_idx - 1);
        rt2 = returns(t_idx - 2);
        
        b(i) = rt;
        if rt1 >= 0.0
            A(i, 1) = 1.0;
            A(i, 2) = rt1;
            A(i, 3) = rt2;
        else
            A(i, 4) = 1.0;
            A(i, 5) = rt1;
            A(i, 6) = rt2;
        end
    end
end

function x = custom_lu_solve(M, rhs)
    % Solver SPL non-singular berukuran kecil dengan LU Partial Pivoting
    n = size(M, 1);
    A = M;
    b = rhs;
    p = 1:n;
    
    for k = 1:(n - 1)
        [~, max_idx] = max(abs(A(k:n, k)));
        piv = k + max_idx - 1;
        if piv ~= k
            A([k, piv], :) = A([piv, k], :);
            p([k, piv]) = p([piv, k]);
        end
        pivot_val = A(k, k);
        if abs(pivot_val) < 1e-15, continue; end
        
        for i = (k + 1):n
            f = A(i, k) / pivot_val;
            A(i, k) = f;
            A(i, (k+1):n) = A(i, (k+1):n) - f * A(k, (k+1):n);
        end
    end
    
    % Forward substitution
    b_p = b(p);
    y = zeros(n, 1);
    for i = 1:n
        y(i) = b_p(i) - A(i, 1:(i-1)) * y(1:(i-1));
    end
    
    % Back substitution
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - A(i, (i+1):n) * x((i+1):n)) / A(i, i);
    end
end

function [H1, v1, alpha, H1A, max_err] = verify_householder_step1(A)
    m = size(A, 1);
    a1 = A(:, 1);
    norm_a1 = norm(a1, 2);
    
    sign_val = sign(a1(1));
    if sign_val == 0, sign_val = 1.0; end
    alpha = -sign_val * norm_a1;
    
    v1 = a1;
    v1(1) = v1(1) - alpha;
    v1 = v1 / norm(v1, 2);
    
    H1 = eye(m) - 2.0 * (v1 * v1');
    H1A = H1 * A;
    max_err = max(abs(H1A(2:end, 1)));
end

function [G1, c, s, G1A, err] = verify_givens_step1(A)
    % Rotasi Givens untuk mengenolkan elemen (2,1) menggunakan pivot baris 1
    a = A(1, 1);
    b_elem = A(2, 1);
    r = hypot(a, b_elem);
    if r == 0
        c = 1.0; s = 0.0;
    else
        c = a / r;
        s = -b_elem / r;
    end
    
    m = size(A, 1);
    G1 = eye(m);
    G1(1, 1) = c;  G1(1, 2) = -s;
    G1(2, 1) = s;  G1(2, 2) = c;
    
    G1A = G1 * A;
    err = abs(G1A(2, 1));
end

function [x_ls, R1, res_norm] = solve_householder_qr(A_in, b_in)
    [m, n] = size(A_in);
    R = A_in;
    c = b_in;
    
    for k = 1:n
        x = R(k:m, k);
        norm_x = norm(x, 2);
        if norm_x < 1e-15, continue; end
        
        sign_val = sign(x(1));
        if sign_val == 0, sign_val = 1.0; end
        alpha = -sign_val * norm_x;
        
        u = x;
        u(1) = u(1) - alpha;
        v = u / norm(u, 2);
        
        % Terapkan transformasi Householder ke submatriks R
        v_dot_R = v' * R(k:m, k:n);
        R(k:m, k:n) = R(k:m, k:n) - 2.0 * (v * v_dot_R);
        
        % Terapkan ke vektor sisi kanan c
        v_dot_c = v' * c(k:m);
        c(k:m) = c(k:m) - 2.0 * v * v_dot_c;
    end
    
    R1 = R(1:n, 1:n);
    c1 = c(1:n);
    
    % Substitusi Mundur
    x_ls = zeros(n, 1);
    for i = n:-1:1
        x_ls(i) = (c1(i) - R1(i, (i+1):n) * x_ls((i+1):n)) / R1(i, i);
    end
    
    res_norm = norm(A_in * x_ls - b_in, 2);
end

function [x_ls, res_norm] = solve_givens_qr(A_in, b_in)
    [m, n] = size(A_in);
    R = A_in;
    c = b_in;
    
    for j = 1:n
        for i = m:-1:(j + 1)
            a_val = R(i-1, j);
            b_val = R(i, j);
            if abs(b_val) < 1e-15, continue; end
            
            hyp = hypot(a_val, b_val);
            c_rot = a_val / hyp;
            s_rot = -b_val / hyp;
            
            % Rotasi baris i-1 dan i
            row_prev = R(i-1, j:n);
            row_curr = R(i, j:n);
            R(i-1, j:n) = c_rot * row_prev - s_rot * row_curr;
            R(i, j:n)   = s_rot * row_prev + c_rot * row_curr;
            
            c_prev = c(i-1);
            c_curr = c(i);
            c(i-1) = c_rot * c_prev - s_rot * c_curr;
            c(i)   = s_rot * c_prev + c_rot * c_curr;
        end
    end
    
    R1 = R(1:n, 1:n);
    c1 = c(1:n);
    
    x_ls = zeros(n, 1);
    for i = n:-1:1
        x_ls(i) = (c1(i) - R1(i, (i+1):n) * x_ls((i+1):n)) / R1(i, i);
    end
    
    res_norm = norm(A_in * x_ls - b_in, 2);
end

function eval_res = evaluate_out_of_sample(train_file, test_file, x_ls)
    T_tr = readtable(train_file);
    T_te = readtable(test_file);
    
    p_tr = T_tr.Close;
    p_te = T_te.Close;
    
    % Train RMSE
    [A_tr, b_tr, ~, ~] = load_and_preprocess_stock(train_file);
    eval_res.rmse_train = sqrt(mean((b_tr - A_tr * x_ls).^2));
    
    % Test Independen (100 obs)
    r_te = (p_te(2:end) - p_te(1:end-1)) ./ p_te(1:end-1);
    m_te = length(r_te) - 2;
    A_te = zeros(m_te, 6);
    b_te = zeros(m_te, 1);
    for i = 1:m_te
        t = i + 2;
        b_te(i) = r_te(t);
        if r_te(t-1) >= 0
            A_te(i, 1) = 1.0; A_te(i, 2) = r_te(t-1); A_te(i, 3) = r_te(t-2);
        else
            A_te(i, 4) = 1.0; A_te(i, 5) = r_te(t-1); A_te(i, 6) = r_te(t-2);
        end
    end
    eval_res.rmse_test_ind = sqrt(mean((b_te - A_te * x_ls).^2));
    
    % Test Kontinu (103 obs dari akhir data latih)
    p_all = [p_tr; p_te];
    r_all = (p_all(2:end) - p_all(1:end-1)) ./ p_all(1:end-1);
    start_idx = length(p_tr) - 1;
    m_cont = length(r_all) - start_idx + 1;
    A_cont = zeros(m_cont, 6);
    b_cont = zeros(m_cont, 1);
    cnt = 1;
    for i = start_idx:length(r_all)
        b_cont(cnt) = r_all(i);
        if r_all(i-1) >= 0
            A_cont(cnt, 1) = 1.0; A_cont(cnt, 2) = r_all(i-1); A_cont(cnt, 3) = r_all(i-2);
        else
            A_cont(cnt, 4) = 1.0; A_cont(cnt, 5) = r_all(i-1); A_cont(cnt, 6) = r_all(i-2);
        end
        cnt = cnt + 1;
    end
    eval_res.rmse_test_cont = sqrt(mean((b_cont - A_cont * x_ls).^2));
end
