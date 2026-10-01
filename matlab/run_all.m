% =========================================================================
% Master Script: run_all.m
% Deskripsi: Menjalankan seluruh pengujian komputasi MATLAB untuk TK 1:
%            - Nomor 1: Analisis Perpindahan Penumpang Antarhalte (Dense vs Banded)
%            - Nomor 2: Model SETAR 2-Rezim (Normal Equations, Householder, Givens)
% =========================================================================

clear; clc; close all;

fprintf('#########################################################################\n');
fprintf('  TUGAS KELOMPOK 1 ANALISIS NUMERIK (CSCM603117) - FASILKOM UI\n');
fprintf('  EKSEKUSI SOLVER KOMPUTASI NUMERIK MANDIRI MENGGUNAKAN MATLAB / OCTAVE\n');
fprintf('#########################################################################\n\n');

% 1. Eksekusi Nomor 1
fprintf('>>> MENJALANKAN EKSPERIMEN NOMOR 1...\n');
no1_solver();

fprintf('\n>>> MENJALANKAN EKSPERIMEN NOMOR 2...\n');
no2_solver();

fprintf('#########################################################################\n');
fprintf('  SELURUH EKSPERIMEN MATLAB BERHASIL DIEKSEKUSI TANPA ERROR (100%% DONE)\n');
fprintf('#########################################################################\n');
