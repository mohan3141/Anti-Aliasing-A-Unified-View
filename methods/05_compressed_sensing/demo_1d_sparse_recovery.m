%% Compressed Sensing Demo: 1D Sparse Spectrum Recovery
% ============================================================
% Scenario: Recover a signal with K sparse frequency components
%           using only M < N measurements (sub-Nyquist)
%
% Key Concepts:
%   - Sparse signal model
%   - Random measurement matrices
%   - OMP (Orthogonal Matching Pursuit) algorithm
%   - L1-minimization (Basis Pursuit)
%   - RIP and coherence
%
% Author: [Your Name]
% Date: 2024
% ============================================================

clear; clc; close all;

fprintf('==========================================================\n');
fprintf('  Compressed Sensing: 1D Sparse Spectrum Recovery\n');
fprintf('==========================================================\n\n');

%% ==================== Parameters ====================
N = 256;          % Signal length (Nyquist samples)
K = 5;            % Sparsity (number of frequency components)
M = 60;           % Number of CS measurements (M << N)
SNR_dB = 40;      % Measurement SNR

fprintf('【Configuration】\n');
fprintf('  Signal length N = %d\n', N);
fprintf('  Sparsity K = %d\n', K);
fprintf('  Measurements M = %d (%.1f%% of N)\n', M, 100*M/N);
fprintf('  Compression ratio: %.1fx\n', N/M);
fprintf('  SNR = %d dB\n\n', SNR_dB);

%% ==================== Generate Sparse Signal ====================
fprintf('【Step 1: Generate K-sparse signal in frequency domain】\n');

% Random sparse frequency support
rng(42);  % For reproducibility
freq_support = sort(randperm(N, K));  % K random frequency indices

% Random amplitudes and phases
amplitudes = 0.5 + rand(K, 1);  % Random amplitudes [0.5, 1.5]
phases = 2*pi * rand(K, 1);     % Random phases [0, 2π]

% Construct sparse frequency vector
s_true = zeros(N, 1);  % Sparse in frequency domain
for i = 1:K
    s_true(freq_support(i)) = amplitudes(i) * exp(1j * phases(i));
end

% Time-domain signal: x = IDFT(s)
x = ifft(s_true) * sqrt(N);  % Normalized IDFT

fprintf('  True frequency support: [%s]\n', num2str(freq_support));
fprintf('  Signal energy: %.4f\n\n', norm(x)^2);

%% ==================== Measurement Matrix ====================
fprintf('【Step 2: Construct measurement matrix Φ】\n');

% Option 1: Gaussian random matrix
Phi = randn(M, N) / sqrt(M);

% Option 2: Bernoulli random matrix (uncomment to use)
% Phi = (2*(rand(M, N) > 0.5) - 1) / sqrt(M);

% Option 3: Subsampled DFT (uncomment to use)
% sample_indices = sort(randperm(N, M));
% F = dftmtx(N) / sqrt(N);
% Phi = real(F(sample_indices, :));

fprintf('  Measurement matrix: %d × %d Gaussian random\n', M, N);
fprintf('  Frobenius norm: %.4f\n\n', norm(Phi, 'fro'));

%% ==================== Sensing Matrix ====================
fprintf('【Step 3: Construct sensing matrix A = Φ·Ψ】\n');

% Sparsity basis: Inverse DFT (signal sparse in frequency)
% Ψ = IDFT matrix, so x = Ψ·s where s is frequency coefficients
Psi = ifft(eye(N)) * sqrt(N);  % IDFT matrix

% Sensing matrix
A = Phi * Psi;

fprintf('  Sparsity basis: Inverse DFT\n');
fprintf('  Sensing matrix A: %d × %d\n', size(A, 1), size(A, 2));

% Compute coherence (for reference)
A_normalized = A ./ vecnorm(A);
coherence = max(max(abs(A_normalized' * A_normalized - eye(N))));
fprintf('  Coherence (approx): %.4f\n\n', coherence);

%% ==================== Take Measurements ====================
fprintf('【Step 4: Acquire M measurements】\n');

% Clean measurements
y_clean = Phi * x;

% Add noise
noise_power = norm(y_clean)^2 / M / (10^(SNR_dB/10));
noise = sqrt(noise_power/2) * (randn(M, 1) + 1j*randn(M, 1));
y = y_clean + noise;

fprintf('  Measurement vector y: %d × 1\n', length(y));
fprintf('  Measurement SNR: %.1f dB\n\n', 10*log10(norm(y_clean)^2/norm(noise)^2));

%% ==================== OMP Recovery Function (inline for Octave compatibility) ====================

% OMP Implementation
function [s, support, residuals] = omp_func(A, y, K)
    [M, N] = size(A);
    r = y;
    support = [];
    s = zeros(N, 1);
    residuals = zeros(K, 1);
    
    for k = 1:K
        correlations = abs(A' * r);
        correlations(support) = 0;
        [~, j_star] = max(correlations);
        support = [support; j_star];
        A_S = A(:, support);
        s_S = A_S \ y;
        r = y - A_S * s_S;
        residuals(k) = norm(r);
    end
    s(support) = s_S;
end

%% ==================== Recovery: OMP ====================
fprintf('【Step 5: OMP Recovery】\n');

tic;
[s_omp, support_omp, residuals_omp] = omp_func(A, y, K);
time_omp = toc;

% Reconstruction error
x_omp = Psi * s_omp;
nmse_omp = norm(x - x_omp)^2 / norm(x)^2;
support_correct_omp = length(intersect(freq_support, support_omp));

fprintf('  Recovered support: [%s]\n', num2str(sort(support_omp)'));
fprintf('  Correct support: %d/%d\n', support_correct_omp, K);
fprintf('  NMSE: %.2e (%.2f dB)\n', nmse_omp, 10*log10(nmse_omp));
fprintf('  Time: %.4f sec\n\n', time_omp);

%% ==================== ISTA Recovery Function (inline for Octave compatibility) ====================

function z = soft_thresh(x, tau)
    z = sign(x) .* max(abs(x) - tau, 0);
end

function [s, cost_history] = ista_func(A, y, lambda, max_iter)
    [~, N] = size(A);
    L = norm(A)^2;
    t = 1 / L;
    s = zeros(N, 1);
    cost_history = zeros(max_iter, 1);
    
    for iter = 1:max_iter
        grad = A' * (A * s - y);
        s_grad = s - t * grad;
        s = soft_thresh(s_grad, lambda * t);
        cost_history(iter) = 0.5 * norm(A*s - y)^2 + lambda * norm(s, 1);
    end
end

%% ==================== Recovery: L1-Minimization (ISTA) ====================
fprintf('【Step 6: L1-Minimization Recovery (ISTA)】\n');

tic;
lambda = 0.01 * norm(A' * y, 'inf');  % Regularization parameter
[s_l1, cost_history] = ista_func(A, y, lambda, 500);
time_l1 = toc;

% Find support by thresholding
threshold = 0.1 * max(abs(s_l1));
support_l1 = find(abs(s_l1) > threshold);

% Reconstruction error
x_l1 = Psi * s_l1;
nmse_l1 = norm(x - x_l1)^2 / norm(x)^2;
support_correct_l1 = length(intersect(freq_support, support_l1));

fprintf('  Detected components: %d\n', length(support_l1));
fprintf('  Correct support: %d/%d\n', support_correct_l1, K);
fprintf('  NMSE: %.2e (%.2f dB)\n', nmse_l1, 10*log10(nmse_l1));
fprintf('  Time: %.4f sec\n\n', time_l1);

%% ==================== Comparison with Pseudo-inverse ====================
fprintf('【Step 7: Comparison with naive methods】\n');

% Pseudo-inverse (no sparsity constraint)
s_pinv = pinv(A) * y;
x_pinv = Psi * s_pinv;
nmse_pinv = norm(x - x_pinv)^2 / norm(x)^2;

fprintf('  Pseudo-inverse NMSE: %.2e (%.2f dB)\n', nmse_pinv, 10*log10(nmse_pinv));
fprintf('  OMP improvement: %.1f dB\n', 10*log10(nmse_pinv) - 10*log10(nmse_omp));
fprintf('  L1 improvement: %.1f dB\n\n', 10*log10(nmse_pinv) - 10*log10(nmse_l1));

%% ==================== Phase Transition Analysis ====================
fprintf('【Step 8: Phase Transition (M vs K)】\n');

% Test different (M, K) combinations
K_range = [2, 5, 10, 15, 20];
M_range = [20, 40, 60, 80, 100, 120];
n_trials = 10;
success_rate = zeros(length(M_range), length(K_range));

fprintf('  Running phase transition analysis...\n');
for ki = 1:length(K_range)
    for mi = 1:length(M_range)
        K_test = K_range(ki);
        M_test = M_range(mi);
        
        n_success = 0;
        for trial = 1:n_trials
            % Generate test signal
            support_test = sort(randperm(N, K_test));
            s_test = zeros(N, 1);
            s_test(support_test) = randn(K_test, 1) + 1j*randn(K_test, 1);
            x_test = Psi * s_test;
            
            % Measurements
            Phi_test = randn(M_test, N) / sqrt(M_test);
            A_test = Phi_test * Psi;
            y_test = Phi_test * x_test;
            
            % OMP recovery
            [s_rec, ~, ~] = omp_func(A_test, y_test, K_test);
            
            % Check success (NMSE < -20dB)
            nmse_test = norm(s_test - s_rec)^2 / norm(s_test)^2;
            if nmse_test < 0.01  % -20 dB
                n_success = n_success + 1;
            end
        end
        success_rate(mi, ki) = n_success / n_trials;
    end
end

fprintf('  Done!\n\n');

%% ==================== Visualization ====================
fprintf('【Generating figures...】\n\n');

figure('Position', [100, 100, 1400, 900]);

% 1. True vs Recovered Spectrum
subplot(2, 3, 1);
stem(1:N, abs(s_true), 'b', 'LineWidth', 1.5, 'MarkerSize', 4);
hold on;
stem(1:N, abs(s_omp), 'r--', 'LineWidth', 1, 'MarkerSize', 3);
xlabel('Frequency Index');
ylabel('Magnitude');
title('Spectrum: True (blue) vs OMP (red)');
legend('True', 'OMP', 'Location', 'best');
xlim([1 N]);
grid on;

% 2. Time-domain signal
subplot(2, 3, 2);
plot(1:N, real(x), 'b-', 'LineWidth', 1.5);
hold on;
plot(1:N, real(x_omp), 'r--', 'LineWidth', 1);
xlabel('Sample Index');
ylabel('Amplitude');
title('Time Domain: True vs OMP');
legend('True', 'OMP', 'Location', 'best');
grid on;

% 3. OMP residual evolution
subplot(2, 3, 3);
semilogy(0:K, [norm(y); residuals_omp(1:K)], 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Iteration');
ylabel('Residual Norm');
title('OMP Residual Decay');
grid on;

% 4. L1 (ISTA) cost function
subplot(2, 3, 4);
semilogy(1:length(cost_history), cost_history, 'g-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Cost Function');
title('ISTA Convergence');
grid on;

% 5. Comparison: True vs OMP vs L1 vs Pinv
subplot(2, 3, 5);
bar_data = [norm(s_true - s_omp)/norm(s_true), ...
            norm(s_true - s_l1)/norm(s_true), ...
            norm(s_true - s_pinv)/norm(s_true)];
bar(bar_data, 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTickLabel', {'OMP', 'L1 (ISTA)', 'Pinv'});
ylabel('Relative Error');
title('Recovery Error Comparison');
grid on;

% 6. Phase Transition Heatmap
subplot(2, 3, 6);
imagesc(K_range, M_range, success_rate);
colorbar;
colormap(flipud(hot));
xlabel('Sparsity K');
ylabel('Measurements M');
title(sprintf('Phase Transition (N=%d)', N));
set(gca, 'YDir', 'normal');
for ki = 1:length(K_range)
    for mi = 1:length(M_range)
        text(K_range(ki), M_range(mi), sprintf('%.0f%%', 100*success_rate(mi,ki)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end

% sgtitle replacement for Octave
try
    sgtitle('Compressed Sensing: 1D Sparse Spectrum Recovery', 'FontSize', 14);
catch
    % For Octave: add title to first subplot instead
    subplot(2, 3, 1);
    title({'Compressed Sensing: 1D Sparse Spectrum Recovery', 'Spectrum: True (blue) vs OMP (red)'});
end

% Save figure
saveas(gcf, 'cs_1d_recovery_results.png');
fprintf('  Figure saved: cs_1d_recovery_results.png\n');

%% ==================== Summary ====================
fprintf('\n==========================================================\n');
fprintf('  SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Signal: N=%d, K=%d sparse in frequency\n', N, K);
fprintf('  Measurements: M=%d (%.1f%% compression)\n', M, 100*M/N);
fprintf('  \n');
fprintf('  Recovery Results:\n');
fprintf('    Method      | NMSE (dB) | Support | Time\n');
fprintf('    ------------|-----------|---------|--------\n');
fprintf('    OMP         | %8.2f  | %d/%d    | %.3fs\n', 10*log10(nmse_omp), support_correct_omp, K, time_omp);
fprintf('    L1 (ISTA)   | %8.2f  | %d/%d    | %.3fs\n', 10*log10(nmse_l1), support_correct_l1, K, time_l1);
fprintf('    Pseudo-inv  | %8.2f  | N/A     | N/A\n', 10*log10(nmse_pinv));
fprintf('==========================================================\n');

%% ==================== End of Script ====================
% Note: Functions omp_func, ista_func, and soft_thresh are defined
% inline above for Octave compatibility
