%% Compressed Sensing Demo: 2D Image Reconstruction
% ============================================================
% Scenario: Reconstruct an image from randomly sampled Fourier
%           measurements (simulating MRI or single-pixel camera)
%
% Key Concepts:
%   - 2D Fourier undersampling
%   - Total Variation (TV) regularization
%   - Variable density sampling patterns
%   - ADMM optimization
%
% Author: [Your Name]
% Date: 2024
% ============================================================

clear; clc; close all;

fprintf('==========================================================\n');
fprintf('  Compressed Sensing: 2D Image Reconstruction\n');
fprintf('==========================================================\n\n');

%% ==================== Parameters ====================
% Image size
N = 128;           % Image dimension (N × N)
sampling_rate = 0.25;  % Fraction of Fourier samples (25%)
SNR_dB = 40;       % Measurement SNR

% TV parameters
lambda_tv = 0.01;  % TV regularization weight
n_iter = 200;      % Number of iterations

fprintf('【Configuration】\n');
fprintf('  Image size: %d × %d\n', N, N);
fprintf('  Sampling rate: %.1f%%\n', 100*sampling_rate);
fprintf('  SNR: %d dB\n', SNR_dB);
fprintf('  TV regularization λ = %.4f\n\n', lambda_tv);

%% ==================== Generate Test Image ====================
fprintf('【Step 1: Generate test image】\n');

% Create a phantom-like test image (piecewise constant)
x_true = create_phantom(N);

fprintf('  Image created: %d × %d phantom\n', N, N);
fprintf('  Dynamic range: [%.2f, %.2f]\n', min(x_true(:)), max(x_true(:)));
fprintf('  Total Variation: %.2f\n\n', compute_tv(x_true));

%% ==================== Sampling Pattern ====================
fprintf('【Step 2: Design sampling pattern】\n');

% Option 1: Uniform random sampling
% mask_uniform = rand(N, N) < sampling_rate;

% Option 2: Variable density sampling (higher density at low frequencies)
% This mimics MRI radial/spiral trajectories
mask = create_variable_density_mask(N, sampling_rate);

M = sum(mask(:));  % Actual number of samples
actual_rate = M / N^2;

fprintf('  Sampling pattern: Variable density\n');
fprintf('  Samples acquired: %d / %d (%.1f%%)\n', M, N^2, 100*actual_rate);
fprintf('  Compression ratio: %.1fx\n\n', N^2/M);

%% ==================== Forward Model ====================
fprintf('【Step 3: Acquire undersampled Fourier data】\n');

% Full Fourier transform
F_full = fft2(x_true) / N;  % Normalized 2D FFT

% Apply sampling mask
F_sampled = F_full .* mask;

% Add noise
noise_power = norm(F_sampled(:))^2 / M / (10^(SNR_dB/10));
noise = sqrt(noise_power/2) * (randn(N,N) + 1j*randn(N,N)) .* mask;
y = F_sampled + noise;

fprintf('  Full k-space: %d × %d\n', N, N);
fprintf('  Measured k-space: %d samples\n', M);
fprintf('  Measurement SNR: %.1f dB\n\n', SNR_dB);

%% ==================== Zero-filled Reconstruction ====================
fprintf('【Step 4: Zero-filled reconstruction (naive)】\n');

% Simply inverse FFT the zero-filled k-space
x_zf = real(ifft2(y) * N);

% Compute error metrics
nmse_zf = norm(x_true(:) - x_zf(:))^2 / norm(x_true(:))^2;
psnr_zf = compute_psnr(x_true, x_zf);

fprintf('  NMSE: %.4f (%.2f dB)\n', nmse_zf, 10*log10(nmse_zf));
fprintf('  PSNR: %.2f dB\n', psnr_zf);
fprintf('  → Severe aliasing artifacts visible!\n\n');

%% ==================== TV-Regularized Reconstruction ====================
fprintf('【Step 5: TV-regularized CS reconstruction】\n');

% Solve: min_x ||F_mask(x) - y||^2 + lambda * TV(x)
% Using ADMM (Alternating Direction Method of Multipliers)

fprintf('  Running ADMM optimization...\n');
tic;
[x_tv, cost_history] = cs_tv_reconstruction(y, mask, lambda_tv, n_iter);
time_tv = toc;

% Compute error metrics
nmse_tv = norm(x_true(:) - x_tv(:))^2 / norm(x_true(:))^2;
psnr_tv = compute_psnr(x_true, x_tv);

fprintf('  Iterations: %d\n', n_iter);
fprintf('  Time: %.2f sec\n', time_tv);
fprintf('  NMSE: %.4f (%.2f dB)\n', nmse_tv, 10*log10(nmse_tv));
fprintf('  PSNR: %.2f dB\n', psnr_tv);
fprintf('  Improvement over zero-fill: %.1f dB\n\n', psnr_tv - psnr_zf);

%% ==================== Wavelet-based Reconstruction ====================
fprintf('【Step 6: Wavelet-based CS reconstruction】\n');

% Use wavelet sparsity instead of TV
fprintf('  Running wavelet-based optimization...\n');
tic;
[x_wav, cost_wav] = cs_wavelet_reconstruction(y, mask, lambda_tv*0.5, n_iter);
time_wav = toc;

nmse_wav = norm(x_true(:) - x_wav(:))^2 / norm(x_true(:))^2;
psnr_wav = compute_psnr(x_true, x_wav);

fprintf('  Time: %.2f sec\n', time_wav);
fprintf('  NMSE: %.4f (%.2f dB)\n', nmse_wav, 10*log10(nmse_wav));
fprintf('  PSNR: %.2f dB\n\n', psnr_wav);

%% ==================== Effect of Sampling Rate ====================
fprintf('【Step 7: Effect of sampling rate】\n');

sampling_rates = [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50];
psnr_vs_rate = zeros(length(sampling_rates), 1);
nmse_vs_rate = zeros(length(sampling_rates), 1);

fprintf('  Testing different sampling rates...\n');
for i = 1:length(sampling_rates)
    rate = sampling_rates(i);
    mask_test = create_variable_density_mask(N, rate);
    y_test = (fft2(x_true)/N + noise*0.1) .* mask_test;
    
    [x_rec, ~] = cs_tv_reconstruction(y_test, mask_test, lambda_tv, 100);
    
    nmse_vs_rate(i) = norm(x_true(:) - x_rec(:))^2 / norm(x_true(:))^2;
    psnr_vs_rate(i) = compute_psnr(x_true, x_rec);
    
    fprintf('    Rate %.0f%%: PSNR = %.1f dB\n', 100*rate, psnr_vs_rate(i));
end
fprintf('\n');

%% ==================== Visualization ====================
fprintf('【Generating figures...】\n\n');

figure('Position', [50, 50, 1600, 900]);

% Row 1: Images
subplot(2, 4, 1);
imagesc(x_true); axis image; colormap(gray); colorbar;
title('Ground Truth');

subplot(2, 4, 2);
imagesc(fftshift(mask)); axis image; colormap(gray); colorbar;
title(sprintf('Sampling Mask (%.0f%%)', 100*actual_rate));

subplot(2, 4, 3);
imagesc(x_zf); axis image; colormap(gray); colorbar;
title(sprintf('Zero-filled (PSNR=%.1fdB)', psnr_zf));

subplot(2, 4, 4);
imagesc(x_tv); axis image; colormap(gray); colorbar;
title(sprintf('TV Recon (PSNR=%.1fdB)', psnr_tv));

% Row 2: Analysis
subplot(2, 4, 5);
imagesc(abs(x_true - x_zf)); axis image; colorbar;
title('Error: Zero-filled');
caxis([0, 0.5]);

subplot(2, 4, 6);
imagesc(abs(x_true - x_tv)); axis image; colorbar;
title('Error: TV Reconstruction');
caxis([0, 0.5]);

subplot(2, 4, 7);
semilogy(1:length(cost_history), cost_history, 'b-', 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('Cost Function');
title('ADMM Convergence');
grid on;

subplot(2, 4, 8);
plot(100*sampling_rates, psnr_vs_rate, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Sampling Rate (%)'); ylabel('PSNR (dB)');
title('PSNR vs Sampling Rate');
grid on;
xlim([5 55]);

sgtitle('Compressed Sensing: 2D Image Reconstruction', 'FontSize', 14);

% Save figure
saveas(gcf, 'cs_2d_image_results.png');
fprintf('  Figure saved: cs_2d_image_results.png\n');

%% ==================== Summary ====================
fprintf('\n==========================================================\n');
fprintf('  SUMMARY\n');
fprintf('==========================================================\n');
fprintf('  Image: %d × %d phantom\n', N, N);
fprintf('  Sampling rate: %.1f%%\n', 100*actual_rate);
fprintf('  \n');
fprintf('  Reconstruction Results:\n');
fprintf('    Method        | PSNR (dB) | NMSE (dB) | Time\n');
fprintf('    --------------|-----------|-----------|------\n');
fprintf('    Zero-filled   | %8.2f  | %8.2f  | instant\n', psnr_zf, 10*log10(nmse_zf));
fprintf('    TV-regularized| %8.2f  | %8.2f  | %.2fs\n', psnr_tv, 10*log10(nmse_tv), time_tv);
fprintf('    Wavelet-based | %8.2f  | %8.2f  | %.2fs\n', psnr_wav, 10*log10(nmse_wav), time_wav);
fprintf('==========================================================\n');

%% ==================== Local Functions ====================

function x = create_phantom(N)
% CREATE_PHANTOM Create a simple piecewise constant phantom image
    x = zeros(N, N);
    
    % Background
    [X, Y] = meshgrid(1:N, 1:N);
    center = N/2;
    
    % Large ellipse (body)
    a1 = N*0.4; b1 = N*0.3;
    mask1 = ((X-center)/a1).^2 + ((Y-center)/b1).^2 <= 1;
    x(mask1) = 0.8;
    
    % Smaller circles (organs)
    r2 = N*0.1;
    mask2 = (X - center*0.7).^2 + (Y - center).^2 <= r2^2;
    x(mask2) = 1.0;
    
    mask3 = (X - center*1.3).^2 + (Y - center).^2 <= r2^2;
    x(mask3) = 0.5;
    
    % Small bright spot
    r4 = N*0.05;
    mask4 = (X - center).^2 + (Y - center*0.7).^2 <= r4^2;
    x(mask4) = 1.2;
    
    % Rectangle
    x(round(N*0.6):round(N*0.7), round(N*0.3):round(N*0.5)) = 0.3;
end

function mask = create_variable_density_mask(N, rate)
% CREATE_VARIABLE_DENSITY_MASK Create sampling mask with higher density at center
    [X, Y] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
    
    % Distance from center (normalized)
    R = sqrt(X.^2 + Y.^2) / (N/2);
    
    % Variable density function (polynomial decay)
    density = max(0, 1 - R.^2);  % Quadratic decay from center
    density = density / sum(density(:));  % Normalize
    
    % Scale to achieve desired rate
    density = density * rate * N^2;
    density = min(density, 1);  % Cap at 1
    
    % Sample according to density
    mask = fftshift(rand(N, N) < density);
    
    % Always include DC component
    mask(1, 1) = 1;
end

function tv = compute_tv(x)
% COMPUTE_TV Compute isotropic total variation
    [gx, gy] = gradient(x);
    tv = sum(sum(sqrt(gx.^2 + gy.^2)));
end

function psnr = compute_psnr(x_true, x_rec)
% COMPUTE_PSNR Compute Peak Signal-to-Noise Ratio
    mse = mean((x_true(:) - x_rec(:)).^2);
    max_val = max(x_true(:));
    psnr = 10 * log10(max_val^2 / mse);
end

function [x, cost_history] = cs_tv_reconstruction(y, mask, lambda, n_iter)
% CS_TV_RECONSTRUCTION Compressed sensing with TV regularization using ADMM
%
%   Solves: min_x ||F_mask(x) - y||^2 + lambda * TV(x)
%
%   Using ADMM with splitting:
%       min ||Ax - y||^2 + lambda * ||z||_1
%       s.t. Dx = z  (D is gradient operator)

    N = size(y, 1);
    
    % ADMM parameters
    rho = 1.0;  % Penalty parameter
    
    % Initialize
    x = real(ifft2(y) * N);  % Zero-filled start
    z = zeros(N, N, 2);      % Auxiliary variable (gradient)
    u = zeros(N, N, 2);      % Dual variable
    
    cost_history = zeros(n_iter, 1);
    
    for iter = 1:n_iter
        % x-update: solve (A'A + rho*D'D)x = A'y + rho*D'(z - u)
        % Using FFT-based solver
        rhs = mask .* y + rho * divergence(z(:,:,1) - u(:,:,1), z(:,:,2) - u(:,:,2));
        x = fft_solve(rhs, mask, rho, N);
        
        % z-update: soft thresholding on gradient
        [gx, gy] = gradient(x);
        z(:,:,1) = soft_threshold_tv(gx + u(:,:,1), lambda/rho);
        z(:,:,2) = soft_threshold_tv(gy + u(:,:,2), lambda/rho);
        
        % u-update: dual ascent
        u(:,:,1) = u(:,:,1) + gx - z(:,:,1);
        u(:,:,2) = u(:,:,2) + gy - z(:,:,2);
        
        % Record cost
        data_fidelity = norm(mask .* (fft2(x)/N - y), 'fro')^2;
        tv_term = lambda * compute_tv(x);
        cost_history(iter) = data_fidelity + tv_term;
    end
end

function x = fft_solve(rhs, mask, rho, N)
% Solve (A'A + rho*D'D)x = rhs in Fourier domain
    % A'A in Fourier: diag(mask)
    % D'D in Fourier: Laplacian ≈ -4 + 2*cos(2πk/N) + 2*cos(2πl/N)
    
    [kx, ky] = meshgrid(0:N-1, 0:N-1);
    laplacian_fft = 4 - 2*cos(2*pi*kx/N) - 2*cos(2*pi*ky/N);
    
    % Solve in Fourier domain
    rhs_fft = fft2(rhs);
    denom = mask + rho * laplacian_fft + 1e-8;  % Regularization
    x_fft = rhs_fft ./ denom;
    x = real(ifft2(x_fft));
end

function div = divergence(gx, gy)
% Compute divergence (negative adjoint of gradient)
    div = -([gx(:,1), diff(gx,1,2)] + [gy(1,:); diff(gy,1,1)]);
end

function z = soft_threshold_tv(x, tau)
% Soft thresholding for TV (isotropic)
    z = sign(x) .* max(abs(x) - tau, 0);
end

function [x, cost_history] = cs_wavelet_reconstruction(y, mask, lambda, n_iter)
% CS_WAVELET_RECONSTRUCTION Compressed sensing with wavelet sparsity
%   Simple ISTA implementation using Haar wavelet

    N = size(y, 1);
    
    % Initialize
    x = real(ifft2(y) * N);
    
    % Step size
    L = 1;
    t = 1 / L;
    
    cost_history = zeros(n_iter, 1);
    
    for iter = 1:n_iter
        % Gradient step (data fidelity term)
        residual = mask .* (fft2(x)/N - y);
        grad = real(ifft2(residual) * N);
        x_grad = x - t * grad;
        
        % Wavelet soft thresholding
        x = wavelet_shrinkage(x_grad, lambda * t);
        
        % Record cost
        data_fidelity = norm(mask .* (fft2(x)/N - y), 'fro')^2;
        cost_history(iter) = data_fidelity + lambda * wavelet_l1(x);
    end
end

function y = wavelet_shrinkage(x, tau)
% Simple 2D Haar wavelet shrinkage
    N = size(x, 1);
    
    % 1-level Haar transform
    [LL, LH, HL, HH] = haar2d(x);
    
    % Shrink detail coefficients
    LH = sign(LH) .* max(abs(LH) - tau, 0);
    HL = sign(HL) .* max(abs(HL) - tau, 0);
    HH = sign(HH) .* max(abs(HH) - tau, 0);
    
    % Inverse transform
    y = ihaar2d(LL, LH, HL, HH);
end

function [LL, LH, HL, HH] = haar2d(x)
% 2D Haar transform (1 level)
    N = size(x, 1);
    M = N/2;
    
    % Row transform
    L = (x(:, 1:2:end) + x(:, 2:2:end)) / sqrt(2);
    H = (x(:, 1:2:end) - x(:, 2:2:end)) / sqrt(2);
    
    % Column transform
    LL = (L(1:2:end, :) + L(2:2:end, :)) / sqrt(2);
    LH = (L(1:2:end, :) - L(2:2:end, :)) / sqrt(2);
    HL = (H(1:2:end, :) + H(2:2:end, :)) / sqrt(2);
    HH = (H(1:2:end, :) - H(2:2:end, :)) / sqrt(2);
end

function x = ihaar2d(LL, LH, HL, HH)
% Inverse 2D Haar transform (1 level)
    M = size(LL, 1);
    N = 2*M;
    
    % Column inverse
    L = zeros(N, M);
    L(1:2:end, :) = (LL + LH) / sqrt(2);
    L(2:2:end, :) = (LL - LH) / sqrt(2);
    
    H = zeros(N, M);
    H(1:2:end, :) = (HL + HH) / sqrt(2);
    H(2:2:end, :) = (HL - HH) / sqrt(2);
    
    % Row inverse
    x = zeros(N, N);
    x(:, 1:2:end) = (L + H) / sqrt(2);
    x(:, 2:2:end) = (L - H) / sqrt(2);
end

function l1 = wavelet_l1(x)
% Compute L1 norm of wavelet coefficients
    [LL, LH, HL, HH] = haar2d(x);
    l1 = sum(abs(LH(:))) + sum(abs(HL(:))) + sum(abs(HH(:)));
end
