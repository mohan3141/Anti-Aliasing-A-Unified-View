%% Compressed Sensing Application: MRI Acceleration
% ============================================================
% Scenario: Accelerate MRI acquisition by undersampling k-space
%           and using CS reconstruction
%
% Key Concepts:
%   - MRI k-space physics
%   - Variable density random sampling
%   - Parallel imaging concepts
%   - Total Variation + Wavelet regularization
%   - Practical considerations (motion, noise)
%
% Author: [Your Name]
% Date: 2024
% ============================================================

clear; clc; close all;

fprintf('==========================================================\n');
fprintf('  Compressed Sensing Application: MRI Acceleration\n');
fprintf('==========================================================\n\n');

%% ==================== MRI Physics Background ====================
fprintf('【MRI Physics Background】\n');
fprintf('  MRI acquires data in k-space (2D Fourier domain)\n');
fprintf('  Each k-space line takes ~10ms to acquire\n');
fprintf('  256×256 image → 256 lines → ~2.5 seconds\n');
fprintf('  \n');
fprintf('  Challenge: Patient motion causes artifacts\n');
fprintf('  Solution: Faster acquisition via undersampling + CS\n\n');

%% ==================== Parameters ====================
N = 256;                    % Image matrix size
acceleration_factors = [2, 4, 6, 8];  % R = 2x, 4x, 6x, 8x acceleration
SNR_dB = 35;               % MRI SNR (typical: 30-40 dB)

% Regularization parameters
lambda_tv = 0.005;
lambda_wav = 0.002;
n_iter = 150;

fprintf('【Configuration】\n');
fprintf('  Image matrix: %d × %d\n', N, N);
fprintf('  Acceleration factors: %s\n', mat2str(acceleration_factors));
fprintf('  SNR: %d dB\n\n', SNR_dB);

%% ==================== Create MRI Brain Phantom ====================
fprintf('【Step 1: Create MRI brain phantom】\n');

x_true = create_brain_phantom(N);

fprintf('  Phantom created: %d × %d\n', N, N);
fprintf('  Tissue contrast simulated\n\n');

%% ==================== Sampling Patterns ====================
fprintf('【Step 2: Design sampling patterns】\n\n');

% Different sampling strategies
sampling_types = {'random_uniform', 'variable_density', 'radial', 'poisson_disk'};

masks = cell(length(acceleration_factors), length(sampling_types));

for ai = 1:length(acceleration_factors)
    R = acceleration_factors(ai);
    rate = 1/R;
    
    fprintf('  R = %dx (%.1f%% sampling):\n', R, 100*rate);
    
    for si = 1:length(sampling_types)
        switch sampling_types{si}
            case 'random_uniform'
                mask = create_random_mask(N, rate);
            case 'variable_density'
                mask = create_variable_density_mask(N, rate, 'polynomial');
            case 'radial'
                mask = create_radial_mask(N, rate);
            case 'poisson_disk'
                mask = create_poisson_disk_mask(N, rate);
        end
        masks{ai, si} = mask;
        
        actual_rate = sum(mask(:)) / N^2;
        fprintf('    %s: %.1f%% actual\n', sampling_types{si}, 100*actual_rate);
    end
    fprintf('\n');
end

%% ==================== Forward Model & Noise ====================
fprintf('【Step 3: Simulate MRI acquisition】\n');

% Full k-space
kspace_full = fft2(x_true) / N;

% Add complex Gaussian noise
noise_power = norm(kspace_full(:))^2 / N^2 / (10^(SNR_dB/10));
kspace_noise = kspace_full + sqrt(noise_power/2) * (randn(N,N) + 1j*randn(N,N));

fprintf('  K-space generated with SNR = %d dB\n\n', SNR_dB);

%% ==================== Reconstruction Comparison ====================
fprintf('【Step 4: Reconstruct from undersampled k-space】\n\n');

% Store results for best sampling pattern (variable density)
si = 2;  % variable_density index

results = struct();
results.psnr = zeros(length(acceleration_factors), 3);  % ZF, TV, TV+Wav
results.ssim = zeros(length(acceleration_factors), 3);
results.time = zeros(length(acceleration_factors), 3);

for ai = 1:length(acceleration_factors)
    R = acceleration_factors(ai);
    mask = masks{ai, si};
    
    fprintf('  R = %dx acceleration:\n', R);
    
    % Undersampled k-space
    kspace_under = kspace_noise .* mask;
    
    % 1. Zero-filled reconstruction
    tic;
    x_zf = abs(ifft2(kspace_under) * N);
    results.time(ai, 1) = toc;
    results.psnr(ai, 1) = compute_psnr(x_true, x_zf);
    results.ssim(ai, 1) = compute_ssim(x_true, x_zf);
    
    % 2. TV reconstruction
    tic;
    [x_tv, ~] = cs_tv_reconstruction(kspace_under, mask, lambda_tv, n_iter);
    results.time(ai, 2) = toc;
    results.psnr(ai, 2) = compute_psnr(x_true, x_tv);
    results.ssim(ai, 2) = compute_ssim(x_true, x_tv);
    
    % 3. TV + Wavelet reconstruction
    tic;
    [x_combined, ~] = cs_combined_reconstruction(kspace_under, mask, lambda_tv, lambda_wav, n_iter);
    results.time(ai, 3) = toc;
    results.psnr(ai, 3) = compute_psnr(x_true, x_combined);
    results.ssim(ai, 3) = compute_ssim(x_true, x_combined);
    
    fprintf('    Zero-fill: PSNR=%.1f dB, SSIM=%.3f\n', results.psnr(ai,1), results.ssim(ai,1));
    fprintf('    TV:        PSNR=%.1f dB, SSIM=%.3f (%.1fs)\n', results.psnr(ai,2), results.ssim(ai,2), results.time(ai,2));
    fprintf('    TV+Wav:    PSNR=%.1f dB, SSIM=%.3f (%.1fs)\n', results.psnr(ai,3), results.ssim(ai,3), results.time(ai,3));
    fprintf('\n');
end

%% ==================== Sampling Pattern Comparison ====================
fprintf('【Step 5: Compare sampling patterns (R=4x)】\n\n');

ai = 2;  % R = 4x
R = acceleration_factors(ai);

pattern_psnr = zeros(length(sampling_types), 1);
pattern_recons = cell(length(sampling_types), 1);

for si = 1:length(sampling_types)
    mask = masks{ai, si};
    kspace_under = kspace_noise .* mask;
    
    [x_rec, ~] = cs_tv_reconstruction(kspace_under, mask, lambda_tv, 100);
    
    pattern_psnr(si) = compute_psnr(x_true, x_rec);
    pattern_recons{si} = x_rec;
    
    fprintf('  %s: PSNR = %.1f dB\n', sampling_types{si}, pattern_psnr(si));
end
fprintf('\n');

%% ==================== Visualization ====================
fprintf('【Generating figures...】\n\n');

% Figure 1: Main results
figure('Position', [50, 50, 1600, 900]);

% Row 1: Ground truth and sampling masks
subplot(2, 4, 1);
imagesc(x_true); axis image; colormap(gray); colorbar;
title('Ground Truth');

subplot(2, 4, 2);
imagesc(fftshift(masks{2, 2})); axis image; colormap(gray);
title(sprintf('Var. Density Mask (R=%dx)', acceleration_factors(2)));

subplot(2, 4, 3);
imagesc(fftshift(masks{3, 2})); axis image; colormap(gray);
title(sprintf('Var. Density Mask (R=%dx)', acceleration_factors(3)));

subplot(2, 4, 4);
imagesc(fftshift(masks{4, 2})); axis image; colormap(gray);
title(sprintf('Var. Density Mask (R=%dx)', acceleration_factors(4)));

% Row 2: Reconstructions at R=4x
subplot(2, 4, 5);
mask_4x = masks{2, 2};
x_zf_4x = abs(ifft2(kspace_noise .* mask_4x) * N);
imagesc(x_zf_4x); axis image; colormap(gray); colorbar;
title(sprintf('R=4x Zero-fill (PSNR=%.1f)', results.psnr(2,1)));

subplot(2, 4, 6);
[x_tv_4x, ~] = cs_tv_reconstruction(kspace_noise .* mask_4x, mask_4x, lambda_tv, n_iter);
imagesc(x_tv_4x); axis image; colormap(gray); colorbar;
title(sprintf('R=4x TV Recon (PSNR=%.1f)', results.psnr(2,2)));

subplot(2, 4, 7);
[x_comb_4x, ~] = cs_combined_reconstruction(kspace_noise .* mask_4x, mask_4x, lambda_tv, lambda_wav, n_iter);
imagesc(x_comb_4x); axis image; colormap(gray); colorbar;
title(sprintf('R=4x TV+Wav (PSNR=%.1f)', results.psnr(2,3)));

subplot(2, 4, 8);
plot(acceleration_factors, results.psnr(:, 1), 'r--o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(acceleration_factors, results.psnr(:, 2), 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
plot(acceleration_factors, results.psnr(:, 3), 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Acceleration Factor R');
ylabel('PSNR (dB)');
legend('Zero-fill', 'TV', 'TV+Wavelet', 'Location', 'best');
title('PSNR vs Acceleration');
grid on;
xlim([1.5 8.5]);

sgtitle('CS-MRI: Effect of Acceleration Factor', 'FontSize', 14);
saveas(gcf, 'app_mri_acceleration_main.png');
fprintf('  Figure saved: app_mri_acceleration_main.png\n');

% Figure 2: Sampling pattern comparison
figure('Position', [100, 100, 1400, 700]);

for si = 1:4
    % Mask
    subplot(2, 4, si);
    imagesc(fftshift(masks{2, si})); axis image; colormap(gray);
    title(sprintf('%s\n(R=4x)', strrep(sampling_types{si}, '_', ' ')));
    
    % Reconstruction
    subplot(2, 4, si + 4);
    imagesc(pattern_recons{si}); axis image; colormap(gray);
    title(sprintf('PSNR = %.1f dB', pattern_psnr(si)));
end

sgtitle('Sampling Pattern Comparison (R=4x)', 'FontSize', 14);
saveas(gcf, 'app_mri_sampling_comparison.png');
fprintf('  Figure saved: app_mri_sampling_comparison.png\n');

%% ==================== Practical Considerations ====================
fprintf('\n【Practical Considerations for Clinical MRI】\n');
fprintf('  \n');
fprintf('  1. Sampling trajectory must be physically realizable\n');
fprintf('     - Gradient slew rate limits\n');
fprintf('     - k-space velocity constraints\n');
fprintf('  \n');
fprintf('  2. Calibration region for SENSE/GRAPPA\n');
fprintf('     - Center of k-space typically fully sampled\n');
fprintf('     - Provides coil sensitivity estimates\n');
fprintf('  \n');
fprintf('  3. Motion artifacts\n');
fprintf('     - Respiratory/cardiac gating\n');
fprintf('     - Navigator echoes\n');
fprintf('  \n');
fprintf('  4. Regularization parameter selection\n');
fprintf('     - Too high: over-smoothing, loss of detail\n');
fprintf('     - Too low: residual aliasing\n');
fprintf('     - Often selected via L-curve or cross-validation\n');
fprintf('  \n');
fprintf('  5. Reconstruction time\n');
fprintf('     - Iterative methods can take seconds to minutes\n');
fprintf('     - GPU acceleration essential for clinical use\n');
fprintf('     - Deep learning approaches now <1 second\n\n');

%% ==================== Summary ====================
fprintf('==========================================================\n');
fprintf('  SUMMARY: CS-MRI Acceleration\n');
fprintf('==========================================================\n');
fprintf('  \n');
fprintf('  Image: %d × %d brain phantom\n', N, N);
fprintf('  Best sampling: Variable density\n');
fprintf('  \n');
fprintf('  Results (Variable Density Sampling):\n');
fprintf('    R  | Zero-fill | TV Recon  | TV+Wavelet\n');
fprintf('    ---|-----------|-----------|----------\n');
for ai = 1:length(acceleration_factors)
    fprintf('    %dx | %5.1f dB  | %5.1f dB  | %5.1f dB\n', ...
        acceleration_factors(ai), results.psnr(ai, 1), results.psnr(ai, 2), results.psnr(ai, 3));
end
fprintf('  \n');
fprintf('  Key findings:\n');
fprintf('  - CS enables 4-8x acceleration with good image quality\n');
fprintf('  - Variable density sampling outperforms uniform random\n');
fprintf('  - Combined TV+Wavelet improves fine structure preservation\n');
fprintf('==========================================================\n');

%% ==================== Local Functions ====================

function x = create_brain_phantom(N)
% CREATE_BRAIN_PHANTOM Create a simplified brain phantom
    x = zeros(N, N);
    [X, Y] = meshgrid(1:N, 1:N);
    cx = N/2; cy = N/2;
    
    % Skull (outer ellipse)
    a1 = N*0.42; b1 = N*0.48;
    mask_skull = ((X-cx)/a1).^2 + ((Y-cy)/b1).^2 <= 1;
    x(mask_skull) = 0.3;  % Bone signal (low in T1)
    
    % Brain (inner ellipse)
    a2 = N*0.38; b2 = N*0.44;
    mask_brain = ((X-cx)/a2).^2 + ((Y-cy)/b2).^2 <= 1;
    x(mask_brain) = 0.7;  % Gray matter
    
    % White matter (central region)
    a3 = N*0.25; b3 = N*0.30;
    mask_wm = ((X-cx)/a3).^2 + ((Y-cy)/b3).^2 <= 1;
    x(mask_wm) = 0.9;  % White matter (bright in T1)
    
    % Ventricles (dark CSF)
    % Lateral ventricles
    v1_cx = cx - N*0.08; v1_cy = cy - N*0.05;
    v1_a = N*0.06; v1_b = N*0.12;
    mask_v1 = ((X-v1_cx)/v1_a).^2 + ((Y-v1_cy)/v1_b).^2 <= 1;
    x(mask_v1) = 0.1;  % CSF (dark in T1)
    
    v2_cx = cx + N*0.08; v2_cy = cy - N*0.05;
    mask_v2 = ((X-v2_cx)/v1_a).^2 + ((Y-v2_cy)/v1_b).^2 <= 1;
    x(mask_v2) = 0.1;
    
    % Third ventricle
    v3_a = N*0.015; v3_b = N*0.08;
    mask_v3 = ((X-cx)/v3_a).^2 + ((Y-cy)/v3_b).^2 <= 1;
    x(mask_v3) = 0.1;
    
    % Some lesions/tumors
    lesion_r = N*0.03;
    lx = cx + N*0.15; ly = cy + N*0.1;
    mask_lesion = (X-lx).^2 + (Y-ly).^2 <= lesion_r^2;
    x(mask_lesion) = 1.0;  % Bright lesion
    
    % Add some texture/noise
    x = x + 0.02 * randn(N, N);
    x = max(0, min(1, x));
    
    % Smooth slightly
    h = fspecial('gaussian', 5, 1);
    x = conv2(x, h, 'same');
end

function mask = create_random_mask(N, rate)
% Uniform random sampling
    mask = rand(N, N) < rate;
    mask(1, 1) = 1;  % Always include DC
end

function mask = create_variable_density_mask(N, rate, type)
% Variable density sampling (higher at center)
    [X, Y] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
    R = sqrt(X.^2 + Y.^2) / (N/2);
    
    switch type
        case 'polynomial'
            density = (1 - R.^3).^2;
        case 'exponential'
            density = exp(-2*R);
        otherwise
            density = 1 - R;
    end
    
    density = max(density, 0);
    density = density / sum(density(:)) * rate * N^2;
    density = min(density, 1);
    
    mask = fftshift(rand(N, N) < density);
    mask(1, 1) = 1;
    
    % Ensure center rows fully sampled (calibration)
    cal_size = round(N * 0.08);
    mask(1:cal_size, :) = 1;
    mask(:, 1:cal_size) = 1;
end

function mask = create_radial_mask(N, rate)
% Radial sampling (like radial MRI)
    mask = zeros(N, N);
    n_spokes = round(rate * N * pi / 2);  % Approximate
    
    [X, Y] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
    
    angles = linspace(0, pi, n_spokes);
    for theta = angles
        for r = 0:N/2
            x_idx = round(N/2 + 1 + r * cos(theta));
            y_idx = round(N/2 + 1 + r * sin(theta));
            if x_idx >= 1 && x_idx <= N && y_idx >= 1 && y_idx <= N
                mask(y_idx, x_idx) = 1;
            end
            % Also negative radius
            x_idx = round(N/2 + 1 - r * cos(theta));
            y_idx = round(N/2 + 1 - r * sin(theta));
            if x_idx >= 1 && x_idx <= N && y_idx >= 1 && y_idx <= N
                mask(y_idx, x_idx) = 1;
            end
        end
    end
    
    mask = fftshift(mask);
end

function mask = create_poisson_disk_mask(N, rate)
% Poisson disk sampling (minimum distance constraint)
    mask = zeros(N, N);
    
    % Target number of samples
    target_samples = round(rate * N^2);
    
    % Minimum distance (approximate)
    min_dist = sqrt(N^2 / target_samples) * 0.7;
    
    % Start with center
    points = [N/2, N/2];
    mask(N/2, N/2) = 1;
    
    % Active list
    active = 1;
    
    max_attempts = 30;
    
    while ~isempty(active) && size(points, 1) < target_samples
        % Pick random active point
        idx = active(randi(length(active)));
        p = points(idx, :);
        
        found = false;
        for attempt = 1:max_attempts
            % Random point in annulus [min_dist, 2*min_dist]
            theta = 2*pi*rand();
            r = min_dist * (1 + rand());
            
            new_p = p + r * [cos(theta), sin(theta)];
            
            % Check bounds
            if new_p(1) < 1 || new_p(1) > N || new_p(2) < 1 || new_p(2) > N
                continue;
            end
            
            % Check distance to all existing points
            dists = sqrt(sum((points - new_p).^2, 2));
            if all(dists >= min_dist * 0.9)
                points = [points; new_p];
                active = [active, size(points, 1)];
                mask(round(new_p(2)), round(new_p(1))) = 1;
                found = true;
                break;
            end
        end
        
        if ~found
            active(active == idx) = [];
        end
    end
    
    mask = fftshift(mask);
    mask(1, 1) = 1;
end

function psnr = compute_psnr(x_true, x_rec)
    mse = mean((x_true(:) - x_rec(:)).^2);
    max_val = max(x_true(:));
    psnr = 10 * log10(max_val^2 / mse);
end

function ssim_val = compute_ssim(x, y)
% Simplified SSIM computation
    C1 = 0.01^2;
    C2 = 0.03^2;
    
    mu_x = mean(x(:));
    mu_y = mean(y(:));
    
    sigma_x = std(x(:));
    sigma_y = std(y(:));
    sigma_xy = mean((x(:) - mu_x) .* (y(:) - mu_y));
    
    ssim_val = (2*mu_x*mu_y + C1) * (2*sigma_xy + C2) / ...
               ((mu_x^2 + mu_y^2 + C1) * (sigma_x^2 + sigma_y^2 + C2));
end

function [x, cost_history] = cs_tv_reconstruction(y, mask, lambda, n_iter)
% TV-regularized CS reconstruction
    N = size(y, 1);
    rho = 1.0;
    
    x = real(ifft2(y) * N);
    z = zeros(N, N, 2);
    u = zeros(N, N, 2);
    
    cost_history = zeros(n_iter, 1);
    
    for iter = 1:n_iter
        rhs = mask .* y + rho * divergence(z(:,:,1) - u(:,:,1), z(:,:,2) - u(:,:,2));
        x = fft_solve(rhs, mask, rho, N);
        
        [gx, gy] = gradient(x);
        z(:,:,1) = soft_threshold_tv(gx + u(:,:,1), lambda/rho);
        z(:,:,2) = soft_threshold_tv(gy + u(:,:,2), lambda/rho);
        
        u(:,:,1) = u(:,:,1) + gx - z(:,:,1);
        u(:,:,2) = u(:,:,2) + gy - z(:,:,2);
        
        cost_history(iter) = norm(mask .* (fft2(x)/N - y), 'fro')^2 + lambda * compute_tv(x);
    end
end

function [x, cost_history] = cs_combined_reconstruction(y, mask, lambda_tv, lambda_wav, n_iter)
% Combined TV + Wavelet reconstruction
    N = size(y, 1);
    
    x = real(ifft2(y) * N);
    cost_history = zeros(n_iter, 1);
    
    for iter = 1:n_iter
        % Gradient step
        residual = mask .* (fft2(x)/N - y);
        grad = real(ifft2(residual) * N);
        
        % TV denoising step
        x_tv = tv_denoise(x - 0.5*grad, lambda_tv);
        
        % Wavelet shrinkage step
        x = wavelet_shrinkage(x_tv, lambda_wav);
        
        cost_history(iter) = norm(mask .* (fft2(x)/N - y), 'fro')^2;
    end
end

function x_out = tv_denoise(x, lambda)
% Simple TV denoising (1 iteration of Chambolle)
    [gx, gy] = gradient(x);
    mag = sqrt(gx.^2 + gy.^2 + 1e-8);
    gx = gx ./ mag;
    gy = gy ./ mag;
    div_p = divergence(gx, gy);
    x_out = x - lambda * div_p;
end

function y = wavelet_shrinkage(x, tau)
    N = size(x, 1);
    [LL, LH, HL, HH] = haar2d(x);
    LH = sign(LH) .* max(abs(LH) - tau, 0);
    HL = sign(HL) .* max(abs(HL) - tau, 0);
    HH = sign(HH) .* max(abs(HH) - tau, 0);
    y = ihaar2d(LL, LH, HL, HH);
end

function x = fft_solve(rhs, mask, rho, N)
    [kx, ky] = meshgrid(0:N-1, 0:N-1);
    laplacian_fft = 4 - 2*cos(2*pi*kx/N) - 2*cos(2*pi*ky/N);
    rhs_fft = fft2(rhs);
    denom = mask + rho * laplacian_fft + 1e-8;
    x_fft = rhs_fft ./ denom;
    x = real(ifft2(x_fft));
end

function div = divergence(gx, gy)
    div = -([gx(:,1), diff(gx,1,2)] + [gy(1,:); diff(gy,1,1)]);
end

function z = soft_threshold_tv(x, tau)
    z = sign(x) .* max(abs(x) - tau, 0);
end

function tv = compute_tv(x)
    [gx, gy] = gradient(x);
    tv = sum(sum(sqrt(gx.^2 + gy.^2)));
end

function [LL, LH, HL, HH] = haar2d(x)
    L = (x(:, 1:2:end) + x(:, 2:2:end)) / sqrt(2);
    H = (x(:, 1:2:end) - x(:, 2:2:end)) / sqrt(2);
    LL = (L(1:2:end, :) + L(2:2:end, :)) / sqrt(2);
    LH = (L(1:2:end, :) - L(2:2:end, :)) / sqrt(2);
    HL = (H(1:2:end, :) + H(2:2:end, :)) / sqrt(2);
    HH = (H(1:2:end, :) - H(2:2:end, :)) / sqrt(2);
end

function x = ihaar2d(LL, LH, HL, HH)
    M = size(LL, 1);
    N = 2*M;
    L = zeros(N, M);
    L(1:2:end, :) = (LL + LH) / sqrt(2);
    L(2:2:end, :) = (LL - LH) / sqrt(2);
    H = zeros(N, M);
    H(1:2:end, :) = (HL + HH) / sqrt(2);
    H(2:2:end, :) = (HL - HH) / sqrt(2);
    x = zeros(N, N);
    x(:, 1:2:end) = (L + H) / sqrt(2);
    x(:, 2:2:end) = (L - H) / sqrt(2);
end
