%% Compressed Sensing Application: Radar Sparse Aperture Imaging
% ============================================================
% Scenario: ISAR/SAR imaging with sparse aperture using CS
%
% Key Concepts:
%   - Radar imaging fundamentals (Range-Doppler)
%   - Sparse aperture problem
%   - Point spread function (PSF) with grating lobes
%   - CS recovery for sparse targets
%   - Comparison with traditional matched filter
%
% Author: [Your Name]
% Date: 2024
% ============================================================

clear; clc; close all;

fprintf('==========================================================\n');
fprintf('  CS Application: Radar Sparse Aperture Imaging\n');
fprintf('==========================================================\n\n');

%% ==================== Radar Basics ====================
fprintf('【Radar Imaging Background】\n');
fprintf('  SAR/ISAR forms images via synthetic aperture\n');
fprintf('  - Range: from pulse compression (chirp bandwidth)\n');
fprintf('  - Cross-range: from Doppler history (aperture length)\n');
fprintf('  \n');
fprintf('  Sparse aperture problem:\n');
fprintf('  - Missing radar pulses → grating lobes\n');
fprintf('  - Causes: sensor failure, jamming, maneuvering\n');
fprintf('  - CS can suppress grating lobes for sparse scenes\n\n');

%% ==================== System Parameters ====================
fprintf('【System Parameters】\n');

% Radar parameters
c = 3e8;                    % Speed of light (m/s)
fc = 10e9;                  % Carrier frequency (X-band, 10 GHz)
lambda = c / fc;            % Wavelength
B = 500e6;                  % Bandwidth (500 MHz)
PRF = 1000;                 % Pulse repetition frequency (Hz)

% Imaging parameters
N_range = 256;              % Range bins
N_azimuth = 256;            % Azimuth bins (pulses)
range_res = c / (2 * B);    % Range resolution
cross_range_res = lambda / (2 * N_azimuth / PRF * 10);  % Simplified

% Sparse aperture
aperture_fill = 0.3;        % 30% of pulses available

fprintf('  Carrier frequency: %.1f GHz\n', fc/1e9);
fprintf('  Bandwidth: %.0f MHz\n', B/1e6);
fprintf('  Range resolution: %.2f m\n', range_res);
fprintf('  Full aperture: %d pulses\n', N_azimuth);
fprintf('  Sparse aperture fill: %.0f%%\n', 100*aperture_fill);
fprintf('\n');

%% ==================== Generate Target Scene ====================
fprintf('【Step 1: Generate sparse target scene】\n');

% Create sparse scene (few point targets)
n_targets = 8;
rng(42);

scene = zeros(N_range, N_azimuth);
target_positions = zeros(n_targets, 2);
target_amplitudes = zeros(n_targets, 1);

for i = 1:n_targets
    % Random position (avoid edges)
    r = randi([30, N_range-30]);
    a = randi([30, N_azimuth-30]);
    
    % Random amplitude
    amp = 0.5 + 0.5 * rand();
    
    scene(r, a) = amp;
    target_positions(i, :) = [r, a];
    target_amplitudes(i) = amp;
end

fprintf('  Number of targets: %d\n', n_targets);
fprintf('  Scene sparsity: %.4f%%\n', 100*n_targets/(N_range*N_azimuth));
fprintf('  Target positions (range, azimuth):\n');
for i = 1:n_targets
    fprintf('    Target %d: (%d, %d), amplitude = %.2f\n', ...
        i, target_positions(i,1), target_positions(i,2), target_amplitudes(i));
end
fprintf('\n');

%% ==================== Full Aperture Imaging ====================
fprintf('【Step 2: Full aperture imaging (reference)】\n');

% Simulate radar data (simplified model)
% Phase history = 2D Fourier transform of scene
full_phase_history = fft2(scene);

% Full aperture image (2D IFFT)
image_full = ifft2(full_phase_history);
image_full_abs = abs(image_full);

% Normalize
image_full_abs = image_full_abs / max(image_full_abs(:));

fprintf('  Full aperture image generated\n');
fprintf('  Peak sidelobe level: %.1f dB\n', ...
    20*log10(max_sidelobe(image_full_abs, target_positions)));
fprintf('\n');

%% ==================== Sparse Aperture (Random Missing Pulses) ====================
fprintf('【Step 3: Sparse aperture acquisition】\n');

% Random pulse selection (simulate missing pulses)
pulse_mask = zeros(1, N_azimuth);
n_pulses = round(aperture_fill * N_azimuth);
selected_pulses = sort(randperm(N_azimuth, n_pulses));
pulse_mask(selected_pulses) = 1;

% 2D mask (all range bins, selected azimuth pulses)
mask_2d = repmat(pulse_mask, N_range, 1);

% Sparse phase history
sparse_phase_history = full_phase_history .* mask_2d;

fprintf('  Pulses acquired: %d / %d (%.0f%%)\n', n_pulses, N_azimuth, 100*aperture_fill);
fprintf('  Missing pulses: %d\n', N_azimuth - n_pulses);
fprintf('\n');

%% ==================== Traditional Imaging (Matched Filter) ====================
fprintf('【Step 4: Traditional imaging (zero-fill + FFT)】\n');

% Zero-filled reconstruction (matched filter)
image_mf = ifft2(sparse_phase_history);
image_mf_abs = abs(image_mf);
image_mf_abs = image_mf_abs / max(image_mf_abs(:));

% Analyze grating lobes
psl_mf = 20*log10(max_sidelobe(image_mf_abs, target_positions));

fprintf('  Matched filter image generated\n');
fprintf('  Peak sidelobe level: %.1f dB (grating lobes!)\n', psl_mf);
fprintf('  → Significant artifacts due to missing pulses\n\n');

%% ==================== CS Reconstruction ====================
fprintf('【Step 5: Compressed sensing reconstruction】\n');

% CS parameters
lambda_cs = 0.01;
n_iter = 200;

fprintf('  Running CS optimization (L1-regularized)...\n');
tic;

% Vectorized problem: y = A*x, recover sparse x
% A = partial Fourier matrix (selected rows)
% Using ADMM for efficiency

[image_cs, cost_history] = cs_radar_reconstruction(...
    sparse_phase_history, mask_2d, lambda_cs, n_iter);

time_cs = toc;

image_cs_abs = abs(image_cs);
image_cs_abs = image_cs_abs / max(image_cs_abs(:));

% Analyze
psl_cs = 20*log10(max_sidelobe(image_cs_abs, target_positions));
target_error = compute_target_localization_error(image_cs_abs, target_positions);

fprintf('  CS reconstruction complete\n');
fprintf('  Time: %.2f seconds\n', time_cs);
fprintf('  Peak sidelobe level: %.1f dB\n', psl_cs);
fprintf('  Sidelobe improvement: %.1f dB\n', psl_mf - psl_cs);
fprintf('  Target localization RMSE: %.2f pixels\n\n', target_error);

%% ==================== OMP Reconstruction ====================
fprintf('【Step 6: OMP (greedy) reconstruction】\n');

fprintf('  Running OMP for comparison...\n');
tic;

[image_omp, ~] = omp_radar_reconstruction(...
    sparse_phase_history, mask_2d, n_targets * 2);  % Allow some extra

time_omp = toc;

image_omp_abs = abs(image_omp);
image_omp_abs = image_omp_abs / max(image_omp_abs(:) + 1e-10);

psl_omp = 20*log10(max_sidelobe(image_omp_abs, target_positions) + 1e-10);

fprintf('  OMP reconstruction complete\n');
fprintf('  Time: %.2f seconds\n', time_omp);
fprintf('  Peak sidelobe level: %.1f dB\n', psl_omp);
fprintf('\n');

%% ==================== Effect of Aperture Fill ====================
fprintf('【Step 7: Effect of aperture fill ratio】\n');

fill_ratios = [0.15, 0.20, 0.30, 0.40, 0.50, 0.70];
psl_vs_fill = zeros(length(fill_ratios), 2);  % MF, CS

for i = 1:length(fill_ratios)
    fill = fill_ratios(i);
    n_p = round(fill * N_azimuth);
    sel_p = sort(randperm(N_azimuth, n_p));
    
    p_mask = zeros(1, N_azimuth);
    p_mask(sel_p) = 1;
    m_2d = repmat(p_mask, N_range, 1);
    
    sparse_ph = full_phase_history .* m_2d;
    
    % Matched filter
    img_mf = abs(ifft2(sparse_ph));
    img_mf = img_mf / max(img_mf(:));
    psl_vs_fill(i, 1) = 20*log10(max_sidelobe(img_mf, target_positions));
    
    % CS
    [img_cs, ~] = cs_radar_reconstruction(sparse_ph, m_2d, lambda_cs, 100);
    img_cs = abs(img_cs);
    img_cs = img_cs / max(img_cs(:));
    psl_vs_fill(i, 2) = 20*log10(max_sidelobe(img_cs, target_positions));
    
    fprintf('  Fill %.0f%%: MF PSL = %.1f dB, CS PSL = %.1f dB\n', ...
        100*fill, psl_vs_fill(i, 1), psl_vs_fill(i, 2));
end
fprintf('\n');

%% ==================== Visualization ====================
fprintf('【Generating figures...】\n\n');

figure('Position', [50, 50, 1600, 900]);

% Row 1: Ground truth, full aperture, sparse aperture mask
subplot(2, 4, 1);
imagesc(scene); axis image; colorbar;
title(sprintf('Ground Truth (%d targets)', n_targets));
xlabel('Azimuth'); ylabel('Range');

subplot(2, 4, 2);
imagesc(20*log10(image_full_abs + 1e-6)); axis image; colorbar;
title('Full Aperture Image');
xlabel('Azimuth'); ylabel('Range');
caxis([-40 0]);

subplot(2, 4, 3);
imagesc(pulse_mask); colormap(gca, gray);
title(sprintf('Aperture Mask (%.0f%%)', 100*aperture_fill));
xlabel('Pulse Index'); ylabel('');
yticks([]);

subplot(2, 4, 4);
plot(cost_history, 'b-', 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('Cost');
title('CS Convergence');
grid on;

% Row 2: Matched filter, CS, OMP, comparison
subplot(2, 4, 5);
imagesc(20*log10(image_mf_abs + 1e-6)); axis image; colorbar;
title(sprintf('Matched Filter (PSL=%.1fdB)', psl_mf));
xlabel('Azimuth'); ylabel('Range');
caxis([-40 0]);

subplot(2, 4, 6);
imagesc(20*log10(image_cs_abs + 1e-6)); axis image; colorbar;
title(sprintf('CS Recon (PSL=%.1fdB)', psl_cs));
xlabel('Azimuth'); ylabel('Range');
caxis([-40 0]);

subplot(2, 4, 7);
imagesc(20*log10(image_omp_abs + 1e-6)); axis image; colorbar;
title(sprintf('OMP Recon (PSL=%.1fdB)', psl_omp));
xlabel('Azimuth'); ylabel('Range');
caxis([-40 0]);

subplot(2, 4, 8);
plot(100*fill_ratios, psl_vs_fill(:,1), 'r--o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(100*fill_ratios, psl_vs_fill(:,2), 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Aperture Fill (%)');
ylabel('Peak Sidelobe Level (dB)');
legend('Matched Filter', 'CS', 'Location', 'best');
title('PSL vs Aperture Fill');
grid on;

sgtitle('Compressed Sensing: Radar Sparse Aperture Imaging', 'FontSize', 14);
saveas(gcf, 'app_radar_imaging_results.png');
fprintf('  Figure saved: app_radar_imaging_results.png\n');

% Figure 2: Cross-sections
figure('Position', [100, 100, 1200, 400]);

% Find strongest target
[~, max_idx] = max(scene(:));
[max_r, max_a] = ind2sub(size(scene), max_idx);

subplot(1, 3, 1);
plot(1:N_azimuth, 20*log10(image_full_abs(max_r, :) + 1e-6), 'k-', 'LineWidth', 1.5);
hold on;
plot(1:N_azimuth, 20*log10(image_mf_abs(max_r, :) + 1e-6), 'r--', 'LineWidth', 1.5);
xlabel('Azimuth'); ylabel('Magnitude (dB)');
title('Azimuth Cut: Full vs MF');
legend('Full Aperture', 'Matched Filter');
ylim([-60 5]); grid on;

subplot(1, 3, 2);
plot(1:N_azimuth, 20*log10(image_full_abs(max_r, :) + 1e-6), 'k-', 'LineWidth', 1.5);
hold on;
plot(1:N_azimuth, 20*log10(image_cs_abs(max_r, :) + 1e-6), 'b-', 'LineWidth', 1.5);
xlabel('Azimuth'); ylabel('Magnitude (dB)');
title('Azimuth Cut: Full vs CS');
legend('Full Aperture', 'CS');
ylim([-60 5]); grid on;

subplot(1, 3, 3);
plot(1:N_range, 20*log10(image_mf_abs(:, max_a) + 1e-6), 'r--', 'LineWidth', 1.5);
hold on;
plot(1:N_range, 20*log10(image_cs_abs(:, max_a) + 1e-6), 'b-', 'LineWidth', 1.5);
xlabel('Range'); ylabel('Magnitude (dB)');
title('Range Cut: MF vs CS');
legend('Matched Filter', 'CS');
ylim([-60 5]); grid on;

sgtitle('Cross-section Analysis', 'FontSize', 12);
saveas(gcf, 'app_radar_imaging_cuts.png');
fprintf('  Figure saved: app_radar_imaging_cuts.png\n');

%% ==================== Summary ====================
fprintf('\n==========================================================\n');
fprintf('  SUMMARY: Radar Sparse Aperture Imaging\n');
fprintf('==========================================================\n');
fprintf('  \n');
fprintf('  Scenario: %d × %d image, %d targets\n', N_range, N_azimuth, n_targets);
fprintf('  Aperture fill: %.0f%%\n', 100*aperture_fill);
fprintf('  \n');
fprintf('  Results:\n');
fprintf('    Method        | PSL (dB) | Time (s)\n');
fprintf('    --------------|----------|----------\n');
fprintf('    Full aperture | %8.1f | N/A\n', 20*log10(max_sidelobe(image_full_abs, target_positions)));
fprintf('    Matched Filter| %8.1f | instant\n', psl_mf);
fprintf('    CS (L1-ADMM)  | %8.1f | %.2f\n', psl_cs, time_cs);
fprintf('    OMP           | %8.1f | %.2f\n', psl_omp, time_omp);
fprintf('  \n');
fprintf('  Key findings:\n');
fprintf('  - Sparse aperture causes severe grating lobes in MF\n');
fprintf('  - CS reduces sidelobes by %.1f dB\n', psl_mf - psl_cs);
fprintf('  - Works well when scene is sparse (few targets)\n');
fprintf('  - Performance degrades for dense scenes\n');
fprintf('==========================================================\n');

%% ==================== Local Functions ====================

function psl = max_sidelobe(image, target_positions)
% Find peak sidelobe level (excluding mainlobes)
    image_copy = image;
    
    % Mask out mainlobes (3x3 region around each target)
    for i = 1:size(target_positions, 1)
        r = target_positions(i, 1);
        a = target_positions(i, 2);
        r_range = max(1, r-2):min(size(image,1), r+2);
        a_range = max(1, a-2):min(size(image,2), a+2);
        image_copy(r_range, a_range) = 0;
    end
    
    psl = max(image_copy(:));
end

function rmse = compute_target_localization_error(image, true_positions)
% Compute RMSE of target localization
    n_targets = size(true_positions, 1);
    errors = zeros(n_targets, 1);
    
    for i = 1:n_targets
        % Find local maximum near true position
        r_true = true_positions(i, 1);
        a_true = true_positions(i, 2);
        
        r_range = max(1, r_true-5):min(size(image,1), r_true+5);
        a_range = max(1, a_true-5):min(size(image,2), a_true+5);
        
        local_patch = image(r_range, a_range);
        [~, max_idx] = max(local_patch(:));
        [dr, da] = ind2sub(size(local_patch), max_idx);
        
        r_est = r_range(1) + dr - 1;
        a_est = a_range(1) + da - 1;
        
        errors(i) = sqrt((r_est - r_true)^2 + (a_est - a_true)^2);
    end
    
    rmse = sqrt(mean(errors.^2));
end

function [image, cost_history] = cs_radar_reconstruction(y, mask, lambda, n_iter)
% CS reconstruction for radar imaging using ADMM
%
%   Solves: min_x ||F_mask(x) - y||^2 + lambda * ||x||_1
%   where F is 2D FFT

    [N_range, N_azimuth] = size(y);
    
    % ADMM parameters
    rho = 1.0;
    
    % Initialize
    x = ifft2(y);  % Zero-filled start
    z = zeros(N_range, N_azimuth);  % Auxiliary
    u = zeros(N_range, N_azimuth);  % Dual
    
    cost_history = zeros(n_iter, 1);
    
    for iter = 1:n_iter
        % x-update: solve (A'A + rho*I)x = A'y + rho*(z - u)
        % In Fourier domain: (mask + rho).*X_fft = mask.*Y + rho*fft2(z-u)
        rhs_fft = mask .* y + rho * fft2(z - u);
        x_fft = rhs_fft ./ (mask + rho);
        x = ifft2(x_fft);
        
        % z-update: soft thresholding
        z = soft_threshold(x + u, lambda/rho);
        
        % u-update
        u = u + x - z;
        
        % Cost
        cost_history(iter) = norm(mask .* (fft2(x) - y), 'fro')^2 + lambda * norm(x(:), 1);
    end
    
    image = x;
end

function [image, support] = omp_radar_reconstruction(y, mask, K)
% OMP for radar imaging (simplified, operates on vectorized problem)
    [N_range, N_azimuth] = size(y);
    N = N_range * N_azimuth;
    
    % Measurement indices
    meas_idx = find(mask(:));
    M = length(meas_idx);
    
    % Measurement vector
    y_vec = y(meas_idx);
    
    % Initialize
    residual = y_vec;
    support = [];
    image = zeros(N_range, N_azimuth);
    
    % Precompute for efficiency (partial Fourier)
    % A(i,j) = exp(-2*pi*1j * (k_r(i)*r(j)/N_range + k_a(i)*a(j)/N_azimuth))
    % This is expensive, so we use FFT tricks
    
    for k = 1:K
        % Correlate residual with all columns (via IFFT)
        temp = zeros(N_range, N_azimuth);
        temp(meas_idx) = residual;
        correlation = ifft2(temp);
        
        % Find max (excluding support)
        corr_abs = abs(correlation);
        corr_abs(support) = 0;
        [~, j_star] = max(corr_abs(:));
        
        % Update support
        support = [support; j_star];
        
        % Solve least squares (simplified)
        % Project y onto selected columns
        image_temp = zeros(N_range, N_azimuth);
        [r_idx, a_idx] = ind2sub([N_range, N_azimuth], support);
        
        % Build small system
        A_small = zeros(M, length(support));
        for si = 1:length(support)
            ei = zeros(N_range, N_azimuth);
            ei(support(si)) = 1;
            Fei = fft2(ei);
            A_small(:, si) = Fei(meas_idx);
        end
        
        % Least squares
        coef = A_small \ y_vec;
        
        % Update residual
        residual = y_vec - A_small * coef;
        
        % Store result
        for si = 1:length(support)
            image(support(si)) = coef(si);
        end
    end
end

function z = soft_threshold(x, tau)
    z = sign(x) .* max(abs(x) - tau, 0);
end
