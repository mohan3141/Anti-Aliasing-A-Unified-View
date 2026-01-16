%% CRT Anti-Aliasing Demo: 2D Phase Unwrapping
% ============================================
% Scenario: InSAR / Multi-baseline Interferometry
%
% Problem: Interferometric phase is wrapped to [-pi, pi]
%   True phase: phi_true = 4*pi*B*h / (lambda*R*sin(theta))
%   Observed:   phi_wrap = mod(phi_true + pi, 2*pi) - pi
%
% Solution: Use multiple baselines (like multiple sample rates!)
%   Baseline B1 -> wrapped phase phi1
%   Baseline B2 -> wrapped phase phi2
%   CRT -> unwrapped true phase
%
% Author: [Your Name]
% Date: 2024
% ============================================

clear; clc; close all;

fprintf('============================================\n');
fprintf('   CRT Anti-Aliasing: 2D Phase Unwrapping\n');
fprintf('============================================\n');

%% ==================== Simulation Parameters ====================
% Image size
Nx = 256;
Ny = 256;
[X, Y] = meshgrid(1:Nx, 1:Ny);

% Create synthetic terrain (elevation map)
% Combination of slopes and features
h_true = 50 * peaks(Nx/4);  % Mountain-like features
h_true = imresize(h_true, [Ny, Nx]);
h_true = h_true + 0.5 * X + 0.3 * Y;  % Add overall slope
h_true = h_true - min(h_true(:));      % Make non-negative

% InSAR parameters
lambda = 0.056;     % Wavelength (C-band SAR, ~5.6 cm)
R = 800e3;          % Slant range (800 km)
theta = 35 * pi/180; % Incidence angle

% Two baselines (like two sample rates)
B1 = 150;           % Baseline 1 (meters)
B2 = 170;           % Baseline 2 (meters) - should be coprime ratio with B1

fprintf('\nSystem Parameters:\n');
fprintf('  Wavelength: %.3f m\n', lambda);
fprintf('  Baseline 1: %d m\n', B1);
fprintf('  Baseline 2: %d m\n', B2);
fprintf('  GCD(B1,B2): %d m\n', gcd(B1, B2));
fprintf('  LCM(B1,B2): %d m\n', lcm(B1, B2));

%% ==================== Generate Wrapped Phases ====================
% Interferometric phase formula:
% phi = 4*pi*B*h / (lambda*R*sin(theta))

% Height-to-phase sensitivity
k1 = 4*pi*B1 / (lambda*R*sin(theta));  % rad/m
k2 = 4*pi*B2 / (lambda*R*sin(theta));  % rad/m

% True (unwrapped) phases
phi1_true = k1 * h_true;
phi2_true = k2 * h_true;

% Wrapped phases (what we actually measure)
phi1_wrap = wrapToPi(phi1_true);
phi2_wrap = wrapToPi(phi2_true);

% Add noise
noise_level = 0.3;  % radians
phi1_wrap = phi1_wrap + noise_level * randn(Ny, Nx);
phi2_wrap = phi2_wrap + noise_level * randn(Ny, Nx);
phi1_wrap = wrapToPi(phi1_wrap);
phi2_wrap = wrapToPi(phi2_wrap);

fprintf('\nPhase Statistics:\n');
fprintf('  True phase range: [%.1f, %.1f] rad\n', min(phi1_true(:)), max(phi1_true(:)));
fprintf('  This spans %.1f wavelengths!\n', (max(phi1_true(:))-min(phi1_true(:)))/(2*pi));

%% ==================== CRT Phase Unwrapping ====================
fprintf('\n--- CRT Unwrapping Algorithm ---\n');

% Ambiguity heights for each baseline
% One 2*pi cycle corresponds to:
h_amb_1 = lambda * R * sin(theta) / (2 * B1);  % meters
h_amb_2 = lambda * R * sin(theta) / (2 * B2);  % meters

fprintf('  Height ambiguity (B1): %.2f m\n', h_amb_1);
fprintf('  Height ambiguity (B2): %.2f m\n', h_amb_2);

% Synthetic ambiguity (from CRT)
% Find the height range that can be uniquely determined
h_amb_synthetic = lcm(round(h_amb_1*100), round(h_amb_2*100)) / 100;
fprintf('  Synthetic ambiguity: %.2f m\n', h_amb_synthetic);

% CRT Unwrapping: Search for height that matches both wrapped phases
h_recovered = zeros(Ny, Nx);

% For efficiency, process in vectorized form
% Search range based on synthetic ambiguity
h_search = linspace(0, max(h_true(:))*1.5, 500);

for i = 1:Ny
    for j = 1:Nx
        best_h = 0;
        min_error = inf;
        
        for h_test = h_search
            % Expected wrapped phases for this height
            phi1_expected = wrapToPi(k1 * h_test);
            phi2_expected = wrapToPi(k2 * h_test);
            
            % Phase difference (handle wrap-around)
            diff1 = abs(wrapToPi(phi1_wrap(i,j) - phi1_expected));
            diff2 = abs(wrapToPi(phi2_wrap(i,j) - phi2_expected));
            
            error = diff1 + diff2;
            
            if error < min_error
                min_error = error;
                best_h = h_test;
            end
        end
        
        h_recovered(i,j) = best_h;
    end
    
    if mod(i, 50) == 0
        fprintf('  Processing row %d/%d...\n', i, Ny);
    end
end

% Calculate errors
error_map = h_recovered - h_true;
rmse = sqrt(mean(error_map(:).^2));
mae = mean(abs(error_map(:)));

fprintf('\nRecovery Results:\n');
fprintf('  RMSE: %.3f m\n', rmse);
fprintf('  MAE:  %.3f m\n', mae);

%% ==================== Single Baseline Unwrapping (for comparison) ====================
% Simple gradient-based unwrapping (will fail for large phase gradients)
fprintf('\n--- Single Baseline Comparison ---\n');

% Naive unwrapping using cumulative sum (very basic)
h_single = phi1_wrap / k1;  % Direct conversion (stays wrapped)

% Try MATLAB's unwrap along each dimension
phi1_unwrap_x = unwrap(phi1_wrap, [], 2);
phi1_unwrap_xy = unwrap(phi1_unwrap_x, [], 1);
h_single_unwrap = phi1_unwrap_xy / k1;

error_single = h_single_unwrap - h_true;
rmse_single = sqrt(mean(error_single(:).^2));

fprintf('  Single baseline RMSE: %.3f m\n', rmse_single);
fprintf('  CRT dual baseline RMSE: %.3f m\n', rmse);
fprintf('  Improvement: %.1fx\n', rmse_single/rmse);

%% ==================== Visualization ====================
figure('Position', [50, 50, 1600, 900]);

% True elevation
subplot(2,4,1);
imagesc(h_true);
colorbar;
colormap(gca, 'parula');
title('True Elevation (m)');
axis image;

% Wrapped phase 1
subplot(2,4,2);
imagesc(phi1_wrap);
colorbar;
colormap(gca, 'hsv');
title(sprintf('Wrapped Phase B1=%dm', B1));
axis image;
caxis([-pi pi]);

% Wrapped phase 2
subplot(2,4,3);
imagesc(phi2_wrap);
colorbar;
colormap(gca, 'hsv');
title(sprintf('Wrapped Phase B2=%dm', B2));
axis image;
caxis([-pi pi]);

% CRT recovered elevation
subplot(2,4,4);
imagesc(h_recovered);
colorbar;
colormap(gca, 'parula');
title('CRT Recovered Elevation (m)');
axis image;

% Single baseline result
subplot(2,4,5);
imagesc(h_single_unwrap);
colorbar;
colormap(gca, 'parula');
title('Single Baseline Unwrap (m)');
axis image;

% CRT error map
subplot(2,4,6);
imagesc(error_map);
colorbar;
colormap(gca, 'jet');
title(sprintf('CRT Error (RMSE=%.2fm)', rmse));
axis image;
caxis([-5 5]);

% Single baseline error
subplot(2,4,7);
imagesc(error_single);
colorbar;
colormap(gca, 'jet');
title(sprintf('Single BL Error (RMSE=%.2fm)', rmse_single));
axis image;
caxis([-50 50]);

% Cross-section comparison
subplot(2,4,8);
mid_row = round(Ny/2);
plot(1:Nx, h_true(mid_row,:), 'k-', 'LineWidth', 2);
hold on;
plot(1:Nx, h_recovered(mid_row,:), 'b--', 'LineWidth', 1.5);
plot(1:Nx, h_single_unwrap(mid_row,:), 'r:', 'LineWidth', 1.5);
xlabel('X Position');
ylabel('Height (m)');
title('Cross-Section Comparison');
legend('True', 'CRT', 'Single BL', 'Location', 'best');
grid on;

sgtitle('CRT Phase Unwrapping: Multi-Baseline InSAR', 'FontSize', 14);

saveas(gcf, 'crt_2d_phase_unwrap.png');
fprintf('\nFigure saved: crt_2d_phase_unwrap.png\n');

%% ==================== Helper Function ====================
function wrapped = wrapToPi(phase)
    % Wrap phase to [-pi, pi]
    wrapped = mod(phase + pi, 2*pi) - pi;
end
