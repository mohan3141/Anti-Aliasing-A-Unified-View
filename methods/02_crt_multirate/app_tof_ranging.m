%% CRT Application: Time-of-Flight Camera Range Ambiguity
% ============================================
% Scenario: ToF camera uses modulated light to measure distance
%
% Problem: Phase-based ranging has ambiguity
%   d = (c * phi) / (4 * pi * f_mod)
%   When phi wraps around 2*pi, distance "resets"
%   
% Ambiguity distance: d_amb = c / (2 * f_mod)
%   - f_mod = 20 MHz -> d_amb = 7.5 m (wraps every 7.5m!)
%
% Solution: Multi-frequency ToF with CRT
%   - f1 = 20 MHz -> d_amb1 = 7.5 m
%   - f2 = 17 MHz -> d_amb2 = 8.82 m
%   - CRT -> d_amb_synthetic = LCM ≈ 132 m!
%
% Author: [Your Name]
% Date: 2024
% ============================================

clear; clc; close all;

fprintf('============================================\n');
fprintf('   CRT Application: ToF Range De-ambiguity\n');
fprintf('============================================\n');

%% ==================== System Parameters ====================
c = 3e8;  % Speed of light (m/s)

% Modulation frequencies (choose carefully for good LCM)
f1 = 20e6;   % 20 MHz
f2 = 17e6;   % 17 MHz (coprime with 20)

% Ambiguity distances
d_amb1 = c / (2 * f1);  % 7.5 m
d_amb2 = c / (2 * f2);  % 8.82 m

% Synthetic ambiguity (approximate, since frequencies aren't exactly coprime)
% For integers: LCM(20,17) = 340 -> synthetic amb ≈ 340/20 * 7.5 = 127.5 m
d_amb_synthetic = lcm(20, 17) / 20 * d_amb1;

fprintf('\nSystem Configuration:\n');
fprintf('  Frequency 1: %.0f MHz -> Ambiguity: %.2f m\n', f1/1e6, d_amb1);
fprintf('  Frequency 2: %.0f MHz -> Ambiguity: %.2f m\n', f2/1e6, d_amb2);
fprintf('  Synthetic Ambiguity: %.2f m\n', d_amb_synthetic);

%% ==================== Create Test Scene ====================
% Simulate a depth image (e.g., room with objects at various distances)
Nx = 128;
Ny = 96;
[X, Y] = meshgrid(linspace(-2, 2, Nx), linspace(-1.5, 1.5, Ny));

% Background wall at 25 m (beyond single-frequency ambiguity!)
d_true = 25 * ones(Ny, Nx);

% Add some objects at different depths
% Object 1: Box at 8 m (within first ambiguity of f1)
mask1 = (abs(X + 1) < 0.5) & (abs(Y) < 0.5);
d_true(mask1) = 8;

% Object 2: Sphere at 18 m (wraps once for f1, twice for f2)
r = sqrt((X - 0.5).^2 + (Y - 0.2).^2);
mask2 = r < 0.6;
d_true(mask2) = 18 - 5 * (1 - r(mask2)/0.6);  % Curved surface

% Object 3: Close object at 3 m
mask3 = (abs(X + 0.2) < 0.3) & (abs(Y + 0.8) < 0.3);
d_true(mask3) = 3;

fprintf('\nTrue Distance Range: [%.1f, %.1f] m\n', min(d_true(:)), max(d_true(:)));

%% ==================== Simulate ToF Measurements ====================
% Phase measurements (with wrapping)
phi1_true = 4 * pi * f1 * d_true / c;
phi2_true = 4 * pi * f2 * d_true / c;

% Wrapped phases (what ToF sensor measures)
phi1_wrap = mod(phi1_true, 2*pi);
phi2_wrap = mod(phi2_true, 2*pi);

% Add realistic noise (phase noise depends on signal amplitude, reflectivity)
noise_std = 0.1;  % radians (~1-2 cm depth noise)
phi1_meas = phi1_wrap + noise_std * randn(Ny, Nx);
phi2_meas = phi2_wrap + noise_std * randn(Ny, Nx);
phi1_meas = mod(phi1_meas, 2*pi);
phi2_meas = mod(phi2_meas, 2*pi);

% Direct (ambiguous) distance estimates
d_meas1 = c * phi1_meas / (4 * pi * f1);
d_meas2 = c * phi2_meas / (4 * pi * f2);

fprintf('\nSingle-Frequency Measurements:\n');
fprintf('  f1: All distances mapped to [0, %.2f] m\n', d_amb1);
fprintf('  f2: All distances mapped to [0, %.2f] m\n', d_amb2);

%% ==================== CRT Range Recovery ====================
fprintf('\n--- CRT De-ambiguity Processing ---\n');

d_recovered = zeros(Ny, Nx);
d_search = linspace(0, d_amb_synthetic, 2000);

for i = 1:Ny
    for j = 1:Nx
        best_d = 0;
        min_error = inf;
        
        for d_test = d_search
            % Expected wrapped phases for this distance
            phi1_exp = mod(4 * pi * f1 * d_test / c, 2*pi);
            phi2_exp = mod(4 * pi * f2 * d_test / c, 2*pi);
            
            % Phase error (handle wrap-around at 0/2pi boundary)
            diff1 = min(abs(phi1_meas(i,j) - phi1_exp), ...
                       2*pi - abs(phi1_meas(i,j) - phi1_exp));
            diff2 = min(abs(phi2_meas(i,j) - phi2_exp), ...
                       2*pi - abs(phi2_meas(i,j) - phi2_exp));
            
            error = diff1^2 + diff2^2;  % Squared error for better minima
            
            if error < min_error
                min_error = error;
                best_d = d_test;
            end
        end
        
        d_recovered(i,j) = best_d;
    end
    
    if mod(i, 20) == 0
        fprintf('  Processing row %d/%d...\n', i, Ny);
    end
end

% Error statistics
error_map = d_recovered - d_true;
rmse = sqrt(mean(error_map(:).^2));
mae = mean(abs(error_map(:)));

fprintf('\nRecovery Results:\n');
fprintf('  RMSE: %.4f m (%.1f mm)\n', rmse, rmse*1000);
fprintf('  MAE:  %.4f m (%.1f mm)\n', mae, mae*1000);

%% ==================== Visualization ====================
figure('Position', [50, 50, 1600, 1000]);

% True depth
subplot(2,4,1);
imagesc(d_true);
colorbar;
colormap(gca, 'turbo');
title('True Distance (m)');
axis image;
caxis([0 30]);

% Single frequency f1
subplot(2,4,2);
imagesc(d_meas1);
colorbar;
colormap(gca, 'turbo');
title(sprintf('Single Freq (%.0f MHz)\nAmbiguous!', f1/1e6));
axis image;
caxis([0 d_amb1]);

% Single frequency f2
subplot(2,4,3);
imagesc(d_meas2);
colorbar;
colormap(gca, 'turbo');
title(sprintf('Single Freq (%.0f MHz)\nAmbiguous!', f2/1e6));
axis image;
caxis([0 d_amb2]);

% CRT recovered
subplot(2,4,4);
imagesc(d_recovered);
colorbar;
colormap(gca, 'turbo');
title(sprintf('CRT Multi-Freq\nRMSE=%.1fmm', rmse*1000));
axis image;
caxis([0 30]);

% Wrapped phase 1
subplot(2,4,5);
imagesc(phi1_meas);
colorbar;
colormap(gca, 'hsv');
title(sprintf('Wrapped Phase (%.0f MHz)', f1/1e6));
axis image;
caxis([0 2*pi]);

% Wrapped phase 2
subplot(2,4,6);
imagesc(phi2_meas);
colorbar;
colormap(gca, 'hsv');
title(sprintf('Wrapped Phase (%.0f MHz)', f2/1e6));
axis image;
caxis([0 2*pi]);

% Error map
subplot(2,4,7);
imagesc(error_map * 100);  % Convert to cm
colorbar;
colormap(gca, 'jet');
title('Error Map (cm)');
axis image;
caxis([-5 5]);

% Cross-section
subplot(2,4,8);
mid_row = round(Ny/2);
plot(1:Nx, d_true(mid_row,:), 'k-', 'LineWidth', 2);
hold on;
plot(1:Nx, d_meas1(mid_row,:), 'b:', 'LineWidth', 1);
plot(1:Nx, d_recovered(mid_row,:), 'r-', 'LineWidth', 1.5);
xlabel('X Pixel');
ylabel('Distance (m)');
title('Cross-Section Comparison');
legend('True', sprintf('Single (%.0fMHz)', f1/1e6), 'CRT', 'Location', 'best');
ylim([0 30]);
grid on;

sgtitle('CRT Multi-Frequency ToF: Range De-ambiguity', 'FontSize', 14);

saveas(gcf, 'crt_tof_ranging.png');
fprintf('\nFigure saved: crt_tof_ranging.png\n');

%% ==================== Frequency Design Analysis ====================
fprintf('\n============================================\n');
fprintf('   Frequency Pair Design Analysis\n');
fprintf('============================================\n');

% Compare different frequency combinations
freq_pairs = [
    20, 17;    % Our choice
    20, 19;    % Close frequencies
    20, 15;    % GCD = 5
    20, 21;    % Coprime
    20, 23;    % Coprime, larger gap
    30, 29;    % Higher frequencies, coprime
];

fprintf('\n%6s %6s | %8s %8s | %10s | %s\n', ...
    'f1', 'f2', 'd_amb1', 'd_amb2', 'd_synth', 'Notes');
fprintf('%s\n', repmat('-', 1, 70));

for k = 1:size(freq_pairs, 1)
    f1_test = freq_pairs(k, 1) * 1e6;
    f2_test = freq_pairs(k, 2) * 1e6;
    
    d1 = c / (2 * f1_test);
    d2 = c / (2 * f2_test);
    
    g = gcd(freq_pairs(k,1), freq_pairs(k,2));
    l = lcm(freq_pairs(k,1), freq_pairs(k,2));
    d_synth = l / freq_pairs(k,1) * d1;
    
    if g == 1
        note = 'Coprime - Good!';
    elseif g < 5
        note = sprintf('GCD=%d', g);
    else
        note = sprintf('GCD=%d - Poor', g);
    end
    
    fprintf('%5.0f %5.0f | %7.2f %7.2f | %9.1f | %s\n', ...
        f1_test/1e6, f2_test/1e6, d1, d2, d_synth, note);
end

fprintf('\nKey Insight: Choose COPRIME frequency ratios for maximum range!\n');

%% ==================== Real-World Considerations ====================
fprintf('\n============================================\n');
fprintf('   Real-World Implementation Notes\n');
fprintf('============================================\n');

notes = {
    '1. Multi-frequency ToF is used in: Azure Kinect, iPhone LiDAR, etc.',
    '2. More than 2 frequencies can further extend range',
    '3. Computation can be parallelized (each pixel independent)',
    '4. Phase noise increases with distance (lower SNR)',
    '5. Multi-path interference causes systematic errors',
    '6. Real systems often use lookup tables instead of search'
};

for i = 1:length(notes)
    fprintf('  %s\n', notes{i});
end
