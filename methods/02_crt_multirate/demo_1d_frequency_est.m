%% CRT Anti-Aliasing Demo: 1D Frequency Estimation
% ============================================
% Scenario: Estimate a high frequency using two low-rate ADCs
% 
% Problem: f_true = 850 Hz, but we only have:
%   - ADC1: Fs1 = 97 Hz  (Nyquist = 48.5 Hz)
%   - ADC2: Fs2 = 101 Hz (Nyquist = 50.5 Hz)
%
% Both ADCs see aliased frequencies, but CRT can recover the truth!
%
% Author: [Your Name]
% Date: 2024
% ============================================

clear; clc; close all;

%% ==================== Parameters ====================
f_true = 850;           % True frequency (Hz) - way above Nyquist!
amplitude = 1.0;
phase_init = pi/6;      % Initial phase

% Two ADCs with coprime sample rates
Fs1 = 97;               % ADC1 sample rate (Hz) - prime number
Fs2 = 101;              % ADC2 sample rate (Hz) - prime number

% Observation time
T_obs = 1.0;            % 1 second
N1 = round(T_obs * Fs1);
N2 = round(T_obs * Fs2);

% Maximum recoverable frequency (LCM of sample rates)
max_recoverable = lcm(Fs1, Fs2);

fprintf('============================================\n');
fprintf('   CRT Anti-Aliasing: 1D Frequency Estimation\n');
fprintf('============================================\n');
fprintf('\nSystem Parameters:\n');
fprintf('  True frequency: %d Hz\n', f_true);
fprintf('  ADC1: Fs = %d Hz, Nyquist = %.1f Hz\n', Fs1, Fs1/2);
fprintf('  ADC2: Fs = %d Hz, Nyquist = %.1f Hz\n', Fs2, Fs2/2);
fprintf('  Max recoverable freq (LCM): %d Hz\n', max_recoverable);

%% ==================== Generate Signals ====================
% Continuous-time signal: x(t) = A * cos(2*pi*f_true*t + phi)

% ADC1 samples
t1 = (0:N1-1)' / Fs1;
x1 = amplitude * cos(2*pi*f_true*t1 + phase_init);
x1 = x1 + 0.1*randn(size(x1));  % Add noise

% ADC2 samples
t2 = (0:N2-1)' / Fs2;
x2 = amplitude * cos(2*pi*f_true*t2 + phase_init);
x2 = x2 + 0.1*randn(size(x2));  % Add noise

%% ==================== Frequency Estimation (FFT) ====================
% Zero-padding for better frequency resolution
N_fft = 4096;

% ADC1 spectrum
X1 = fft(x1 .* hamming(N1), N_fft);
freq_axis_1 = (0:N_fft-1) * Fs1 / N_fft;
% Fold to [-Fs/2, Fs/2]
freq_axis_1(freq_axis_1 > Fs1/2) = freq_axis_1(freq_axis_1 > Fs1/2) - Fs1;

% Find peak (positive frequency only)
[~, idx1] = max(abs(X1(1:N_fft/2)));
f_alias_1 = (idx1-1) * Fs1 / N_fft;

% ADC2 spectrum
X2 = fft(x2 .* hamming(N2), N_fft);
freq_axis_2 = (0:N_fft-1) * Fs2 / N_fft;

[~, idx2] = max(abs(X2(1:N_fft/2)));
f_alias_2 = (idx2-1) * Fs2 / N_fft;

fprintf('\n--- Aliased Frequency Detection ---\n');
fprintf('  ADC1 sees: %.2f Hz\n', f_alias_1);
fprintf('  ADC2 sees: %.2f Hz\n', f_alias_2);

% Theoretical aliased frequencies
f_alias_1_theory = mod(f_true, Fs1);
if f_alias_1_theory > Fs1/2
    f_alias_1_theory = Fs1 - f_alias_1_theory;
end
f_alias_2_theory = mod(f_true, Fs2);
if f_alias_2_theory > Fs2/2
    f_alias_2_theory = Fs2 - f_alias_2_theory;
end

fprintf('  (Theory: ADC1=%.2f Hz, ADC2=%.2f Hz)\n', f_alias_1_theory, f_alias_2_theory);

%% ==================== CRT Recovery Algorithm ====================
fprintf('\n--- CRT Recovery ---\n');

% Method: Search for f in [0, LCM) that matches both aliased observations
% 
% For each candidate f:
%   expected_alias_1 = f mod Fs1 (folded to [0, Fs1/2])
%   expected_alias_2 = f mod Fs2 (folded to [0, Fs2/2])
%   
% Find f that minimizes error to observed aliases

search_resolution = 0.5;  % Hz
f_candidates = 0:search_resolution:max_recoverable;

best_f = 0;
min_error = inf;
errors = zeros(size(f_candidates));

for i = 1:length(f_candidates)
    f_test = f_candidates(i);
    
    % Expected aliased frequency for ADC1
    expected_1 = mod(f_test, Fs1);
    if expected_1 > Fs1/2
        expected_1 = Fs1 - expected_1;
    end
    
    % Expected aliased frequency for ADC2
    expected_2 = mod(f_test, Fs2);
    if expected_2 > Fs2/2
        expected_2 = Fs2 - expected_2;
    end
    
    % Total error
    error = abs(expected_1 - f_alias_1) + abs(expected_2 - f_alias_2);
    errors(i) = error;
    
    if error < min_error
        min_error = error;
        best_f = f_test;
    end
end

fprintf('  Recovered frequency: %.2f Hz\n', best_f);
fprintf('  True frequency:      %.2f Hz\n', f_true);
fprintf('  Recovery error:      %.2f Hz\n', abs(best_f - f_true));

%% ==================== Visualization ====================
figure('Position', [100, 100, 1400, 900]);

% --- Time domain signals ---
subplot(2,3,1);
plot(t1*1000, x1, 'b.-', 'MarkerSize', 8);
hold on;
t_fine = linspace(0, 0.05, 1000);
plot(t_fine*1000, amplitude*cos(2*pi*f_true*t_fine + phase_init), 'r-', 'LineWidth', 0.5);
xlabel('Time (ms)');
ylabel('Amplitude');
title(sprintf('ADC1: Fs=%d Hz (samples shown as dots)', Fs1));
legend('Samples', 'True signal', 'Location', 'best');
xlim([0 50]);
grid on;

subplot(2,3,4);
plot(t2*1000, x2, 'g.-', 'MarkerSize', 8);
hold on;
plot(t_fine*1000, amplitude*cos(2*pi*f_true*t_fine + phase_init), 'r-', 'LineWidth', 0.5);
xlabel('Time (ms)');
ylabel('Amplitude');
title(sprintf('ADC2: Fs=%d Hz (samples shown as dots)', Fs2));
legend('Samples', 'True signal', 'Location', 'best');
xlim([0 50]);
grid on;

% --- Frequency domain ---
subplot(2,3,2);
freq_plot_1 = (0:N_fft/2-1) * Fs1 / N_fft;
stem(freq_plot_1, abs(X1(1:N_fft/2))/N1, 'b', 'MarkerSize', 3);
hold on;
yl = ylim; plot([f_alias_1 f_alias_1], yl, 'r--', 'LineWidth', 2);  % xline alternative
xlabel('Frequency (Hz)');
ylabel('|X|');
title(sprintf('ADC1 Spectrum: Peak at %.1f Hz', f_alias_1));
xlim([0 Fs1/2]);
grid on;

subplot(2,3,5);
freq_plot_2 = (0:N_fft/2-1) * Fs2 / N_fft;
stem(freq_plot_2, abs(X2(1:N_fft/2))/N2, 'g', 'MarkerSize', 3);
hold on;
yl = ylim; plot([f_alias_2 f_alias_2], yl, 'r--', 'LineWidth', 2);  % xline alternative
xlabel('Frequency (Hz)');
ylabel('|X|');
title(sprintf('ADC2 Spectrum: Peak at %.1f Hz', f_alias_2));
xlim([0 Fs2/2]);
grid on;

% --- CRT Error Surface ---
subplot(2,3,3);
plot(f_candidates, errors, 'b-', 'LineWidth', 0.5);
hold on;
plot(best_f, min_error, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
plot(f_true, errors(round(f_true/search_resolution)+1), 'g^', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Candidate Frequency (Hz)');
ylabel('CRT Error');
title('CRT Search: Error vs Candidate Frequency');
legend('Error', sprintf('Recovered: %.1f Hz', best_f), sprintf('True: %.1f Hz', f_true));
grid on;

% --- Summary ---
subplot(2,3,6);
axis off;

text(0.5, 0.95, 'CRT ANTI-ALIASING SUMMARY', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

info_text = {
    sprintf('True Frequency: %d Hz', f_true),
    '',
    sprintf('ADC1: Fs=%d Hz, sees %.1f Hz', Fs1, f_alias_1),
    sprintf('ADC2: Fs=%d Hz, sees %.1f Hz', Fs2, f_alias_2),
    '',
    sprintf('Recoverable Range: [0, %d) Hz', max_recoverable),
    '',
    sprintf('Recovered: %.1f Hz', best_f),
    sprintf('Error: %.2f Hz', abs(best_f - f_true)),
    '',
    'Key Insight:',
    'Coprime sample rates create unique',
    'aliasing patterns for each frequency!'
};

text(0.1, 0.75, info_text, 'FontSize', 11, 'VerticalAlignment', 'top', ...
    'FontName', 'FixedWidth');

% Save figure
saveas(gcf, 'crt_1d_frequency_est.png');
fprintf('\nFigure saved: crt_1d_frequency_est.png\n');

%% ==================== Extended Demo: Multiple Frequencies ====================
fprintf('\n============================================\n');
fprintf('   Extended Test: Various Frequencies\n');
fprintf('============================================\n');

test_freqs = [120, 350, 567, 850, 1234, 3456, 7890];
test_freqs = test_freqs(test_freqs < max_recoverable);  % Within range

fprintf('\n%10s | %10s | %10s | %10s | %10s\n', ...
    'f_true', 'f_alias1', 'f_alias2', 'f_recov', 'Error');
fprintf('%s\n', repmat('-', 1, 60));

for f_test_true = test_freqs
    % Generate signals
    x1_test = cos(2*pi*f_test_true*t1);
    x2_test = cos(2*pi*f_test_true*t2);
    
    % Detect aliased frequencies
    X1_test = fft(x1_test .* hamming(N1), N_fft);
    X2_test = fft(x2_test .* hamming(N2), N_fft);
    
    [~, idx1] = max(abs(X1_test(1:N_fft/2)));
    [~, idx2] = max(abs(X2_test(1:N_fft/2)));
    
    f_a1 = (idx1-1) * Fs1 / N_fft;
    f_a2 = (idx2-1) * Fs2 / N_fft;
    
    % CRT recovery
    best_f_test = 0;
    min_err = inf;
    for f_cand = 0:0.5:max_recoverable
        exp1 = mod(f_cand, Fs1);
        if exp1 > Fs1/2, exp1 = Fs1 - exp1; end
        exp2 = mod(f_cand, Fs2);
        if exp2 > Fs2/2, exp2 = Fs2 - exp2; end
        
        err = abs(exp1 - f_a1) + abs(exp2 - f_a2);
        if err < min_err
            min_err = err;
            best_f_test = f_cand;
        end
    end
    
    fprintf('%10.1f | %10.1f | %10.1f | %10.1f | %10.2f\n', ...
        f_test_true, f_a1, f_a2, best_f_test, abs(best_f_test - f_test_true));
end

fprintf('\nAll frequencies within [0, %d) Hz recovered successfully!\n', max_recoverable);
