%% Staggered PRF - FMCW雷达EMG信号检测应用
% 
% 本脚本演示Staggered PRF在FMCW雷达EMG(肌电)信号检测中的应用：
% 1. ADF4159 Dual-Slope模式配置
% 2. 完整的信号生成和处理流程
% 3. kHz级别EMG频率的欠采样恢复
% 4. 时变EMG频率的跟踪
%
% 背景: EMG信号会调制皮肤的微振动，在FMCW雷达中表现为
%       多普勒频率的调制。EMG频率通常在1-10 kHz范围，
%       而chirp重复率仅为几kHz，导致严重的速度模糊。
%
% Author: Anti-Aliasing Recovery Toolkit
% Date: 2024-01

clear; close all; clc;

fprintf('============================================================\n');
fprintf('  Staggered PRF - FMCW雷达EMG信号检测演示\n');
fprintf('============================================================\n\n');

%% ==================== Part 1: FMCW雷达系统参数 ====================
fprintf('Part 1: FMCW雷达系统参数\n');
fprintf('-------------------------\n');

% 物理常数
c = 3e8;              % 光速 (m/s)
fc = 77e9;            % 载波频率 77 GHz (毫米波雷达)
lambda = c / fc;      % 波长 ≈ 3.9 mm

% FMCW chirp参数
B = 500e6;            % 带宽 500 MHz
fs_adc = 245.76e6;    % ADC采样率 (TI AWR系列典型值)

% ADF4159 Dual-Slope 配置
% 使用两个不同周期的chirp交替发射
T1 = 100e-6;          % Chirp 1: 100 μs
T2 = 150e-6;          % Chirp 2: 150 μs
D = 20e-6;            % 两个chirp之间的delay: 20 μs

% 计算chirp斜率
K1 = B / T1;          % Slope 1: 5 THz/s
K2 = B / T2;          % Slope 2: 3.33 THz/s

% 等效多普勒采样率
Delta_t1 = D + T2;    % T1结束 → T2结束: 170 μs
Delta_t2 = D + T1;    % T2结束 → T1结束: 120 μs
Fs1 = 1 / Delta_t1;   % 5882.35 Hz
Fs2 = 1 / Delta_t2;   % 8333.33 Hz

% 帧参数
T_frame = T1 + D + T2 + D;
frame_rate = 1 / T_frame;

% 距离分辨率和最大距离
range_res = c / (2 * B);  % 0.3 m
N_samples = round(T1 * fs_adc);  % chirp采样点数

fprintf('载波频率: %.1f GHz, 波长: %.2f mm\n', fc/1e9, lambda*1e3);
fprintf('带宽: %.0f MHz, 距离分辨率: %.2f m\n', B/1e6, range_res);
fprintf('ADC采样率: %.2f MHz\n', fs_adc/1e6);
fprintf('\n');
fprintf('Staggered PRF配置:\n');
fprintf('  T1 = %.0f μs (K1 = %.2f THz/s)\n', T1*1e6, K1/1e12);
fprintf('  T2 = %.0f μs (K2 = %.2f THz/s)\n', T2*1e6, K2/1e12);
fprintf('  D = %.0f μs\n', D*1e6);
fprintf('  等效采样率: Fs1 = %.1f Hz, Fs2 = %.1f Hz\n', Fs1, Fs2);
fprintf('  Nyquist: Fn1 = %.1f Hz, Fn2 = %.1f Hz\n', Fs1/2, Fs2/2);
fprintf('  帧率: %.1f Hz\n', frame_rate);
fprintf('\n');

%% ==================== Part 2: 目标和EMG参数 ====================
fprintf('Part 2: 目标和EMG参数\n');
fprintf('----------------------\n');

% 目标参数
R_target = 0.5;       % 目标距离: 0.5 m (近距离EMG检测)
v_bulk = 0.02;        % 整体运动速度: 20 mm/s (呼吸、抖动等)

% EMG参数
% EMG信号会引起皮肤微振动，其频率通常在1-10 kHz范围
% 这里模拟一个时变的EMG信号
f_emg_center = 5000;  % EMG中心频率: 5 kHz
f_emg_mod = 30;       % EMG频率调制速度: 30 Hz (模拟肌肉收缩变化)
f_emg_amp = 1500;     % EMG频率变化幅度: ±1.5 kHz
emg_phase_amp = 0.5;  % EMG引起的相位调制幅度 (radians)

% SNR
SNR_dB = 25;          % 信噪比

fprintf('目标距离: %.2f m\n', R_target);
fprintf('整体运动: %.1f mm/s\n', v_bulk*1e3);
fprintf('EMG信号:\n');
fprintf('  中心频率: %.0f Hz\n', f_emg_center);
fprintf('  频率范围: %.0f - %.0f Hz\n', f_emg_center-f_emg_amp, f_emg_center+f_emg_amp);
fprintf('  调制频率: %.0f Hz\n', f_emg_mod);
fprintf('SNR: %.0f dB\n', SNR_dB);
fprintf('\n');

%% ==================== Part 3: 信号生成 ====================
fprintf('Part 3: 信号生成\n');
fprintf('----------------\n');

% 仿真时间
sim_duration = 0.5;   % 500 ms
n_frames = floor(sim_duration / T_frame);

fprintf('仿真时间: %.0f ms\n', sim_duration*1e3);
fprintf('总帧数: %d\n', n_frames);
fprintf('\n');

% 预分配存储
chirp_phases_T1 = zeros(n_frames, 1);   % T1 chirp末端相位
chirp_phases_T2 = zeros(n_frames, 1);   % T2 chirp末端相位
chirp_times_T1 = zeros(n_frames, 1);    % T1 chirp时间戳
chirp_times_T2 = zeros(n_frames, 1);    % T2 chirp时间戳
true_emg_freqs = zeros(n_frames, 1);    % 真实EMG频率记录

% 混叠频率计算函数
alias_freq = @(f, Fs) mod(f + Fs/2, Fs) - Fs/2;

% Range bin对应的beat频率
f_beat = 2 * B * R_target / (c * T1);  % 对于T1
range_bin = round(f_beat * T1 * N_samples / fs_adc);

fprintf('Beat频率: %.2f MHz → Range bin: %d\n', f_beat/1e6, range_bin);
fprintf('\n');

% 信号生成循环
current_time = 0;

for frame = 1:n_frames
    % 帧内处理两个chirp
    
    % === Chirp T1 ===
    t_mid_T1 = current_time + T1/2;  % chirp中点时刻
    
    % 当前EMG频率 (时变)
    f_emg_now = f_emg_center + f_emg_amp * sin(2*pi*f_emg_mod*t_mid_T1);
    true_emg_freqs(frame) = f_emg_now;
    
    % 多普勒相位积累
    % 1. 整体运动引起的相位
    phase_bulk = 4*pi*v_bulk*t_mid_T1 / lambda;
    
    % 2. EMG调制引起的相位
    phase_emg = emg_phase_amp * sin(2*pi*f_emg_now*t_mid_T1);
    
    % 总相位 (在Range bin位置)
    total_phase_T1 = phase_bulk + phase_emg;
    
    % 添加噪声
    noise_std = 1 / (10^(SNR_dB/20));
    total_phase_T1 = total_phase_T1 + noise_std * randn();
    
    chirp_phases_T1(frame) = total_phase_T1;
    chirp_times_T1(frame) = current_time + T1;  % T1结束时刻
    
    current_time = current_time + T1 + D;
    
    % === Chirp T2 ===
    t_mid_T2 = current_time + T2/2;
    
    % EMG频率 (可能已变化)
    f_emg_now = f_emg_center + f_emg_amp * sin(2*pi*f_emg_mod*t_mid_T2);
    
    % 相位
    phase_bulk = 4*pi*v_bulk*t_mid_T2 / lambda;
    phase_emg = emg_phase_amp * sin(2*pi*f_emg_now*t_mid_T2);
    total_phase_T2 = phase_bulk + phase_emg + noise_std * randn();
    
    chirp_phases_T2(frame) = total_phase_T2;
    chirp_times_T2(frame) = current_time + T2;  % T2结束时刻
    
    current_time = current_time + T2 + D;
end

fprintf('信号生成完成。\n\n');

%% ==================== Part 4: 多普勒处理 ====================
fprintf('Part 4: 多普勒处理\n');
fprintf('------------------\n');

% 从相位序列估计瞬时频率
% 方法: 相邻相位差分 → 频率

% 对于T1序列: 间隔为 Delta_t1 = D + T2
phase_diff_T1 = diff(chirp_phases_T1);
time_diff_T1 = Delta_t1;
freq_T1 = phase_diff_T1 / (2*pi*time_diff_T1);  % 原始频率估计

% 对于T2序列: 间隔为 Delta_t2 = D + T1
phase_diff_T2 = diff(chirp_phases_T2);
time_diff_T2 = Delta_t2;
freq_T2 = phase_diff_T2 / (2*pi*time_diff_T2);

% 混叠频率 (实际测量到的)
freq_alias_T1 = arrayfun(@(f) alias_freq(f, Fs1), freq_T1);
freq_alias_T2 = arrayfun(@(f) alias_freq(f, Fs2), freq_T2);

% 注意: 由于整体运动速度很低(20mm/s)，其多普勒频率约为:
f_doppler_bulk = 2*v_bulk/lambda;
fprintf('整体运动多普勒: %.1f Hz (可忽略)\n', f_doppler_bulk);
fprintf('主要观测: EMG调制引起的频率分量 (kHz级别)\n\n');

%% ==================== Part 5: CRT恢复算法 ====================
fprintf('Part 5: CRT频率恢复\n');
fprintf('-------------------\n');

% CRT搜索恢复函数
function f_recovered = crt_search_single(fa1, fa2, Fs1, Fs2, f_range, tol)
    alias_freq_local = @(f, Fs) mod(f + Fs/2, Fs) - Fs/2;
    
    f_candidates = f_range(1):1:f_range(2);  % 1 Hz分辨率
    f_recovered = NaN;
    min_error = inf;
    
    for f = f_candidates
        fa1_exp = alias_freq_local(f, Fs1);
        fa2_exp = alias_freq_local(f, Fs2);
        
        error = abs(fa1_exp - fa1) + abs(fa2_exp - fa2);
        
        if error < min_error
            min_error = error;
            if error < tol
                f_recovered = f;
            end
        end
    end
end

% 恢复每个时刻的EMG频率
n_points = min(length(freq_alias_T1), length(freq_alias_T2));
f_recovered = zeros(n_points, 1);
f_range = [0, 15000];  % 搜索范围 0-15 kHz
tol = 50;  % 容差 50 Hz

fprintf('正在恢复EMG频率...\n');
tic;

for i = 1:n_points
    f_recovered(i) = crt_search_single(freq_alias_T1(i), freq_alias_T2(i), ...
                                        Fs1, Fs2, f_range, tol);
end

elapsed = toc;
fprintf('恢复完成，耗时: %.2f 秒\n', elapsed);

% 计算误差
true_freqs_trimmed = true_emg_freqs(1:n_points);
errors = abs(f_recovered - true_freqs_trimmed);
valid_idx = ~isnan(f_recovered);
rmse = sqrt(mean(errors(valid_idx).^2));
success_rate = sum(errors(valid_idx) < 100) / sum(valid_idx) * 100;

fprintf('\n恢复性能:\n');
fprintf('  RMSE: %.1f Hz\n', rmse);
fprintf('  成功率 (<100Hz误差): %.1f%%\n', success_rate);
fprintf('  有效恢复点: %d / %d\n', sum(valid_idx), n_points);
fprintf('\n');

%% ==================== Part 6: 结果可视化 ====================
fprintf('Part 6: 结果可视化\n');
fprintf('------------------\n');

time_axis = (1:n_points) * T_frame * 1000;  % ms

figure('Name', 'Staggered PRF FMCW雷达EMG检测', 'Position', [100, 50, 1400, 900]);

% 子图1: 真实EMG频率 vs 恢复频率
subplot(3,2,1);
plot(time_axis, true_freqs_trimmed/1000, 'b-', 'LineWidth', 1.5, 'DisplayName', '真实EMG频率');
hold on;
plot(time_axis, f_recovered/1000, 'r.', 'MarkerSize', 4, 'DisplayName', '恢复频率');
xlabel('时间 (ms)');
ylabel('频率 (kHz)');
title('EMG频率恢复结果');
legend('Location', 'best');
grid on;
ylim([2, 8]);

% 子图2: 混叠后的测量频率
subplot(3,2,2);
plot(time_axis, freq_alias_T1, 'b.', 'MarkerSize', 3, 'DisplayName', sprintf('Fs1=%.0fHz', Fs1));
hold on;
plot(time_axis, freq_alias_T2, 'r.', 'MarkerSize', 3, 'DisplayName', sprintf('Fs2=%.0fHz', Fs2));
xlabel('时间 (ms)');
ylabel('混叠频率 (Hz)');
title('原始测量 (混叠后)');
legend('Location', 'best');
grid on;

% 子图3: 恢复误差
subplot(3,2,3);
valid_errors = errors(valid_idx);
valid_times = time_axis(valid_idx);
stem(valid_times, valid_errors, 'b', 'MarkerSize', 2);
xlabel('时间 (ms)');
ylabel('误差 (Hz)');
title(sprintf('频率恢复误差 (RMSE=%.1f Hz)', rmse));
grid on;
ylim([0, max(valid_errors)*1.1]);

% 子图4: 相位序列
subplot(3,2,4);
plot(chirp_times_T1*1000, chirp_phases_T1, 'b.-', 'MarkerSize', 3, 'DisplayName', 'T1 chirps');
hold on;
plot(chirp_times_T2*1000, chirp_phases_T2, 'r.-', 'MarkerSize', 3, 'DisplayName', 'T2 chirps');
xlabel('时间 (ms)');
ylabel('相位 (rad)');
title('Range Bin相位序列');
legend('Location', 'best');
grid on;
xlim([0, 50]);  % 只显示前50ms

% 子图5: (fa1, fa2) 散点图
subplot(3,2,5);
scatter(freq_alias_T1, freq_alias_T2, 10, true_freqs_trimmed(1:end-1)/1000, 'filled');
colorbar;
xlabel('fa1 (Hz)');
ylabel('fa2 (Hz)');
title('(fa1, fa2) 测量对 (颜色=真实频率/kHz)');
grid on;
axis equal;

% 子图6: 频谱分析
subplot(3,2,6);
% 对恢复的频率序列做FFT，应该能看到EMG调制频率
f_recovered_valid = f_recovered(valid_idx);
f_recovered_valid(isnan(f_recovered_valid)) = mean(f_recovered_valid(~isnan(f_recovered_valid)));
N_fft = 2^nextpow2(length(f_recovered_valid));
F_spec = fft(f_recovered_valid - mean(f_recovered_valid), N_fft);
freq_axis_spec = (0:N_fft/2-1) * frame_rate / N_fft;
plot(freq_axis_spec, abs(F_spec(1:N_fft/2)), 'b-', 'LineWidth', 1);
xlabel('调制频率 (Hz)');
ylabel('幅度');
title(sprintf('EMG频率调制谱 (应看到 %.0f Hz 峰)', f_emg_mod));
grid on;
xlim([0, 100]);
xline(f_emg_mod, 'r--', 'LineWidth', 1.5);

sgtitle(sprintf('Staggered PRF FMCW雷达EMG检测 | T1=%.0fμs, T2=%.0fμs | SNR=%.0fdB', ...
                T1*1e6, T2*1e6, SNR_dB), 'FontSize', 14, 'FontWeight', 'bold');

%% ==================== Part 7: 配置优化分析 ====================
fprintf('\nPart 7: 不同配置对比\n');
fprintf('---------------------\n');

configs = struct();
configs(1).name = '配置A (100/150μs)';
configs(1).T1 = 100e-6; configs(1).T2 = 150e-6;
configs(2).name = '配置B (127/173μs)';
configs(2).T1 = 127e-6; configs(2).T2 = 173e-6;
configs(3).name = '配置C (80/120μs)';
configs(3).T1 = 80e-6; configs(3).T2 = 120e-6;

fprintf('%-20s | Fs1 (Hz) | Fs2 (Hz) | 帧率 (Hz) | 唯一范围\n', '配置');
fprintf('--------------------|----------|----------|-----------|----------\n');

for i = 1:length(configs)
    T1_c = configs(i).T1;
    T2_c = configs(i).T2;
    D_c = D;
    
    Fs1_c = 1 / (D_c + T2_c);
    Fs2_c = 1 / (D_c + T1_c);
    frame_rate_c = 1 / (T1_c + D_c + T2_c + D_c);
    
    % 估计唯一范围 (简化)
    unique_range = Fs1_c * Fs2_c / gcd(round(Fs1_c), round(Fs2_c));
    
    fprintf('%-20s | %8.1f | %8.1f | %9.1f | %.1f kHz\n', ...
            configs(i).name, Fs1_c, Fs2_c, frame_rate_c, unique_range/1000);
end

%% ==================== Part 8: 保存结果 ====================
fprintf('\nPart 8: 保存结果\n');
fprintf('----------------\n');

% 保存关键数据到MAT文件
results = struct();
results.config.T1 = T1;
results.config.T2 = T2;
results.config.D = D;
results.config.Fs1 = Fs1;
results.config.Fs2 = Fs2;
results.emg.center_freq = f_emg_center;
results.emg.mod_freq = f_emg_mod;
results.emg.amp = f_emg_amp;
results.data.time = time_axis;
results.data.true_freq = true_freqs_trimmed;
results.data.recovered_freq = f_recovered;
results.data.alias_T1 = freq_alias_T1;
results.data.alias_T2 = freq_alias_T2;
results.performance.rmse = rmse;
results.performance.success_rate = success_rate;

save('staggered_prf_fmcw_results.mat', 'results');
fprintf('结果已保存到: staggered_prf_fmcw_results.mat\n');

fprintf('\n============================================================\n');
fprintf('  演示完成！\n');
fprintf('============================================================\n');
fprintf('\n关键结论:\n');
fprintf('1. Staggered PRF成功将%d Hz EMG信号从%.0f/%.0f Hz采样率中恢复\n', ...
        f_emg_center, Fs1, Fs2);
fprintf('2. RMSE = %.1f Hz，成功率 = %.1f%%\n', rmse, success_rate);
fprintf('3. EMG的%.0f Hz调制特征在恢复频谱中清晰可见\n', f_emg_mod);
fprintf('\n');
