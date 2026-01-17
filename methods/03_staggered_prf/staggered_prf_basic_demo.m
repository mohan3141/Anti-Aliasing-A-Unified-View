%% Staggered PRF 基础原理演示
% 
% 本脚本演示Staggered PRF（交错脉冲重复频率）方法的核心原理：
% 1. 混叠现象的可视化
% 2. 双采样率下的混叠模式差异
% 3. 利用CRT搜索恢复真实频率
% 4. 性能分析与可视化
%
% Author: Anti-Aliasing Recovery Toolkit
% Date: 2024-01

clear; close all; clc;

fprintf('========================================\n');
fprintf('  Staggered PRF 基础原理演示\n');
fprintf('========================================\n\n');

%% ==================== Part 1: 系统参数配置 ====================
fprintf('Part 1: 系统参数配置\n');
fprintf('--------------------\n');

% Staggered PRF 配置
T1 = 100e-6;      % Chirp 1 周期: 100 μs
T2 = 150e-6;      % Chirp 2 周期: 150 μs
D = 20e-6;        % 两个chirp之间的delay: 20 μs

% 计算等效采样率
% 时序: [T1][D][T2][D][T1][D][T2]...
% Δt1 = T1结束到T2结束的时间 = D + T2
% Δt2 = T2结束到T1结束的时间 = D + T1
Delta_t1 = D + T2;   % 170 μs
Delta_t2 = D + T1;   % 120 μs

Fs1 = 1 / Delta_t1;  % 等效采样率1: 5882.35 Hz
Fs2 = 1 / Delta_t2;  % 等效采样率2: 8333.33 Hz

% 帧率
T_frame = T1 + D + T2 + D;  % 一个完整帧
frame_rate = 1 / T_frame;

fprintf('Chirp周期: T1 = %.0f μs, T2 = %.0f μs\n', T1*1e6, T2*1e6);
fprintf('Delay: D = %.0f μs\n', D*1e6);
fprintf('等效采样间隔: Δt1 = %.0f μs, Δt2 = %.0f μs\n', Delta_t1*1e6, Delta_t2*1e6);
fprintf('等效采样率: Fs1 = %.2f Hz, Fs2 = %.2f Hz\n', Fs1, Fs2);
fprintf('Nyquist频率: Fn1 = %.2f Hz, Fn2 = %.2f Hz\n', Fs1/2, Fs2/2);
fprintf('帧周期: %.0f μs, 帧率: %.2f Hz\n', T_frame*1e6, frame_rate);
fprintf('\n');

%% ==================== Part 2: 混叠频率计算 ====================
fprintf('Part 2: 混叠频率可视化\n');
fprintf('----------------------\n');

% 混叠频率计算函数
alias_freq = @(f, Fs) mod(f + Fs/2, Fs) - Fs/2;

% 测试频率范围
f_test = 0:10:15000;  % 0 - 15 kHz

% 计算混叠后的频率
f_alias1 = arrayfun(@(f) alias_freq(f, Fs1), f_test);
f_alias2 = arrayfun(@(f) alias_freq(f, Fs2), f_test);

% 绘图
figure('Name', 'Staggered PRF - 混叠映射', 'Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(f_test/1000, f_alias1, 'b-', 'LineWidth', 1.5);
hold on;
yline(0, 'k--', 'LineWidth', 0.5);
xline(Fs1/2/1000, 'r--', 'Fn1', 'LineWidth', 1);
xlabel('真实频率 (kHz)');
ylabel('混叠频率 (Hz)');
title(sprintf('Fs1 = %.1f Hz 的混叠映射', Fs1));
grid on;
ylim([-Fs1/2*1.1, Fs1/2*1.1]);
legend('f_{alias}', '', 'Nyquist', 'Location', 'best');

subplot(2,2,2);
plot(f_test/1000, f_alias2, 'r-', 'LineWidth', 1.5);
hold on;
yline(0, 'k--', 'LineWidth', 0.5);
xline(Fs2/2/1000, 'b--', 'Fn2', 'LineWidth', 1);
xlabel('真实频率 (kHz)');
ylabel('混叠频率 (Hz)');
title(sprintf('Fs2 = %.1f Hz 的混叠映射', Fs2));
grid on;
ylim([-Fs2/2*1.1, Fs2/2*1.1]);
legend('f_{alias}', '', 'Nyquist', 'Location', 'best');

subplot(2,2,[3,4]);
plot(f_test/1000, f_alias1, 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Fs1=%.0fHz', Fs1));
hold on;
plot(f_test/1000, f_alias2, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('Fs2=%.0fHz', Fs2));
xlabel('真实频率 (kHz)');
ylabel('混叠频率 (Hz)');
title('双采样率混叠映射对比');
grid on;
legend('Location', 'best');

% 标注几个关键频率点
key_freqs = [3000, 5000, 7000, 10000];
for f = key_freqs
    fa1 = alias_freq(f, Fs1);
    fa2 = alias_freq(f, Fs2);
    plot(f/1000, fa1, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
    plot(f/1000, fa2, 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    text(f/1000+0.2, max(fa1,fa2)+200, sprintf('%dHz', f), 'FontSize', 9);
end

fprintf('关键频率的混叠结果:\n');
fprintf('  f(Hz)   |  fa1(Fs1)  |  fa2(Fs2)\n');
fprintf('  --------|------------|----------\n');
for f = key_freqs
    fa1 = alias_freq(f, Fs1);
    fa2 = alias_freq(f, Fs2);
    fprintf('  %5d   |  %+7.1f   |  %+7.1f\n', f, fa1, fa2);
end
fprintf('\n');

%% ==================== Part 3: CRT搜索恢复算法 ====================
fprintf('Part 3: CRT搜索恢复算法\n');
fprintf('------------------------\n');

% 目标：给定混叠后的测量值 (fa1, fa2)，恢复真实频率 f_true

function f_recovered = crt_search(fa1_measured, fa2_measured, Fs1, Fs2, f_range, tol)
    % CRT搜索算法恢复真实频率
    %
    % 输入:
    %   fa1_measured - Fs1采样率下的混叠频率测量值
    %   fa2_measured - Fs2采样率下的混叠频率测量值
    %   Fs1, Fs2     - 两个采样率
    %   f_range      - 搜索范围 [f_min, f_max]
    %   tol          - 匹配容差
    %
    % 输出:
    %   f_recovered  - 恢复的真实频率（若未找到则为 NaN）
    
    alias_freq_local = @(f, Fs) mod(f + Fs/2, Fs) - Fs/2;
    
    f_candidates = f_range(1):0.5:f_range(2);  % 0.5 Hz 分辨率
    f_recovered = NaN;
    min_error = inf;
    
    for f = f_candidates
        fa1_expected = alias_freq_local(f, Fs1);
        fa2_expected = alias_freq_local(f, Fs2);
        
        error = abs(fa1_expected - fa1_measured) + abs(fa2_expected - fa2_measured);
        
        if error < min_error
            min_error = error;
            if error < tol
                f_recovered = f;
            end
        end
    end
end

% 测试恢复算法
test_freqs = [1000, 3500, 5000, 7200, 10000, 12500];
f_range = [0, 25000];  % 搜索范围
tol = 1;  % 容差 1 Hz

fprintf('频率恢复测试 (无噪声):\n');
fprintf('  真实频率  |  fa1      |  fa2      |  恢复频率  |  误差\n');
fprintf('  ----------|-----------|-----------|-----------|------\n');

for f_true = test_freqs
    fa1 = alias_freq(f_true, Fs1);
    fa2 = alias_freq(f_true, Fs2);
    f_rec = crt_search(fa1, fa2, Fs1, Fs2, f_range, tol);
    error = abs(f_rec - f_true);
    fprintf('  %7d   | %+8.1f  | %+8.1f  | %9.1f  | %.1f Hz\n', ...
            f_true, fa1, fa2, f_rec, error);
end
fprintf('\n');

%% ==================== Part 4: 噪声鲁棒性分析 ====================
fprintf('Part 4: 噪声鲁棒性分析\n');
fprintf('----------------------\n');

% 在有噪声的情况下测试恢复性能
SNR_values = [30, 20, 15, 10, 5];  % dB
N_trials = 100;  % 每个SNR的试验次数
f_true = 5000;   % 测试频率

% 估计噪声方差 (假设信号幅度为1)
signal_power = 1;

results = struct();

for snr_idx = 1:length(SNR_values)
    SNR_dB = SNR_values(snr_idx);
    noise_power = signal_power / (10^(SNR_dB/10));
    noise_std = sqrt(noise_power);
    
    errors = zeros(N_trials, 1);
    success_count = 0;
    
    for trial = 1:N_trials
        % 添加噪声到混叠频率测量
        fa1_true = alias_freq(f_true, Fs1);
        fa2_true = alias_freq(f_true, Fs2);
        
        % 相位噪声转换为频率噪声 (简化模型)
        freq_noise_std = noise_std * 10;  % Hz
        fa1_noisy = fa1_true + randn() * freq_noise_std;
        fa2_noisy = fa2_true + randn() * freq_noise_std;
        
        % 恢复
        f_rec = crt_search(fa1_noisy, fa2_noisy, Fs1, Fs2, f_range, tol + freq_noise_std*2);
        
        if ~isnan(f_rec)
            errors(trial) = abs(f_rec - f_true);
            if errors(trial) < 50  % 50 Hz以内算成功
                success_count = success_count + 1;
            end
        else
            errors(trial) = NaN;
        end
    end
    
    results(snr_idx).SNR = SNR_dB;
    results(snr_idx).RMSE = sqrt(nanmean(errors.^2));
    results(snr_idx).success_rate = success_count / N_trials * 100;
    
    fprintf('SNR = %2d dB: RMSE = %6.1f Hz, 成功率 = %.1f%%\n', ...
            SNR_dB, results(snr_idx).RMSE, results(snr_idx).success_rate);
end

%% ==================== Part 5: 唯一恢复范围分析 ====================
fprintf('\nPart 5: 唯一恢复范围分析\n');
fprintf('------------------------\n');

% 分析 (fa1, fa2) 对在不同频率下是否唯一
f_analysis = 0:1:50000;  % 分析0-50kHz
fa_pairs = zeros(length(f_analysis), 2);

for i = 1:length(f_analysis)
    fa_pairs(i,1) = alias_freq(f_analysis(i), Fs1);
    fa_pairs(i,2) = alias_freq(f_analysis(i), Fs2);
end

% 四舍五入到1Hz精度，检查唯一性
fa_pairs_rounded = round(fa_pairs);
[unique_pairs, ia, ic] = unique(fa_pairs_rounded, 'rows');

fprintf('分析范围: 0 - 50 kHz\n');
fprintf('总频率点数: %d\n', length(f_analysis));
fprintf('唯一 (fa1, fa2) 对数: %d\n', size(unique_pairs, 1));

% 找到第一个重复的频率
for i = 1:length(f_analysis)
    count = sum(ic == ic(i));
    if count > 1
        % 找到所有映射到相同 (fa1, fa2) 的频率
        same_indices = find(ic == ic(i));
        if f_analysis(same_indices(1)) ~= f_analysis(i) && f_analysis(i) > 0
            fprintf('\n首个重复点:\n');
            fprintf('  频率 %d Hz 和 %d Hz 映射到相同的 (fa1, fa2) = (%.1f, %.1f)\n', ...
                    f_analysis(same_indices(1)), f_analysis(i), ...
                    fa_pairs(i,1), fa_pairs(i,2));
            fprintf('  → 唯一恢复范围约为: 0 - %d Hz\n', f_analysis(i));
            break;
        end
    end
end

% 理论唯一范围 (近似)
% 对于非整数采样率，唯一范围约为 lcm(Fs1, Fs2)
% 这里用数值方法估计
unique_range = f_analysis(i);
fprintf('  → 理论估计: Fs1 × Fs2 / gcd(Fs1, Fs2) ≈ %.0f Hz\n', ...
        Fs1 * Fs2 / gcd(round(Fs1), round(Fs2)));

%% ==================== Part 6: 结果可视化 ====================
fprintf('\nPart 6: 结果可视化\n');
fprintf('------------------\n');

figure('Name', 'Staggered PRF - 恢复性能', 'Position', [150, 150, 1000, 600]);

% 子图1: (fa1, fa2) 唯一性映射
subplot(2,2,1);
scatter(fa_pairs(1:100:end,1), fa_pairs(1:100:end,2), 5, f_analysis(1:100:end)/1000, 'filled');
colorbar;
xlabel('fa1 (Hz)');
ylabel('fa2 (Hz)');
title('(fa1, fa2) 对的唯一性 (颜色=真实频率/kHz)');
grid on;
axis equal;

% 子图2: 时序图
subplot(2,2,2);
t_plot = [0, T1, T1, T1+D, T1+D, T1+D+T2, T1+D+T2, T_frame]*1e6;
y_plot = [1, 1, 0, 0, 1, 1, 0, 0];
plot(t_plot, y_plot, 'b-', 'LineWidth', 2);
hold on;
% 标注
text(T1/2*1e6, 1.1, 'T1', 'HorizontalAlignment', 'center', 'FontSize', 11);
text((T1+D/2)*1e6, 0.1, 'D', 'HorizontalAlignment', 'center', 'FontSize', 11);
text((T1+D+T2/2)*1e6, 1.1, 'T2', 'HorizontalAlignment', 'center', 'FontSize', 11);
text((T1+D+T2+D/2)*1e6, 0.1, 'D', 'HorizontalAlignment', 'center', 'FontSize', 11);
xlabel('时间 (μs)');
ylabel('Chirp 活动');
title('Staggered PRF 时序图');
ylim([-0.2, 1.5]);
grid on;

% 子图3: SNR vs 性能
subplot(2,2,3);
yyaxis left;
plot([results.SNR], [results.RMSE], 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
ylabel('RMSE (Hz)');
yyaxis right;
plot([results.SNR], [results.success_rate], 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
ylabel('成功率 (%)');
xlabel('SNR (dB)');
title('频率恢复性能 vs SNR');
grid on;
legend('RMSE', '成功率', 'Location', 'best');

% 子图4: 单个频率的恢复演示
subplot(2,2,4);
f_demo = 7200;  % 示例频率
fa1_demo = alias_freq(f_demo, Fs1);
fa2_demo = alias_freq(f_demo, Fs2);

f_search = 0:10:15000;
fa1_search = arrayfun(@(f) alias_freq(f, Fs1), f_search);
fa2_search = arrayfun(@(f) alias_freq(f, Fs2), f_search);

plot(f_search/1000, fa1_search, 'b-', 'LineWidth', 1, 'DisplayName', 'fa1 曲线');
hold on;
plot(f_search/1000, fa2_search, 'r-', 'LineWidth', 1, 'DisplayName', 'fa2 曲线');
yline(fa1_demo, 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf('fa1测量=%.0f', fa1_demo));
yline(fa2_demo, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('fa2测量=%.0f', fa2_demo));
xline(f_demo/1000, 'g-', 'LineWidth', 2, 'DisplayName', sprintf('真实频率=%dHz', f_demo));
xlabel('频率 (kHz)');
ylabel('混叠频率 (Hz)');
title('CRT搜索原理示意');
legend('Location', 'best');
grid on;

sgtitle('Staggered PRF 抗混叠方法演示', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\n演示完成！\n');
fprintf('图形已生成，展示了Staggered PRF方法的核心原理。\n');

%% ==================== 辅助函数定义 ====================
% 注意: 在MATLAB中，局部函数需要放在脚本末尾或单独的文件中
% 这里的 crt_search 函数已经在上面定义为嵌套函数
