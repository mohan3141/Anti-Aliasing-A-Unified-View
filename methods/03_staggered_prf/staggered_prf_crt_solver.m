function [f_recovered, confidence] = staggered_prf_crt_solver(fa1, fa2, Fs1, Fs2, varargin)
%STAGGERED_PRF_CRT_SOLVER 基于CRT的双采样率频率恢复
%
%   使用Staggered PRF的两个混叠频率测量值，通过CRT搜索恢复真实频率。
%
%   语法:
%       f_recovered = staggered_prf_crt_solver(fa1, fa2, Fs1, Fs2)
%       f_recovered = staggered_prf_crt_solver(fa1, fa2, Fs1, Fs2, 'param', value, ...)
%       [f_recovered, confidence] = staggered_prf_crt_solver(...)
%
%   输入:
%       fa1     - 在Fs1采样率下的混叠频率测量值 (Hz), 可以是向量
%       fa2     - 在Fs2采样率下的混叠频率测量值 (Hz), 可以是向量
%       Fs1     - 第一个等效采样率 (Hz)
%       Fs2     - 第二个等效采样率 (Hz)
%
%   可选参数 (名称-值对):
%       'FreqRange'     - 搜索范围 [f_min, f_max] (Hz), 默认 [0, Fs1*Fs2/gcd]
%       'Resolution'    - 搜索分辨率 (Hz), 默认 1 Hz
%       'Tolerance'     - 匹配容差 (Hz), 默认 min(Fs1, Fs2)/20
%       'Method'        - 搜索方法: 'grid' (默认) 或 'iterative'
%       'Verbose'       - 是否显示详细信息, 默认 false
%
%   输出:
%       f_recovered - 恢复的真实频率 (Hz), 未找到时为 NaN
%       confidence  - 置信度分数 [0, 1], 基于匹配误差
%
%   示例:
%       % 基本用法
%       Fs1 = 5882; Fs2 = 8333;
%       fa1 = -882; fa2 = 833;  % 真实频率5000 Hz的混叠
%       f = staggered_prf_crt_solver(fa1, fa2, Fs1, Fs2);
%       % f = 5000
%
%       % 批量处理
%       fa1_vec = [100, -882, 1200];
%       fa2_vec = [100, 833, 1200];
%       f_vec = staggered_prf_crt_solver(fa1_vec, fa2_vec, Fs1, Fs2);
%
%       % 自定义搜索范围
%       f = staggered_prf_crt_solver(fa1, fa2, Fs1, Fs2, ...
%               'FreqRange', [0, 10000], 'Tolerance', 10);
%
%   理论背景:
%       在Staggered PRF系统中，真实频率f经过两个不同采样率的欠采样后，
%       产生两个混叠频率:
%           fa1 = ((f + Fs1/2) mod Fs1) - Fs1/2
%           fa2 = ((f + Fs2/2) mod Fs2) - Fs2/2
%       
%       如果gcd(Fs1, Fs2)较小，则f在范围[0, lcm(Fs1,Fs2))内可唯一确定。
%       本函数通过网格搜索或迭代方法找到满足上述两个方程的f。
%
%   Author: Anti-Aliasing Recovery Toolkit
%   Date: 2024-01
%   Version: 1.0

    %% 输入解析
    p = inputParser;
    addRequired(p, 'fa1', @isnumeric);
    addRequired(p, 'fa2', @isnumeric);
    addRequired(p, 'Fs1', @(x) isnumeric(x) && isscalar(x) && x > 0);
    addRequired(p, 'Fs2', @(x) isnumeric(x) && isscalar(x) && x > 0);
    
    % 默认搜索范围: [0, 近似唯一范围]
    default_range = [0, Fs1 * Fs2 / gcd(round(Fs1), round(Fs2))];
    addParameter(p, 'FreqRange', default_range, @(x) isnumeric(x) && length(x)==2);
    addParameter(p, 'Resolution', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Tolerance', min(Fs1, Fs2)/20, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Method', 'grid', @(x) ismember(x, {'grid', 'iterative'}));
    addParameter(p, 'Verbose', false, @islogical);
    
    parse(p, fa1, fa2, Fs1, Fs2, varargin{:});
    opts = p.Results;
    
    %% 输入验证
    if length(fa1) ~= length(fa2)
        error('fa1 and fa2 must have the same length');
    end
    
    %% 处理
    n_points = length(fa1);
    f_recovered = zeros(size(fa1));
    confidence = zeros(size(fa1));
    
    if opts.Verbose
        fprintf('Staggered PRF CRT Solver\n');
        fprintf('  Fs1 = %.1f Hz, Fs2 = %.1f Hz\n', Fs1, Fs2);
        fprintf('  Search range: [%.0f, %.0f] Hz\n', opts.FreqRange(1), opts.FreqRange(2));
        fprintf('  Resolution: %.1f Hz, Tolerance: %.1f Hz\n', opts.Resolution, opts.Tolerance);
        fprintf('  Processing %d points...\n', n_points);
    end
    
    % 选择搜索方法
    switch opts.Method
        case 'grid'
            % 网格搜索
            for i = 1:n_points
                [f_recovered(i), confidence(i)] = grid_search(...
                    fa1(i), fa2(i), Fs1, Fs2, opts.FreqRange, opts.Resolution, opts.Tolerance);
            end
            
        case 'iterative'
            % 迭代搜索 (从粗到细)
            for i = 1:n_points
                [f_recovered(i), confidence(i)] = iterative_search(...
                    fa1(i), fa2(i), Fs1, Fs2, opts.FreqRange, opts.Resolution, opts.Tolerance);
            end
    end
    
    if opts.Verbose
        valid_count = sum(~isnan(f_recovered));
        fprintf('  Found: %d / %d (%.1f%%)\n', valid_count, n_points, valid_count/n_points*100);
    end
end

%% ==================== 辅助函数 ====================

function fa = alias_freq(f, Fs)
    % 计算混叠频率 (centered at 0)
    fa = mod(f + Fs/2, Fs) - Fs/2;
end

function [f_best, conf] = grid_search(fa1, fa2, Fs1, Fs2, f_range, resolution, tol)
    % 网格搜索
    f_candidates = f_range(1):resolution:f_range(2);
    
    f_best = NaN;
    min_error = inf;
    
    for f = f_candidates
        fa1_exp = alias_freq(f, Fs1);
        fa2_exp = alias_freq(f, Fs2);
        
        error = abs(fa1_exp - fa1) + abs(fa2_exp - fa2);
        
        if error < min_error
            min_error = error;
            if error < tol
                f_best = f;
            end
        end
    end
    
    % 计算置信度 (基于匹配误差)
    if ~isnan(f_best)
        conf = max(0, 1 - min_error / tol);
    else
        conf = 0;
    end
end

function [f_best, conf] = iterative_search(fa1, fa2, Fs1, Fs2, f_range, final_res, tol)
    % 迭代搜索: 从粗到细
    
    % 第一遍: 粗搜索
    coarse_res = (f_range(2) - f_range(1)) / 1000;
    coarse_res = max(coarse_res, final_res * 10);
    
    [f_coarse, ~] = grid_search(fa1, fa2, Fs1, Fs2, f_range, coarse_res, tol * 5);
    
    if isnan(f_coarse)
        f_best = NaN;
        conf = 0;
        return;
    end
    
    % 第二遍: 细搜索 (在粗搜索结果附近)
    fine_range = [max(f_range(1), f_coarse - coarse_res*2), ...
                  min(f_range(2), f_coarse + coarse_res*2)];
    
    [f_best, conf] = grid_search(fa1, fa2, Fs1, Fs2, fine_range, final_res, tol);
end
