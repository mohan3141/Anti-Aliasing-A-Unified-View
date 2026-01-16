function f_alias = alias_freq(f_true, Fs, mode)
%ALIAS_FREQ Calculate the aliased frequency given true frequency and sample rate
%
%   f_alias = alias_freq(f_true, Fs)
%   f_alias = alias_freq(f_true, Fs, mode)
%
%   Inputs:
%       f_true - True frequency (can be vector/matrix)
%       Fs     - Sample rate (scalar)
%       mode   - 'centered' (default): output in [-Fs/2, Fs/2]
%                'positive': output in [0, Fs/2]
%
%   Output:
%       f_alias - Aliased frequency
%
%   Examples:
%       % 850 Hz sampled at 97 Hz
%       alias_freq(850, 97)           % Returns -47 Hz (centered)
%       alias_freq(850, 97, 'positive')  % Returns 47 Hz
%
%       % Vectorized
%       alias_freq([100, 200, 300], 97)  % Returns [3, -9, 9] Hz
%
%   Theory:
%       Aliasing occurs when f > Fs/2 (Nyquist frequency)
%       The aliased frequency is: f_alias = ((f + Fs/2) mod Fs) - Fs/2
%
%   Author: [Your Name]
%   Date: 2024

    if nargin < 3
        mode = 'centered';
    end
    
    % Fold to [-Fs/2, Fs/2]
    f_alias = mod(f_true + Fs/2, Fs) - Fs/2;
    
    % If positive mode requested, take absolute value
    if strcmpi(mode, 'positive')
        f_alias = abs(f_alias);
    end
end
