function [f_recovered, error, all_candidates] = crt_solver(f_alias, Fs, f_max, resolution)
%CRT_SOLVER Recover true frequency from aliased observations using CRT
%
%   [f_recovered, error] = crt_solver(f_alias, Fs)
%   [f_recovered, error] = crt_solver(f_alias, Fs, f_max)
%   [f_recovered, error] = crt_solver(f_alias, Fs, f_max, resolution)
%   [f_recovered, error, all_candidates] = crt_solver(...)
%
%   Inputs:
%       f_alias    - Vector of observed aliased frequencies [f1, f2, ...]
%       Fs         - Vector of sample rates [Fs1, Fs2, ...]
%       f_max      - Maximum frequency to search (default: LCM of Fs)
%       resolution - Search resolution in Hz (default: 0.1)
%
%   Outputs:
%       f_recovered    - Recovered true frequency
%       error          - Minimum error achieved
%       all_candidates - Structure with all search results
%
%   Example:
%       % True frequency 850 Hz, sampled at 97 Hz and 101 Hz
%       f_alias = [47.0, 41.5];  % Observed aliased frequencies
%       Fs = [97, 101];
%       f_true = crt_solver(f_alias, Fs);
%       % Returns f_true ≈ 850 Hz
%
%   Theory:
%       The Chinese Remainder Theorem guarantees a unique solution
%       in [0, LCM(Fs1, Fs2, ...)) when the sample rates are coprime.
%
%   Author: [Your Name]
%   Date: 2024

    % Input validation
    if length(f_alias) ~= length(Fs)
        error('f_alias and Fs must have the same length');
    end
    
    if length(f_alias) < 2
        error('Need at least 2 measurements for CRT');
    end
    
    % Default parameters
    if nargin < 4
        resolution = 0.1;
    end
    
    if nargin < 3 || isempty(f_max)
        % Calculate LCM of all sample rates
        f_max = Fs(1);
        for i = 2:length(Fs)
            f_max = lcm(round(f_max), round(Fs(i)));
        end
    end
    
    % Search grid
    f_candidates = 0:resolution:f_max;
    n_candidates = length(f_candidates);
    errors = zeros(n_candidates, 1);
    
    % For each candidate frequency
    for i = 1:n_candidates
        f_test = f_candidates(i);
        total_error = 0;
        
        % Check against each aliased observation
        for j = 1:length(f_alias)
            % Expected aliased frequency
            expected = mod(f_test, Fs(j));
            if expected > Fs(j)/2
                expected = Fs(j) - expected;
            end
            
            % Handle sign ambiguity
            observed = abs(f_alias(j));
            
            % Circular error (for frequencies near 0 or Fs/2)
            err = min(abs(expected - observed), ...
                     min(abs(expected - observed + Fs(j)/2), ...
                         abs(expected - observed - Fs(j)/2)));
            
            total_error = total_error + err^2;  % Squared error
        end
        
        errors(i) = sqrt(total_error);
    end
    
    % Find minimum
    [error, idx] = min(errors);
    f_recovered = f_candidates(idx);
    
    % Return all candidates if requested
    if nargout > 2
        all_candidates.frequencies = f_candidates;
        all_candidates.errors = errors;
        all_candidates.f_max = f_max;
        all_candidates.Fs = Fs;
    end
end
