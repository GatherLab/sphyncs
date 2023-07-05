function [p, peaks_fit] = fit_peaks(measured_peaks, n_ext,d_bead,modes,n_bead, Modestart,np)
%[p, peaks_fit] = fit_peaks(measured_peaks)
% Finds properties from peak positions in the spectra
%   measured_peaks is a list of spectral peak locations
%   p is a struct where p{1} is the list of mode numbers, p{2} is n_ext and
%   p{3} is radius of bead

% Assumptions:
% 1. There are even number of peaks correspoding to TM/TE pairs
if mod(length(measured_peaks), 2) == 0
    num_modes = length(measured_peaks) / 2;
else
    error("number of peaks must be even")
end

% 2. Bead refractive index is:
n_bead;

% 3. Search bounds are:
lb = [modes(1), n_ext(1), d_bead(1)/2];
ub = [modes(2), n_ext(2), d_bead(2)/2];

%4.3 mode numbers or 2?

if Modestart ==1;% 2 mode pairs
    
    % Objective is MSE of peak locations
    estimated_peaks = @(x) spectral_peaks(round(x(1):x(1) + num_modes - 1), ... % modes
        x(2), ... % n_ext
        x(3), ... % r
        n_bead,Modestart);
    objective = @(x) sum((measured_peaks - estimated_peaks(x)).^2);
    
end


if Modestart == 0; %1 extra mode number for the 1st TE
    estimated_peaks = @(x) spectral_peaks(round(x(1):x(1) + num_modes-1 ), ... % modes
        x(2), ... % n_ext
        x(3), ... % r
        n_bead,Modestart);
    objective = @(x) sum((measured_peaks - estimated_peaks(x)).^2);
end

% Global constrained minimisation required for this problem
ms = MultiStart('StartPointsToRun', 'bounds', 'UseParallel', true, 'Display', 'off');
opts = optimoptions(@fmincon, 'Algorithm', 'sqp');

% Starting parameters - need a list of points for global minimisation
n_ext_guess = @(x) rand([x, 1]) * (ub(2) - lb(2)) + lb(2);
r_guess = @(x) rand([x, 1]) * (ub(3) - lb(3)) + lb(3);

mode_search = (lb(1):ub(1))';
nr = length(mode_search);
mode_guess = @(x) mode_search(mod(1:x, nr) + 1);

 % number of search points np
start_points = cat(2, mode_guess(np), n_ext_guess(np), r_guess(np));
x_rng = start_points(round(np / 2), :);

problem = createOptimProblem('fmincon', 'objective', objective, ...
    'x0', x_rng, 'lb', lb, 'ub', ub, 'options', opts);
rs = CustomStartPointSet(start_points);

[x, ~] = run(ms, problem, rs);
x(1) = round(x(1));


% % Refinement - second pass improves mode guesses a bit
% mode_search = (x(1) + (-3:3))';
% nr = length(mode_search);
% mode_guess = @(x) mode_search(mod(1:x, nr) + 1);
% 
% np = 51; % number of search points
% start_points = cat(2, mode_guess(np), n_ext_guess(np), r_guess(np));
% x_rng = start_points(round(np / 2), :);
% 
% problem = createOptimProblem('fmincon', 'objective', objective, ...
%     'x0', x_rng, 'lb', lb, 'ub', ub, 'options', opts);
% rs = CustomStartPointSet(start_points);
% 
% [x, ~] = run(ms, problem, rs);
% x(1) = round(x(1));


% Output struct
p = {};
p{1} = round(x(1):x(1) + num_modes - 1);
p{2} = x(2);
p{3} = x(3);
peaks_fit = estimated_peaks(x);

end

