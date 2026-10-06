function out = generate_trial_null(data, n_perm, n1)
% GENERATE_TRIAL_NULL helper function to run permutations on gpu
%
% INPUTS:  data - double, array of size observations by timepoints
%
%          n_perm - double, number of permutations to run
%
%          n1 - double, number of observations in condition 1, assumes
%               indices 1:n1 in data correspond to this condition and
%               indices n1+1:end correspond to condition 2
%
% OUTPUTS: out - double, data for generating null distribution/permutation
%                statistics at the group level

% Copyright (C) 2021  James Kragel

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

n_tp = size(data,2);
out = nan(n_tp, n_perm);

perm_count = 0;
n_subperm = 100;

while perm_count < n_perm
    
    if rem(perm_count, 1000) == 0
        disp(perm_count)
    end

    % preallocate permuted data - can play with n_subperm to optimize memory,
    % currently doing 100 permutations per computation
    b = gpuArray(zeros(size(data,1), size(data,2), n_subperm));
    
    for i = 1:n_subperm
        b(:,:,i) = data(randperm(size(data,1)), :);
    end
    
    out(:, (perm_count+1):(perm_count+n_subperm)) = gather(mean(b(1:n1,:,:),'omitnan')- ...
        mean(b(n1+1:end,:,:),'omitnan'));

    perm_count = perm_count+n_subperm;
end
