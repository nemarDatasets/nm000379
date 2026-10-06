function run_plv_stats()
% RUN_PLV_STATS computes phase-locking value connectivity metric at the
% subject level and saves to an intermediate file
%
% INPUTS:  no inputs  - simple wrapper to run analysis over subjects  
%
% OUTPUTS: no outputs - saves psi data to an intermediate data directory

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

subjects = {'S1', 'S2', 'S4', 'S5', 'S6'};
for s = 1:length(subjects)

    % compare plv at encoding vs. permutation (across trials) control
    res = run_subj_plv_baseline(subjects{s});

	if ~exist('..\..\scratch_data\connectivity\', 'dir')
		mkdir('..\..\scratch_data\connectivity\') 
	end

    save(['..\..\scratch_data\connectivity\' subjects{s} '_plv_base_wvis.mat'], 'res');
    
end

