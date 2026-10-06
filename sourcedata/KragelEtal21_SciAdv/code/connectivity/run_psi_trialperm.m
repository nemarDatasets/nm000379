function run_psi_trialperm()
% RUN_PSI_TRIALPERM computes phase-slope index and permutation stats at the
% fixation level prior to group analysis
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
    
    res = run_subj_psi_trialperm(subjects{s});

    if ~exist('..\..\scratch_data\connectivity\', 'dir')
        mkdir('..\..\scratch_data\connectivity\') 
    end

    save(['..\..\scratch_data\connectivity\' subjects{s} '_psi_wvis_trialperm.mat'], 'res');    
    
end

