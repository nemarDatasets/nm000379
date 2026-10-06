function stats = run_group_psi_stats(data, c1_idx, c2_idx)
% RUN_GROUP_PSI_STATS compares phase slope index across the group
%
% INPUTS: data - fieldtrip structure containing plv data per condition
%		  
%		  c1_idx - logical, index of condition 1 data
%
%		  c2_idx - logical, index of condition 2 data
%
% OUTPUTS: stats - fieldtrip structure containing permutations stats

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

cfg                  = [];
cfg.parameter        = 'avg';
cfg.method           = 'montecarlo';

cfg.design           = nan(2, sum(c1_idx)+sum(c2_idx));
cfg.design(1,c1_idx) = 1;
cfg.design(1,c2_idx) = 2;
cfg.design(2,:)      = [1:sum(c1_idx), 1:sum(c2_idx)]; 
cfg.ivar             = 1;
cfg.uvar             = 2;

cfg.numrandomization = 'all'; 
cfg.correctm         = 'fdr';
cfg.correcttail      = 'alpha';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.channel          = 'all';
cfg.latency          = 'all'; % 

stats = ft_timelockstatistics(cfg, data);

end

