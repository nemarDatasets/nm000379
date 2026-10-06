function stats = run_group_plv_stats(data, c1_idx, c2_idx, f_avg, t_avg)
% RUN_GROUP_PLV_STATS computes group level plv stats
%
% INPUTS: data - fieldtrip structure containing plv data per condition
%		  
%		  c1_idx - logical, index of condition 1 data
%
%		  c2_idx - logical, index of condition 2 data
%
%		  f_avg  - logical, whether to average over frequencies before stats
%
%		  t_avg  - logical, whether to average over time before stats
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


if ~exist('f_avg', 'var')
    f_avg = false;
end

if ~exist('t_avg', 'var')
    t_avg = false;
end

cfg                  = [];
cfg.parameter        = 'plv';
cfg.method           = 'montecarlo';

cfg.design           = nan(2, sum(c1_idx)+sum(c2_idx));
cfg.design(1,c1_idx) = 1;
cfg.design(1,c2_idx) = 2;
cfg.design(2,:)      = [1:sum(c1_idx), 1:sum(c2_idx)]; 
cfg.ivar             = 1;
cfg.uvar             = 2;

cfg.numrandomization = 'all'; 
cfg.correctm         = 'fdr';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.channel          = 'all';
cfg.latency          = 'all';
cfg.tail			 = 1; % greater than zero 

if f_avg
   cfg.avgoverfreq = 'yes';
end

if t_avg
   cfg.avgovertime = 'yes';
end

stats = ft_freqstatistics(cfg, data);

end

