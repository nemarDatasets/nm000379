function stats = run_power_stats(data, c1_idx, c2_idx, f_avg, t_avg)
% RUN_POWER_STATS runs permutation statistics between conditions
%
% INPUTS: data - fieldtrip data structure containing time-frequency data
%				 to compare
%
%		  c1_idx - logical index for condition 1 in the data structure
%
%		  c2_idx - logical index for condition 2 in the data structure
%
%		  f_avg  - logical flag specifying if analysis averages over frequencies
%
%	      t_avg  - logical flag specifying if analysis averages over time
%
%
% OUTPUTS: stats - fieldtrip structure containing permutations stats
%
% NOTE: this function does not correct for multiple comparisons

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
cfg.parameter        = 'powspctrm';
cfg.neighbours       = [];
cfg.method           = 'montecarlo';

cfg.design           = nan(2, sum(c1_idx)+sum(c2_idx));
cfg.design(1,c1_idx) = 1;
cfg.design(1,c2_idx) = 2;
cfg.design(2,:)      = [1:sum(c1_idx), 1:sum(c2_idx)]; 
cfg.ivar             = 1;
cfg.uvar             = 2;

cfg.numrandomization = 'all'; 
cfg.correctm         = 'none';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.channel          = 'all';
cfg.latency          = 'all';

if f_avg
   cfg.avgoverfreq = 'yes';
end

if t_avg
   cfg.avgovertime = 'yes';
end

stats = ft_freqstatistics(cfg, data);

end

