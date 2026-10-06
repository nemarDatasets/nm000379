function run_revisitation_healthy_controls()
% RUN_REVISITATION_HEALTHY_CONTROLS computes effects of revisitation on subsequent
% fixation-sequence replay in three independent datasets
%
% INPUTS:  no inputs, loads necessary data from disk
%
% OUTPUTS: no outputs, generates figures          

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

% find revisitation fixations and their relation to subsequent fixation-sequence replay
[kk_res, kk_null_res] = compute_kaspkoni_revisitation();
[figrim_res, figrim_null_res] = compute_figrim_revisitation();

% organize data into format for plotting
revisit_res_mem1 = run_kk_revisitation_memory(kk_res, kk_null_res, 'mem1');
revisit_res_mem2 = run_kk_revisitation_memory(kk_res, kk_null_res, 'mem2');
revisit_res_figrim = run_figrim_revisitation_memory(figrim_res, figrim_null_res);

% overall effect
dat.m_viewed = cat(1, revisit_res_mem1.m_viewed, ...
    revisit_res_mem2.m_viewed, ...
    revisit_res_figrim.m_viewed);

dat.m_reins = cat(1, revisit_res_mem1.m_reins, ...
    revisit_res_mem2.m_reins, ...
    revisit_res_figrim.m_reins);

dat.m_rev_reins = cat(1, revisit_res_mem1.m_rev_reins, ...
    revisit_res_mem2.m_rev_reins, ...
    revisit_res_figrim.m_rev_reins);

dataset = [ones(1,size(revisit_res_mem1.m_viewed,1)) ...
     2*ones(1,size(revisit_res_mem2.m_viewed,1)) ...
     3*ones(1,size(revisit_res_figrim.m_viewed,1))]';

% and generate figures/stats
plot_revisitation_healthy_controls(dat, dataset, 1);

% stimulus-driven effects and comparison
dat.null_viewed = cat(1, revisit_res_mem1.null_viewed, ...
    revisit_res_mem2.null_viewed, ...
    revisit_res_figrim.null_viewed);

dat.null_reins = cat(1, revisit_res_mem1.null_reins, ...
    revisit_res_mem2.null_reins, ...
    revisit_res_figrim.null_reins);

dat.null_rev_reins = cat(1, revisit_res_mem1.null_rev_reins, ...
    revisit_res_mem2.null_rev_reins, ...
    revisit_res_figrim.null_rev_reins);

dataset = [ones(1,size(revisit_res_mem1.null_viewed,1)) ...
     2*ones(1,size(revisit_res_mem2.null_viewed,1)) ...
     3*ones(1,size(revisit_res_figrim.null_viewed,1))]';

% and generate figures/stats
plot_revisitation_healthy_controls_stimulus(dat, dataset, 1);


end

