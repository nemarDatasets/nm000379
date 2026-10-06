function [low_stats, high_stats, stats, diff_stats, stats_pre, stats_post] = run_encoding_stats(res)
% RUN_ENCODING_STATS runs stats comparing p_episode for encoding conditions
%
% INPUTS:  res - structure, contains epoched saccade, fixation, and trial level data from
%                run_bosc_epoch.m
%
% OUTPUTS: low_stats - structure, stats structure for revisit vs. other fixations for low theta
%
%          high_stats - structure, stats structure for revisit vs. other fixations for high theta
%
%          stats - structure, stats structure for revisit vs. other fixations for theta
%
%          diff_stats - structure, stats structure interactions of high vs. low theta for
%                       revisitation vs. other fixations
%
%          stats_pre - structure, effect size estimates for significant differences preceding
%                      fixations of interest
%
%          stats_post - structure, effect size estimates for significant differences following
%                      fixations of interest

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

[low_nr_elec_z, low_nr_subj_z, low_nr_trial_z] = get_stats(res, 'not_repeated_low');
[low_r_elec_z, low_r_subj_z, low_r_trial_z] = get_stats(res, 'repeated_low');

[high_nr_elec_z, high_nr_subj_z, high_nr_trial_z] = get_stats(res, 'not_repeated_high');
[high_r_elec_z, high_r_subj_z, high_r_trial_z] = get_stats(res, 'repeated_high');

% low theta

low_stats = run_ft_group_stats(low_r_subj_z, ...
    low_nr_subj_z, ...
    [-.7 .7]);

[hl, hp] = plot_ft_res(low_r_subj_z, low_nr_subj_z, low_stats, [-1.5 1.5] );
set(gcf,'Position',[7.0644   15.4781    7.1438    7.4348]);

% high theta

high_stats = run_ft_group_stats(high_r_subj_z, ...
    high_nr_subj_z, ...
    [-.7 .7]);

[hl, hp] = plot_ft_res(high_r_subj_z, high_nr_subj_z, high_stats, [-2 2] );
set(gcf,'Position',[7.0644   15.4781    7.1438    7.4348]);

% high and low
for s = 1:size(high_r_subj_z,1)
    r(s,:) = nanmean(cat(1,high_r_subj_z(s,:),low_r_subj_z(s,:)),1);
    o(s,:) = nanmean(cat(1,high_nr_subj_z(s,:),low_nr_subj_z(s,:)),1);
end

stats = run_ft_group_stats(r, ...
    o, ...
    [-.70 .70]);

[hl, hp] = plot_ft_res(r, o, stats, [-1.5 1.5] );
set(gcf,'Position',[7.0644   15.4781    7.1438    7.4348]);

%
idx_pre = stats.time < 0 & stats.mask == 1;
idx_post = stats.time > 0 & stats.mask == 1;

time = (-750:1:750)/1000;
t_idx = find(time >= -.700 & time <= .700);

if any(idx_pre)
    stats_pre = mes(nanmean(r(:,t_idx(idx_pre)),2), ...
        nanmean(o(:,t_idx(idx_pre)),2), ...
        'hedgesg', ...
        'isDep', 1);
else
    stats_pre = [];
end

if any(idx_post)
    stats_post = mes(nanmean(r(:,t_idx(idx_post)),2), ...
        nanmean(o(:,t_idx(idx_post)),2), ...
        'hedgesg', ...
        'isDep', 1);
else
    stats_post = [];
end

diff_stats = run_ft_group_stats(low_r_subj_z - low_nr_subj_z, ...
    high_r_subj_z - low_nr_subj_z, ...
    [-.70 .70]);

colors(1,:) = [0 102 204]/255;
colors(2,:) = [209 119 35]/255;
labels = {'Low Theta','High Theta'};


[hl, hp] = plot_ft_res(low_r_subj_z - low_nr_subj_z, ...
    high_r_subj_z - low_nr_subj_z, ...
    diff_stats, [-1.5 1.5], ...
    colors, labels);

set(gcf,'Position',[7.0644   15.4781    7.1438    7.4348]);

end

function stat = run_ft_group_stats(dat1, dat2, toi)

cond1.dimord = 'subj_chan_time';
cond1.avg = reshape(dat1, size(dat1,1), 1, size(dat1,2));
cond1.time = (-750:1:750)/1000;
cond1.label = {'x'};

cond2.dimord = 'subj_chan_time';
cond2.avg = reshape(dat2, size(dat2,1), 1, size(dat2,2));
cond2.time = (-750:1:750)/1000;
cond2.label = {'x'};

cfg = [];

cfg.latency = toi;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'fdr';

cfg.tail = 0;
cfg.correctail = 'alpha';
cfg.alpha = 0.05;               % alpha level of the permutation test
cfg.numrandomization = 'all';      % number of draws from the permutation distribution

design = zeros(2,2*size(dat1,1));
design(1,1:size(dat1,1)) = 1;
design(1,size(dat1,1)+1:(2*size(dat1,1)))= 2;
design(2,:) = [1:size(dat1,1) 1:size(dat1,1)];

cfg.design = design;             % design matrix
cfg.ivar  = 1;
cfg.uvar  = 2;

stat = ft_timelockstatistics(cfg, cond1, cond2);

end

function [elec_z, subj_z, trial_z] = get_stats(res, fieldname)

stats = run_stats(res, fieldname);
elec_z = cat(1,stats(:).z);

for s = 1:length(res)

    tmp = stats(s).z;
    subj_z(s,:) = nanmean(tmp);
    trial_z(s) = run_trial_stats(res(s), fieldname);

end

end

function stats = run_stats(res, fieldname)


for s = 1:length(res)
    
    stats(s).p = nan(size(res(s).initial_low,2), size(res(s).initial_low,3));
    stats(s).z = stats(s).p;
    for e = 1:size(res(s).initial_low,2)
        
        post_idx = 1:1501;
        post = squeeze(nanmean(res(s).(fieldname)(:, e, post_idx),3));
        
        for t = 1:1501
            
            pre_idx = t;
            
            pre = squeeze(nanmean(res(s).(fieldname)(:, e, pre_idx),3));
            
            
            if all(isnan(pre)) || all(isnan(post)) || all(isnan(pre) | isnan(post))
                continue
            end
            
            [stats(s).p(e,t), h, st] = signrank(pre, post, 'method', 'approximate');

            if isfield(st, 'zval')
                stats(s).z(e,t) = st.zval;
            else
                stats(s).z(e,t) = nan;
            end
        end
        
    end
end



end

function stats = run_trial_stats(res, fieldname)


for s = 1:length(res)
    
    dat = res(s).(fieldname);
    trial_m = squeeze(nanmean(dat,3));
    mu = nanmean(trial_m); sd = sqrt(mu.*(1-mu)./size(dat,1));
        
    z_dat = nan(size(dat));
        for e = 1:size(dat,2)
            z_dat(:,e,:) = (dat(:,e,:) - mu(e))./sd(e);
        end
     
        stats(s).z = z_dat;
        
end

end

