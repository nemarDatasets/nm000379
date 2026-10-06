function [stats_high_revisit, stats_high_other, stats_pre, stats_post] = run_retrieval_stats(res, enc_res)
% RUN_RETRIEVAL_STATS runs stats comparing p_episode for retrieval conditions
%
% INPUTS:  res - structure, contains epoched saccade, fixation, and trial level data from
%                run_bosc_epoch_ret.m
%
%          enc_res - structure, contains epoched saccade, fixation, and trial level data from
%                run_bosc_epoch.m
%
% OUTPUTS: stats_high_revisit - structure, stats structure for revisit vs.
%                               zero
%
%          stats_high_other - structure, stats structure for other vs. zero
%
%          stats - structure, stats structure for revisit vs. other fixations for theta
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

% For null results with low theta, you can run this commented code
% [~, low_f_subj_z] = get_stats(res, 'forward');
% [~, low_b_subj_z] = get_stats(res, 'backward');
% [~, low_nr_subj_z] = get_stats(res, 'no_replay');
% [~, low_n_subj_z] = get_stats(res, 'novel');
% [~, low_new_subj_z] = get_stats(res, 'new');

% retrieval, high theta
[~, high_f_subj_z] = get_stats(res, 'forward_high');
[~, high_b_subj_z] = get_stats(res, 'backward_high');
[~, high_nr_subj_z] = get_stats(res, 'no_replay_high');

% encoding

% For low theta, you can execute the commented code below
% [~, enc_low_nr_subj_z] = get_stats(enc_res, 'not_repeated_low');
% [~, enc_low_r_subj_z] = get_stats(enc_res, 'repeated_low');

[~, enc_high_nr_subj_z] = get_stats(enc_res, 'not_repeated_high');
[~, enc_high_r_subj_z] = get_stats(enc_res, 'repeated_high');

% enc_all_r_subj_z = nanmean(cat(3, enc_low_r_subj_z, enc_high_r_subj_z),3);
% enc_all_nr_subj_z = nanmean(cat(3, enc_low_nr_subj_z, enc_high_nr_subj_z),3);

% average lag one
% low_abs1_subj_z = nanmean(cat(3, low_f_subj_z, low_b_subj_z),3);
% high_abs1_subj_z = nanmean(cat(3, high_f_subj_z, high_b_subj_z),3);

% average old
% low_old_subj_z = nanmean(cat(3, low_f_subj_z, low_b_subj_z, low_nr_subj_z),3);
high_old_subj_z = nanmean(cat(3, high_f_subj_z, high_b_subj_z, high_nr_subj_z),3);

fig('FontSize', 14, ...
'Font', 'Arial', ...
'border','on', ...
'units', 'centimeters', ...
'width', 12, ...
'height', 12);

old_color = [70 70 70]/255;

ax = subplot(1, 2, 1);
stats_high_revisit = run_comparison(high_old_subj_z, enc_high_r_subj_z, [old_color; 0 1 0], {'Old','Revisit'}, ax);

axis square
ax = subplot(1, 2, 2);
stats_high_other = run_comparison(high_old_subj_z, enc_high_nr_subj_z, [old_color; .7 .7 .7], {'Old','Other'}, ax);
axis square

[stats_pre, stats_post] = compute_es(stats, high_old_subj_z, enc_high_nr_subj_z);

set(gcf, 'Position', [7.0644   15.4781   23.2304    7.4348]);

end

function [stats_pre, stats_post] = compute_es(stats, dat1, dat2)


idx_pre = stats.time < 0 & stats.mask == 1;
idx_post = stats.time > 0 & stats.mask == 1;

time = (-750:1:750)/1000;
t_idx = find(time >= -.700 & time <= .500);

if any(idx_pre)
    stats_pre = mes(nanmean(dat1(:,t_idx(idx_pre)),2), ...
        nanmean(dat2(:,t_idx(idx_pre)),2), ...
        'hedgesg', ...
        'isDep', 1);
else
    idx_pre = abs(stats.stat) == max(abs(stats.stat(stats.time < 0)));
    
    stats_pre = mes(nanmean(dat1(:,t_idx(idx_pre)),2), ...
        nanmean(dat2(:,t_idx(idx_pre)),2), ...
        'hedgesg', ...
        'isDep', 1);    
end

if any(idx_post)
    stats_post = mes(nanmean(dat1(:,t_idx(idx_post)),2), ...
        nanmean(dat2(:,t_idx(idx_post)),2), ...
        'hedgesg', ...
        'isDep', 1);
else
    idx_post = abs(stats.stat) == max(abs(stats.stat(stats.time > 0)));
    
    stats_post = mes(nanmean(dat1(:,t_idx(idx_post)),2), ...
        nanmean(dat2(:,t_idx(idx_post)),2), ...
        'hedgesg', ...
        'isDep', 1); 
end

end

function stats = run_comparison(dat1, dat2, colors, labels, ax)

toi = [-.7 .5];
ylim = [-1.5 1.5];

stats = run_ft_group_stats(dat1, ...
    dat2, ...
    toi);

if ~exist('ax','var')
    [hl1, hp1, hl2, hp2] = plot_ft_res(dat1, dat2, stats, ylim, colors, labels);
else
    [hl1, hp1, hl2, hp2] = plot_ft_res(dat1, dat2, stats, ylim, colors, labels, ax);
    set(gcf,'Position',[7.0644   15.4781    7.1438    7.4348]);
end

end

function stat = run_ft_group_stats(dat1, dat2, toi)

% stats, before the fixation
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

function [elec_z, subj_z] = get_stats(res, fieldname)

stats = run_stats(res, fieldname);
elec_z = cat(1,stats(:).z);

for s = 1:length(res)

    tmp = stats(s).z;   
    subj_z(s,:) = nanmean(tmp,1);
    
end

end

function stats = run_stats(res, fieldname)

fieldstr = fieldnames(res);

for s = 1:length(res)
    
    stats(s).p = nan(size(res(s).(fieldstr{1}),2), size(res(s).(fieldstr{1}),3));
    stats(s).z = stats(s).p;
    for e = 1:size(res(s).(fieldstr{1}),2)
        
        dat = res(s).(fieldname);
        if length(size(dat)) == 2
            dat = reshape(dat, 1, size(dat,1), size(dat,2));
        end
        
        post_idx = 1:1501; % exclude edges
        post = squeeze(nanmean(dat(:, e, post_idx),3));
        
        for t = 1:1501
            
            pre_idx = t;
            
            pre = squeeze(nanmean(dat(:, e, pre_idx),3));
            
            
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


