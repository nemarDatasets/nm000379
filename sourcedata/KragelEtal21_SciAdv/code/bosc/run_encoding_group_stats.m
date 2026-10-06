function [low_stats, high_stats, stats, diff_stats] = run_encoding_group_stats(vis_flag, ltc_flag)
% RUN_ENCODING_GROUP_STATS runs group stats
%
% INPUTS:  vis_flag - logical, flag specifying whether analysis runs on DAN\VN contacts
%
%          ltc_flag - logical, flag specifying whether analysis runs on LTC contacts
%
% OUTPUTS: low_stats - structure, stats structure for revisit vs. other fixations for low theta
%
%          high_stats - structure, stats structure for revisit vs. other fixations for high theta
%
%          stats - structure, stats structure for revisit vs. other fixations for theta
%
%          diff_stats - structure, stats structure interactions of high vs. low theta for
%                       revisitation vs. other fixations

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

if ~exist('vis_flag','var')
    vis_flag = false;
end

if ~exist('ltc_flag','var')
    ltc_flag = false;
end

subjects = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6'};

t_idx = 1:1501;

enc_subview = [];
enc_reins = [];

rec_viewed = [];
rec_reins = [];

base_dir = '..\..\scratch_data\bosc\';
if ltc_flag
    load([base_dir 'ltc_theta_res.mat'],'theta_res');
elseif vis_flag
    load([base_dir 'vis_theta_res.mat'],'theta_res');
else
    load([base_dir 'theta_res.mat'],'theta_res');
end

for s = 1:length(subjects)
    
    [enc_res, ~, events] = run_bosc_epoch(subjects{s}, vis_flag, ltc_flag, theta_res(s));
    
    enc_idx = strcmp({events.phase}, 'encode');
    
    fix_idx = strcmp({events.type}, 'fixation');
    
    enc_fix_events = events(fix_idx & enc_idx);

    enc_fix_info = ev2fix(enc_fix_events);
    
    % set up condition labels for selecting data
    sel_labels = ones(size(enc_fix_info.repeated));
    sel_labels(enc_fix_info.repeated == 1) = 1;
    sel_labels(enc_fix_info.repeated == 0) = 2;
    sel_labels(isnan(enc_fix_info.repeated)) = 3;
    
    res(s).initial_low = get_condition(enc_res, 'p_ep_fix_low', sel_labels==1, false, t_idx);
    res(s).repeated_low = get_condition(enc_res, 'p_ep_fix_low', sel_labels==2, false, t_idx);
    res(s).not_repeated_low = get_condition(enc_res, 'p_ep_fix_low', sel_labels==3, false, t_idx);
    
    res(s).initial_high = get_condition(enc_res, 'p_ep_fix_high', sel_labels==1, false, t_idx);
    res(s).repeated_high = get_condition(enc_res, 'p_ep_fix_high', sel_labels==2, false, t_idx);
    res(s).not_repeated_high = get_condition(enc_res, 'p_ep_fix_high', sel_labels==3, false, t_idx);

    
end

[low_stats, high_stats, stats, diff_stats] = run_encoding_stats(res);

end

function m = get_condition(res, osc_name, idx, mean_flag, t_idx)

target_samples = 3751; % 1 kHz
n_samples = size(res.(osc_name), 3);

dat = res.(osc_name);
dat = dat(idx==1,:,:);
if target_samples ~= n_samples
    dat_rs = nan(size(dat,1), size(dat,2), target_samples);
    for c = 1:size(dat,2)
        dat_rs(:,c,:) = resample(squeeze(dat(:,c,:))', target_samples, n_samples)';
    end
    dat_rs(dat_rs < 0) = 0; dat_rs(dat_rs > 1) = 1; % bounds from resampling
    dat = dat_rs;
end

if mean_flag
    m = squeeze(nanmean(dat(:,:,t_idx)));
else
    m = squeeze(dat(:,:,t_idx));
end

end
