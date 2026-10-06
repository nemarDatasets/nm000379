function stats = run_encoding_trialperm_stats(res, ltc_flag)
% RUN_TRIALPERM_STATS computes stats based on permutation at trial level
%
% INPUTS:  res - structure, contains results from epoched p_episode data
%
%          ltc_flag - logical, flag specifying whether analysis runs on LTC contacts
%
% OUTPUTS: stats - structure, contains p-values replicating analyses with trial-level
%                  permutation

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

for i = 1:6
    low_nr_trial_z(i).z = low_nr_trial_z(i).z - nanmean(low_nr_trial_z(i).z,3);
    low_r_trial_z(i).z = low_r_trial_z(i).z - nanmean(low_r_trial_z(i).z,3);
    
    high_nr_trial_z(i).z = high_nr_trial_z(i).z - nanmean(high_nr_trial_z(i).z,3);
    high_r_trial_z(i).z = high_r_trial_z(i).z - nanmean(high_r_trial_z(i).z,3);
end
    
n_perm = 1000;

null_diff_low = nan(length(res), size(low_nr_trial_z(i).z,3), n_perm);
for s = 1:length(res)
    % low
    data = cat(1, low_r_trial_z(s).z, low_nr_trial_z(s).z);
    data = gpuArray(squeeze(nanmean(data,2)));
    
    n1 = size(low_r_trial_z(s).z,1);
    null_diff_low(s,:,:) = generate_trial_null(data, n_perm, n1);
    
end

[~,~,~,null_st_low] = ttest(null_diff_low);
null_stat_low = squeeze(null_st_low.tstat);

null_diff_high = nan(length(res), size(high_nr_trial_z(i).z,3), n_perm);
for s = 1:length(res)
    % high
    data = cat(1, high_r_trial_z(s).z, high_nr_trial_z(s).z);
    data = gpuArray(squeeze(nanmean(data,2)));
    
    n1 = size(high_r_trial_z(s).z,1);
    null_diff_high(s,:,:) = generate_trial_null(data, n_perm, n1);
    
end

[~,~,~,null_st_high] = ttest(null_diff_high);
null_stat_high = squeeze(null_st_high.tstat);

[~,~,~,null_st_diff] = ttest(null_diff_low, null_diff_high);
null_stat_lowvhigh = squeeze(null_st_diff.tstat);

% for both

lnr = nan(6,1501);
lr = lnr;
hr = lr;
hnr = lnr;

for i = 1:6
    lnr(i,:) = nanmean(nanmean(low_nr_trial_z(i).z,2));
    lr(i,:)=  nanmean(nanmean(low_r_trial_z(i).z,2));
    
    hnr(i,:) = nanmean(nanmean(high_nr_trial_z(i).z,2));
    hr(i,:)=  nanmean(nanmean(high_r_trial_z(i).z,2));
end

if ~ltc_flag
    load('hc_theta.mat', 'high_stats', 'low_stats', 'diff_stats', 'stats');
else
    load('ltc_theta.mat', 'high_stats', 'low_stats', 'diff_stats', 'stats');
end

time = (-750:1:750)/1000;
t_idx = find(time >= -.700 & time <= .700);

% observed

low_idx = low_stats.mask == 1 & low_stats.stat < 0;
[~, ~, ~, st] = ttest(lr(:,t_idx(low_idx)), lnr(:,t_idx(low_idx)));
low_st = nanmean(st.tstat);
p_pos = (sum(squeeze(nanmean(null_stat_low(t_idx(low_idx),:),1))' > low_st)+1)/(n_perm+1);
p_neg = (sum(squeeze(nanmean(null_stat_low(t_idx(low_idx),:),1))' < low_st)+1)/(n_perm+1);
stats.p_low_pre = min([p_pos; p_neg]);

low_idx = low_stats.mask == 1 & low_stats.stat > 0;
[~, ~, ~, st] = ttest(lr(:,t_idx(low_idx)), lnr(:,t_idx(low_idx)));
low_st = nanmean(st.tstat);
p_pos = (sum(squeeze(nanmean(null_stat_low(t_idx(low_idx),:),1))' > low_st)+1)/(n_perm+1);
p_neg = (sum(squeeze(nanmean(null_stat_low(t_idx(low_idx),:),1))' < low_st)+1)/(n_perm+1);
stats.p_low_post = min([p_pos; p_neg]);

% high
high_idx = high_stats.mask == 1 & high_stats.stat < 0;
[~, ~, ~, st] = ttest(hr(:,t_idx(high_idx)), hnr(:,t_idx(high_idx)));
high_st = nanmean(st.tstat);
p_pos = (sum(squeeze(nanmean(null_stat_high(t_idx(high_idx),:),1))' > high_st)+1)/(n_perm+1);
p_neg = (sum(squeeze(nanmean(null_stat_high(t_idx(high_idx),:),1))' < high_st)+1)/(n_perm+1);
stats.p_high_pre = min([p_pos; p_neg]);

high_idx = high_stats.mask == 1 & high_stats.stat > 0;
[~, ~, ~, st] = ttest(hr(:,t_idx(high_idx)), hnr(:,t_idx(high_idx)));
high_st = nanmean(st.tstat);
p_pos = (sum(squeeze(nanmean(null_stat_high(t_idx(high_idx),:),1))' > high_st)+1)/(n_perm+1);
p_neg = (sum(squeeze(nanmean(null_stat_high(t_idx(high_idx),:),1))' < high_st)+1)/(n_perm+1);
stats.p_high_post = min([p_pos; p_neg]);

% diff
hvl= lr - lnr - hr + hnr;
null_stat_diff = null_diff_low - null_diff_high;

d_idx = stats.mask == 1 & stats.stat < 0;
[~,~,~, st] = ttest(null_stat_diff);
null_stat_diff = nanmean(squeeze(st.tstat(1,t_idx(d_idx),:)));
[~, ~, ~, st] = ttest(hvl(:, t_idx(d_idx)));
d_st = nanmean(st.tstat);

p_pos = (sum(squeeze(null_stat_diff)' > d_st)+1)/(n_perm+1);
p_neg = (sum(squeeze(null_stat_diff)' < d_st)+1)/(n_perm+1);
stats.p_diff_pre = min([p_pos; p_neg]);

hvl= lr - lnr - hr + hnr;
null_stat_diff = null_diff_low - null_diff_high;

d_idx = stats.mask == 1 & stats.stat > 0;
[~,~,~, st] = ttest(null_stat_diff);
null_stat_diff = nanmean(squeeze(st.tstat(1,t_idx(d_idx),:)));
[~, ~, ~, st] = ttest(hvl(:, t_idx(d_idx)));
d_st = nanmean(st.tstat);

p_pos = (sum(squeeze(null_stat_diff)' > d_st)+1)/(n_perm+1);
p_neg = (sum(squeeze(null_stat_diff)' < d_st)+1)/(n_perm+1);
stats.p_diff_post = min([p_pos; p_neg]);

% and both
b_diff = nanmean(cat(3, lr-lnr, hr-hnr),3);
null_stat_both = nanmean(cat(4, null_diff_low, null_diff_high),4);

b_idx = stats.mask == 1 & stats.stat < 0;
[~,~,~, st] = ttest(null_stat_both);
null_stat_both = nanmean(squeeze(st.tstat(1,t_idx(b_idx),:)));
[~, ~, ~, st] = ttest(b_diff(:, t_idx(b_idx)));
b_st = nanmean(st.tstat);

p_pos = (sum(squeeze(null_stat_both)' > b_st)+1)/(n_perm+1);
p_neg = (sum(squeeze(null_stat_both)' < b_st)+1)/(n_perm+1);
stats.p_both_pre = min([p_pos; p_neg]);

b_diff = nanmean(cat(3, lr-lnr, hr-hnr),3);
null_stat_both = nanmean(cat(4, null_diff_low, null_diff_high),4);

b_idx = stats.mask == 1 & stats.stat > 0;
[~,~,~, st] = ttest(null_stat_both);
null_stat_both = nanmean(squeeze(st.tstat(1,t_idx(b_idx),:)));
[~, ~, ~, st] = ttest(b_diff(:, t_idx(b_idx)));
b_st = nanmean(st.tstat);

p_pos = (sum(squeeze(null_stat_both)' > b_st)+1)/(n_perm+1);
p_neg = (sum(squeeze(null_stat_both)' < b_st)+1)/(n_perm+1);
stats.p_both_post = min([p_pos; p_neg]);

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
        
        post_idx = 1:1501; % exclude edges
        post = squeeze(nanmean(res(s).(fieldname)(:, e, post_idx),3));
        
        for t = 1:1501
            
            pre_idx = t;
            
            pre = squeeze(nanmean(res(s).(fieldname)(:, e, pre_idx),3));
            
            
            if all(isnan(pre)) || all(isnan(post)) || all(isnan(pre) | isnan(post))
                continue
            end
            
            stats(s).z(e,t) = mysignrank(pre, post);

        end
        
    end
end

end

function stats = run_trial_stats(res, fieldname)


for s = 1:length(res)
        
    dat = res(s).(fieldname);
    trial_m = squeeze(nanmean(dat,3));
    mu = nanmean(trial_m);
    sd = ones(size(mu)); % just mean center
    
    z_dat = nan(size(dat));
    for e = 1:size(dat,2)
        z_dat(:,e,:) = (dat(:,e,:) - mu(e))./sd(e);
    end
    
    stats(s).z = z_dat;
    
end


end

