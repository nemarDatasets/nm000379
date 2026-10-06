function stats = run_retrieval_trialperm_stats(res, enc_res, ltc_flag)
% RUN_RETRIEVAL_TRIALPERM_STATS computes stats based on permutation at trial level
%
% INPUTS:  res - structure, contains retrieval results from epoched p_episode data
%
%          enc_res - structure, contains encoding results from epoched p_episode data
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

[~, ~, low_f_trial_z] = get_stats(res, 'forward');
[~, ~, low_b_trial_z] = get_stats(res, 'backward');
[~, ~, low_nr_trial_z] = get_stats(res, 'no_replay');

[~, ~, high_f_trial_z] = get_stats(res, 'forward_high');
[~, ~, high_b_trial_z] = get_stats(res, 'backward_high');
[~, ~, high_nr_trial_z] = get_stats(res, 'no_replay_high');

% encoding
[~, ~, enc_low_nr_trial_z] = get_stats(enc_res, 'not_repeated_low');
[~, ~, enc_low_r_trial_z] = get_stats(enc_res, 'repeated_low');

[~, ~, enc_high_nr_trial_z] = get_stats(enc_res, 'not_repeated_high');
[~, ~, enc_high_r_trial_z] = get_stats(enc_res, 'repeated_high');

% and now for some stats

% compare 'other' vs. retrieval

n_perm = 1000;

% old (test) vs. revisit (study): ovr 
null_ovr_high = nan(length(res), size(enc_high_r_trial_z(1).z,3), n_perm);
null_ovo_high = nan(length(res), size(enc_high_r_trial_z(1).z,3), n_perm);
for s = 1:length(res)
    % high -- old vs. revisit
    high_old_z = squeeze(nanmean(cat(1, high_f_trial_z(s).z, ...
                                        high_b_trial_z(s).z, ...
                                        high_nr_trial_z(s).z),2));
    
    data = cat(1, high_old_z, squeeze(nanmean(enc_high_r_trial_z(s).z,2)));
    data = gpuArray(data);
    
    n1 = size(high_old_z,1);
    null_ovr_high(s,:,:) = generate_trial_null(data, n_perm, n1);
    
    % high -- old vs. other
    data = cat(1, high_old_z, squeeze(nanmean(enc_high_nr_trial_z(s).z,2)));
    data = gpuArray(data);
    
    n1 = size(high_old_z,1);
    null_ovo_high(s,:,:) = generate_trial_null(data, n_perm, n1);
end

% old (test) vs. revisit (study): ovr 
null_ovr_low = nan(length(res), size(enc_low_r_trial_z(1).z,3), n_perm);
null_ovo_low = nan(length(res), size(enc_low_r_trial_z(1).z,3), n_perm);
for s = 1:length(res)
    % low -- old vs. revisit
    low_old_z = squeeze(nanmean(cat(1, low_f_trial_z(s).z, ...
                                       low_b_trial_z(s).z, ...
                                       low_nr_trial_z(s).z),2));
    
    data = cat(1, low_old_z, squeeze(nanmean(enc_low_r_trial_z(s).z,2)));
    data = gpuArray(data);
    
    n1 = size(low_old_z,1);
    null_ovr_low(s,:,:) = generate_trial_null(data, n_perm, n1);
    
    % low -- old vs. other
    data = cat(1, low_old_z, squeeze(nanmean(enc_low_nr_trial_z(s).z,2)));
    data = gpuArray(data);
    
    n1 = size(low_old_z,1);
    null_ovo_low(s,:,:) = generate_trial_null(data, n_perm, n1);
end

% high
[~,~,~,null_st_high] = ttest(null_ovr_high);
null_stat_ovr_high = squeeze(null_st_high.tstat);

[~,~,~,null_st_high] = ttest(null_ovo_high);
null_stat_ovo_high = squeeze(null_st_high.tstat);

% low
[~,~,~,null_st_low] = ttest(null_ovo_low);
null_stat_ovo_low = squeeze(null_st_low.tstat);

[~,~,~,null_st_low] = ttest(null_ovr_low);
null_stat_ovr_low = squeeze(null_st_low.tstat);

% high
for s = 1:6
    h_old_z(s,:) = squeeze(nanmean(nanmean(cat(1, high_f_trial_z(s).z, ...
        high_b_trial_z(s).z, ...
        high_nr_trial_z(s).z),2)));
      
    h_encr_z(s,:) = squeeze(nanmean(nanmean(enc_high_r_trial_z(s).z, 2)));
    h_encnr_z(s,:) = squeeze(nanmean(nanmean(enc_high_nr_trial_z(s).z, 2)));
end

% low
for s = 1:6
    l_old_z(s,:) = squeeze(nanmean(nanmean(cat(1, low_f_trial_z(s).z, ...
        low_b_trial_z(s).z, ...
        low_nr_trial_z(s).z),2)));
    
    l_encr_z(s,:) = squeeze(nanmean(nanmean(enc_low_r_trial_z(s).z, 2)));
    l_encnr_z(s,:) = squeeze(nanmean(nanmean(enc_low_nr_trial_z(s).z, 2)));
end


time = (-750:1:750)/1000;
t_idx = time >= -.700 & time <= .500;

pre_idx = time < 0;
post_idx = time > 0;

pre_idx = pre_idx(t_idx);
post_idx = post_idx(t_idx);

t_idx = find(t_idx);

% old vs. revisit

h_idx = stats_high_revisit.mask == 1 & pre_idx;

if any(h_idx)
    [~, ~, ~, st] = ttest(h_old_z(:,t_idx(h_idx)), h_encr_z(:,t_idx(h_idx)));
    h_st = nanmean(st.tstat);
    p_pos = (sum(squeeze(nanmean(null_stat_ovr_high(t_idx(h_idx),:),1))' > h_st)+1)/(n_perm+1);
    p_neg = (sum(squeeze(nanmean(null_stat_ovr_high(t_idx(h_idx),:),1))' < h_st)+1)/(n_perm+1);
    p_h_pre = min([p_pos; p_neg]);
else
    h_st = nan;
    p_h_pre = nan;
end

stats.high_revisit.t_pre = h_st;
stats.high_revisit.p_pre = p_h_pre;

h_idx = stats_high_revisit.mask == 1 & post_idx;
if any(h_idx)
    [~, ~, ~, st] = ttest(h_old_z(:,t_idx(h_idx)), h_encr_z(:,t_idx(h_idx)));
    h_st = nanmean(st.tstat);
    p_pos = (sum(squeeze(nanmean(null_stat_ovr_high(t_idx(h_idx),:),1))' > h_st)+1)/(n_perm+1);
    p_neg = (sum(squeeze(nanmean(null_stat_ovr_high(t_idx(h_idx),:),1))' < h_st)+1)/(n_perm+1);
    p_h_post = min([p_pos; p_neg]);
else
    h_st = nan;
    p_h_post = nan;
end

stats.high_revisit.t_post = h_st;
stats.high_revisit.p_post = p_h_post;

% old vs. other

h_idx = stats_high_other.mask == 1 & pre_idx;

if any(h_idx)
    [~, ~, ~, st] = ttest(h_old_z(:,t_idx(h_idx)), h_encnr_z(:,t_idx(h_idx)));
    h_st = nanmean(st.tstat);
    p_pos = (sum(squeeze(nanmean(null_stat_ovo_high(t_idx(h_idx),:),1))' > h_st)+1)/(n_perm+1);
    p_neg = (sum(squeeze(nanmean(null_stat_ovo_high(t_idx(h_idx),:),1))' < h_st)+1)/(n_perm+1);
    p_h_pre = min([p_pos; p_neg]);
else
    h_st = nan;
    p_h_pre = nan;
end

stats.high_other.t_pre = h_st;
stats.high_other.p_pre = p_h_pre;

h_idx = stats_high_other.mask == 1 & post_idx;

[~, ~, ~, st] = ttest(h_old_z(:,t_idx(h_idx)), h_encnr_z(:,t_idx(h_idx)));
h_st = nanmean(st.tstat);
p_pos = (sum(squeeze(nanmean(null_stat_ovo_high(t_idx(h_idx),:),1))' > h_st)+1)/(n_perm+1);
p_neg = (sum(squeeze(nanmean(null_stat_ovo_high(t_idx(h_idx),:),1))' < h_st)+1)/(n_perm+1);
p_h_post = min([p_pos; p_neg]);

stats.high_other.t_post = h_st;
stats.high_other.p_post = p_h_post;

end

function [elec_z, subj_z, trial_z] = get_stats(res, fieldname)

stats = run_stats(res, fieldname);
elec_z = cat(1,stats(:).z);

for s = 1:length(res)
    tmp = stats(s).z;    
    subj_z(s,:) = nanmean(tmp,1);
    trial_z(s) = run_trial_stats(res(s), fieldname);
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

function stats = run_trial_stats(res, fieldname)

for s = 1:length(res)
        
    dat = res(s).(fieldname);
    trial_m = squeeze(nanmean(dat,3));
    mu = nanmean(trial_m); sd = ones(size(mu)); % mean center only
    
    z_dat = nan(size(dat));
    for e = 1:size(dat,2)
        z_dat(:,e,:) = (dat(:,e,:) - mu(e))./sd(e);
    end
    
    stats(s).z = z_dat;
    
end


end
